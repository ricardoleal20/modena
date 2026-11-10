//! Full velocity-Verlet experiment. JSON schema version 1 is intended for an
//! external campaign driver, which freezes binaries and supplies SHA-256,
//! process CPU/RSS/load and AB/BA repetition provenance. No timing acceptance
//! decision is made here. Example:
//! cargo run --release -p modena-bench --example layout_step --
//!   --mode cluster4 --n 864 --threads 1 --steps 10 --rebuild-policy natural --seed 17
use std::{
    path::PathBuf,
    sync::atomic::{AtomicU64, Ordering},
    time::Instant,
};

use clap::{Parser, ValueEnum};
use modena_bench::packed_pair::{validate_domain, Layout, PackedForceModel};
use modena_core::{
    Energies, ForceModel, Integrator, ModenaError, NeighborList, Result, SimBox, System, Vec3,
    VelocityVerlet, Virial,
};
use serde_json::json;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Interleaved,
    InterleavedSerial,
    Soa,
    Aosoa4,
    Cluster4,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
enum RebuildPolicy {
    Natural,
    NoRebuild,
    Forced,
}

#[derive(Debug, Parser)]
#[command(about = "B03 scalar packed layout full-step experiment; no SIMD claim")]
struct Args {
    #[arg(long, value_enum)]
    mode: Mode,
    #[arg(long)]
    n: usize,
    /// Zero uses Rayon's actual default count (reported in output).
    #[arg(long, default_value_t = 0)]
    threads: usize,
    #[arg(long)]
    steps: usize,
    #[arg(long, value_enum, default_value_t = RebuildPolicy::Natural)]
    rebuild_policy: RebuildPolicy,
    #[arg(long, default_value_t = 17)]
    seed: u64,
    #[arg(long, default_value_t = 0.0001)]
    dt: f64,
    /// Optional full initial/final JSON state for independent auditing; written
    /// outside elapsed time. Parent directory must already exist.
    #[arg(long)]
    snapshot: Option<PathBuf>,
}

struct CountingModel {
    inner: PackedForceModel,
    forced: bool,
    calls: AtomicU64,
    neighbor_pair_visits: AtomicU64,
}

impl ForceModel for CountingModel {
    fn conserves_total_momentum(&self, _: &System) -> bool {
        true
    }
    fn neighbor_cutoff(&self) -> Option<f64> {
        Some(2.5)
    }
    fn compute(&self, s: &mut System, v: &mut Virial, nl: &mut NeighborList) -> Result<f64> {
        // The integrator calls here AFTER drift, so forced construction belongs
        // to this force evaluation and is inside the full-step timer.
        if self.forced {
            nl.force_rebuild(s);
        }
        let energy = self.inner.compute(s, v, nl)?;
        self.calls.fetch_add(1, Ordering::Relaxed);
        self.neighbor_pair_visits
            .fetch_add(nl.n_pairs() as u64, Ordering::Relaxed);
        Ok(energy)
    }
}

// Fixed, explicit SplitMix64 stream avoids changing RNG dependencies. It is not
// a thermostat or a claim of a canonical velocity distribution.
fn sample(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *seed;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^= z >> 31;
    ((z >> 11) as f64) * (1.0 / ((1_u64 << 53) as f64)) - 0.5
}

fn initial_state(n: usize, mut seed: u64) -> System {
    // Fill an FCC lattice at a=1.6796, truncating only the final cells for
    // non-4*cube sizes (notably N10000). This exact workload is reported.
    let cells = ((n as f64 / 4.).cbrt().ceil() as usize).max(1);
    let a = 1.6796;
    let mut s = System::new(SimBox::cubic(a * cells as f64));
    s.resize(n);
    let basis = [[0., 0., 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]];
    for i in 0..n {
        let cell = i / 4;
        let xyz = [cell / (cells * cells), (cell / cells) % cells, cell % cells];
        let mut p = [0.; 3];
        let mut v = [0.; 3];
        for k in 0..3 {
            p[k] = a * (xyz[k] as f64 + basis[i % 4][k]) + 0.002 * sample(&mut seed);
            v[k] = 0.1 * sample(&mut seed);
        }
        s.set_pos(i, Vec3(p));
        s.set_vel(i, Vec3(v));
    }
    s.wrap_positions();
    let com = s.total_momentum() * (1. / n as f64);
    for i in 0..n {
        s.set_vel(i, s.vel(i) - com);
    }
    s
}

/// Reproducibility checksum, explicitly noncryptographic FNV-1a over LE bits.
fn checksum(s: &System) -> String {
    let mut hash = 0xcbf29ce484222325_u64;
    for value in s
        .box_size
        .size
        .0
        .iter()
        .chain(&s.positions)
        .chain(&s.velocities)
        .chain(&s.forces)
        .chain(&s.masses)
    {
        for byte in value.to_bits().to_le_bytes() {
            hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
        }
    }
    format!("{hash:016x}")
}

fn state_snapshot(
    s: &System,
    energies: Option<Energies>,
    steps: usize,
    nl: Option<&NeighborList>,
) -> serde_json::Value {
    json!({
        "n":s.n(),"logical_ids":(0..s.n()).collect::<Vec<_>>(),
        "cell_vectors":s.box_size.basis().map(|v| v.0),"periodic":[true,true,true],
        "masses":s.masses,"species":s.species,"charges":s.charges,"atomic_numbers":s.atomic_numbers,
        "positions":s.positions,"velocities":s.velocities,"forces":s.forces,
        "energies":energies,"steps_completed":steps,"checksum_fnv1a64":checksum(s),
        "neighbor_pairs":nl.map(|n| n.pairs()),
        "pair_scales":nl.map(|n| n.pair_scales().iter().map(|p| json!({"i":p.i,"j":p.j,"lj":p.lj,"coulomb":p.coulomb})).collect::<Vec<_>>())
    })
}

fn run(args: &Args) -> Result<serde_json::Value> {
    if args.n == 0 || args.steps == 0 || !args.dt.is_finite() || args.dt <= 0. {
        return Err(ModenaError::Input(
            "n and steps must be positive; dt must be finite positive".into(),
        ));
    }
    let setup = Instant::now();
    let mut s = initial_state(args.n, args.seed);
    if s.box_size.size.0.iter().any(|&l| l < 5.6) {
        return Err(ModenaError::Input(
            "experiment box must span at least twice cutoff+skin (use n>=109)".into(),
        ));
    }
    let initial_checksum = checksum(&s);
    let initial_input = args
        .snapshot
        .as_ref()
        .map(|_| state_snapshot(&s, None, 0, None));
    let inner = match args.mode {
        Mode::Interleaved => PackedForceModel::new(Layout::Interleaved),
        Mode::InterleavedSerial => PackedForceModel::serial_baseline(),
        Mode::Soa => PackedForceModel::new(Layout::Soa),
        Mode::Aosoa4 => PackedForceModel::new(Layout::Aosoa4),
        Mode::Cluster4 => PackedForceModel::new(Layout::Cluster4),
    };
    let model = CountingModel {
        inner,
        forced: args.rebuild_policy == RebuildPolicy::Forced,
        calls: AtomicU64::new(0),
        neighbor_pair_visits: AtomicU64::new(0),
    };
    let mut nl = NeighborList::with_skin(2.5, 0.3);
    validate_domain(&s, &nl)?;
    let mut virial = Virial::default();
    let pe = model.compute(&mut s, &mut virial, &mut nl)?;
    let mut energies = Energies::from(&s, pe, &virial);
    let initial_force_ready = args
        .snapshot
        .as_ref()
        .map(|_| state_snapshot(&s, Some(energies), 0, Some(&nl)));
    let initial_pairs = nl.n_pairs();
    let setup_rebuilds = nl.rebuild_count();
    let setup_stats = model.inner.stats();
    let setup_calls = model.calls.load(Ordering::Relaxed);
    let setup_pair_visits = model.neighbor_pair_visits.load(Ordering::Relaxed);
    let setup_seconds = setup.elapsed().as_secs_f64();
    let mut integrator = VelocityVerlet::new();
    let start = Instant::now();
    for _ in 0..args.steps {
        energies = integrator.step(&mut s, &model, &mut nl, args.dt)?;
    }
    let elapsed_seconds = start.elapsed().as_secs_f64();
    let rebuilds = nl.rebuild_count() - setup_rebuilds;
    if args.rebuild_policy == RebuildPolicy::NoRebuild && rebuilds != 0 {
        return Err(ModenaError::Input(format!(
            "no-rebuild contract failed: {rebuilds} actual rebuilds; reduce steps/dt"
        )));
    }
    if !energies.potential.is_finite()
        || !energies.virial.is_finite()
        || !energies.total.is_finite()
        || s.positions
            .iter()
            .chain(&s.velocities)
            .chain(&s.forces)
            .any(|x| !x.is_finite())
    {
        return Err(ModenaError::Input("nonfinite final state/E/F/W".into()));
    }
    let stats = model.inner.stats();
    // Final active count is OUTSIDE elapsed time, identically for every mode.
    let final_active_pairs = nl
        .pairs()
        .iter()
        .filter(|&&(i, j)| {
            let r2 = s.box_size.min_image(s.pos(i) - s.pos(j)).norm2();
            (1e-14..6.25).contains(&r2) && nl.pair_scale(i, j).0 != 0.
        })
        .count();
    let baseline = matches!(args.mode, Mode::Interleaved | Mode::InterleavedSerial);
    let production_auto = matches!(args.mode, Mode::Interleaved);
    let row = json!({
        "schema_version":1,"experiment":"b03-packed-layout","mode":args.mode.to_possible_value().unwrap().get_name(),
        "pair_execution":if production_auto { "production-auto" } else { "serial-scalar" },
        "comparison_class":if production_auto && rayon::current_num_threads()>1 { "production-threading-control" } else { "serial-layout" },
        "simd_claim":false,"n":args.n,"requested_threads":args.threads,"actual_threads":rayon::current_num_threads(),
        "steps_requested":args.steps,"steps_completed":model.calls.load(Ordering::Relaxed)-setup_calls,
        "seed":args.seed,"dt":args.dt,"rebuild_policy":args.rebuild_policy.to_possible_value().unwrap().get_name(),
        "workload":"truncated-fcc-a1.6796-splitmix64-jitter0.002-velocity0.1-com-removed-v1",
        "box_lengths":s.box_size.size.0,"potential":"shifted-force-lj-sigma1-epsilon1-mass1-cutoff2.5-skin0.3",
        "setup_seconds":setup_seconds,"elapsed_seconds":elapsed_seconds,
        "setup_rebuilds":setup_rebuilds,"rebuilds":rebuilds,"initial_neighbor_pairs":initial_pairs,
        "final_neighbor_pairs":nl.n_pairs(),"neighbor_pair_visits":model.neighbor_pair_visits.load(Ordering::Relaxed)-setup_pair_visits,
        "final_active_pairs":final_active_pairs,
        "evaluated_pairs":if baseline { None } else { Some(stats.evaluated_pairs-setup_stats.evaluated_pairs) },
        "setup_cluster_reconstructions":setup_stats.cluster_reconstructions,
        "cluster_reconstructions":stats.cluster_reconstructions-setup_stats.cluster_reconstructions,
        "retained_packed_bytes":stats.retained_packed_bytes,
        "memory_scope":"packed Vec capacities only; excludes core neighbor/baseline workspace and transient reconstruction map; process RSS required externally",
        "initial_checksum_fnv1a64":initial_checksum,"final_checksum_fnv1a64":checksum(&s),
        "energies":energies,"force_max_abs":s.forces.iter().fold(0_f64,|m,f| m.max(f.abs())),"finite_efw":true
    });
    if let Some(path) = &args.snapshot {
        let snapshot = json!({"schema_version":1,"row":row,"initial_input":initial_input,
            "initial_force_ready":initial_force_ready,
            "final":state_snapshot(&s,Some(energies),args.steps,Some(&nl))});
        let file = std::fs::File::create(path)
            .map_err(|e| ModenaError::Input(format!("cannot create snapshot: {e}")))?;
        serde_json::to_writer(file, &snapshot)
            .map_err(|e| ModenaError::Input(format!("cannot write snapshot: {e}")))?;
    }
    Ok(row)
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let pool_start = Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()?;
    let pool_setup_seconds = pool_start.elapsed().as_secs_f64();
    let mut row = pool.install(|| run(&args))?;
    row["pool_setup_seconds"] = json!(pool_setup_seconds);
    println!("{}", serde_json::to_string(&row)?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_rejects_unknown_modes_potentials_and_policies() {
        for extra in [
            vec!["--mode", "simd"],
            vec!["--potential", "morse"],
            vec!["--rebuild-policy", "never"],
        ] {
            let mut args = vec!["layout_step", "--mode", "soa", "--n", "864", "--steps", "2"];
            args.extend(extra);
            assert!(Args::try_parse_from(args).is_err());
        }
    }

    #[test]
    fn snapshot_contains_real_states_and_neighbor_pairs_for_independent_audit() {
        let path =
            std::env::temp_dir().join(format!("modena-b03-snapshot-{}.json", std::process::id()));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let row = pool
            .install(|| {
                run(&Args {
                    mode: Mode::Cluster4,
                    n: 256,
                    threads: 1,
                    steps: 3,
                    rebuild_policy: RebuildPolicy::Forced,
                    seed: 17,
                    dt: 0.0001,
                    snapshot: Some(path.clone()),
                })
            })
            .unwrap();
        let snapshot: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        std::fs::remove_file(path).unwrap();
        assert_eq!(snapshot["row"], row);
        for key in ["initial_input", "initial_force_ready", "final"] {
            let state = &snapshot[key];
            for field in ["positions", "velocities", "forces"] {
                assert_eq!(state[field].as_array().unwrap().len(), 768);
            }
            for field in ["masses", "species", "logical_ids"] {
                assert_eq!(state[field].as_array().unwrap().len(), 256);
            }
            assert_eq!(state["cell_vectors"].as_array().unwrap().len(), 3);
        }
        assert_eq!(
            snapshot["initial_input"]["checksum_fnv1a64"],
            row["initial_checksum_fnv1a64"]
        );
        assert_eq!(
            snapshot["final"]["checksum_fnv1a64"],
            row["final_checksum_fnv1a64"]
        );
        assert_eq!(snapshot["final"]["steps_completed"], 3);
        assert_eq!(
            snapshot["final"]["neighbor_pairs"]
                .as_array()
                .unwrap()
                .len(),
            row["final_neighbor_pairs"].as_u64().unwrap() as usize
        );
        assert_eq!(snapshot["final"]["energies"], row["energies"]);
        assert_ne!(
            snapshot["initial_input"]["positions"],
            snapshot["final"]["positions"]
        );
    }

    #[test]
    fn no_rebuild_rejects_observed_rebuild_instead_of_mislabeling_row() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let error = pool
            .install(|| {
                run(&Args {
                    mode: Mode::Cluster4,
                    n: 256,
                    threads: 1,
                    steps: 2,
                    rebuild_policy: RebuildPolicy::NoRebuild,
                    seed: 17,
                    dt: 1.,
                    snapshot: None,
                })
            })
            .unwrap_err();
        assert!(
            error.to_string().contains("no-rebuild contract failed"),
            "{error}"
        );
    }

    #[test]
    fn full_step_rows_agree_and_forced_rebuilds_are_counted() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        pool.install(|| {
            for policy in [
                RebuildPolicy::Natural,
                RebuildPolicy::NoRebuild,
                RebuildPolicy::Forced,
            ] {
                let mut reference: Option<serde_json::Value> = None;
                for mode in [
                    Mode::Interleaved,
                    Mode::InterleavedSerial,
                    Mode::Soa,
                    Mode::Aosoa4,
                    Mode::Cluster4,
                ] {
                    let row = run(&Args {
                        mode,
                        n: 256,
                        threads: 1,
                        steps: 3,
                        rebuild_policy: policy,
                        seed: 17,
                        dt: 0.0001,
                        snapshot: None,
                    })
                    .unwrap();
                    assert_eq!(row["steps_completed"], 3);
                    assert_eq!(
                        row["rebuilds"],
                        if policy == RebuildPolicy::Forced {
                            3
                        } else {
                            0
                        }
                    );
                    assert_eq!(row["finite_efw"], true);
                    if let Some(expected) = &reference {
                        assert_eq!(
                            row["initial_checksum_fnv1a64"],
                            expected["initial_checksum_fnv1a64"]
                        );
                        assert_eq!(
                            row["neighbor_pair_visits"],
                            expected["neighbor_pair_visits"]
                        );
                        for field in ["potential", "virial", "total"] {
                            let a = row["energies"][field].as_f64().unwrap();
                            let b = expected["energies"][field].as_f64().unwrap();
                            assert!((a - b).abs() <= 1e-8 + 1e-12 * b.abs());
                        }
                    } else {
                        reference = Some(row);
                    }
                }
            }
        });
    }
}
