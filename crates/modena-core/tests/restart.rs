//! Restart/checkpoint gates: exact state round-trip, seamless NVE resume
//! across a checkpoint seam, and the water equilibrate-then-measure workflow
//! (NVT phase writes a checkpoint, NVE phase restores it and continues the
//! step counter).

use std::path::{Path, PathBuf};

use modena_core::io::restart;
use modena_core::{InputConfig, Simulation};

/// Minimal tempdir (avoid pulling a dev-dependency just for this).
fn tempdir() -> PathBuf {
    static SEQUENCE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let mut p = std::env::temp_dir();
    p.push(format!(
        "modena_restart_{}_{}_{}",
        std::process::id(),
        SEQUENCE.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&p).unwrap();
    p
}

/// Write `toml` into `dir` and build a `Simulation` with `dir` as base dir
/// (relative log/restart paths resolve there, mirroring the CLI).
fn sim_from_toml(dir: &Path, toml: &str) -> Simulation {
    let path = dir.join(format!(
        "deck_{}.toml",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&path, toml).unwrap();
    let config = InputConfig::from_path(&path).unwrap();
    Simulation::from_config(&config, dir).unwrap()
}

fn build_err_from_toml(dir: &Path, toml: &str) -> String {
    let path = dir.join(format!(
        "deck_{}.toml",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::write(&path, toml).unwrap();
    let config = InputConfig::from_path(&path).unwrap();
    match Simulation::from_config(&config, dir) {
        Ok(_) => panic!("expected a build error"),
        Err(e) => e.to_string(),
    }
}

/// One `f64` per logged row: `HashMap`-free `(step, total)` pairs from a CSV
/// log written by the runner.
fn log_rows(dir: &Path, file: &str) -> Vec<(usize, f64)> {
    let text = std::fs::read_to_string(dir.join(file)).unwrap();
    text.lines()
        .skip(1) // header
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let cols: Vec<&str> = l.split(',').collect();
            (
                cols[0].parse::<usize>().unwrap(),
                cols[4].parse::<f64>().unwrap(),
            )
        })
        .collect()
}

/// A checkpoint written after step K (directly, and by the run loop's
/// end-of-run write) restores bit-identical positions, velocities, box, and
/// step counter.
#[test]
fn checkpoint_roundtrip_is_bit_identical() {
    let dir = tempdir();
    let toml = r#"
units = "reduced"

[structure]
type = "fcc"
lattice_constant = 1.5444
replicas = [3, 3, 3]
element = "Ar"

[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

[integrator]
type = "velocity_verlet"
dt = 0.005

[run]
steps = 137
initial_temperature = 1.2
seed = 42
restart_file = "run.chk"
restart_every = 500
"#;
    let mut sim = sim_from_toml(&dir, toml);
    sim.run().unwrap();
    assert_eq!(sim.current_step(), 137);

    // Direct save after the run (the "call the save after step K" path)…
    let manual = dir.join("manual.chk");
    sim.write_checkpoint(&manual).unwrap();
    // …and the run loop's own end-of-run write.
    let run_written = dir.join("run.chk");
    for path in [&manual, &run_written] {
        let cp = restart::load(path).unwrap();
        let sys = sim.system();
        assert_eq!(cp.atoms, sys.n(), "atom count");
        assert_eq!(cp.step, 137, "step counter");
        assert_eq!(cp.positions, sys.positions, "positions (bits)");
        assert_eq!(cp.velocities, sys.velocities, "velocities (bits)");
        assert_eq!(cp.box_size, sys.box_size, "box (bits)");
        // The two files carry the same state byte-for-byte.
        assert_eq!(
            std::fs::read_to_string(&manual).unwrap(),
            std::fs::read_to_string(path).unwrap(),
            "manual and run-loop checkpoints differ"
        );
    }
    std::fs::remove_dir_all(&dir).ok();
}

/// An NVE run split at a checkpoint seam is the same run: the seam row is
/// continuous, every logged step matches the one-go run, and the combined
/// drift equals the one-go drift.
#[test]
fn nve_seam_continuity_and_combined_drift() {
    let dir = tempdir();
    const BASE: &str = r#"
units = "reduced"

[structure]
type = "fcc"
lattice_constant = 1.5444
replicas = [3, 3, 3]
element = "Ar"

[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5

[integrator]
type = "velocity_verlet"
dt = 0.005

[run]
initial_temperature = 1.2
seed = 42
log_every = 10
restart_every = 500
"#;
    // Phase A: L = 120 NVE steps, checkpoint at the end.
    let a = format!("{BASE}\nsteps = 120\nrestart_file = \"a.chk\"\nlog_file = \"a.csv\"\n");
    let mut sim_a = sim_from_toml(&dir, &a);
    let sum_a = sim_a.run().unwrap();
    assert_eq!(sum_a.steps, 120);

    // Phase B: restore + M = 80 more NVE steps.
    let b = format!("{BASE}\nsteps = 80\nrestore_file = \"a.chk\"\nlog_file = \"b.csv\"\n");
    let mut sim_b = sim_from_toml(&dir, &b);
    assert_eq!(
        sim_b.current_step(),
        120,
        "restored simulation must continue from the checkpoint step"
    );
    let sum_b = sim_b.run().unwrap();
    assert_eq!(sum_b.steps, 80, "RunSummary counts this invocation's steps");

    // One-go reference: L + M = 200 steps from the same initial state.
    let c = format!("{BASE}\nsteps = 200\nlog_file = \"c.csv\"\n");
    let mut sim_c = sim_from_toml(&dir, &c);
    let sum_c = sim_c.run().unwrap();

    // Seam continuity: A's last logged row and B's first logged row are both
    // the checkpoint state — the totals must agree to ~1e-12 relative.
    let rows_a = log_rows(&dir, "a.csv");
    let rows_b = log_rows(&dir, "b.csv");
    let (seam_a, e_a) = *rows_a.last().unwrap();
    let (seam_b, e_b) = rows_b[0];
    assert_eq!(seam_a, 120, "A's last row");
    assert_eq!(
        seam_b, 120,
        "B's first row must be labeled with the seam step"
    );
    let seam_rel = (e_a - e_b).abs() / e_a.abs().max(1e-30);
    assert!(
        seam_rel < 1e-12,
        "seam total discontinuous: {e_a} vs {e_b} (rel {seam_rel:.3e})"
    );

    // Every logged step of the split run (A's rows + B's rows past the seam)
    // matches the one-go run at the same step.
    let rows_c = log_rows(&dir, "c.csv");
    let by_step_c: std::collections::HashMap<usize, f64> = rows_c.iter().copied().collect();
    let mut split = rows_a.clone();
    split.extend(rows_b.iter().skip(1).copied());
    for (step, total) in &split {
        let one_go = by_step_c
            .get(step)
            .unwrap_or_else(|| panic!("one-go log missing step {step}"));
        let rel = (total - one_go).abs() / one_go.abs().max(1e-30);
        assert!(
            rel < 1e-12,
            "step {step}: split total {total} vs one-go {one_go} (rel {rel:.3e})"
        );
    }

    // Combined drift over L + M equals the one-go drift (tight), and the
    // final full-precision energies agree.
    let drift = |rows: &[(usize, f64)]| {
        let e0 = rows[0].1;
        rows.iter()
            .map(|(_, e)| (e - e0).abs())
            .fold(0.0_f64, f64::max)
            / e0.abs().max(1e-30)
    };
    let drift_split = drift(&split);
    let drift_one = drift(&rows_c);
    assert!(
        (drift_split - drift_one).abs() < 1e-12,
        "combined drift {drift_split:.3e} vs one-go {drift_one:.3e}"
    );
    let rel_final = (sum_b.final_energies.total - sum_c.final_energies.total).abs()
        / sum_c.final_energies.total.abs().max(1e-30);
    assert!(
        rel_final < 1e-12,
        "final total: split {} vs one-go {} (rel {rel_final:.3e})",
        sum_b.final_energies.total,
        sum_c.final_energies.total
    );
    std::fs::remove_dir_all(&dir).ok();
}

/// The water two-phase workflow: NVT equilibration writes a checkpoint, an
/// NVE run restores it and continues the step counter; the NVE segment is
/// conservative (loose bound — the initial-state relaxation caveat is
/// documented in docs/WATER_MODEL.md), and a wrong-atom-count restore is a
/// config error.
#[test]
fn water_two_phase_nvt_then_nve_restore() {
    let dir = tempdir();
    const WATER: &str = r#"
units = "real"

[structure]
type = "water"
n_waters = 64
density = 0.997
seed = 7

[model]
type = "lennard_jones"
cutoff = 9.0
[[model.types]]
sigma = 3.166
epsilon = 0.1553
[[model.types]]
sigma = 3.166
epsilon = 0.0

[electrostatics]
method = "spme"
kappa = 0.34
real_cutoff = 9.0
kspace_cutoff = 4.0

[integrator]
type = "velocity_verlet"
dt = 0.0005

[run]
initial_temperature = 298.15
seed = 87287
log_every = 50
restart_every = 500
"#;
    // Phase 1: NVT (stochastic velocity rescaling at 298.15 K). 1000 steps
    // at dt = 0.0005 ps with strong coupling (tau = 0.005 ps — the
    // cold-lattice collapse dumps heat far faster than tau ~ 0.1 ps can
    // drain, which leaves the melt at thousands of K): measured to leave
    // the subsequent 500-step NVE segment at ~1.5e-2 drift. 300 steps
    // (the deck's smoke budget) still leaves ~6.6e-2 — see
    // bench/validation/water/modena_spce_nvt.toml.
    let nvt = format!(
        "{WATER}\nsteps = 1000\nrestart_file = \"nvt.chk\"\nlog_file = \"nvt.csv\"\n\
         \n[thermostat]\ntype = \"velocity_rescaling\"\ntarget_temperature = 298.15\ntau = 0.005\n"
    );
    let mut sim1 = sim_from_toml(&dir, &nvt);
    let sum1 = sim1.run().unwrap();
    assert_eq!(sum1.steps, 1000);
    assert!(sim1.current_step() == 1000);

    // Phase 2: NVE restore + 500 steps (no thermostat block).
    let nve = format!("{WATER}\nsteps = 500\nrestart_mode = \"state_only\"\nrestore_file = \"nvt.chk\"\nlog_file = \"nve.csv\"\n");
    let mut sim2 = sim_from_toml(&dir, &nve);
    assert_eq!(
        sim2.current_step(),
        1000,
        "restored water run must continue from step 1000"
    );
    let sum2 = sim2.run().unwrap();
    assert_eq!(sum2.steps, 500);

    // The restored run's log labels continue logically: first row is the
    // seam (step 1000), last row is step 1500.
    let rows = log_rows(&dir, "nve.csv");
    assert_eq!(rows[0].0, 1000, "first NVE log row");
    assert_eq!(rows.last().unwrap().0, 1500, "last NVE log row");

    // The restored state really is phase 1's end state (not a fresh draw):
    // the seam row of the NVE log equals phase 1's last logged total.
    let rows_nvt = log_rows(&dir, "nvt.csv");
    let seam_rel = (rows[0].1 - rows_nvt.last().unwrap().1).abs()
        / rows_nvt.last().unwrap().1.abs().max(1e-30);
    assert!(
        seam_rel < 1e-9,
        "seam mismatch: NVT end {} vs NVE start {} (rel {seam_rel:.3e})",
        rows_nvt.last().unwrap().1,
        rows[0].1
    );

    // Conservative NVE segment — deliberately loose (|drift| < 5e-2): the
    // restored state comes from a short equilibration off a cold lattice, so
    // the relaxation caveat in docs/WATER_MODEL.md applies; the exact-restore
    // and seam gates above are the strong checks. (The absolute value
    // matters — `relative_drift` is signed.)
    let drift = sum2.conservation.relative_drift;
    assert!(
        drift.abs() < 5e-2,
        "restored NVE water segment drifts {drift:.3e} (> 5e-2)"
    );

    // Wrong atom count → config error (125 waters = 375 atoms vs the
    // 192-atom checkpoint).
    let toml = WATER.replace("n_waters = 64", "n_waters = 125").replace(
        "restart_every = 500",
        "restart_every = 500\nrestore_file = \"nvt.chk\"",
    ) + "\nsteps = 10\n";
    let err = build_err_from_toml(&dir, &toml);
    assert!(
        err.contains("atom"),
        "expected an atom-count mismatch error, got: {err}"
    );
    std::fs::remove_dir_all(&dir).ok();
}
