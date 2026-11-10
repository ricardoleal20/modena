//! Instrumented B04 cost harness. Setup/validation/serialization are untimed.
//! No selector promotion or uninstrumented-throughput claim is made here.
use clap::Parser;
use modena_core::{
    ForceModel, Integrator, LennardJones, NeighborBuildStrategy, NeighborList, PairForceModel,
    SimBox, System, Vec3, VelocityVerlet, Virial,
};
use serde_json::{json, Value};
use std::alloc::{GlobalAlloc, Layout, System as SystemAllocator};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

struct CountingAllocator;
static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static LIVE_BYTES: AtomicU64 = AtomicU64::new(0);
// Counts requested Rust heap bytes, never allocator usable size or OS RSS.
// Relaxed atomics are present in BOTH strategies' timed allocation paths.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        let p = SystemAllocator.alloc(layout);
        if !p.is_null() {
            LIVE_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        p
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        let p = SystemAllocator.alloc_zeroed(layout);
        if !p.is_null() {
            LIVE_BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, layout: Layout) {
        LIVE_BYTES.fetch_sub(layout.size() as u64, Ordering::Relaxed);
        SystemAllocator.dealloc(p, layout);
    }
    unsafe fn realloc(&self, p: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        let new = SystemAllocator.realloc(p, layout, size);
        if !new.is_null() {
            if size >= layout.size() {
                LIVE_BYTES.fetch_add((size - layout.size()) as u64, Ordering::Relaxed);
            } else {
                LIVE_BYTES.fetch_sub((layout.size() - size) as u64, Ordering::Relaxed);
            }
        }
        new
    }
}
#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

type Result<T> = std::result::Result<T, String>;

#[derive(Parser, Debug)]
struct Args {
    #[arg(long)]
    snapshot: PathBuf,
    #[arg(long, value_parser = ["scalar", "fractional", "auto", "paired"])]
    strategy: String,
    #[arg(long, default_value = "fractional", value_parser = ["fractional", "auto"])]
    candidate: String,
    #[arg(long, value_parser = ["rebuild", "step"])]
    mode: String,
    #[arg(long)]
    skin: f64,
    #[arg(long)]
    threads: usize,
    #[arg(long, default_value_t = 5)]
    warmup_blocks: usize,
    #[arg(long, default_value_t = 20)]
    sample_blocks: usize,
    #[arg(long, default_value_t = 10)]
    block_steps: usize,
    #[arg(long, default_value_t = 0.0005)]
    dt: f64,
    #[arg(long)]
    output: PathBuf,
}

fn check(ok: bool, message: &str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(message.to_owned())
    }
}
fn keys(v: &Value, names: &[&str]) -> Result<()> {
    let object = v.as_object().ok_or("expected object")?;
    check(
        object.len() == names.len() && names.iter().all(|n| object.contains_key(*n)),
        "missing/unknown fields",
    )
}
fn real(v: &Value) -> Result<f64> {
    v.as_f64()
        .filter(|x| x.is_finite())
        .ok_or_else(|| "expected finite number".into())
}
fn matrix(v: &Value) -> Result<[Vec3; 3]> {
    let rows = v
        .as_array()
        .filter(|a| a.len() == 3)
        .ok_or("expected 3x3 matrix")?;
    let mut out = [Vec3::new(0., 0., 0.); 3];
    for (i, row) in rows.iter().enumerate() {
        let r = row
            .as_array()
            .filter(|a| a.len() == 3)
            .ok_or("expected length-three vector")?;
        out[i] = Vec3::new(real(&r[0])?, real(&r[1])?, real(&r[2])?);
    }
    Ok(out)
}

// serde_json parses all value syntax first. This second lexical walk rejects
// duplicate object keys (including escaped spellings) instead of last-wins data.
fn reject_duplicate_keys(raw: &[u8]) -> Result<()> {
    fn ws(b: &[u8], at: &mut usize) {
        while *at < b.len() && b[*at].is_ascii_whitespace() {
            *at += 1;
        }
    }
    fn string(b: &[u8], at: &mut usize) -> Result<String> {
        let start = *at;
        *at += 1;
        while *at < b.len() {
            match b[*at] {
                b'\\' => *at += 2,
                b'"' => {
                    *at += 1;
                    return serde_json::from_slice(&b[start..*at]).map_err(|e| e.to_string());
                }
                _ => *at += 1,
            }
        }
        Err("unterminated JSON string".into())
    }
    fn value(b: &[u8], at: &mut usize) -> Result<()> {
        ws(b, at);
        match b.get(*at) {
            Some(b'{') => {
                *at += 1;
                ws(b, at);
                let mut seen = BTreeSet::new();
                if b.get(*at) == Some(&b'}') {
                    *at += 1;
                    return Ok(());
                }
                loop {
                    ws(b, at);
                    let key = string(b, at)?;
                    check(seen.insert(key), "duplicate JSON key")?;
                    ws(b, at);
                    *at += 1;
                    value(b, at)?;
                    ws(b, at);
                    if b.get(*at) == Some(&b'}') {
                        *at += 1;
                        break;
                    }
                    *at += 1;
                }
            }
            Some(b'[') => {
                *at += 1;
                ws(b, at);
                if b.get(*at) == Some(&b']') {
                    *at += 1;
                    return Ok(());
                }
                loop {
                    value(b, at)?;
                    ws(b, at);
                    if b.get(*at) == Some(&b']') {
                        *at += 1;
                        break;
                    }
                    *at += 1;
                }
            }
            Some(b'"') => {
                string(b, at)?;
            }
            Some(_) => {
                while *at < b.len() && !b",]} \n\r\t".contains(&b[*at]) {
                    *at += 1;
                }
            }
            None => return Err("missing JSON value".into()),
        }
        Ok(())
    }
    value(raw, &mut 0)
}

struct Snapshot {
    system: System,
    id: String,
    family: String,
    cell_after: Option<SimBox>,
}
fn snapshot(raw: &[u8]) -> Result<Snapshot> {
    let v: Value = serde_json::from_slice(raw).map_err(|e| e.to_string())?;
    reject_duplicate_keys(raw)?;
    keys(
        &v,
        &[
            "schema_version",
            "case_id",
            "family",
            "n",
            "basis",
            "positions",
            "velocities",
            "physical",
            "protocol",
            "cell_after",
            "generation",
        ],
    )?;
    check(v["schema_version"].as_u64() == Some(1), "snapshot schema")?;
    let n = v["n"]
        .as_u64()
        .filter(|n| (2..=8192).contains(n))
        .ok_or("snapshot N")? as usize;
    let mut system = System::new(SimBox::triclinic(matrix(&v["basis"])?));
    check(
        system.box_size.volume().is_finite() && system.box_size.volume() > 0.,
        "invalid cell",
    )?;
    check(
        v["physical"]
            == json!({"model":"LennardJones", "precision":"f64", "sigma":1.0,"epsilon":1.0,"mass":1.0,"cutoff":2.5,"cutoff_treatment":"shifted_force","tail_correction":false}),
        "physical contract",
    )?;
    check(
        v["protocol"] == json!({"dt":0.0005,"warmup_blocks":5,"sample_blocks":20,"block_steps":10}),
        "snapshot protocol",
    )?;
    keys(
        &v["generation"],
        &["seed", "modulo", "integer_arithmetic", "density", "side"],
    )?;
    check(
        v["generation"]["seed"].as_u64() == Some(20261106)
            && v["generation"]["modulo"].as_u64() == Some(65536)
            && v["generation"]["integer_arithmetic"] == "u64_no_overflow_then_modulo"
            && real(&v["generation"]["density"])? > 0.
            && v["generation"]["side"].as_u64().is_some(),
        "generation contract",
    )?;
    system.resize(n);
    for (field, velocity) in [("positions", false), ("velocities", true)] {
        let points = v[field]
            .as_array()
            .filter(|p| p.len() == n)
            .ok_or("physical array length")?;
        for (i, point) in points.iter().enumerate() {
            let xyz = point
                .as_array()
                .filter(|p| p.len() == 3)
                .ok_or("physical vector length")?;
            let vector = Vec3::new(real(&xyz[0])?, real(&xyz[1])?, real(&xyz[2])?);
            if velocity {
                system.set_vel(i, vector);
            } else {
                system.set_pos(i, vector);
            }
        }
    }
    let id = v["case_id"]
        .as_str()
        .filter(|s| !s.is_empty())
        .ok_or("case_id")?
        .to_owned();
    let family = v["family"]
        .as_str()
        .filter(|s| {
            [
                "moderate",
                "substantial",
                "unreduced",
                "elongated",
                "grid_cap",
                "boundary",
                "basis_change",
            ]
            .contains(s)
        })
        .ok_or("family")?
        .to_owned();
    let cell_after = if v["cell_after"].is_null() {
        None
    } else {
        Some(SimBox::triclinic(matrix(&v["cell_after"])?))
    };
    check(
        (family == "basis_change") == cell_after.is_some(),
        "basis-change contract",
    )?;
    Ok(Snapshot {
        system,
        id,
        family,
        cell_after,
    })
}

fn strategy(s: &str) -> NeighborBuildStrategy {
    match s {
        "scalar" => NeighborBuildStrategy::Scalar,
        "fractional" => NeighborBuildStrategy::FractionalStencil,
        _ => NeighborBuildStrategy::Auto,
    }
}
fn strategy_name(s: NeighborBuildStrategy) -> &'static str {
    match s {
        NeighborBuildStrategy::Scalar => "scalar",
        NeighborBuildStrategy::FractionalStencil => "fractional",
        NeighborBuildStrategy::Auto => "auto",
    }
}
fn fallback_name(s: modena_core::NeighborBuildFallbackReason) -> &'static str {
    use modena_core::NeighborBuildFallbackReason::*;
    match s {
        UncertifiedGeometry => "uncertified_geometry",
        FractionalBoundary => "fractional_boundary",
        ImageDomain => "image_domain",
        GridCap => "grid_cap",
        StencilCap => "stencil_cap",
        WorkCap => "work_cap",
        MemoryCap => "memory_cap",
        AllBins => "all_bins",
        CandidateCost => "candidate_cost",
    }
}

struct Replica {
    system: System,
    neighbors: NeighborList,
    model: PairForceModel,
    integrator: VelocityVerlet,
    energy: f64,
    virial: f64,
}
impl Replica {
    fn new(system: &System, requested: &str, skin: f64) -> Result<Self> {
        let mut r = Self {
            system: system.clone(),
            neighbors: NeighborList::with_skin(2.5, skin)
                .with_triclinic_strategy_for_diagnostics(strategy(requested)),
            model: PairForceModel::new(LennardJones::new(1., 1., 2.5).map_err(|e| e.to_string())?),
            integrator: VelocityVerlet::new(),
            energy: 0.,
            virial: 0.,
        };
        r.evaluate()?;
        Ok(r)
    }
    fn evaluate(&mut self) -> Result<()> {
        self.system.zero_forces();
        let mut w = Virial::default();
        self.energy = self
            .model
            .compute(&mut self.system, &mut w, &mut self.neighbors)
            .map_err(|e| e.to_string())?;
        self.virial = w.w;
        Ok(())
    }
    fn advance(&mut self, dt: f64) -> Result<()> {
        let e = self
            .integrator
            .step(&mut self.system, &self.model, &mut self.neighbors, dt)
            .map_err(|e| e.to_string())?;
        self.energy = e.potential;
        self.virial = e.virial;
        Ok(())
    }
}

#[derive(Clone, Copy, Default)]
struct Errors {
    energy: f64,
    force: f64,
    virial: f64,
    position: f64,
    velocity: f64,
}
fn compare(a: &Replica, b: &Replica) -> Result<Errors> {
    fn diff(a: f64, b: f64) -> Result<f64> {
        check(a.is_finite() && b.is_finite(), "nonfinite physical output")?;
        let e = (a - b).abs() / (1. + a.abs().max(b.abs()));
        check(e <= 1e-10, "physical parity failure (atol=rtol=1e-10)")?;
        Ok(e)
    }
    fn maximum(a: &[f64], b: &[f64]) -> Result<f64> {
        check(a.len() == b.len(), "state length mismatch")?;
        a.iter()
            .zip(b)
            .try_fold(0_f64, |m, (&a, &b)| Ok(m.max(diff(a, b)?)))
    }
    check(
        a.neighbors.pairs() == b.neighbors.pairs(),
        "sorted pair equality failure",
    )?;
    check(
        a.neighbors.pairs().windows(2).all(|w| w[0] < w[1]),
        "unsorted/duplicate pairs",
    )?;
    Ok(Errors {
        energy: diff(a.energy, b.energy)?,
        force: maximum(&a.system.forces, &b.system.forces)?,
        virial: diff(a.virial, b.virial)?,
        position: maximum(&a.system.positions, &b.system.positions)?,
        velocity: maximum(&a.system.velocities, &b.system.velocities)?,
    })
}

fn expected_control_reason(family: &str) -> Option<modena_core::NeighborBuildFallbackReason> {
    use modena_core::NeighborBuildFallbackReason::{FractionalBoundary, GridCap};
    match family {
        "boundary" => Some(FractionalBoundary),
        "grid_cap" => Some(GridCap),
        _ => None,
    }
}

fn initial_diagnostic(r: &Replica, requested: &str) -> Result<Value> {
    let d = r
        .neighbors
        .triclinic_build_diagnostics()
        .ok_or("missing initial build diagnostics")?;
    Ok(
        json!({"requested":requested,"selected":strategy_name(d.selected),"fallback":d.fallback.map(fallback_name)}),
    )
}

fn validate_control_initial(r: &Replica, family: &str, requested: &str) -> Result<()> {
    if let Some(expected) = expected_control_reason(family) {
        let d = r
            .neighbors
            .triclinic_build_diagnostics()
            .ok_or("missing control initial diagnostics")?;
        check(
            d.selected == NeighborBuildStrategy::Scalar,
            "control initial selection must be scalar",
        )?;
        if requested == "fractional" {
            check(
                d.fallback == Some(expected),
                "control initial fallback reason mismatch",
            )?;
        }
    }
    Ok(())
}

#[derive(Default)]
struct Counters {
    elapsed: u128,
    rebuilds: usize,
    checks: u64,
    workspace: usize,
    allocations: u64,
    retained: u64,
    selected: [usize; 2],
    fallback: BTreeMap<&'static str, usize>,
    last_selected: Option<NeighborBuildStrategy>,
    last_fallback: Option<&'static str>,
    n_pairs: usize,
}
fn validate_control_counts(c: &Counters, family: &str, requested: &str) -> Result<()> {
    if let Some(expected) = expected_control_reason(family) {
        // Zero relevant builds is analyzed as inconclusive, never coverage.
        if c.rebuilds == 0 {
            return Ok(());
        }
        check(
            c.last_selected == Some(NeighborBuildStrategy::Scalar) && c.selected == [c.rebuilds, 0],
            "every control new build must select scalar",
        )?;
        if requested == "fractional" {
            let expected = fallback_name(expected);
            check(
                c.last_fallback == Some(expected)
                    && c.fallback.get(expected) == Some(&c.rebuilds)
                    && c.fallback
                        .iter()
                        .all(|(&reason, &count)| count == 0 || reason == expected),
                "every control new build must have expected fallback reason",
            )?;
        }
    }
    Ok(())
}
fn account(r: &Replica, before: usize, c: &mut Counters) -> Result<()> {
    let delta = r
        .neighbors
        .rebuild_count()
        .checked_sub(before)
        .ok_or("rebuild counter decreased")?;
    check(
        delta <= 1,
        "more than one build in a kernel; diagnostics insufficient",
    )?;
    if delta == 1 {
        let d = r
            .neighbors
            .triclinic_build_diagnostics()
            .ok_or("missing actual build diagnostics")?;
        c.rebuilds += 1;
        c.checks = c
            .checks
            .checked_add(d.candidate_checks)
            .ok_or("check overflow")?;
        c.workspace = c.workspace.max(d.requested_workspace_bytes);
        c.selected[usize::from(d.selected == NeighborBuildStrategy::FractionalStencil)] += 1;
        c.last_selected = Some(d.selected);
        c.last_fallback = d.fallback.map(fallback_name);
        c.n_pairs = r.neighbors.n_pairs();
        if let Some(f) = d.fallback {
            *c.fallback.entry(fallback_name(f)).or_default() += 1;
        }
    }
    Ok(())
}
fn block(
    r: &mut Replica,
    args: &Args,
    cell_after: Option<SimBox>,
    start: Instant,
) -> Result<Counters> {
    let mut c = Counters {
        n_pairs: r.neighbors.n_pairs(),
        ..Counters::default()
    };
    // Reserve every possible reason outside timers: account updates never allocate.
    for name in [
        "uncertified_geometry",
        "fractional_boundary",
        "image_domain",
        "grid_cap",
        "stencil_cap",
        "work_cap",
        "memory_cap",
        "all_bins",
        "candidate_cost",
    ] {
        c.fallback.insert(name, 0);
    }
    if args.mode == "step" {
        let allocations = ALLOCATIONS.load(Ordering::Relaxed);
        let clock = Instant::now();
        for _ in 0..args.block_steps {
            let before = r.neighbors.rebuild_count();
            r.advance(args.dt)?;
            account(r, before, &mut c)?;
        }
        c.elapsed = clock.elapsed().as_nanos();
        c.allocations = ALLOCATIONS.load(Ordering::Relaxed) - allocations;
        c.retained = LIVE_BYTES.load(Ordering::Relaxed);
    } else {
        let original = r.system.box_size;
        for k in 0..args.block_steps {
            if let Some(after) = cell_after {
                r.system.box_size = if k % 2 == 0 { after } else { original };
            }
            let before = r.neighbors.rebuild_count();
            let allocations = ALLOCATIONS.load(Ordering::Relaxed);
            let clock = Instant::now();
            r.neighbors.force_rebuild(&r.system);
            c.elapsed += clock.elapsed().as_nanos();
            c.allocations += ALLOCATIONS.load(Ordering::Relaxed) - allocations;
            c.retained = LIVE_BYTES.load(Ordering::Relaxed);
            account(r, before, &mut c)?;
            // Advance the real continuous replica outside the build-only timer.
            // Recompute F(t) after a cell change before the next VV step.
            if cell_after.is_some() {
                r.evaluate()?;
            }
            r.advance(args.dt)?;
            check(
                start.elapsed() < Duration::from_secs(600),
                "case_cap: 600 seconds",
            )?;
        }
    }
    check(c.elapsed > 0, "zero elapsed time")?;
    check(
        start.elapsed() < Duration::from_secs(600),
        "case_cap: 600 seconds",
    )?;
    c.fallback.retain(|_, count| *count != 0);
    Ok(c)
}

fn record(
    r: &Replica,
    c: &Counters,
    e: Errors,
    phase: &str,
    i: usize,
    order: usize,
    requested: &str,
) -> Value {
    let d = r
        .neighbors
        .triclinic_build_diagnostics()
        .expect("validated nonempty diagnostic");
    let mut selections = BTreeMap::new();
    if c.selected[0] > 0 {
        selections.insert("scalar", c.selected[0]);
    }
    if c.selected[1] > 0 {
        selections.insert("fractional", c.selected[1]);
    }
    json!({"phase":phase,"sample_index":i,"order":order,"strategy":requested,"requested":requested,
        "selected":strategy_name(c.last_selected.unwrap_or(d.selected)),"fallback":if c.last_selected.is_some() {c.last_fallback} else {d.fallback.map(fallback_name)},
        "selected_counts":selections,"fallback_counts":c.fallback,"block_steps":10,
        "elapsed_ns":c.elapsed as u64,"rebuild_count":c.rebuilds,"candidate_checks":c.checks,"n_pairs":c.n_pairs,
        "requested_bytes":c.workspace,"retained_bytes":c.retained,"allocation_count":c.allocations,
        "max_energy_error":e.energy,"max_force_error":e.force,"virial_error":e.virial,
        "max_position_error":e.position,"max_velocity_error":e.velocity,"pair_equal":true,"physics_passed":true})
}

fn output(command: &mut Command) -> Result<String> {
    let p = command.output().map_err(|e| e.to_string())?;
    check(
        p.status.success(),
        &format!(
            "metadata command failed: {}",
            String::from_utf8_lossy(&p.stderr)
        ),
    )?;
    String::from_utf8(p.stdout)
        .map(|s| s.trim().to_owned())
        .map_err(|e| e.to_string())
}
fn hash(path: &Path) -> Result<String> {
    Ok(
        output(Command::new("shasum").args(["-a", "256"]).arg(path))?
            .split_whitespace()
            .next()
            .ok_or("hash output")?
            .to_owned(),
    )
}
fn metadata(args: &Args, s: &Snapshot) -> Result<Value> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .ok_or("root")?;
    let identity = output(Command::new("python3").current_dir(root).args(["-c", "import runpy,json; m=runpy.run_path('bench/triclinic_neighbors.py'); print(json.dumps(m['source_identity']()))"]))?;
    let identity: Value = serde_json::from_str(&identity).map_err(|e| e.to_string())?;
    let build_rev = option_env!("B04_BUILD_REV")
        .ok_or("build missing B04_BUILD_REV; freeze source before compilation")?;
    let build_source =
        option_env!("B04_BUILD_SOURCE_SHA256").ok_or("build missing B04_BUILD_SOURCE_SHA256")?;
    let compiler = option_env!("B04_RUSTC_VERSION").ok_or("build missing B04_RUSTC_VERSION")?;
    check(
        identity[0] == build_rev && identity[1] == build_source,
        "binary/source runtime epoch differs from declared build epoch",
    )?;
    Ok(
        json!({"schema_version":1,"base_rev":identity[0],"source_sha256":identity[1],"build_rev":build_rev,"build_source_sha256":build_source,
        "binary_sha256":hash(&std::env::current_exe().map_err(|e|e.to_string())?)?,"snapshot_sha256":hash(&args.snapshot)?,
        "case_id":s.id,"family":s.family,"n":s.system.n(),"model":"LennardJones","precision":"f64","cutoff_treatment":"shifted_force",
        "epsilon":1.0,"sigma":1.0,"mass":1.0,"cutoff":2.5,"tail_correction":false,
        "host":output(Command::new("uname").args(["-a"]))?,"compiler":compiler,"profile":if cfg!(debug_assertions) {"debug"} else {"release"},"features":["neighbor-diagnostics"],"threads":args.threads,"mode":args.mode,"skin":args.skin,"dt":args.dt,
        "warmup_blocks":args.warmup_blocks,"sample_blocks":args.sample_blocks,"block_steps":args.block_steps,"candidate":args.candidate,"paired":args.strategy=="paired",
        "memory_scope":"requested_neighbor_workspace; retained_process_requested_live_bytes",
        "allocation_scope":"process_alloc_and_realloc_calls_in_timed_kernel",
        "instrumentation":"counting_global_allocator_relaxed_atomics",
        "timing_scope":if args.mode=="step" {"complete_steps_with_per_step_diagnostic_accounting"} else {"forced_neighbor_builds_only_continuous_VV_advance_outside_timers"}}),
    )
}

fn preflight(s: &System, skin: f64, start: Instant) -> Result<(f64, f64)> {
    let mut distance = f64::INFINITY;
    let mut margin = f64::INFINITY;
    for i in 0..s.n() {
        for j in i + 1..s.n() {
            let r = s.box_size.min_image(s.pos(j) - s.pos(i)).norm();
            check(
                r.is_finite() && r > 0.,
                "initial overlap/nonfinite distance",
            )?;
            distance = distance.min(r);
            margin = margin.min((r - 2.5).abs()).min((r - (2.5 + skin)).abs());
        }
        check(
            start.elapsed() < Duration::from_secs(600),
            "case_cap during preflight",
        )?;
    }
    check(margin > 1e-8, "unresolved initial cutoff fixture")?;
    Ok((distance, margin))
}

fn execute(args: &Args, document: &mut Value) -> Result<()> {
    check(
        [0., 0.3].contains(&args.skin) && [1, 4].contains(&args.threads),
        "skin/threads contract",
    )?;
    check(
        args.warmup_blocks == 5
            && args.sample_blocks == 20
            && args.block_steps == 10
            && args.dt == 0.0005,
        "fixed 5/20/10/dt protocol",
    )?;
    let start = Instant::now();
    let raw = fs::read(&args.snapshot).map_err(|e| e.to_string())?;
    let s = snapshot(&raw)?;
    check(
        s.cell_after.is_none() || args.mode == "rebuild",
        "cell-change is a separate rebuild case",
    )?;
    document["metadata"] = metadata(args, &s)?;
    let (mut minimum, mut margin) = preflight(&s.system, args.skin, start)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .map_err(|e| e.to_string())?;
    pool.install(|| -> Result<()> {
        let requested = if args.strategy=="paired" { args.candidate.as_str() } else { args.strategy.as_str() };
        let mut a = Replica::new(&s.system,"scalar",args.skin)?;
        let mut b = Replica::new(&s.system,requested,args.skin)?;
        document["validation"] = json!({"initial":false,"final":false,"minimum_distance":minimum,"minimum_cutoff_margin":margin,
            "initial_diagnostics":{"scalar":initial_diagnostic(&a,"scalar")?,"candidate":initial_diagnostic(&b,requested)?}});
        compare(&a,&b)?;
        validate_control_initial(&b,&s.family,requested)?;
        if let Some(after) = s.cell_after {
            let mut changed = s.system.clone(); changed.box_size=after;
            let (d,m)=preflight(&changed,args.skin,start)?; minimum=minimum.min(d); margin=margin.min(m);
            let ca=Replica::new(&changed,"scalar",args.skin)?;
            let cb=Replica::new(&changed,requested,args.skin)?;
            compare(&ca,&cb)?;
        }
        document["validation"]["initial"]=json!(true);
        document["validation"]["minimum_distance"]=json!(minimum);
        document["validation"]["minimum_cutoff_margin"]=json!(margin);
        let records = document["records"].as_array_mut().ok_or("records")?;
        records.reserve_exact(50);
        for (phase,count) in [("warmup",args.warmup_blocks),("measured",args.sample_blocks)] {
            for i in 0..count {
                let (ca,cb) = if i%2==0 { (block(&mut a,args,s.cell_after,start)?,block(&mut b,args,s.cell_after,start)?) }
                              else { let cb=block(&mut b,args,s.cell_after,start)?; (block(&mut a,args,s.cell_after,start)?,cb) };
                let comparison = compare(&a,&b);
                let e = comparison.as_ref().copied().unwrap_or_default();
                let mut ra=record(&a,&ca,e,phase,i,usize::from(i%2!=0),"scalar");
                let mut rb=record(&b,&cb,e,phase,i,usize::from(i%2==0),requested);
                if comparison.is_err() {
                    for r in [&mut ra,&mut rb] {
                        r["pair_equal"]=json!(a.neighbors.pairs()==b.neighbors.pairs()); r["physics_passed"]=json!(false);
                        for k in ["max_energy_error","max_force_error","virial_error","max_position_error","max_velocity_error"] { r[k]=Value::Null; }
                    }
                }
                if args.strategy=="paired" {
                    if i%2==0 { records.push(ra); records.push(rb); } else { records.push(rb); records.push(ra); }
                } else { records.push(rb); }
                comparison?;
                validate_control_counts(&cb,&s.family,requested)?;
            }
        }
        compare(&a,&b)?;
        document["validation"]["final"] = json!(true);
        Ok(())
    })?;
    document["status"] = json!("complete");
    Ok(())
}

fn main() {
    let args = Args::parse();
    // Exclusive creation prevents accidental stale-result reuse, even on error.
    let mut file = match OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args.output)
    {
        Ok(f) => f,
        Err(e) => {
            eprintln!("output: {e}");
            std::process::exit(1);
        }
    };
    let mut document = json!({"schema_version":1,"status":"invalid","metadata":null,"validation":null,"records":[],"error":null});
    let result = execute(&args, &mut document);
    if let Err(ref e) = result {
        document["error"] = json!(e);
        eprintln!("{e}");
    }
    let written = serde_json::to_writer_pretty(&mut file, &document)
        .and_then(|_| file.write_all(b"\n").map_err(serde_json::Error::io));
    if let Err(e) = written {
        eprintln!("write: {e}");
        std::process::exit(1);
    }
    if result.is_err() {
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn duplicates_include_escaped_keys() {
        for raw in [br#"{"a":1,"a":2}"#.as_slice(), br#"{"a":1,"\u0061":2}"#] {
            assert!(reject_duplicate_keys(raw).is_err());
        }
        assert!(reject_duplicate_keys(br#"{"a":[{"a":2}],"b":3}"#).is_ok());
    }
    #[test]
    fn cache_report_does_not_add_old_build_checks() {
        let mut s = System::new(SimBox::triclinic([
            Vec3::new(10., 0., 0.),
            Vec3::new(2., 10., 0.),
            Vec3::new(1., 1.5, 10.),
        ]));
        s.resize(2);
        s.set_pos(0, Vec3::new(0.31, 0.41, 0.51));
        s.set_pos(1, Vec3::new(1.31, 0.41, 0.51));
        let mut r = Replica::new(&s, "scalar", 0.3).unwrap();
        let before = r.neighbors.rebuild_count();
        r.advance(0.0005).unwrap();
        let mut c = Counters::default();
        account(&r, before, &mut c).unwrap();
        assert_eq!(c.rebuilds, 0);
        assert_eq!(c.checks, 0);
        assert_eq!(c.workspace, 0);
    }
    #[test]
    fn parity_guard_checks_actual_force_not_checksum() {
        let mut s = System::new(SimBox::triclinic([
            Vec3::new(10., 0., 0.),
            Vec3::new(2., 10., 0.),
            Vec3::new(1., 1.5, 10.),
        ]));
        s.resize(2);
        s.set_pos(1, Vec3::new(1.1, 0., 0.));
        let a = Replica::new(&s, "scalar", 0.).unwrap();
        let mut b = Replica::new(&s, "fractional", 0.).unwrap();
        assert!(compare(&a, &b).is_ok());
        b.system.forces[0] += 1e-3;
        assert!(compare(&a, &b).is_err());
    }

    fn control_fixture(family: &str) -> Snapshot {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap();
        let code = format!("import runpy,json; m=runpy.run_path('bench/triclinic_neighbors.py'); print(json.dumps(m['make_control_snapshot']({family:?})))");
        let raw = Command::new("python3")
            .current_dir(root)
            .args(["-c", &code])
            .output()
            .unwrap();
        assert!(
            raw.status.success(),
            "{}",
            String::from_utf8_lossy(&raw.stderr)
        );
        snapshot(&raw.stdout).unwrap()
    }

    #[test]
    fn wrong_real_fallback_is_rejected_for_named_control() {
        let mut s = System::new(SimBox::triclinic([
            Vec3::new(3., 0., 0.),
            Vec3::new(0.6, 3., 0.),
            Vec3::new(0.3, 0.45, 3.),
        ]));
        s.resize(2);
        s.set_pos(1, Vec3::new(1.1, 0., 0.));
        let r = Replica::new(&s, "fractional", 0.).unwrap();
        assert_eq!(
            r.neighbors.triclinic_build_diagnostics().unwrap().fallback,
            Some(modena_core::NeighborBuildFallbackReason::AllBins)
        );
        assert!(validate_control_initial(&r, "boundary", "fractional").is_err());
        assert!(validate_control_initial(&r, "grid_cap", "fractional").is_err());
    }

    #[test]
    fn real_controls_keep_expected_reason_through_continuous_steps() {
        for family in ["boundary", "grid_cap"] {
            let s = control_fixture(family);
            for skin in [0., 0.3] {
                let mut a = Replica::new(&s.system, "scalar", skin).unwrap();
                let mut b = Replica::new(&s.system, "fractional", skin).unwrap();
                validate_control_initial(&b, family, "fractional").unwrap();
                let expected = if family == "boundary" {
                    modena_core::NeighborBuildFallbackReason::FractionalBoundary
                } else {
                    modena_core::NeighborBuildFallbackReason::GridCap
                };
                assert_eq!(
                    b.neighbors.triclinic_build_diagnostics().unwrap().fallback,
                    Some(expected)
                );
                assert_eq!(
                    b.neighbors.triclinic_build_diagnostics().unwrap().selected,
                    NeighborBuildStrategy::Scalar
                );
                compare(&a, &b).unwrap();
                let mut new_builds = 0;
                for _ in 0..250 {
                    let before = b.neighbors.rebuild_count();
                    a.advance(0.0005).unwrap();
                    b.advance(0.0005).unwrap();
                    compare(&a, &b).unwrap();
                    if b.neighbors.rebuild_count() > before {
                        new_builds += 1;
                        let d = b.neighbors.triclinic_build_diagnostics().unwrap();
                        assert_eq!(d.selected, NeighborBuildStrategy::Scalar);
                        assert_eq!(d.fallback, Some(expected));
                    }
                    if family == "boundary" {
                        assert_eq!(&b.system.positions[..3], &[0., 0., 0.]);
                        assert_eq!(&b.system.velocities[..3], &[0., 0., 0.]);
                        assert_eq!(&b.system.forces[..3], &[0., 0., 0.]);
                    }
                }
                if skin == 0. {
                    assert_eq!(new_builds, 250);
                }
            }
        }
    }
}
