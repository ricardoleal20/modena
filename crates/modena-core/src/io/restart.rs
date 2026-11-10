//! Versioned, validated and atomic scientific continuation checkpoints.
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::io::Write;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::constraint_stress::{BathImpulse, ProjectionReport, ProjectionStage};
use crate::error::{ModenaError, Result};
use crate::stress::{
    pressure_total, ConstraintSolverOptions, IntervalPressure, MissingTensorReason,
    ModelEvaluation, PressureObservation, SymTensor6, TensorAvailability,
};
use crate::system::SimBox;
use crate::vec3::Vec3;

/// Format tag (first token of the header line).
const MAGIC: &str = "MODENA_RESTART";
/// Current format version (second token of the header line).
const VERSION: u32 = 2;
static TEMP_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// A loaded/saved checkpoint. All values in reduced units.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    /// Atom count (sanity-checked against the restoring config's structure).
    pub atoms: usize,
    /// Global step counter at save time.
    pub step: usize,
    /// Reduced time `step · dt` at save time (informational; the step is
    /// authoritative — `dt` belongs to the config, not the checkpoint).
    pub time: f64,
    /// Simulation box (orthorhombic edges, or a triclinic cell).
    pub box_size: SimBox,
    /// Flat SoA positions, length `3 · atoms`.
    pub positions: Vec<f64>,
    /// Flat SoA velocities, length `3 · atoms`.
    pub velocities: Vec<f64>,
    /// None for legacy v1 or mechanical-only imports. Full continuation needs v2.
    pub continuation: Option<ContinuationState>,
}

/// All identity, cached mechanics, stochastic and extended variables needed
/// for exact continuation. Physical signatures are normalized configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ContinuationState {
    pub forces: Vec<f64>,
    pub masses: Vec<f64>,
    pub species: Vec<u32>,
    /// Explicit atomic identity; absent in older v2 checkpoints. Empty may
    /// match only an unspecified native system during continuation.
    #[serde(default)]
    pub atomic_numbers: Vec<u8>,
    pub charges: Vec<f64>,
    pub elements: Vec<String>,
    pub constraint_count: usize,
    pub conserved_com_dof: usize,
    pub potential_energy: f64,
    pub virial: f64,
    pub physical_signature: Value,
    pub integrator: Value,
    pub thermostat: Value,
    pub barostat: Value,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Version2 {
    atoms: usize,
    step: usize,
    time: f64,
    box_edges: [f64; 3],
    box_vectors: Option<[[f64; 3]; 3]>,
    positions: Vec<f64>,
    velocities: Vec<f64>,
    continuation: Option<ContinuationState>,
}

/// Mechanical continuation and its versioned pressure-estimator sibling.
#[derive(Debug, Clone)]
pub struct StressCheckpoint {
    pub mechanical: Checkpoint,
    pub stress: Option<RestartStress>,
}

/// Reduced-unit endpoint after controllers and history of the last pre-controller interval.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RestartStress {
    pub estimator_version: u32,
    pub solver: ConstraintSolverOptions,
    pub model: ModelEvaluation,
    pub endpoint: PressureObservation,
    pub interval: Option<IntervalPressure>,
    pub last_stages: Vec<ProjectionReport>,
    pub last_bath: Option<BathImpulse>,
    pub step: usize,
    pub dt: f64,
    pub physical_signature: Value,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Version3 {
    mechanical: Version2,
    #[serde(deserialize_with = "Deserialize::deserialize")]
    stress: Option<RestartStress>,
}

fn mechanical_payload(cp: &Checkpoint) -> Version2 {
    Version2 {
        atoms: cp.atoms,
        step: cp.step,
        time: cp.time,
        box_edges: cp.box_size.size.0,
        box_vectors: cp.box_size.vectors().map(|v| v.map(|x| x.0)),
        positions: cp.positions.clone(),
        velocities: cp.velocities.clone(),
        continuation: cp.continuation.clone(),
    }
}

/// Save the original v2 mechanical format.
pub fn save<P: AsRef<Path>>(path: P, cp: &Checkpoint) -> Result<()> {
    validate(cp)?;
    atomic_save(path.as_ref(), VERSION, &mechanical_payload(cp))
}

/// Save mechanical and pressure continuation together in a v3 file.
pub fn save_with_stress<P: AsRef<Path>>(path: P, cp: &StressCheckpoint) -> Result<()> {
    validate(&cp.mechanical)?;
    if let Some(stress) = &cp.stress {
        validate_stress(&cp.mechanical, stress)?;
    }
    atomic_save(
        path.as_ref(),
        3,
        &Version3 {
            mechanical: mechanical_payload(&cp.mechanical),
            stress: cp.stress.clone(),
        },
    )
}

/// Sync a unique sibling, atomically rename, then sync the containing directory.
fn atomic_save(path: &Path, version: u32, payload: &impl Serialize) -> Result<()> {
    let json = serde_json::to_vec(payload)
        .map_err(|e| ModenaError::Input(format!("checkpoint serialization: {e}")))?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    std::fs::create_dir_all(parent)?;
    let sequence = TEMP_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let file_name = path
        .file_name()
        .ok_or_else(|| ModenaError::Input("checkpoint path has no filename".into()))?;
    let temp = parent.join(format!(
        ".{}.{}.{}.tmp",
        file_name.to_string_lossy(),
        std::process::id(),
        sequence
    ));
    let mut created = false;
    let result = (|| -> Result<()> {
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp)?;
        created = true;
        writeln!(file, "{MAGIC} {version}")?;
        file.write_all(&json)?;
        file.write_all(b"\n")?;
        file.sync_all()?;
        drop(file);
        std::fs::rename(&temp, path)?;
        std::fs::File::open(parent)?.sync_all()?;
        Ok(())
    })();
    if result.is_err() && created {
        let _ = std::fs::remove_file(&temp);
    }
    result
}

fn validate(cp: &Checkpoint) -> Result<()> {
    let bad = |s: &str| ModenaError::Input(format!("checkpoint: {s}"));
    let count = cp
        .atoms
        .checked_mul(3)
        .ok_or_else(|| bad("atom buffer size overflow"))?;
    if cp.atoms == 0 || cp.positions.len() != count || cp.velocities.len() != count {
        return Err(bad("invalid atom count or buffer lengths"));
    }
    if !cp.time.is_finite()
        || cp.time < 0.0
        || cp
            .positions
            .iter()
            .chain(&cp.velocities)
            .any(|x| !x.is_finite())
    {
        return Err(bad("nonfinite or invalid mechanical state"));
    }
    if cp
        .box_size
        .size
        .0
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !cp.box_size.volume().is_finite()
        || cp.box_size.volume() <= 0.0
        || cp
            .box_size
            .inv_basis_rows()
            .iter()
            .flat_map(|x| x.0)
            .any(|x| !x.is_finite())
    {
        return Err(bad("box must have finite positive edges and volume"));
    }
    if let Some(c) = &cp.continuation {
        if c.forces.len() != count
            || c.masses.len() != cp.atoms
            || c.species.len() != cp.atoms
            || (!c.atomic_numbers.is_empty() && c.atomic_numbers.len() != cp.atoms)
            || c.charges.len() != cp.atoms
            || c.elements.len() != cp.atoms
        {
            return Err(bad("continuation atom buffers have inconsistent lengths"));
        }
        if c.forces.iter().chain(&c.charges).any(|x| !x.is_finite())
            || c.masses.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || !c.potential_energy.is_finite()
            || !c.virial.is_finite()
        {
            return Err(bad(
                "invalid continuation masses, energy, forces or charges",
            ));
        }
        if c.constraint_count
            .checked_add(c.conserved_com_dof)
            .is_none_or(|x| x > count)
            || ![0, 3].contains(&c.conserved_com_dof)
            || !c.physical_signature.is_object()
        {
            return Err(bad("invalid physical DOF or signature"));
        }
    }
    Ok(())
}

pub(crate) fn stress_signature(mechanical: &Value, solver: ConstraintSolverOptions) -> Value {
    let mut signature = mechanical.clone();
    signature["pressure_estimator"] = serde_json::json!({
        "version": 1, "solver": solver,
        "endpoint_timing": "endpoint_after_controllers", "interval_timing": "interval_pre_controllers",
    });
    signature
}

fn invalid_stress(message: &str) -> ModenaError {
    ModenaError::Input(format!("checkpoint stress: {message}"))
}

fn validate_stress(cp: &Checkpoint, s: &RestartStress) -> Result<()> {
    let state = cp
        .continuation
        .as_ref()
        .ok_or_else(|| invalid_stress("missing mechanical continuation"))?;
    if s.estimator_version != 1 {
        return Err(invalid_stress("unsupported estimator version"));
    }
    s.solver.validate()?;
    let signature = &state.physical_signature;
    if s.step != cp.step
        || !s.dt.is_finite()
        || s.dt <= 0.
        || cp.time != s.step as f64 * s.dt
        || signature["integrator"]["dt"].as_f64() != Some(s.dt)
        || s.physical_signature != stress_signature(signature, s.solver)
    {
        return Err(invalid_stress(
            "step, timestep or estimator physical signature mismatch",
        ));
    }
    let end = &s.endpoint;
    // NHC has no staged estimator. An accepted box move still refreshes the
    // full conservative cache, while its legacy endpoint remains unavailable.
    let legacy_endpoint = signature["integrator"]["type"] == "nose_hoover_chain"
        && !signature["barostat"].is_null()
        && s.step > 0
        && s.interval.is_none()
        && matches!(&end.model,TensorAvailability::Missing(r) if r.len()==1 && r[0].component=="integrator without step stress" && r[0].reason==MissingTensorReason::LegacyForceModel);
    if s.model.energy != state.potential_energy
        || s.model.scalar_virial != state.virial
        || (s.model.tensor != end.model && !legacy_endpoint)
        || end.volume != cp.box_size.volume()
    {
        return Err(invalid_stress(
            "endpoint model, scalar cache or volume mismatch",
        ));
    }
    if let TensorAvailability::Complete(t) = s.model.tensor {
        let trace = t.trace();
        let tolerance = 5e-12 * (1. + t.0[0].abs() + t.0[1].abs() + t.0[2].abs());
        if !t.finite()
            || !trace.is_finite()
            || !tolerance.is_finite()
            || (trace - state.virial).abs() > tolerance
        {
            return Err(invalid_stress(
                "model tensor trace differs from scalar virial",
            ));
        }
    }
    let kinetic =
        cp.velocities
            .chunks_exact(3)
            .zip(&state.masses)
            .fold(SymTensor6::ZERO, |sum, (v, m)| {
                let v = Vec3::new(v[0], v[1], v[2]);
                sum.plus(SymTensor6::dyad(v, v).scaled(*m))
            });
    if kinetic != end.kinetic || (state.constraint_count == 0 && end.constraint != SymTensor6::ZERO)
    {
        return Err(invalid_stress(
            "endpoint kinetic or constraint tensor inconsistent with mechanics",
        ));
    }
    check_pressure(
        end.kinetic,
        &end.model,
        end.constraint,
        end.volume,
        &end.total,
        end.scalar,
    )?;
    let Some(interval) = &s.interval else {
        if !s.last_stages.is_empty() || s.last_bath.is_some() {
            return Err(invalid_stress(
                "stages or bath without a completed interval",
            ));
        }
        return Ok(());
    };
    let transformed = matches!(&interval.model_average,TensorAvailability::Missing(why) if why.iter().any(|r|r.reason==MissingTensorReason::PostStepTransform));
    if s.step == 0
        || interval.dt != s.dt
        || (!transformed && interval.volume != end.volume)
        || (transformed && signature["barostat"].is_null())
    {
        return Err(invalid_stress(
            "interval step, volume, timestep or transform inconsistent",
        ));
    }
    check_pressure(
        interval.kinetic_average,
        &interval.model_average,
        interval.constraint_average,
        interval.volume,
        &interval.total,
        interval.scalar,
    )?;
    validate_history(s, state, interval)
}

fn check_pressure(
    kinetic: SymTensor6,
    model: &TensorAvailability,
    constraint: SymTensor6,
    volume: f64,
    total: &TensorAvailability,
    scalar: Option<f64>,
) -> Result<()> {
    if matches!(model,TensorAvailability::Missing(why) if why.is_empty() || why.iter().any(|r|r.component.is_empty()))
    {
        return Err(invalid_stress(
            "unavailable tensor needs component provenance",
        ));
    }
    let (expected, p) = pressure_total(kinetic, model, constraint, volume)?;
    if &expected != total || p != scalar {
        return Err(invalid_stress(
            "pressure components, total or scalar trace disagree",
        ));
    }
    Ok(())
}

fn validate_history(
    s: &RestartStress,
    state: &ContinuationState,
    interval: &IntervalPressure,
) -> Result<()> {
    use ProjectionStage::*;
    let signature = &state.physical_signature;
    let expected: &[ProjectionStage] = match signature["integrator"]["type"].as_str() {
        Some("velocity_verlet")
            if state.constraint_count > 0 || !signature["constraints"].is_null() =>
        {
            &[VvPosition, VvFinal]
        }
        Some("velocity_verlet") => &[],
        Some("langevin_baoab") => &[
            BaoabB0,
            BaoabA1Pre,
            BaoabA1Post,
            BaoabO,
            BaoabA2Pre,
            BaoabA2Post,
            BaoabB1,
        ],
        _ => return Err(invalid_stress("unsupported integrator for staged interval")),
    };
    if !s.last_stages.iter().map(|r| &r.stage).eq(expected.iter()) {
        return Err(invalid_stress("incomplete or out-of-order stages"));
    }
    let mut mechanical = SymTensor6::ZERO;
    for report in &s.last_stages {
        let mut tensor = SymTensor6::ZERO;
        for impulse in &report.impulses {
            if impulse.edge >= state.constraint_count
                || impulse
                    .lever
                    .0
                    .iter()
                    .chain(&impulse.impulse_on_j.0)
                    .any(|x| !x.is_finite())
            {
                return Err(invalid_stress("invalid constraint edge or impulse"));
            }
            tensor = tensor.plus(SymTensor6::dyad(impulse.lever, impulse.impulse_on_j));
        }
        let target_residual = match report.stage {
            VvPosition | BaoabA1Pre | BaoabA2Pre => report.max_position_residual,
            VvFinal | BaoabB0 | BaoabA1Post | BaoabO | BaoabA2Post | BaoabB1 => {
                report.max_velocity_residual
            }
        };
        if !tensor.finite()
            || target_residual > s.solver.tolerance
            || tensor != report.tensor_impulse
            || [report.max_position_residual, report.max_velocity_residual]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.)
        {
            return Err(invalid_stress(
                "invalid projection residual or ledger tensor",
            ));
        }
        if report.stage != BaoabO {
            mechanical = mechanical.plus(tensor);
        }
    }
    if SymTensor6(mechanical.0.map(|x| x / s.dt)) != interval.constraint_average {
        return Err(invalid_stress(
            "interval reaction does not match mechanical stages",
        ));
    }
    let ou = s.last_stages.iter().find(|r| r.stage == BaoabO);
    match (&s.last_bath, ou) {
        (Some(b), Some(o))
            if b.unprojected.len() == state.masses.len()
                && b.reaction.len() == state.masses.len()
                && (state.constraint_count != 0
                    || b.reaction.iter().all(|j| *j == Vec3::zero()))
                && b.unprojected
                    .iter()
                    .chain(&b.reaction)
                    .flat_map(|v| v.0)
                    .all(f64::is_finite)
                && b.unprojected_moment.finite()
                && b.reaction_moment == o.tensor_impulse =>
        {
            Ok(())
        }
        (None, None) => Ok(()),
        _ => Err(invalid_stress("missing, unexpected or invalid OU bath")),
    }
}

/// Parse a checkpoint from `path`.
pub fn load<P: AsRef<Path>>(path: P) -> Result<Checkpoint> {
    parse(path.as_ref(), false).map(|(cp, _)| cp)
}

/// Load v3 pressure continuation, or import v1/v2 without estimator history.
pub fn load_with_stress<P: AsRef<Path>>(path: P) -> Result<StressCheckpoint> {
    let (mechanical, stress) = parse(path.as_ref(), true)?;
    Ok(StressCheckpoint { mechanical, stress })
}

fn unpack_mechanical(payload: Version2) -> Result<Checkpoint> {
    let bad = |s: &str| ModenaError::Input(format!("checkpoint: {s}"));
    let box_size = match payload.box_vectors {
        Some(v) => {
            let vectors = v.map(|a| Vec3::new(a[0], a[1], a[2]));
            if v.iter().flatten().any(|x| !x.is_finite())
                || !vectors[0].dot(vectors[1].cross(vectors[2])).is_finite()
                || vectors[0].dot(vectors[1].cross(vectors[2])).abs() == 0.0
            {
                return Err(bad("invalid triclinic basis"));
            }
            let box_size = SimBox::triclinic(vectors);
            if box_size.size.0 != payload.box_edges {
                return Err(bad("inconsistent triclinic box edges"));
            }
            box_size
        }
        None => SimBox::from_array(payload.box_edges),
    };
    let cp = Checkpoint {
        atoms: payload.atoms,
        step: payload.step,
        time: payload.time,
        box_size,
        positions: payload.positions,
        velocities: payload.velocities,
        continuation: payload.continuation,
    };
    validate(&cp)?;
    Ok(cp)
}

fn parse(path: &Path, read_stress: bool) -> Result<(Checkpoint, Option<RestartStress>)> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text.lines();

    let bad = |msg: &str| ModenaError::Input(format!("restart file `{}`: {msg}", path.display()));

    // Header line: MAGIC VERSION.
    let header = lines.next().ok_or_else(|| bad("file is empty"))?;
    let mut header_parts = header.split_whitespace();
    if header_parts.next() != Some(MAGIC) {
        return Err(bad("not a MODENA restart file (bad magic line)"));
    }
    let version: u32 = header_parts
        .next()
        .ok_or_else(|| bad("missing format version"))?
        .parse()
        .map_err(|e| bad(&format!("bad format version: {e}")))?;
    if header_parts.next().is_some() {
        return Err(bad("trailing header tokens"));
    }
    if version == VERSION {
        let payload: Version2 = serde_json::from_str(&lines.collect::<Vec<_>>().join("\n"))
            .map_err(|e| bad(&format!("malformed v2 payload: {e}")))?;
        return Ok((unpack_mechanical(payload)?, None));
    }
    if version == 3 {
        if !read_stress {
            return Err(bad(
                "version 3 includes pressure continuation; use load_with_stress",
            ));
        }
        let payload: Version3 = serde_json::from_str(&lines.collect::<Vec<_>>().join("\n"))
            .map_err(|e| bad(&format!("malformed v3 payload: {e}")))?;
        let mechanical = unpack_mechanical(payload.mechanical)?;
        if let Some(stress) = &payload.stress {
            validate_stress(&mechanical, stress)?;
        }
        return Ok((mechanical, payload.stress));
    }
    if version != 1 {
        return Err(bad(&format!(
            "unsupported version {version} (expected {VERSION})"
        )));
    }

    // Key/value header lines until the `positions` marker.
    let mut atoms: Option<usize> = None;
    let mut step = None;
    let mut time = None;
    let mut keys = std::collections::HashSet::new();
    let mut box_edges: Option<[f64; 3]> = None;
    let mut triclinic: Option<[Vec3; 3]> = None;
    loop {
        let line = lines
            .next()
            .ok_or_else(|| bad("missing `positions` block"))?;
        let mut it = line.split_whitespace();
        let key = it.next();
        if !keys.insert(key) {
            return Err(bad("duplicate header key"));
        }
        match key {
            Some("atoms") => {
                atoms = Some(parse_usize(it.next(), &bad)?);
            }
            Some("step") => step = Some(parse_usize(it.next(), &bad)?),
            Some("time") => time = Some(parse_f64_tokens(&mut it, 1, &bad)?[0]),
            Some("box") => {
                let e = parse_f64_tokens(&mut it, 3, &bad)?;
                box_edges = Some([e[0], e[1], e[2]]);
            }
            Some("triclinic") => {
                let v = parse_f64_tokens(&mut it, 9, &bad)?;
                triclinic = Some([
                    Vec3::new(v[0], v[1], v[2]),
                    Vec3::new(v[3], v[4], v[5]),
                    Vec3::new(v[6], v[7], v[8]),
                ]);
            }
            Some("positions") => {
                if it.next().is_some() {
                    return Err(bad("trailing marker tokens"));
                }
                break;
            }
            _ => return Err(bad(&format!("malformed header line: {line}"))),
        }
        if it.next().is_some() {
            return Err(bad("trailing header values"));
        }
    }
    let atoms = atoms.ok_or_else(|| bad("missing `atoms` header"))?;
    let edges = box_edges.ok_or_else(|| bad("missing `box` header"))?;
    let box_size = match triclinic {
        Some(v) => {
            if v.iter().flat_map(|x| x.0).any(|x| !x.is_finite())
                || !v[0].dot(v[1].cross(v[2])).is_finite()
                || v[0].dot(v[1].cross(v[2])).abs() == 0.0
            {
                return Err(bad("invalid triclinic basis"));
            }
            let b = SimBox::triclinic(v);
            if b.size.0 != edges {
                return Err(bad("inconsistent triclinic box edges"));
            }
            b
        }
        None => SimBox::from_array(edges),
    };

    let count = atoms
        .checked_mul(3)
        .ok_or_else(|| bad("atom buffer size overflow"))?;
    if count > lines.clone().count() {
        return Err(bad("atom count exceeds available data"));
    }
    let positions = parse_block(&mut lines, count, &bad)?;
    expect_marker(&mut lines, "velocities", &bad)?;
    let velocities = parse_block(&mut lines, count, &bad)?;
    if lines.any(|line| !line.trim().is_empty()) {
        return Err(bad("trailing corrupt data"));
    }

    let cp = Checkpoint {
        atoms,
        step: step.ok_or_else(|| bad("missing step"))?,
        time: time.ok_or_else(|| bad("missing time"))?,
        box_size,
        positions,
        velocities,
        continuation: None,
    };
    validate(&cp)?;
    Ok((cp, None))
}

/// Parse one integer (usize) header token.
fn parse_usize(tok: Option<&str>, bad: &impl Fn(&str) -> ModenaError) -> Result<usize> {
    tok.and_then(|t| t.parse().ok())
        .ok_or_else(|| bad("expected an unsigned integer header value"))
}

/// Parse exactly `n` whitespace-separated f64 tokens from a header line.
fn parse_f64_tokens(
    it: &mut std::str::SplitWhitespace<'_>,
    n: usize,
    bad: &impl Fn(&str) -> ModenaError,
) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let tok = it.next().ok_or_else(|| bad("truncated header line"))?;
        out.push(
            tok.parse()
                .map_err(|e| bad(&format!("bad float `{tok}`: {e}")))?,
        );
    }
    Ok(out)
}

/// Parse `n` one-float-per-line values.
fn parse_block<'a>(
    lines: &mut std::str::Lines<'a>,
    n: usize,
    bad: &impl Fn(&str) -> ModenaError,
) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let line = lines.next().ok_or_else(|| bad("truncated float block"))?;
        out.push(
            line.trim()
                .parse()
                .map_err(|e| bad(&format!("bad float line `{line}`: {e}")))?,
        );
    }
    Ok(out)
}

fn expect_marker<'a>(
    lines: &mut std::str::Lines<'a>,
    marker: &str,
    bad: &impl Fn(&str) -> ModenaError,
) -> Result<()> {
    match lines.next() {
        Some(l) if l.trim() == marker => Ok(()),
        _ => Err(bad(&format!("missing `{marker}` marker"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write as _;

    /// 17 significant digits round-trip exactly, including the awkward ones.
    #[test]
    fn float_roundtrip_is_exact() {
        let awkward = [
            0.1_f64 + 0.2,
            f64::MIN_POSITIVE,
            -1.0e-300,
            1234567890123456.7,
        ];
        for &v in &awkward {
            let mut s = String::new();
            writeln!(s, "{v:0.16e}").unwrap();
            let back: f64 = s.trim().parse().unwrap();
            assert_eq!(back.to_bits(), v.to_bits(), "{v} -> {s}");
        }
    }

    #[test]
    fn save_load_roundtrip() {
        let dir = std::env::temp_dir().join(format!(
            "modena_restart_ut_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("cp.chk");
        let cp = Checkpoint {
            continuation: None,
            atoms: 2,
            step: 91,
            time: 0.455,
            box_size: SimBox::from_array([3.0, 4.0, 5.0]),
            positions: vec![0.1 + 0.2, 1.0, 2.0, 3.0, 4.0, 5.0],
            velocities: vec![-0.5, 0.25, 0.125, 1e-300, 6.0, 7.0],
        };
        save(&path, &cp).unwrap();
        let back = load(&path).unwrap();
        assert_eq!(back.atoms, cp.atoms);
        assert_eq!(back.step, cp.step);
        assert_eq!(back.time.to_bits(), cp.time.to_bits());
        assert_eq!(back.box_size, cp.box_size);
        assert_eq!(back.positions, cp.positions);
        assert_eq!(back.velocities, cp.velocities);
        // Triclinic cell round-trips through its vectors.
        let tri = Checkpoint {
            box_size: SimBox::triclinic([
                Vec3::new(2.0, 0.1, 0.0),
                Vec3::new(0.0, 2.0, 0.2),
                Vec3::new(0.05, 0.0, 2.0),
            ]),
            ..cp
        };
        save(&path, &tri).unwrap();
        let back = load(&path).unwrap();
        assert!(back.box_size.is_triclinic());
        assert_eq!(back.box_size, tri.box_size);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn load_rejects_garbage() {
        let dir = std::env::temp_dir().join(format!(
            "modena_restart_bad_{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bad.chk");
        std::fs::write(&path, "hello world\n").unwrap();
        assert!(load(&path).is_err());
        std::fs::write(&path, "MODENA_RESTART 2\natoms 1\n").unwrap();
        assert!(load(&path).is_err());
        std::fs::remove_dir_all(&dir).ok();
    }
}
