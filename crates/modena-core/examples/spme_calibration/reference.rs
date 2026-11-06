//! Independent all-image real-space / full ±k Ewald, protocol §3.
//! No production force, neighbor-list, pair-policy or minimum-image implementation is used.
use super::{
    invalid, CalibrationFailure, Evaluation, Limits, ReferenceSettings, Snapshot, WorkCounts,
};
use std::{
    f64::consts::PI,
    mem::size_of,
    time::{Duration, Instant},
};

fn cap(message: impl Into<String>) -> CalibrationFailure {
    CalibrationFailure::BudgetExceeded(message.into())
}
fn precision() -> CalibrationFailure {
    CalibrationFailure::ReferencePrecisionLimited("nonfinite independent Ewald intermediate".into())
}
fn finite(x: f64) -> Result<f64, CalibrationFailure> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(precision())
    }
}
fn add(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_add(b)
        .ok_or_else(|| cap("reference work/heap addition overflow"))
}
fn mul(a: u64, b: u64) -> Result<u64, CalibrationFailure> {
    a.checked_mul(b)
        .ok_or_else(|| cap("reference work/heap product overflow"))
}
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    // Neumaier compensation, including force components and charge/structure sums.
    fn push(&mut self, x: f64) -> Result<(), CalibrationFailure> {
        finite(x)?;
        let t = finite(self.value + x)?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= x.abs() {
                    (self.value - t) + x
                } else {
                    (x - t) + self.value
                },
        )?;
        self.value = t;
        Ok(())
    }
    fn get(self) -> Result<f64, CalibrationFailure> {
        finite(self.value + self.correction)
    }
}

/// Serializable per-call diagnostics. Upper bounds and actual work are distinct.
/// `real_terms` counts examined ordered nonzero-self image tuples (including rejected radii).
/// `atom_k_terms` counts actual atom visits in the rho and force passes for included modes.
/// Preflight uses twice N times the entire nonzero reciprocal bounding rectangle.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceAudit {
    pub real_bounds: [i32; 3],
    pub reciprocal_bounds: [i32; 3],
    pub prechecked_work: WorkCounts,
    pub work: WorkCounts,
    pub real_included_terms: u64,
    pub reciprocal_vectors: u64,
    pub reciprocal_enumerated_vectors: u64,
    pub real_membership_hash: u64,
    pub reciprocal_membership_hash: u64,
    pub owned_heap_bytes: u64,
    pub deadline_checks: u64,
    #[serde(
        serialize_with = "super::serialize_finite",
        deserialize_with = "super::deserialize_finite"
    )]
    pub background_energy: f64,
    #[serde(
        serialize_with = "super::serialize_finite",
        deserialize_with = "super::deserialize_finite"
    )]
    pub background_w: f64,
}

/// Reuse one budget across a family. Resume from retained actual counters/elapsed time;
/// retained snapshots/results outside this evaluator belong in `owned_heap_baseline`.
/// The convenience `reference` has a fresh budget; family orchestration must use this context.
pub struct ReferenceBudget {
    pub work: WorkCounts,
    pub calls: u64,
    pub owned_heap_baseline: u64,
    limits: Limits,
    started: Instant,
    elapsed_before: Duration,
    deadline_checks: u64,
    ticks: u16,
}
impl ReferenceBudget {
    pub fn new(limits: &Limits) -> Self {
        Self {
            work: WorkCounts {
                real_terms: 0,
                atom_k_terms: 0,
            },
            calls: 0,
            owned_heap_baseline: 0,
            limits: limits.clone(),
            started: Instant::now(),
            elapsed_before: Duration::ZERO,
            deadline_checks: 0,
            ticks: 0,
        }
    }
    pub fn resume(
        limits: &Limits,
        work: WorkCounts,
        calls: u64,
        elapsed_before: Duration,
        owned_heap_baseline: u64,
    ) -> Result<Self, CalibrationFailure> {
        let mut result = Self {
            work,
            calls,
            owned_heap_baseline,
            elapsed_before,
            ..Self::new(limits)
        };
        result.checkpoint()?;
        if result.work.real_terms > limits.max_real_terms_per_family
            || result.work.atom_k_terms > limits.max_atom_k_terms_per_family
            || result.calls > limits.max_accuracy_calls
            || result.owned_heap_baseline > limits.max_owned_heap_bytes
        {
            return Err(cap("resumed reference family exceeds retained caps"));
        }
        Ok(result)
    }
    pub fn elapsed(&self) -> Duration {
        self.elapsed_before.saturating_add(self.started.elapsed())
    }
    pub fn checkpoint(&mut self) -> Result<(), CalibrationFailure> {
        self.deadline_checks = add(self.deadline_checks, 1)?;
        if self.elapsed() >= Duration::from_secs(self.limits.max_family_wall_seconds) {
            return Err(cap("cumulative reference family wall deadline"));
        }
        Ok(())
    }
    fn tick(&mut self) -> Result<(), CalibrationFailure> {
        self.ticks += 1;
        if self.ticks == 1024 {
            self.ticks = 0;
            self.checkpoint()?;
        }
        Ok(())
    }
}

fn bounds(values: [f64; 3]) -> Result<[i32; 3], CalibrationFailure> {
    let mut out = [0; 3];
    for (i, v) in values.into_iter().enumerate() {
        if !v.is_finite() || v < 0. || v > f64::from(i32::MAX - 1) {
            return Err(cap("reference enumeration bound overflow"));
        }
        out[i] = v as i32;
    }
    Ok(out)
}
fn rectangle(b: [i32; 3]) -> Result<u64, CalibrationFailure> {
    b.into_iter().try_fold(1, |product, axis| {
        mul(product, add(mul(axis as u64, 2)?, 1)?)
    })
}
fn buffer<T: Clone>(n: usize, value: T) -> Result<Vec<T>, CalibrationFailure> {
    let mut out = Vec::new();
    out.try_reserve_exact(n)
        .map_err(|e| cap(format!("reference allocation: {e}")))?;
    out.resize(n, value);
    Ok(out)
}
fn hash_word(hash: &mut u64, value: u64) {
    for byte in value.to_le_bytes() {
        *hash = (*hash ^ u64::from(byte)).wrapping_mul(1_099_511_628_211);
    }
}
fn input_heap(s: &Snapshot) -> Result<u64, CalibrationFailure> {
    let mut bytes = 0;
    for (capacity, width) in [
        (s.positions.capacity(), size_of::<[f64; 3]>()),
        (s.charges.capacity(), size_of::<f64>()),
        (s.scales.capacity(), size_of::<super::Scale>()),
        (s.id.capacity(), 1),
        (s.family.capacity(), 1),
        (s.units.capacity(), 1),
    ] {
        bytes = add(bytes, mul(capacity as u64, width as u64)?)?;
    }
    Ok(bytes)
}
/// Basic physics validation deliberately accepts internal FD, strain, and raw wraps.
/// Calibration eligibility remains enforced by the separate exact frozen/hash boundary.
fn physical_validation(
    s: &Snapshot,
    limits: &Limits,
    budget: &mut ReferenceBudget,
) -> Result<(), CalibrationFailure> {
    if s.units != "reduced" || s.pbc != [true; 3] {
        return Err(CalibrationFailure::UnsupportedDomain(format!(
            "{}: reference requires reduced units and full orthorhombic PBC",
            s.id
        )));
    }
    let n = s.positions.len();
    if n == 0 || n != s.charges.len() {
        return Err(invalid(format!(
            "{}: nonempty matching coordinates/charges required",
            s.id
        )));
    }
    if n > limits.max_atoms {
        return Err(cap(format!("{}: reference atom cap", s.id)));
    }
    if s.cell.iter().any(|x| !x.is_finite() || *x <= 0.)
        || s.positions
            .iter()
            .flatten()
            .chain(&s.charges)
            .any(|x| !x.is_finite())
        || s.charges.iter().all(|q| *q == 0.)
    {
        return Err(invalid(format!(
            "{}: invalid finite reference physical inputs",
            s.id
        )));
    }
    let max_pairs = mul(n as u64, n as u64 - 1)? / 2;
    if s.scales.len() as u64 > max_pairs {
        return Err(invalid(format!(
            "{}: too many reference pair policies",
            s.id
        )));
    }
    let mut last = None;
    for c in &s.scales {
        budget.tick()?;
        if c.i >= c.j
            || c.j >= n
            || !c.coulomb.is_finite()
            || !(0.0..=1.0).contains(&c.coulomb)
            || last.is_some_and(|key| key >= (c.i, c.j))
        {
            return Err(invalid(format!(
                "{}: reference requires canonical unique valid pair policies",
                s.id
            )));
        }
        last = Some((c.i, c.j));
    }
    Ok(())
}

pub fn reference(
    s: &Snapshot,
    p: ReferenceSettings,
    limits: &Limits,
) -> Result<Evaluation, CalibrationFailure> {
    reference_with_budget(s, p, limits, &mut ReferenceBudget::new(limits)).map(|(value, _)| value)
}

pub fn reference_with_budget(
    s: &Snapshot,
    p: ReferenceSettings,
    limits: &Limits,
    budget: &mut ReferenceBudget,
) -> Result<(Evaluation, ReferenceAudit), CalibrationFailure> {
    if limits != &budget.limits {
        return Err(invalid("reference budget limits changed within family"));
    }
    let checks_before = budget.deadline_checks;
    budget.checkpoint()?;
    physical_validation(s, limits, budget)?;
    super::validate_reference_settings(&p)?;
    let n = s.positions.len();
    let real_bounds = bounds(s.cell.map(|l| (p.real_cutoff / l).ceil() + 1.))?;
    let reciprocal_bounds = bounds(s.cell.map(|l| (p.kspace_cutoff * l / (2. * PI)).ceil()))?;
    let estimated = WorkCounts {
        real_terms: mul(mul(n as u64, n as u64)?, rectangle(real_bounds)?)? - n as u64,
        atom_k_terms: mul(mul(2, n as u64)?, rectangle(reciprocal_bounds)? - 1)?,
    };
    let mut heap = add(input_heap(s)?, budget.owned_heap_baseline)?;
    let baseline = heap;
    // Simultaneously live wrapped coordinates, compensated forces, and output forces.
    heap = add(
        heap,
        mul(
            n as u64,
            (size_of::<[f64; 3]>() * 2 + size_of::<[Sum; 3]>()) as u64,
        )?,
    )?;
    super::check_work([1; 3], &estimated, heap, limits)?;
    if add(budget.work.real_terms, estimated.real_terms)? > limits.max_real_terms_per_family
        || add(budget.work.atom_k_terms, estimated.atom_k_terms)?
            > limits.max_atom_k_terms_per_family
        || add(budget.calls, 1)? > limits.max_accuracy_calls
    {
        return Err(cap("reference prechecked cumulative family work/call cap"));
    }
    let mut positions = buffer(n, [0.; 3])?;
    let mut forces = buffer(n, [Sum::default(); 3])?;
    let mut output = buffer(n, [0.; 3])?;
    heap = add(
        baseline,
        add(
            mul(
                add(positions.capacity() as u64, output.capacity() as u64)?,
                size_of::<[f64; 3]>() as u64,
            )?,
            mul(forces.capacity() as u64, size_of::<[Sum; 3]>() as u64)?,
        )?,
    )?;
    if heap > limits.max_owned_heap_bytes {
        return Err(cap(
            "reference actual vector capacities exceed owned heap cap",
        ));
    }
    let before = budget.work.clone();
    budget.calls = add(budget.calls, 1)?;
    for (dst, src) in positions.iter_mut().zip(&s.positions) {
        for a in 0..3 {
            dst[a] = finite(src[a].rem_euclid(s.cell[a]))?;
            if dst[a] == s.cell[a] {
                dst[a] = 0.;
            }
        }
    }
    for i in 0..n {
        for j in 0..i {
            budget.tick()?;
            let mut r2 = 0.;
            for (a, edge) in s.cell.iter().enumerate() {
                let d = positions[i][a] - positions[j][a];
                let d = d - edge * (d / edge).round();
                r2 += d * d;
            }
            if finite(r2)? == 0. {
                return Err(invalid(format!(
                    "{}: overlapping reference particles",
                    s.id
                )));
            }
        }
    }
    let mut audit = ReferenceAudit {
        real_bounds,
        reciprocal_bounds,
        prechecked_work: estimated,
        work: WorkCounts {
            real_terms: 0,
            atom_k_terms: 0,
        },
        real_included_terms: 0,
        reciprocal_vectors: 0,
        reciprocal_enumerated_vectors: 0,
        real_membership_hash: 14_695_981_039_346_656_037,
        reciprocal_membership_hash: 14_695_981_039_346_656_037,
        owned_heap_bytes: heap,
        deadline_checks: 0,
        background_energy: 0.,
        background_w: 0.,
    };
    let mut energy = Sum::default();
    let mut w = Sum::default();
    let kappa2 = finite(p.kappa * p.kappa)?;
    let rcut2 = finite(p.real_cutoff * p.real_cutoff)?;
    let kcut2 = finite(p.kspace_cutoff * p.kspace_cutoff)?;
    let screening = 2. * p.kappa / PI.sqrt();
    // Ordered i,j,n: energy/W have 1/2, force on i has no 1/2. Include self images.
    for i in 0..n {
        for j in 0..n {
            let qq = finite(s.charges[i] * s.charges[j])?;
            for x in -real_bounds[0]..=real_bounds[0] {
                for y in -real_bounds[1]..=real_bounds[1] {
                    for z in -real_bounds[2]..=real_bounds[2] {
                        budget.tick()?;
                        if i == j && [x, y, z] == [0; 3] {
                            continue;
                        }
                        budget.work.real_terms += 1;
                        let image = [x, y, z];
                        let d: [f64; 3] = std::array::from_fn(|a| {
                            positions[i][a] - positions[j][a] + f64::from(image[a]) * s.cell[a]
                        });
                        let r2 = finite(d.iter().map(|v| v * v).sum())?;
                        if r2 >= rcut2 {
                            continue;
                        }
                        if r2 == 0. {
                            return Err(precision());
                        }
                        audit.real_included_terms += 1;
                        for value in [i as u64, j as u64, x as u64, y as u64, z as u64] {
                            hash_word(&mut audit.real_membership_hash, value);
                        }
                        let r = r2.sqrt();
                        let erfc = libm::erfc(p.kappa * r);
                        let gaussian = (-kappa2 * r2).exp();
                        energy.push(0.5 * qq * erfc / r)?;
                        w.push(0.5 * qq * (erfc / r + screening * gaussian))?;
                        let factor = qq * (erfc / (r2 * r) + screening * gaussian / r2);
                        for (a, da) in d.iter().enumerate() {
                            forces[i][a].push(factor * da)?;
                        }
                    }
                }
            }
        }
    }
    budget.checkpoint()?;
    let volume = finite(s.cell.iter().product())?;
    if volume <= 0. {
        return Err(precision());
    }
    let c = finite(2. * PI / volume)?;
    for x in -reciprocal_bounds[0]..=reciprocal_bounds[0] {
        for y in -reciprocal_bounds[1]..=reciprocal_bounds[1] {
            for z in -reciprocal_bounds[2]..=reciprocal_bounds[2] {
                budget.tick()?;
                if [x, y, z] == [0; 3] {
                    continue;
                }
                audit.reciprocal_enumerated_vectors += 1;
                let mode = [x, y, z];
                let k: [f64; 3] = std::array::from_fn(|a| 2. * PI * f64::from(mode[a]) / s.cell[a]);
                let k2 = finite(k.iter().map(|v| v * v).sum())?;
                if k2 >= kcut2 {
                    continue;
                }
                if k2 == 0. {
                    return Err(precision());
                }
                audit.reciprocal_vectors += 1;
                for m in mode {
                    hash_word(&mut audit.reciprocal_membership_hash, m as u64);
                }
                let mut re = Sum::default();
                let mut im = Sum::default();
                for (r, q) in positions.iter().zip(&s.charges) {
                    budget.tick()?;
                    budget.work.atom_k_terms += 1;
                    let phase = finite((0..3).map(|a| k[a] * r[a]).sum())?;
                    re.push(q * phase.cos())?;
                    im.push(q * phase.sin())?;
                }
                let re = re.get()?;
                let im = im.get()?;
                let weight = finite(c * (-k2 / (4. * kappa2)).exp() / k2)?;
                let mode_energy = finite(weight * (re * re + im * im))?;
                energy.push(mode_energy)?;
                w.push(mode_energy * (1. - k2 / (2. * kappa2)))?;
                for (i, (r, q)) in positions.iter().zip(&s.charges).enumerate() {
                    budget.tick()?;
                    budget.work.atom_k_terms += 1;
                    let phase = finite((0..3).map(|a| k[a] * r[a]).sum())?;
                    let factor = 2. * weight * q * (re * phase.sin() - im * phase.cos());
                    for (a, ka) in k.iter().enumerate() {
                        forces[i][a].push(factor * ka)?;
                    }
                }
            }
        }
    }
    let mut q = Sum::default();
    let mut q2 = Sum::default();
    for charge in &s.charges {
        budget.tick()?;
        q.push(*charge)?;
        q2.push(charge * charge)?;
    }
    energy.push(-p.kappa / PI.sqrt() * q2.get()?)?;
    let q = q.get()?;
    audit.background_energy = finite(-PI * q * q / (2. * volume * kappa2))?;
    audit.background_w = finite(3. * audit.background_energy)?;
    energy.push(audit.background_energy)?;
    w.push(audit.background_w)?;
    // Bare minimum-image correction once per canonical policy, including outside R.
    for scale in &s.scales {
        budget.tick()?;
        let d: [f64; 3] = std::array::from_fn(|a| {
            let v = positions[scale.i][a] - positions[scale.j][a];
            v - s.cell[a] * (v / s.cell[a]).round()
        });
        let r2 = finite(d.iter().map(|v| v * v).sum())?;
        if r2 == 0. {
            return Err(precision());
        }
        let removed = (1. - scale.coulomb) * s.charges[scale.i] * s.charges[scale.j] / r2.sqrt();
        energy.push(-removed)?;
        w.push(-removed)?;
        for (a, da) in d.iter().enumerate() {
            let f = removed * da / r2;
            forces[scale.i][a].push(-f)?;
            forces[scale.j][a].push(f)?;
        }
    }
    for (dst, src) in output.iter_mut().zip(forces) {
        for a in 0..3 {
            dst[a] = src[a].get()?;
        }
    }
    budget.checkpoint()?;
    audit.work = WorkCounts {
        real_terms: budget.work.real_terms - before.real_terms,
        atom_k_terms: budget.work.atom_k_terms - before.atom_k_terms,
    };
    audit.deadline_checks = budget.deadline_checks - checks_before;
    Ok((
        Evaluation {
            energy: energy.get()?,
            strain_w: w.get()?,
            forces: output,
        },
        audit,
    ))
}
