//! Periodic constraint geometry, projection impulses, and endpoint mechanical reactions.
//! Geometry is runtime-only: a stage freezes its cell and bonded image shifts.
use crate::constraints::Constraints;
use crate::error::{ModenaError, Result};
use crate::periodic_geometry::{
    certified_cell, edge_displacement, global_shift_bound, image_bounds, CertifiedCell,
    GeometryOptions,
};
use crate::stress::SymTensor6;
use crate::{SimBox, System, Vec3};
use serde::{Deserialize, Deserializer, Serialize};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    content = "value",
    rename_all = "snake_case",
    deny_unknown_fields
)]
pub enum ProjectionStage {
    VvPosition,
    VvFinal,
    BaoabB0,
    BaoabA1Pre,
    BaoabA1Post,
    BaoabO,
    BaoabA2Pre,
    BaoabA2Post,
    BaoabB1,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConstraintEdge {
    pub i: usize,
    pub j: usize,
    #[serde(deserialize_with = "finite_scalar")]
    pub distance: f64,
    /// `r = q_j - q_i + H * image_shift`.
    pub image_shift: [i64; 3],
    pub reference: Vec3,
}
#[derive(Debug, Clone)]
pub struct ConstraintGeometry {
    pub edges: Vec<ConstraintEdge>,
    pub cell: SimBox,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConstraintImpulse {
    pub edge: usize,
    pub lever: Vec3,
    pub impulse_on_j: Vec3,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProjectionReport {
    pub stage: ProjectionStage,
    pub impulses: Vec<ConstraintImpulse>,
    pub tensor_impulse: SymTensor6,
    /// Correction passes: each water Newton solve and each nonempty loose-bond
    /// GS sweep counts once. Convergence-only checks do not count.
    pub iterations: usize,
    #[serde(deserialize_with = "finite_scalar")]
    pub max_position_residual: f64,
    #[serde(deserialize_with = "finite_scalar")]
    pub max_velocity_residual: f64,
}
fn finite_scalar<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<f64, D::Error> {
    let x = f64::deserialize(d)?;
    if x.is_finite() {
        Ok(x)
    } else {
        Err(serde::de::Error::custom("expected finite scalar"))
    }
}
pub(crate) fn geometry_error(message: impl std::fmt::Display) -> ModenaError {
    ModenaError::Physics(format!("constraint geometry: {message}"))
}
impl ProjectionReport {
    pub(crate) fn new(stage: ProjectionStage) -> Self {
        Self {
            stage,
            impulses: Vec::new(),
            tensor_impulse: SymTensor6::ZERO,
            iterations: 0,
            max_position_residual: 0.,
            max_velocity_residual: 0.,
        }
    }
    /// Finite coordinates or velocities can lose a correction to rounding even
    /// when the final constraint residual is zero. Compare each atom component
    /// to the §9 impulse tolerance using only its incident ledger increments.
    pub(crate) fn check_reconstruction(
        &self,
        geometry: &ConstraintGeometry,
        masses: &[f64],
        before: &[f64],
        after: &[f64],
        interval: f64,
    ) -> Result<()> {
        const RECONSTRUCTION_TOLERANCE: f64 = 5e-10;
        let mut impulses = vec![(Vec3::zero(), Vec3::zero()); masses.len()];
        for update in &self.impulses {
            let edge = &geometry.edges[update.edge];
            let magnitude = Vec3(update.impulse_on_j.0.map(f64::abs));
            for (atom, contribution) in [
                (edge.i, -update.impulse_on_j),
                (edge.j, update.impulse_on_j),
            ] {
                let (sum, absolute_sum) = &mut impulses[atom];
                *sum += contribution;
                *absolute_sum += magnitude;
            }
        }
        for (i, (sum, absolute_sum)) in impulses.iter().enumerate() {
            for axis in 0..3 {
                let k = 3 * i + axis;
                let represented = masses[i] * (after[k] - before[k]) / interval;
                let difference = (represented - sum.0[axis]).abs();
                let tolerance = RECONSTRUCTION_TOLERANCE * (1. + absolute_sum.0[axis]);
                if !difference.is_finite() || !tolerance.is_finite() || difference > tolerance {
                    return Err(geometry_error(format!("stage {:?}, atom {i}, axis {axis}: represented impulse {represented:e} differs from ledger {:e} (tolerance {tolerance:e})",self.stage,sum.0[axis])));
                }
            }
        }
        Ok(())
    }
    pub(crate) fn record(&mut self, edge: usize, lever: Vec3, impulse_on_j: Vec3) -> Result<()> {
        let tensor = self
            .tensor_impulse
            .plus(SymTensor6::dyad(lever, impulse_on_j));
        if !tensor.finite() || impulse_on_j.0.iter().any(|x| !x.is_finite()) {
            return Err(geometry_error("nonfinite projection impulse"));
        }
        self.impulses.push(ConstraintImpulse {
            edge,
            lever,
            impulse_on_j,
        });
        self.tensor_impulse = tensor;
        Ok(())
    }
}
/// Stable order: loose bonds, then OH1/OH2/H1H2 for each water.
pub(crate) fn topology(c: &Constraints) -> Vec<(usize, usize, f64)> {
    c.bonds
        .iter()
        .map(|e| (e.i, e.j, e.d))
        .chain(c.waters.iter().flat_map(|w| {
            [
                (w.first, w.first + 1, w.d_oh),
                (w.first, w.first + 2, w.d_oh),
                (w.first + 1, w.first + 2, w.d_hh),
            ]
        }))
        .collect()
}
/// B04 inverse-error bounds certify ALL images inside a Cartesian radius,
/// including skew cells. Reject excess work rather than truncate the search.
fn images(
    cell: &CertifiedCell,
    a: Vec3,
    b: Vec3,
    radius: f64,
    mut visit: impl FnMut([i64; 3], Vec3) -> Result<()>,
) -> Result<()> {
    let mut options = GeometryOptions {
        cutoff: radius,
        pbc: [true; 3],
        max_abs_shift: 1_i64 << 52,
    };
    options.max_abs_shift = global_shift_bound(cell, a, b, options)
        .map_err(|e| geometry_error(format!("image bound {e:?}")))?;
    let bounds = image_bounds(cell, a, b, options)
        .map_err(|e| geometry_error(format!("image bounds {e:?}")))?;
    let count = bounds.iter().try_fold(1_u64, |n, (lo, hi)| {
        n.checked_mul(u64::try_from(hi.checked_sub(*lo)?.checked_add(1)?).ok()?)
    });
    if count.is_none_or(|n| n > 1_000_000) {
        return Err(geometry_error("certified image search exceeds work limit"));
    }
    for x in bounds[0].0..=bounds[0].1 {
        for y in bounds[1].0..=bounds[1].1 {
            for z in bounds[2].0..=bounds[2].1 {
                let shift = [x, y, z];
                let r = edge_displacement(cell, a, b, shift)
                    .map_err(|e| geometry_error(format!("image displacement {e:?}")))?;
                visit(shift, r)?;
            }
        }
    }
    Ok(())
}
/// Orthorhombic minimization separates by axis. The next closest image
/// differs by exactly one adjacent integer on one axis; the gap in squared
/// length on axis k is L_k (L_k - 2 |r_k|). For long unwrapped coordinates
/// use the certified Cartesian search instead of assuming safe range reduction.
fn closest(
    sim_box: SimBox,
    certified: Option<&CertifiedCell>,
    a: Vec3,
    b: Vec3,
    distance: f64,
) -> Result<([i64; 3], Vec3)> {
    let basis = sim_box.basis();
    let raw = b - a;
    let scale = basis
        .iter()
        .map(|v| v.norm2())
        .sum::<f64>()
        .max(1.)
        .max(distance * distance);
    let tie = 64. * f64::EPSILON * scale;
    if !sim_box.is_triclinic()
        && (0..3).all(|k| a.0[k].abs().max(b.0[k].abs()) <= 8. * sim_box.size.0[k])
    {
        let reference = sim_box.min_image(raw);
        let mut shift = [0; 3];
        for (k, value) in shift.iter_mut().enumerate() {
            let length = sim_box.size.0[k];
            *value = -(raw.0[k] / length).round() as i64;
            if length * (length - 2. * reference.0[k].abs()) <= tie {
                return Err(geometry_error("ambiguous nearest constraint image"));
            }
        }
        if reference.norm2() <= 0. || !reference.norm2().is_finite() {
            return Err(geometry_error("invalid reference length"));
        }
        return Ok((shift, reference));
    }
    let fallback;
    let cell = match certified {
        Some(c) => c,
        None => {
            fallback = certify(sim_box)?;
            &fallback
        }
    };
    // Seed the B04 search with a real integer shift. Even if rounding the
    // fractional coordinates gives a suboptimal seed, its actual Cartesian
    // length bounds the true minimum. Bound the search BEFORE enumeration.
    let mut seed = [0; 3];
    for (k, slot) in seed.iter_mut().enumerate() {
        let x = -(cell.inverse[k].dot(raw)).round();
        *slot = crate::periodic_geometry::checked_shift_endpoint(x, 1_i64 << 52)
            .map_err(|e| geometry_error(format!("image seed {e:?}")))?;
    }
    let candidate = edge_displacement(cell, a, b, seed)
        .map_err(|e| geometry_error(format!("image seed {e:?}")))?;
    let radius = (candidate.norm2() + 2. * tie).sqrt().next_up();
    let mut best = f64::INFINITY;
    let mut second = f64::INFINITY;
    let mut result = ([0; 3], Vec3::zero());
    images(cell, a, b, radius, |shift, r| {
        let d2 = r.norm2();
        if !d2.is_finite() {
            return Err(geometry_error("nonfinite image length"));
        }
        if d2 < best {
            second = best;
            best = d2;
            result = (shift, r);
        } else {
            second = second.min(d2);
        }
        Ok(())
    })?;
    if !best.is_finite() || best <= 0. || second - best <= tie {
        return Err(geometry_error(
            "zero-length or ambiguous nearest constraint image",
        ));
    }
    Ok(result)
}
fn shortest_lattice(sim_box: SimBox) -> Result<f64> {
    if !sim_box.is_triclinic() {
        validate_orthorhombic(sim_box)?;
        return Ok(sim_box.size.0.into_iter().fold(f64::INFINITY, f64::min));
    }
    let cell = certify(sim_box)?;
    // A basis vector bounds the shortest vector from above. B04 enumerates
    // every vector inside that radius; no 27-image assumption is made.
    let mut shortest = cell
        .basis
        .iter()
        .map(|v| v.norm())
        .fold(f64::INFINITY, f64::min);
    images(
        &cell,
        Vec3::zero(),
        Vec3::zero(),
        shortest.next_up(),
        |shift, r| {
            if shift != [0; 3] {
                shortest = shortest.min(r.norm());
            }
            Ok(())
        },
    )?;
    if !shortest.is_finite() || shortest <= 0. {
        return Err(geometry_error("invalid shortest lattice vector"));
    }
    Ok(shortest)
}
/// Select unique closest images before prediction and validate graph closure
/// and independence. Networks with nonzero periodic winding are unsupported.
pub fn constraint_geometry(
    system: &System,
    constraints: &Constraints,
) -> Result<ConstraintGeometry> {
    geometry_from_positions(system, constraints, &system.positions)
}
fn validate_orthorhombic(cell: SimBox) -> Result<()> {
    if cell.size.0.iter().all(|x| x.is_finite() && *x > 0.)
        && cell.volume().is_finite()
        && cell.volume() > 0.
    {
        Ok(())
    } else {
        Err(geometry_error("invalid orthorhombic cell"))
    }
}
fn certify(sim_box: SimBox) -> Result<CertifiedCell> {
    certified_cell(sim_box).map_err(|e| geometry_error(format!("invalid cell {e:?}")))
}
pub(crate) fn geometry_from_positions(
    system: &System,
    constraints: &Constraints,
    positions: &[f64],
) -> Result<ConstraintGeometry> {
    let geometry = select_geometry(system, constraints, positions)?;
    check_graph_rank(&geometry, system)?;
    Ok(geometry)
}
fn select_geometry(
    system: &System,
    constraints: &Constraints,
    positions: &[f64],
) -> Result<ConstraintGeometry> {
    constraints.validate(system)?;
    if positions.len() != system.positions.len() || positions.iter().any(|x| !x.is_finite()) {
        return Err(geometry_error("invalid reference positions"));
    }
    let cell = if system.box_size.is_triclinic() {
        Some(certify(system.box_size)?)
    } else {
        validate_orthorhombic(system.box_size)?;
        None
    };
    let pos = |i: usize| Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
    let mut edges = Vec::new();
    for (i, j, distance) in topology(constraints) {
        let (image_shift, reference) =
            closest(system.box_size, cell.as_ref(), pos(i), pos(j), distance)?;
        edges.push(ConstraintEdge {
            i,
            j,
            distance,
            image_shift,
            reference,
        });
    }
    let geometry = ConstraintGeometry {
        edges,
        cell: system.box_size,
    };
    validate_connectivity(&geometry, system)?;
    Ok(geometry)
}
/// Connected blocks keep disjoint-water rank storage to 3x3 per block.
pub(crate) fn connected_blocks(edges: &[ConstraintEdge], n: usize) -> Vec<Vec<usize>> {
    let mut adjacent = vec![Vec::new(); n];
    for (k, e) in edges.iter().enumerate() {
        adjacent[e.i].push(k);
        adjacent[e.j].push(k);
    }
    let mut seen = vec![false; edges.len()];
    let mut blocks = Vec::new();
    for first in 0..edges.len() {
        if seen[first] {
            continue;
        }
        let mut block = Vec::new();
        let mut queue = VecDeque::from([first]);
        seen[first] = true;
        while let Some(k) = queue.pop_front() {
            block.push(k);
            for atom in [edges[k].i, edges[k].j] {
                for &l in &adjacent[atom] {
                    if !seen[l] {
                        seen[l] = true;
                        queue.push_back(l);
                    }
                }
            }
        }
        blocks.push(block);
    }
    blocks
}
fn validate_connectivity(geometry: &ConstraintGeometry, s: &System) -> Result<()> {
    let mut pairs = BTreeSet::new();
    let mut adjacent = vec![Vec::new(); s.n()];
    for (k, e) in geometry.edges.iter().enumerate() {
        if e.i >= s.n()
            || e.j >= s.n()
            || e.i == e.j
            || !e.distance.is_finite()
            || e.distance <= 0.
            || !(e.distance * e.distance).is_finite()
            || e.reference.0.iter().any(|x| !x.is_finite())
            || !e.reference.norm2().is_finite()
            || e.reference.norm2() <= 0.
            || e.image_shift
                .iter()
                .any(|x| x.unsigned_abs() > (1_u64 << 52))
        {
            return Err(geometry_error(format!("invalid edge {k}")));
        }
        if !pairs.insert((e.i.min(e.j), e.i.max(e.j))) {
            return Err(geometry_error(format!("duplicate edge {k}")));
        }
        adjacent[e.i].push((e.j, k, 1_i64));
        adjacent[e.j].push((e.i, k, -1_i64));
    }
    let mut shifts = vec![None; s.n()];
    let mut points = vec![Vec3::zero(); s.n()];
    for root in 0..s.n() {
        if shifts[root].is_some() {
            continue;
        }
        shifts[root] = Some([0_i64; 3]);
        let mut queue = VecDeque::from([root]);
        while let Some(i) = queue.pop_front() {
            let shift = shifts[i].unwrap();
            for &(j, k, sign) in &adjacent[i] {
                let e = &geometry.edges[k];
                let mut next = [0; 3];
                for a in 0..3 {
                    next[a] = shift[a]
                        .checked_add(sign * e.image_shift[a])
                        .ok_or_else(|| geometry_error("graph shift overflow"))?;
                }
                let point = points[i] + e.reference * (sign as f64);
                if let Some(old) = shifts[j] {
                    let scale = 1_f64.max(point.norm()).max(points[j].norm());
                    if old != next || (point - points[j]).norm() > 256. * f64::EPSILON * scale {
                        return Err(geometry_error(format!(
                            "winding or inconsistent cycle at edge {k}"
                        )));
                    }
                } else {
                    shifts[j] = Some(next);
                    points[j] = point;
                    queue.push_back(j);
                }
            }
        }
    }
    Ok(())
}
fn check_graph_rank(geometry: &ConstraintGeometry, s: &System) -> Result<()> {
    for block in connected_blocks(&geometry.edges, s.n()) {
        let mut matrix = gram_matrix(&geometry.edges, &s.masses, &block)?;
        factor_gram(&mut matrix, &block)?;
    }
    Ok(())
}
fn gram_matrix(edges: &[ConstraintEdge], masses: &[f64], block: &[usize]) -> Result<Vec<f64>> {
    let n = block.len();
    let size = n
        .checked_mul(n)
        .ok_or_else(|| geometry_error("constraint rank allocation overflow"))?;
    let mut a = Vec::new();
    a.try_reserve_exact(size)
        .map_err(|_| geometry_error("constraint rank allocation failed"))?;
    a.resize(size, 0.);
    let mut rows = BTreeMap::<usize, Vec<(usize, Vec3)>>::new();
    for (row, &k) in block.iter().enumerate() {
        let e = &edges[k];
        rows.entry(e.i).or_default().push((row, -e.reference));
        rows.entry(e.j).or_default().push((row, e.reference));
    }
    for (atom, row) in rows {
        for &(i, gi) in &row {
            for &(j, gj) in &row {
                a[i * n + j] += gi.dot(gj) / masses[atom];
            }
        }
    }
    Ok(a)
}

/// Form P D⁻¹ A D⁻¹ Pᵀ = L Lᵀ in place. The lower triangle contains L.
/// Scales stay in original order; permutation[k] is the original row at pivot k.
fn factor_gram(a: &mut [f64], block: &[usize]) -> Result<(Vec<f64>, Vec<usize>)> {
    let n = block.len();
    let mut permutation: Vec<_> = (0..n).collect();
    let scale: Vec<_> = (0..n).map(|i| a[i * n + i].sqrt()).collect();
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = (a[i * n + j] / scale[i]) / scale[j];
        }
    }
    if a.iter().any(|x| !x.is_finite()) {
        return Err(geometry_error(format!("nonfinite rank block {block:?}")));
    }
    let max_diagonal = (0..n).map(|i| a[i * n + i]).fold(0., f64::max);
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| a[i * n + i].total_cmp(&a[j * n + j]))
            .unwrap();
        let pivot = a[p * n + p];
        if pivot <= 1e-12 * max_diagonal {
            return Err(geometry_error(format!(
                "singular constraint block {block:?}, scaled pivot {pivot:e}, max diagonal {max_diagonal:e}"
            )));
        }
        for j in 0..n {
            a.swap(k * n + j, p * n + j);
        }
        for i in 0..n {
            a.swap(i * n + k, i * n + p);
        }
        permutation.swap(k, p);
        a[k * n + k] = pivot.sqrt();
        for i in k + 1..n {
            a[i * n + k] /= a[k * n + k];
        }
        for i in k + 1..n {
            for j in i..n {
                let value = a[i * n + j] - a[i * n + k] * a[j * n + k];
                if !value.is_finite() {
                    return Err(geometry_error(format!(
                        "nonfinite factorization of block {block:?}"
                    )));
                }
                a[i * n + j] = value;
                a[j * n + i] = value;
            }
        }
    }
    Ok((scale, permutation))
}
impl ConstraintGeometry {
    pub(crate) fn validate_stage(
        &self,
        s: &System,
        c: &Constraints,
        velocity: bool,
        caller_supplied: bool,
    ) -> Result<f64> {
        if caller_supplied {
            c.validate(s)?;
            if self.cell.basis() != s.box_size.basis() || self.edges.len() != c.count() {
                return Err(geometry_error("stage cell or topology mismatch"));
            }
            for (e, (i, j, d)) in self.edges.iter().zip(topology(c)) {
                if (e.i, e.j, e.distance) != (i, j, d) {
                    return Err(geometry_error("stage topology mismatch"));
                }
            }
            validate_connectivity(self, s)?;
            check_graph_rank(self, s)?;
            // Public geometry can be constructed or edited by callers. Its frozen
            // reference must itself be a unique nearest image, even for SHAKE
            // when only the predictor (not the starting coordinates) is available.
            let cell = if self.cell.is_triclinic() {
                Some(certify(self.cell)?)
            } else {
                validate_orthorhombic(self.cell)?;
                None
            };
            for e in &self.edges {
                if closest(
                    self.cell,
                    cell.as_ref(),
                    Vec3::zero(),
                    e.reference,
                    e.distance,
                )?
                .0 != [0; 3]
                {
                    return Err(geometry_error("nonminimal bonded reference image"));
                }
            }
        }
        let half = shortest_lattice(self.cell)? * 0.5;
        self.check_branch(s, half)?;
        if velocity {
            for e in &self.edges {
                let r = self.vector(s, e)?;
                if (r - e.reference).norm() > 64. * f64::EPSILON * 1_f64.max(r.norm()) {
                    return Err(geometry_error(
                        "velocity stage has stale reference geometry",
                    ));
                }
            }
        }
        Ok(half)
    }
    pub(crate) fn vector(&self, s: &System, e: &ConstraintEdge) -> Result<Vec3> {
        // Never select a new minimum image during projection.
        let mut r = s.pos(e.j) - s.pos(e.i);
        for (h, shift) in self.cell.basis().iter().zip(e.image_shift) {
            r += *h * (shift as f64);
        }
        if r.0.iter().any(|x| !x.is_finite()) || !r.norm2().is_finite() {
            return Err(geometry_error("nonfinite continuous stage vector"));
        }
        Ok(r)
    }
    pub(crate) fn check_branch(&self, s: &System, half: f64) -> Result<()> {
        for k in 0..self.edges.len() {
            self.check_edge(s, half, k)?;
        }
        Ok(())
    }
    pub(crate) fn check_edge(&self, s: &System, half: f64, k: usize) -> Result<()> {
        let e = &self.edges[k];
        let r = self.vector(s, e)?;
        let dot = r.dot(e.reference);
        if dot <= 0. || !dot.is_finite() || (r - e.reference).norm() >= half {
            return Err(geometry_error(format!(
                "stage loses bonded image at edge {k}"
            )));
        }
        Ok(())
    }
    pub(crate) fn residuals(&self, s: &System) -> Result<(f64, f64)> {
        let (mut position, mut velocity) = (0_f64, 0_f64);
        for e in &self.edges {
            let r = self.vector(s, e)?;
            let p = (r.norm2() - e.distance * e.distance).abs();
            let v = r.dot(s.vel(e.j) - s.vel(e.i)).abs();
            if !p.is_finite() || !v.is_finite() {
                return Err(geometry_error("nonfinite stage residual"));
            }
            position = position.max(p);
            velocity = velocity.max(v);
        }
        Ok((position, velocity))
    }
}

/// Instantaneous deterministic constraint forces and their unnormalized virial.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EndpointReaction {
    #[serde(deserialize_with = "finite_multipliers")]
    pub multipliers: Vec<f64>,
    pub forces: Vec<Vec3>,
    pub tensor: SymTensor6,
    #[serde(deserialize_with = "finite_scalar")]
    pub relative_linear_residual: f64,
}
fn finite_multipliers<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<Vec<f64>, D::Error> {
    let values = Vec::<f64>::deserialize(d)?;
    if values.iter().all(|x| x.is_finite()) {
        Ok(values)
    } else {
        Err(serde::de::Error::custom(
            "expected finite constraint multipliers",
        ))
    }
}

/// Solve the acceleration constraints at the supplied final positions. Geometry
/// is public, so certify its image choices against this snapshot before use.
pub fn endpoint_reaction(
    system: &System,
    geometry: &ConstraintGeometry,
) -> Result<EndpointReaction> {
    if geometry.cell.basis() != system.box_size.basis() {
        return Err(geometry_error("endpoint cell mismatch"));
    }
    let topology = Constraints {
        bonds: geometry
            .edges
            .iter()
            .map(|e| crate::constraints::Constraint {
                i: e.i,
                j: e.j,
                d: e.distance,
            })
            .collect(),
        waters: vec![],
    };
    let current = select_geometry(system, &topology, &system.positions)?;
    for (id, (e, now)) in geometry.edges.iter().zip(&current.edges).enumerate() {
        let difference = (e.reference - now.reference).norm();
        if e.image_shift != now.image_shift
            || !difference.is_finite()
            || difference > 64. * f64::EPSILON * now.reference.norm().max(1.)
        {
            return Err(geometry_error(format!(
                "stale endpoint geometry at edge {id}, reference difference {difference:e}"
            )));
        }
    }
    solve_endpoint(system, &current)
}

pub(crate) fn reaction_from_constraints(
    system: &System,
    constraints: &Constraints,
) -> Result<EndpointReaction> {
    let geometry = select_geometry(system, constraints, &system.positions)?;
    solve_endpoint(system, &geometry)
}

fn solve_endpoint(system: &System, g: &ConstraintGeometry) -> Result<EndpointReaction> {
    if system.forces.iter().any(|f| !f.is_finite()) {
        return Err(geometry_error("nonfinite conservative endpoint force"));
    }
    let mut reaction = EndpointReaction {
        multipliers: vec![0.; g.edges.len()],
        forces: vec![Vec3::zero(); system.n()],
        tensor: SymTensor6::ZERO,
        relative_linear_residual: 0.,
    };
    if g.edges.is_empty() {
        return Ok(reaction);
    }
    for block in connected_blocks(&g.edges, system.n()) {
        let mut rhs = Vec::with_capacity(block.len());
        for &id in &block {
            let e = &g.edges[id];
            let r = e.reference;
            let relative_v = system.vel(e.j) - system.vel(e.i);
            let position = (r.norm2() - e.distance * e.distance).abs();
            let tangent = r.dot(relative_v).abs();
            if !position.is_finite() || !tangent.is_finite() || position > 1e-10 || tangent > 1e-10
            {
                return Err(geometry_error(format!("endpoint block {block:?} edge {id} position residual {position:e}, tangent residual {tangent:e}")));
            }
            let acceleration = |i: usize| Vec3(system.force(i).0.map(|f| f / system.masses[i]));
            let b = relative_v.norm2() + r.dot(acceleration(e.j) - acceleration(e.i));
            if !b.is_finite() {
                return Err(geometry_error(format!(
                    "nonfinite acceleration in endpoint block {block:?} edge {id}"
                )));
            }
            rhs.push(b);
        }
        let a = gram_matrix(&g.edges, &system.masses, &block)?;
        let (lambda, residual) = solve_block(&a, &rhs, &block)?;
        reaction.relative_linear_residual = reaction.relative_linear_residual.max(residual);
        for (&id, multiplier) in block.iter().zip(lambda) {
            let e = &g.edges[id];
            reaction.multipliers[id] = multiplier;
            let on_j = e.reference * -multiplier;
            reaction.forces[e.i] -= on_j;
            reaction.forces[e.j] += on_j;
            reaction.tensor = reaction.tensor.plus(SymTensor6::dyad(e.reference, on_j));
            if !reaction.tensor.finite()
                || [e.i, e.j]
                    .iter()
                    .any(|&i| reaction.forces[i].0.iter().any(|f| !f.is_finite()))
            {
                return Err(geometry_error(format!(
                    "nonfinite reaction in block {block:?} at edge {id}"
                )));
            }
        }
    }
    Ok(reaction)
}

fn solve_block(a: &[f64], rhs: &[f64], block: &[usize]) -> Result<(Vec<f64>, f64)> {
    let n = block.len();
    // Keep the unscaled A and b for an independent residual after undoing
    // both the permutation and the scaling.
    let mut l = Vec::new();
    l.try_reserve_exact(a.len())
        .map_err(|_| geometry_error("constraint factor allocation failed"))?;
    l.extend_from_slice(a);
    let (scale, perm) = factor_gram(&mut l, block)?;
    let mut x = vec![0.; n];
    for i in 0..n {
        let inner: f64 = (0..i).map(|j| l[i * n + j] * x[j]).sum();
        x[i] = (rhs[perm[i]] / scale[perm[i]] - inner) / l[i * n + i];
    }
    for i in (0..n).rev() {
        let inner: f64 = (i + 1..n).map(|j| l[j * n + i] * x[j]).sum();
        x[i] = (x[i] - inner) / l[i * n + i];
    }
    let mut lambda = vec![0.; n];
    for (i, &row) in perm.iter().enumerate() {
        lambda[row] = x[i] / scale[row];
    }
    let residual = linear_residual(a, rhs, &lambda, block)?;
    Ok((lambda, residual))
}

/// Normwise backward error ||Aλ-b||∞ / (||A||∞ ||λ||∞ + ||b||∞).
fn linear_residual(a: &[f64], rhs: &[f64], lambda: &[f64], block: &[usize]) -> Result<f64> {
    let mut numerator = 0_f64;
    let mut matrix_norm = 0_f64;
    for (row, b) in a.chunks_exact(block.len()).zip(rhs) {
        let error = (row.iter().zip(lambda).map(|(a, x)| a * x).sum::<f64>() - b).abs();
        let row_sum = row.iter().map(|a| a.abs()).sum::<f64>();
        if !error.is_finite() || !row_sum.is_finite() {
            return Err(geometry_error(format!(
                "nonfinite linear residual in endpoint block {block:?}"
            )));
        }
        numerator = numerator.max(error);
        matrix_norm = matrix_norm.max(row_sum);
    }
    let norm = |x: &[f64]| x.iter().map(|v| v.abs()).fold(0., f64::max);
    let denominator = matrix_norm * norm(lambda) + norm(rhs);
    let relative = if denominator == 0. && numerator == 0. {
        0.
    } else {
        numerator / denominator
    };
    if !denominator.is_finite() || !relative.is_finite() || relative > 1e-11 {
        return Err(geometry_error(format!("endpoint block {block:?} linear residual {relative:e} (numerator {numerator:e}, denominator {denominator:e})")));
    }
    Ok(relative)
}

#[cfg(test)]
mod endpoint_tests {
    use super::linear_residual;

    #[test]
    fn backward_error_guard_rejects_wrong_or_nonfinite_solution() {
        let a = [2., 0., 0., 1.];
        let rhs = [4., 1.];
        let block = [7, 9];
        assert_eq!(linear_residual(&a, &rhs, &[2., 1.], &block).unwrap(), 0.);
        let err = linear_residual(&a, &rhs, &[2., 1.001], &block).unwrap_err();
        assert!(err.to_string().contains("[7, 9] linear residual"));
        assert!(linear_residual(&a, &rhs, &[f64::NAN, 1.], &block).is_err());
        assert!(linear_residual(&[f64::MAX], &[f64::MAX], &[1.], &[1]).is_err());
    }
}
