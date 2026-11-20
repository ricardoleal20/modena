//! SHAKE (position) and RATTLE (velocity) constraints for fixed bond lengths.
//!
//! A `Constraints` set holds rigid distance constraints `|r_ij| = d` (typically
//! the X–H bonds, or all bonds in a rigid-molecule model). The integrator
//! applies them after the unconstrained position drift ([`Constraints::shake`])
//! and after the final velocity half-kick ([`Constraints::rattle`]).
//!
//! SHAKE solves the position constraints iteratively (Gauss–Seidel over the
//! constraint list, `max_iter` sweeps) using the linearized Lagrange-multiplier
//! update
//!
//! ```text
//!   g = (|r_ij'|² − d²) / (2 (1/m_i + 1/m_j) (r_ij' · r_ij(t)))
//!   Δx_i = +(g/m_i) r_ij(t),   Δx_j = −(g/m_j) r_ij(t)
//! ```
//!
//! where `r_ij'` is the predicted (post-drift) bond vector and `r_ij(t)` the
//! reference vector at the start of the step. RATTLE then projects the
//! velocities onto the constraint tangent (`r_ij · v_ij = 0`) with
//!
//! ```text
//!   μ = (r_ij · v_ij) / ((1/m_i + 1/m_j) |r_ij|²)
//!   Δv_i = +(μ/m_i) r_ij,   Δv_j = −(μ/m_j) r_ij
//! ```
//!
//! For a single, non-overlapping set of constraints this converges in a handful
//! of iterations; the tolerance and iteration cap are configurable. The
//! combination SHAKE+RATTLE on velocity-Verlet is the symplectic RATTLE scheme
//! (Andersen, 1983) — it conserves a shadow Hamiltonian on the constraint
//! manifold, verified by the constrained-NVE test in `tests/bonded.rs`.

use crate::constraint_stress::{
    constraint_geometry, geometry_from_positions, ConstraintGeometry, ProjectionReport,
    ProjectionStage,
};
use crate::error::{ModenaError, Result};
use serde::{Deserialize, Serialize};

use crate::system::{SimBox, System};
use crate::vec3::Vec3;

/// One rigid distance constraint: atoms `i`, `j` held at fixed length `d`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Constraint {
    pub i: usize,
    pub j: usize,
    pub d: f64,
}

/// One rigid water molecule: the atoms at `first`, `first + 1`, `first + 2`
/// (O, H, H — a selected contiguous triple; unrelated atoms may coexist) held at
/// O–H length `d_oh` and H–H length `d_hh` (SPC/E: 1.0 Å / 1.633 Å). Solved
/// per molecule with a Newton solve of the three coupled constraints
/// (M-SHAKE class). Drift corrections follow saved reference directions to
/// form a RATTLE step; standalone initialization uses current directions.
/// Singular or unconverged projections return errors; large drifts can fail.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RigidWater {
    pub first: usize,
    pub d_oh: f64,
    pub d_hh: f64,
}

/// Newton cap for the per-molecule rigid-water projection. From a violated
/// state one iteration lands at O(violation²); three reach machine precision
/// for any physically attainable distortion — the cap is a runaway backstop.
const WATER_NEWTON_ITER: usize = 20;

/// A set of distance constraints applied by the integrator: loose SHAKE/RATTLE
/// `bonds` plus rigid `waters` solved exactly per molecule.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Constraints {
    pub bonds: Vec<Constraint>,
    pub waters: Vec<RigidWater>,
}

pub(crate) fn validate_state_buffers(system: &System) -> Result<()> {
    let n = system.n();
    if system.positions.len() != 3 * n
        || system.velocities.len() != 3 * n
        || system.forces.len() != 3 * n
    {
        return Err(constraint_error(
            "system buffer lengths do not match atom count",
        ));
    }
    if system.masses.iter().any(|m| !m.is_finite() || *m <= 0.0) {
        return Err(constraint_error("masses must be finite and positive"));
    }
    if system
        .positions
        .iter()
        .chain(&system.velocities)
        .any(|x| !x.is_finite())
    {
        return Err(constraint_error("positions/velocities must be finite"));
    }
    Ok(())
}

impl Constraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.bonds.is_empty() && self.waters.is_empty()
    }

    /// Number of independent constraints, assuming the provided distances are
    /// nonredundant. A nonlinear rigid water removes three kinetic DOF.
    pub fn count(&self) -> usize {
        self.bonds.len() + 3 * self.waters.len()
    }

    pub(crate) fn validate_parameters(&self) -> Result<()> {
        if self
            .bonds
            .iter()
            .any(|c| c.i == c.j || !c.d.is_finite() || c.d <= 0.0 || !(c.d * c.d).is_finite())
            || self.waters.iter().any(|w| {
                !w.d_oh.is_finite()
                    || !w.d_hh.is_finite()
                    || w.d_oh <= 0.0
                    || w.d_hh <= 0.0
                    || w.d_hh >= 2.0 * w.d_oh
                    || !(w.d_oh * w.d_oh).is_finite()
                    || !(w.d_hh * w.d_hh).is_finite()
            })
        {
            return Err(constraint_error("invalid constraint geometry parameters"));
        }
        Ok(())
    }

    /// Validate buffers, atom indices, positive masses and target geometry.
    pub fn validate(&self, system: &System) -> Result<()> {
        self.validate_parameters()?;
        validate_state_buffers(system)?;
        let n = system.n();
        for c in &self.bonds {
            if c.i >= n || c.j >= n || c.i == c.j || !c.d.is_finite() || c.d <= 0.0 {
                return Err(constraint_error("invalid bond indices or target length"));
            }
        }
        let mut occupied = vec![false; n];
        for w in &self.waters {
            if w.first.checked_add(2).is_none_or(|end| end >= n)
                || !w.d_oh.is_finite()
                || !w.d_hh.is_finite()
                || w.d_oh <= 0.0
                || w.d_hh <= 0.0
                || w.d_hh >= 2.0 * w.d_oh
            {
                return Err(constraint_error(
                    "invalid water indices or triangle geometry",
                ));
            }
            for used in &mut occupied[w.first..w.first + 3] {
                if *used {
                    return Err(constraint_error("overlapping rigid waters"));
                }
                *used = true;
            }
        }
        Ok(())
    }

    /// Geometric initialization projection of selected water triples,
    /// preserving their mass-weighted COM. Dynamics must use `shake` with its
    /// saved pre-drift reference directions instead of this current-normal solve.
    /// Water indices need not cover all atoms or start at atom zero.
    pub fn settle(&self, system: &mut System, tol: f64) -> Result<()> {
        self.validate(system)?;
        validate_tolerance(tol)?;
        let box_size = system.box_size;
        for &w in &self.waters {
            settle_one(
                &box_size,
                &mut system.positions[3 * w.first..3 * w.first + 9],
                &system.masses,
                w,
                tol,
            )?;
        }
        check_residual(
            self.water_residual(system),
            tol,
            "water position projection",
        )
    }

    /// Project velocities onto each rigid water's tangent space.
    pub fn settle_velocities(&self, system: &mut System, tol: f64) -> Result<()> {
        self.validate(system)?;
        validate_tolerance(tol)?;
        let box_size = system.box_size;
        for &w in &self.waters {
            settle_vel_one(
                &box_size,
                &system.positions,
                &system.masses,
                &mut system.velocities[3 * w.first..3 * w.first + 9],
                w,
                tol,
            )?;
        }
        check_residual(
            self.water_velocity_residual(system),
            tol,
            "water velocity projection",
        )
    }

    /// SHAKE positions, returning an error on a singular solve or any residual
    /// above the absolute length-squared tolerance after `max_iter` sweeps.
    pub fn shake(
        &self,
        system: &mut System,
        x_ref: &[f64],
        tol: f64,
        max_iter: usize,
    ) -> Result<()> {
        let geometry = geometry_from_positions(system, self, x_ref)?;
        self.shake_at_geometry(system, &geometry, tol, max_iter)
    }

    /// Use a pre-drift geometry without recording a ledger on the legacy path.
    pub(crate) fn shake_at_geometry(
        &self,
        system: &mut System,
        geometry: &ConstraintGeometry,
        tol: f64,
        max_iter: usize,
    ) -> Result<()> {
        self.project(
            system,
            geometry,
            Projection::Position { interval: 1. },
            tol,
            max_iter,
            None,
        )
    }

    /// Velocity projection using the same physical solver as `rattle_report`,
    /// without allocating or accumulating an impulse ledger.
    pub fn rattle(&self, system: &mut System, tol: f64, max_iter: usize) -> Result<()> {
        let geometry = constraint_geometry(system, self)?;
        self.project(system, &geometry, Projection::Velocity, tol, max_iter, None)
    }

    /// Project an unwrapped predictor using pre-drift geometry. Leaves the
    /// correction unwrapped so the caller can apply delta-q/interval to velocity
    /// before wrapping once. On failure restores the input position buffer.
    pub fn shake_report(
        &self,
        system: &mut System,
        geometry: &ConstraintGeometry,
        interval: f64,
        stage: ProjectionStage,
        tol: f64,
        max_iter: usize,
    ) -> Result<ProjectionReport> {
        self.project_report(
            system,
            geometry,
            Projection::Position { interval },
            stage,
            tol,
            max_iter,
        )
    }

    /// Project only velocities at fixed geometry. A stale lever is rejected.
    /// On failure restores input velocities and discards the incomplete ledger.
    pub fn rattle_report(
        &self,
        system: &mut System,
        geometry: &ConstraintGeometry,
        stage: ProjectionStage,
        tol: f64,
        max_iter: usize,
    ) -> Result<ProjectionReport> {
        self.project_report(system, geometry, Projection::Velocity, stage, tol, max_iter)
    }

    fn project_report(
        &self,
        system: &mut System,
        geometry: &ConstraintGeometry,
        mode: Projection,
        stage: ProjectionStage,
        tol: f64,
        max_iter: usize,
    ) -> Result<ProjectionReport> {
        mode.validate_stage(stage)?;
        let before = mode.buffer(system).to_vec();
        let mut report = ProjectionReport::new(stage);
        let result = self
            .project(system, geometry, mode, tol, max_iter, Some(&mut report))
            .and_then(|()| {
                report.check_reconstruction(
                    geometry,
                    &system.masses,
                    &before,
                    mode.buffer(system),
                    mode.interval(),
                )
            });
        if let Err(error) = result {
            match mode {
                Projection::Position { .. } => system.positions = before,
                Projection::Velocity => system.velocities = before,
            }
            return Err(error);
        }
        Ok(report)
    }

    fn validate_projection(
        &self,
        system: &System,
        geometry: &ConstraintGeometry,
        mode: Projection,
        tol: f64,
        max_iter: usize,
        caller_supplied: bool,
    ) -> Result<f64> {
        validate_tolerance(tol)?;
        if max_iter == 0 {
            return Err(constraint_error("iteration cap must be positive"));
        }
        if let Projection::Position { interval } = mode {
            if !interval.is_finite() || interval <= 0. {
                return Err(constraint_error(
                    "drift interval must be finite and positive",
                ));
            }
        }
        geometry.validate_stage(
            system,
            self,
            matches!(mode, Projection::Velocity),
            caller_supplied,
        )
    }

    fn project(
        &self,
        system: &mut System,
        geometry: &ConstraintGeometry,
        mode: Projection,
        tol: f64,
        max_iter: usize,
        mut report: Option<&mut ProjectionReport>,
    ) -> Result<()> {
        let half =
            self.validate_projection(system, geometry, mode, tol, max_iter, report.is_some())?;
        // Only edges incident on a moved atom can lose their continuous image.
        let guard = ProjectionGuard::new(geometry, mode, half, system.n());
        for _ in 0..max_iter {
            // Keep the existing water-before-loose-bonds update order.
            for (water_index, &water) in self.waters.iter().enumerate() {
                let first_edge = self.bonds.len() + 3 * water_index;
                self.project_water(
                    system,
                    geometry,
                    first_edge,
                    water,
                    mode,
                    tol,
                    &guard,
                    &mut report,
                )?;
            }
            let mut changed = false;
            for (index, e) in geometry.edges[..self.bonds.len()].iter().enumerate() {
                let r = geometry.vector(system, e)?;
                let reference = match mode {
                    Projection::Position { .. } => e.reference,
                    Projection::Velocity => r,
                };
                let (residual, denom) = match mode {
                    Projection::Position { .. } => (
                        r.norm2() - e.distance * e.distance,
                        2. * (1. / system.masses[e.i] + 1. / system.masses[e.j]) * r.dot(reference),
                    ),
                    Projection::Velocity => (
                        r.dot(system.vel(e.j) - system.vel(e.i)),
                        (1. / system.masses[e.i] + 1. / system.masses[e.j]) * r.norm2(),
                    ),
                };
                if !denom.is_finite() || denom.abs() < 1e-20 || !residual.is_finite() {
                    return Err(constraint_error("singular or nonfinite projection solve"));
                }
                if residual.abs() <= tol {
                    continue;
                }
                let multiplier = residual / denom;
                let on_i = reference * (multiplier / system.masses[e.i]);
                let on_j = reference * (-multiplier / system.masses[e.j]);
                match mode {
                    Projection::Position { .. } => {
                        system.add_pos(e.i, on_i);
                        system.add_pos(e.j, on_j);
                    }
                    Projection::Velocity => {
                        system.add_vel(e.i, on_i);
                        system.add_vel(e.j, on_j);
                    }
                }
                if let Some(ledger) = report.as_deref_mut() {
                    // Loose solver: atom j receives MINUS multiplier*r.
                    ledger.record(
                        index,
                        reference,
                        reference * (-multiplier / mode.interval()),
                    )?;
                }
                changed = true;
                guard.check(system, &[e.i, e.j])?;
            }
            if changed {
                if let Some(ledger) = report.as_deref_mut() {
                    ledger.iterations += 1;
                }
            }
            let (position, velocity) = geometry.residuals(system)?;
            let residual = match mode {
                Projection::Position { .. } => position,
                Projection::Velocity => velocity,
            };
            if residual <= tol {
                if let Some(ledger) = report {
                    ledger.max_position_residual = position;
                    ledger.max_velocity_residual = velocity;
                }
                return Ok(());
            }
        }
        Err(constraint_error("projection iteration limit"))
    }

    #[allow(clippy::too_many_arguments)]
    fn project_water(
        &self,
        system: &mut System,
        geometry: &ConstraintGeometry,
        first_edge: usize,
        water: RigidWater,
        mode: Projection,
        tol: f64,
        guard: &ProjectionGuard,
        report: &mut Option<&mut ProjectionReport>,
    ) -> Result<()> {
        let edges = &geometry.edges[first_edge..first_edge + 3];
        let masses = [
            system.masses[water.first],
            system.masses[water.first + 1],
            system.masses[water.first + 2],
        ];
        // Preserve the existing inner Newton cap; velocity uses one linear solve.
        for iteration in 0..=WATER_NEWTON_ITER {
            let current = [
                geometry.vector(system, &edges[0])?,
                geometry.vector(system, &edges[1])?,
                geometry.vector(system, &edges[2])?,
            ];
            let mut residual = [0.; 3];
            for k in 0..3 {
                residual[k] = match mode {
                    Projection::Position { .. } => {
                        current[k].norm2() - edges[k].distance * edges[k].distance
                    }
                    Projection::Velocity => {
                        current[k].dot(system.vel(edges[k].j) - system.vel(edges[k].i))
                    }
                };
            }
            if residual.iter().any(|x| !x.is_finite()) {
                return Err(constraint_error("nonfinite water residual"));
            }
            if residual.iter().all(|x| x.abs() < tol) {
                return Ok(());
            }
            if iteration == WATER_NEWTON_ITER {
                return Err(constraint_error("water Newton iteration limit"));
            }
            let reference = match mode {
                Projection::Position { .. } => {
                    [edges[0].reference, edges[1].reference, edges[2].reference]
                }
                Projection::Velocity => current,
            };
            let factor = match mode {
                Projection::Position { .. } => 2.,
                Projection::Velocity => 1.,
            };
            let matrix = water_matrix(current, reference, masses, factor);
            let Some(lam) = solve3(&matrix, &residual.map(|x| -x)) else {
                return Err(constraint_error("singular water projection solve"));
            };
            // Preserve simultaneous mass-weighted updates and arithmetic order.
            let delta = water_deltas(reference, lam, masses);
            for (k, d) in delta.into_iter().enumerate() {
                match mode {
                    Projection::Position { .. } => system.add_pos(water.first + k, d),
                    Projection::Velocity => system.add_vel(water.first + k, d),
                }
            }
            if let Some(ledger) = report.as_deref_mut() {
                for k in 0..3 {
                    // Water Newton: atom j receives PLUS lam*r.
                    ledger.record(
                        first_edge + k,
                        reference[k],
                        reference[k] * (lam[k] / mode.interval()),
                    )?;
                }
                ledger.iterations += 1;
            }
            guard.check(system, &[water.first, water.first + 1, water.first + 2])?;
            if matches!(mode, Projection::Velocity) {
                let residual = edges
                    .iter()
                    .zip(current)
                    .map(|(e, r)| r.dot(system.vel(e.j) - system.vel(e.i)).abs())
                    .fold(0., f64::max);
                return check_residual(residual, tol, "water velocity projection");
            }
        }
        Err(constraint_error("water Newton iteration limit"))
    }

    fn water_pairs(&self) -> impl Iterator<Item = Constraint> + '_ {
        self.waters.iter().flat_map(|w| {
            [
                Constraint {
                    i: w.first,
                    j: w.first + 1,
                    d: w.d_oh,
                },
                Constraint {
                    i: w.first,
                    j: w.first + 2,
                    d: w.d_oh,
                },
                Constraint {
                    i: w.first + 1,
                    j: w.first + 2,
                    d: w.d_hh,
                },
            ]
        })
    }
    fn position_residual(system: &System, c: Constraint) -> f64 {
        let d = system.box_size.min_image(system.pos(c.j) - system.pos(c.i));
        let residual = (d.norm2() - c.d * c.d).abs();
        if residual.is_finite() {
            residual
        } else {
            f64::INFINITY
        }
    }
    fn velocity_residual(system: &System, c: Constraint) -> f64 {
        let d = system.box_size.min_image(system.pos(c.j) - system.pos(c.i));
        let residual = d.dot(system.vel(c.j) - system.vel(c.i)).abs();
        if residual.is_finite() {
            residual
        } else {
            f64::INFINITY
        }
    }
    fn water_residual(&self, system: &System) -> f64 {
        self.water_pairs()
            .map(|c| Self::position_residual(system, c))
            .fold(0.0, f64::max)
    }
    fn water_velocity_residual(&self, system: &System) -> f64 {
        self.water_pairs()
            .map(|c| Self::velocity_residual(system, c))
            .fold(0.0, f64::max)
    }
    /// Worst distance-squared residual, including all water distances.
    pub fn max_residual(&self, system: &System) -> f64 {
        self.bonds
            .iter()
            .copied()
            .chain(self.water_pairs())
            .map(|c| Self::position_residual(system, c))
            .fold(0.0, f64::max)
    }
    /// Worst tangent residual, including water distances.
    pub fn max_velocity_residual(&self, system: &System) -> f64 {
        self.bonds
            .iter()
            .copied()
            .chain(self.water_pairs())
            .map(|c| Self::velocity_residual(system, c))
            .fold(0.0, f64::max)
    }
}

#[derive(Clone, Copy)]
enum Projection {
    Position { interval: f64 },
    Velocity,
}
impl Projection {
    fn validate_stage(self, stage: ProjectionStage) -> Result<()> {
        use ProjectionStage::*;
        let position_stage = match stage {
            VvPosition | BaoabA1Pre | BaoabA2Pre => true,
            VvFinal | BaoabB0 | BaoabA1Post | BaoabO | BaoabA2Post | BaoabB1 => false,
        };
        if position_stage != matches!(self, Self::Position { .. }) {
            return Err(constraint_error(&format!(
                "stage {stage:?} does not match projection operation"
            )));
        }
        Ok(())
    }
    fn buffer(self, s: &System) -> &[f64] {
        match self {
            Self::Position { .. } => &s.positions,
            Self::Velocity => &s.velocities,
        }
    }
    fn interval(self) -> f64 {
        match self {
            Self::Position { interval } => interval,
            Self::Velocity => 1.,
        }
    }
}
struct ProjectionGuard<'a> {
    geometry: &'a ConstraintGeometry,
    mode: Projection,
    half: f64,
    adjacent: Vec<Vec<usize>>,
}
impl<'a> ProjectionGuard<'a> {
    fn new(geometry: &'a ConstraintGeometry, mode: Projection, half: f64, n: usize) -> Self {
        let mut adjacent = Vec::new();
        if matches!(mode, Projection::Position { .. }) {
            adjacent.resize_with(n, Vec::new);
            for (k, e) in geometry.edges.iter().enumerate() {
                adjacent[e.i].push(k);
                adjacent[e.j].push(k);
            }
        }
        Self {
            geometry,
            mode,
            half,
            adjacent,
        }
    }
    fn check(&self, s: &System, moved: &[usize]) -> Result<()> {
        for &i in moved {
            let value = match self.mode {
                Projection::Position { .. } => s.pos(i),
                Projection::Velocity => s.vel(i),
            };
            if value.0.iter().any(|x| !x.is_finite()) {
                return Err(constraint_error("nonfinite projection update"));
            }
            if let Some(incident) = self.adjacent.get(i) {
                for &k in incident {
                    self.geometry.check_edge(s, self.half, k)?;
                }
            }
        }
        Ok(())
    }
}

fn constraint_error(message: &str) -> ModenaError {
    ModenaError::Physics(format!("constraint projection: {message}"))
}
fn validate_tolerance(tol: f64) -> Result<()> {
    if tol.is_finite() && tol > 0.0 {
        Ok(())
    } else {
        Err(constraint_error("tolerance must be finite and positive"))
    }
}
fn check_residual(residual: f64, tol: f64, method: &str) -> Result<()> {
    if residual.is_finite() && residual <= tol {
        Ok(())
    } else {
        Err(constraint_error(&format!(
            "{method}: residual {residual:e} exceeds tolerance {tol:e}"
        )))
    }
}

/// Standalone geometry repair uses current normals at each Newton iteration.
/// It has no elapsed physical interval and produces no impulse ledger.
fn settle_one(
    box_size: &SimBox,
    chunk: &mut [f64],
    masses: &[f64],
    w: RigidWater,
    tol: f64,
) -> Result<()> {
    let (mo, m1, m2) = (masses[w.first], masses[w.first + 1], masses[w.first + 2]);
    let doh2 = w.d_oh * w.d_oh;
    let dhh2 = w.d_hh * w.d_hh;
    let p = |chunk: &[f64], k: usize| Vec3::new(chunk[3 * k], chunk[3 * k + 1], chunk[3 * k + 2]);
    for iteration in 0..=WATER_NEWTON_ITER {
        let (o, h1, h2) = (p(chunk, 0), p(chunk, 1), p(chunk, 2));
        let v1 = box_size.min_image(h1 - o); // O–H1
        let v2 = box_size.min_image(h2 - o); // O–H2
        let v3 = box_size.min_image(h2 - h1); // H1–H2
        let g = [v1.norm2() - doh2, v2.norm2() - doh2, v3.norm2() - dhh2];
        if g[0].abs() < tol && g[1].abs() < tol && g[2].abs() < tol {
            return Ok(());
        }
        if iteration == WATER_NEWTON_ITER {
            return Err(constraint_error("water Newton iteration limit"));
        }
        let r = [v1, v2, v3];
        let j = water_matrix(r, r, [mo, m1, m2], 2.);
        let Some(lam) = solve3(&j, &[-g[0], -g[1], -g[2]]) else {
            return Err(constraint_error("singular water position solve"));
        };
        for (atom, delta) in water_deltas(r, lam, [mo, m1, m2]).into_iter().enumerate() {
            add3(chunk, atom, delta);
        }
    }
    Err(constraint_error("water Newton iteration limit"))
}

/// Project one water's velocities so `r_k · v_jk = 0` for all three
/// constraints — the same 3×3 system as the position Jacobian (with the
/// position vectors in place of the constraint vectors), linear in the
/// multipliers: one solve.
fn settle_vel_one(
    box_size: &SimBox,
    positions: &[f64],
    masses: &[f64],
    velocities: &mut [f64],
    w: RigidWater,
    tol: f64,
) -> Result<()> {
    let (mo, m1, m2) = (masses[w.first], masses[w.first + 1], masses[w.first + 2]);
    let f = &positions[w.first * 3..w.first * 3 + 9];
    let pos = |k: usize| Vec3::new(f[3 * k], f[3 * k + 1], f[3 * k + 2]);
    let vel = |k: usize| {
        Vec3::new(
            velocities[3 * k],
            velocities[3 * k + 1],
            velocities[3 * k + 2],
        )
    };
    let r1 = box_size.min_image(pos(1) - pos(0));
    let r2 = box_size.min_image(pos(2) - pos(0));
    let r3 = box_size.min_image(pos(2) - pos(1));
    let c = [
        r1.dot(vel(1) - vel(0)),
        r2.dot(vel(2) - vel(0)),
        r3.dot(vel(2) - vel(1)),
    ];
    if c[0].abs() < tol && c[1].abs() < tol && c[2].abs() < tol {
        if r1.norm2() < 1e-20 || r2.norm2() < 1e-20 || r3.norm2() < 1e-20 {
            return Err(constraint_error("degenerate water velocity geometry"));
        }
        return Ok(());
    }
    let r = [r1, r2, r3];
    let m = water_matrix(r, r, [mo, m1, m2], 1.);
    let Some(mu) = solve3(&m, &[-c[0], -c[1], -c[2]]) else {
        return Err(constraint_error("singular water velocity solve"));
    };
    for (atom, delta) in water_deltas(r, mu, [mo, m1, m2]).into_iter().enumerate() {
        add3(velocities, atom, delta);
    }
    Ok(())
}

/// Water Newton Jacobian: rows are current edges, columns the actual impulse
/// directions. Velocity projection uses current edges for both and factor 1.
fn water_matrix(
    current: [Vec3; 3],
    reference: [Vec3; 3],
    masses: [f64; 3],
    factor: f64,
) -> [[f64; 3]; 3] {
    let [v1, v2, v3] = current;
    let [r1, r2, r3] = reference;
    let [mo, m1, m2] = masses;
    [
        [
            factor * v1.dot(r1) * (1. / mo + 1. / m1),
            factor * v1.dot(r2) / mo,
            -factor * v1.dot(r3) / m1,
        ],
        [
            factor * v2.dot(r1) / mo,
            factor * v2.dot(r2) * (1. / mo + 1. / m2),
            factor * v2.dot(r3) / m2,
        ],
        [
            -factor * v3.dot(r1) / m1,
            factor * v3.dot(r2) / m2,
            factor * v3.dot(r3) * (1. / m1 + 1. / m2),
        ],
    ]
}
fn water_deltas(reference: [Vec3; 3], lam: [f64; 3], masses: [f64; 3]) -> [Vec3; 3] {
    let [r1, r2, r3] = reference;
    let [mo, m1, m2] = masses;
    [
        (r1 * lam[0] + r2 * lam[1]) * (-1. / mo),
        (r1 * lam[0] - r3 * lam[2]) * (1. / m1),
        (r2 * lam[1] + r3 * lam[2]) * (1. / m2),
    ]
}

/// Add `d` to atom `k` of a contiguous 3-atom chunk buffer.
#[inline]
fn add3(buf: &mut [f64], k: usize, d: Vec3) {
    buf[3 * k] += d.0[0];
    buf[3 * k + 1] += d.0[1];
    buf[3 * k + 2] += d.0[2];
}

/// General 3×3 solve via Cramer's rule. `None` on a singular/nonfinite system.
fn solve3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if !det.is_finite() || det.abs() < 1e-300 {
        return None;
    }
    let x0 = b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2])
        + a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2]);
    let x1 = a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2])
        - b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0]);
    let x2 = a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1])
        - a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0])
        + b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let solution = [x0 / det, x1 / det, x2 / det];
    solution.iter().all(|x| x.is_finite()).then_some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::SimBox;

    /// SHAKE/RATTLE on a single bond: two atoms, stretched predictor, fixed to
    /// the target length, and the velocity projected onto the tangent.
    #[test]
    fn shake_single_bond_satisfies_constraint() {
        let mut sys = System::new(SimBox::cubic(100.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.5, 0.0, 0.0)); // stretched (target 1.0)
                                                  // Reference (start-of-step) bond at equilibrium along x.
        let x_ref = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let cons = Constraints {
            bonds: vec![Constraint { i: 0, j: 1, d: 1.0 }],
            waters: vec![],
        };
        cons.shake(&mut sys, &x_ref, 1e-12, 50).unwrap();
        let r = sys.box_size.min_image(sys.pos(1) - sys.pos(0));
        assert!(
            (r.norm() - 1.0).abs() < 1e-9,
            "bond length after SHAKE: {}",
            r.norm()
        );
        // COM preserved (equal masses): the predicted midpoint (0.75) is unchanged.
        let com = (sys.pos(0) + sys.pos(1)) * 0.5;
        assert!(
            (com - Vec3::new(0.75, 0.0, 0.0)).norm() < 1e-9,
            "COM drifted: {com:?}"
        );
    }

    #[test]
    fn rattle_removes_radial_velocity() {
        let mut sys = System::new(SimBox::cubic(100.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.0, 0.0, 0.0));
        // Pure radial velocity (stretches the bond) — must be removed.
        sys.set_vel(0, Vec3::new(-1.0, 0.0, 0.0));
        sys.set_vel(1, Vec3::new(2.0, 0.5, 0.0));
        let cons = Constraints {
            bonds: vec![Constraint { i: 0, j: 1, d: 1.0 }],
            waters: vec![],
        };
        cons.rattle(&mut sys, 1e-12, 50).unwrap();
        let rij = sys.box_size.min_image(sys.pos(1) - sys.pos(0));
        let vij = sys.vel(1) - sys.vel(0);
        assert!(
            rij.dot(vij).abs() < 1e-9,
            "radial velocity after RATTLE: {}",
            rij.dot(vij)
        );
    }

    /// SPC/E-geometry water, badly distorted, restored to machine precision
    /// by the exact per-molecule Newton projection — with COM preserved
    /// (the impulses are internal and mass-weighted).
    #[test]
    fn settle_restores_distorted_water() {
        let d_oh = 1.0_f64;
        let d_hh = 2.0 * d_oh * (109.47_f64.to_radians() / 2.0).sin();
        let mut sys = System::new(SimBox::cubic(20.0));
        sys.resize(3);
        let half = 109.47_f64.to_radians() / 2.0;
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(d_oh * half.cos(), d_oh * half.sin(), 0.0));
        sys.set_pos(2, Vec3::new(d_oh * half.cos(), -d_oh * half.sin(), 0.0));
        sys.masses = vec![15.9994, 1.008, 1.008];
        // Distort far beyond any physical dt (15% bond stretch).
        sys.set_pos(1, sys.pos(1) + Vec3::new(0.15, -0.1, 0.07));
        sys.set_pos(2, sys.pos(2) + Vec3::new(-0.12, 0.09, 0.11));
        let com_before = mass_weighted_com(&sys);
        let cons = Constraints {
            bonds: vec![],
            waters: vec![RigidWater {
                first: 0,
                d_oh,
                d_hh,
            }],
        };
        cons.settle(&mut sys, 1e-12).unwrap();
        let oh1 = (sys.pos(1) - sys.pos(0)).norm();
        let oh2 = (sys.pos(2) - sys.pos(0)).norm();
        let hh = (sys.pos(2) - sys.pos(1)).norm();
        assert!((oh1 - d_oh).abs() < 1e-10, "OH1 {oh1}");
        assert!((oh2 - d_oh).abs() < 1e-10, "OH2 {oh2}");
        assert!((hh - d_hh).abs() < 1e-10, "HH {hh}");
        let com_after = mass_weighted_com(&sys);
        assert!(
            (com_after - com_before).norm() < 1e-12,
            "COM moved: {com_before:?} -> {com_after:?}"
        );
    }

    /// The velocity projection removes the radial component along all three
    /// constraint directions at once (one exact 3×3 solve).
    #[test]
    fn settle_projects_water_velocities() {
        let d_oh = 1.0_f64;
        let d_hh = 2.0 * d_oh * (109.47_f64.to_radians() / 2.0).sin();
        let mut sys = System::new(SimBox::cubic(20.0));
        sys.resize(3);
        let half = 109.47_f64.to_radians() / 2.0;
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(d_oh * half.cos(), d_oh * half.sin(), 0.0));
        sys.set_pos(2, Vec3::new(d_oh * half.cos(), -d_oh * half.sin(), 0.0));
        sys.masses = vec![15.9994, 1.008, 1.008];
        // Radial (stretching) + tangential velocities on every bond.
        sys.set_vel(0, Vec3::new(0.3, -0.2, 0.1));
        sys.set_vel(1, Vec3::new(-0.5, 0.4, 0.2));
        sys.set_vel(2, Vec3::new(0.2, 0.1, -0.3));
        let cons = Constraints {
            bonds: vec![],
            waters: vec![RigidWater {
                first: 0,
                d_oh,
                d_hh,
            }],
        };
        cons.settle_velocities(&mut sys, 1e-12).unwrap();
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            let r = sys.pos(j) - sys.pos(i);
            let v = sys.vel(j) - sys.vel(i);
            assert!(
                r.dot(v).abs() < 1e-10,
                "radial velocity on bond ({i},{j}): {}",
                r.dot(v)
            );
        }
    }

    fn mass_weighted_com(sys: &System) -> Vec3 {
        let mut m = 0.0;
        let mut c = Vec3::zero();
        for i in 0..sys.n() {
            c += sys.pos(i) * sys.masses[i];
            m += sys.masses[i];
        }
        c * (1.0 / m)
    }
}
