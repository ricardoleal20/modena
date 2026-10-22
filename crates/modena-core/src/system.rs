//! The particle system: positions, velocities, forces, masses, species, box.
//!
//! M7 layout: the per-component position/velocity/force arrays are stored as
//! flat `Vec<f64>` of length `3N` (structure-of-arrays) for cache-line and
//! SIMD friendliness in the hot force kernel, and so downstream milestones
//! (bonded topology, Ewald FFT grid, GPU upload) operate on contiguous buffers.
//! [`Vec3`] remains the per-atom math helper; accessors [`pos`](Self::pos)/
//! [`set_pos`](Self::set_pos)/[`vel`](Self::vel)/... bridge the two views.
//! Masses, species, charges, and element labels stay per-atom (`N`).

use crate::error::{ModenaError, Result};
use crate::vec3::Vec3;

/// Simulation box.
///
/// Orthorhombic by default (stored as the per-axis edge `size`). M4 adds
/// **triclinic** support: a general cell defined by three Bravais vectors
/// `a`, `b`, `c` (the cell matrix `h = [a b c]`, columns). All minimum-image
/// dispatch routes through [`min_image`](Self::min_image); triclinic images
/// use a bounded nearest-lattice search in Cartesian distance.
///
/// `size` is always the **orthorhombic bounding box** of the cell — for an
/// orthorhombic cell it *is* the cell; for a triclinic cell it is the
/// enclosing box (edge_k = |a_k|+|b_k|+|c_k|). Cartesian periodic cell stencils
/// are insufficient for general oblique faces; neighbor providers must use a
/// geometry-safe search. `volume` returns `|det h|` for a triclinic cell.
///
/// The M5 barostats scale `size` and positions isotropically and are
/// orthorhombic-only; combining a triclinic box with a barostat is a
/// configuration error (the cell matrix would not stay consistent).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimBox {
    /// Orthorhombic bounding-box edges. For an orthorhombic cell, the cell
    /// edges themselves.
    pub size: Vec3,
    /// Triclinic cell vectors `[a, b, c]` (the columns of `h`); `None` for an
    /// orthorhombic cell.
    vectors: Option<[Vec3; 3]>,
    /// Cached inverse `h⁻¹` rows: fractional `s_k = inv_rows[k] · r`. `None`
    /// for an orthorhombic cell.
    inv_rows: Option<[Vec3; 3]>,
}

impl SimBox {
    pub fn cubic(edge: f64) -> Self {
        Self {
            size: Vec3::splat(edge),
            vectors: None,
            inv_rows: None,
        }
    }

    pub fn new(size: Vec3) -> Self {
        Self {
            size,
            vectors: None,
            inv_rows: None,
        }
    }

    pub fn from_array(a: [f64; 3]) -> Self {
        Self {
            size: Vec3::new(a[0], a[1], a[2]),
            vectors: None,
            inv_rows: None,
        }
    }

    /// Triclinic cell from the three Bravais vectors `[a, b, c]`. The
    /// orthorhombic bounding box is computed and cached in `size`; the inverse
    /// `h⁻¹` (as rows) is cached for fractional-coordinate wrapping.
    pub fn triclinic(vectors: [Vec3; 3]) -> Self {
        let [a, b, c] = vectors;
        // Bounding-box edges: the cell vertices range over ±a ±b ±c, so the
        // extent along axis k is |a_k| + |b_k| + |c_k|.
        let size = Vec3::new(
            a.0[0].abs() + b.0[0].abs() + c.0[0].abs(),
            a.0[1].abs() + b.0[1].abs() + c.0[1].abs(),
            a.0[2].abs() + b.0[2].abs() + c.0[2].abs(),
        );
        let inv_rows = inverse_rows(a, b, c);
        Self {
            size,
            vectors: Some(vectors),
            inv_rows: Some(inv_rows),
        }
    }

    /// Whether this is a (non-orthorhombic) triclinic cell.
    pub fn is_triclinic(&self) -> bool {
        self.vectors.is_some()
    }

    /// The cell vectors `[a, b, c]` if triclinic.
    pub fn vectors(&self) -> Option<&[Vec3; 3]> {
        self.vectors.as_ref()
    }

    /// The real-space Bravais basis `[a, b, c]` (columns of `h`). For an
    /// orthorhombic cell these are the axis-aligned box edges; for a triclinic
    /// cell the stored vectors. Used by Ewald to bound the k-vector sphere.
    pub fn basis(&self) -> [Vec3; 3] {
        match self.vectors {
            Some(v) => v,
            None => [
                Vec3::new(self.size.0[0], 0.0, 0.0),
                Vec3::new(0.0, self.size.0[1], 0.0),
                Vec3::new(0.0, 0.0, self.size.0[2]),
            ],
        }
    }

    /// The inverse-cell rows `h⁻¹` (fractional `s_k = row_k · r`). For an
    /// orthorhombic cell, `1/L_k` along each axis. These are also the
    /// reciprocal-lattice directions: a k-vector is `2π Σ n_k · row_k`.
    pub fn inv_basis_rows(&self) -> [Vec3; 3] {
        match self.inv_rows {
            Some(v) => v,
            None => [
                Vec3::new(1.0 / self.size.0[0], 0.0, 0.0),
                Vec3::new(0.0, 1.0 / self.size.0[1], 0.0),
                Vec3::new(0.0, 0.0, 1.0 / self.size.0[2]),
            ],
        }
    }

    pub fn volume(&self) -> f64 {
        match self.vectors {
            Some([a, b, c]) => (a.cross(b).dot(c)).abs(),
            None => self.size.0[0] * self.size.0[1] * self.size.0[2],
        }
    }

    /// Minimum-image of a displacement vector. The single dispatch point for
    /// periodic wrapping in the engine. Orthorhombic fast path: per-axis
    /// `dx -= L·round(dx/L)`. Triclinic: bounded lattice enumeration.
    #[inline]
    pub fn min_image(&self, d: Vec3) -> Vec3 {
        match self.inv_rows {
            None => {
                let mut r = d.0;
                let b = self.size.0;
                for k in 0..3 {
                    r[k] -= b[k] * (r[k] / b[k]).round();
                }
                Vec3(r)
            }
            Some(inv) => {
                let s = [
                    wrap_half(inv[0].dot(d)),
                    wrap_half(inv[1].dot(d)),
                    wrap_half(inv[2].dot(d)),
                ];
                let [a, b, c] = self.vectors.unwrap();
                let mut best = a * s[0] + b * s[1] + c * s[2];
                let mut best2 = best.norm2();
                // Any closer image q has |row_k·q| <= ||row_k|| R,
                // so |s_k-n_k| <= ||row_k|| R. This finite bounding box
                // contains every improvement, even for highly skew bases.
                let radius = best2.sqrt();
                let bounds = std::array::from_fn::<_, 3, _>(|k| {
                    let width = inv[k].norm() * radius;
                    let margin = 32.0 * f64::EPSILON * (1.0 + width);
                    (
                        (s[k] - width - margin).ceil() as i64,
                        (s[k] + width + margin).floor() as i64,
                    )
                });
                for i in bounds[0].0..=bounds[0].1 {
                    for j in bounds[1].0..=bounds[1].1 {
                        for k in bounds[2].0..=bounds[2].1 {
                            let candidate = a * (s[0] - i as f64)
                                + b * (s[1] - j as f64)
                                + c * (s[2] - k as f64);
                            let r2 = candidate.norm2();
                            if r2 < best2 {
                                best = candidate;
                                best2 = r2;
                            }
                        }
                    }
                }
                best
            }
        }
    }

    /// Wrap a position into the primary periodic image.
    pub fn wrap(&self, r: Vec3) -> Vec3 {
        match self.inv_rows {
            None => {
                let mut out = r.0;
                let b = self.size.0;
                for k in 0..3 {
                    out[k] -= b[k] * (out[k] / b[k]).floor();
                }
                Vec3(out)
            }
            Some(inv) => {
                // Fractional position, wrapped to [0, 1).
                let s0 = wrap_unit(inv[0].dot(r));
                let s1 = wrap_unit(inv[1].dot(r));
                let s2 = wrap_unit(inv[2].dot(r));
                let [a, b, c] = self.vectors.unwrap();
                a * s0 + b * s1 + c * s2
            }
        }
    }

    /// Scale one box axis by `factor` (used by M5 barostats; orthorhombic).
    pub fn scale_axis(&mut self, axis: usize, factor: f64) {
        if let Some(mut vectors) = self.vectors {
            for v in &mut vectors {
                v.0[axis] *= factor;
            }
            *self = Self::triclinic(vectors);
        } else {
            self.size.0[axis] *= factor;
        }
    }

    /// Isotropic scale of all three axes (used by M5 barostats).
    pub fn scale_uniform(&mut self, factor: f64) {
        if let Some(vectors) = self.vectors {
            *self = Self::triclinic(vectors.map(|v| v * factor));
        } else {
            self.size *= factor;
        }
    }
}

/// Wrap `x` to the range `[-0.5, 0.5)` (minimum-image convention).
#[inline]
fn wrap_half(x: f64) -> f64 {
    x - x.round()
}

/// Wrap `x` to the range `[0, 1)` (primary-cell convention).
#[inline]
fn wrap_unit(x: f64) -> f64 {
    x - x.floor()
}

/// Inverse of the cell matrix `h = [a b c]` (columns), returned as its **rows**:
/// `inv_rows[0]` is the first row of `h⁻¹`, so `s_k = inv_rows[k] · r`. Computed
/// via the adjugate: `h⁻¹ = adj(h) / det(h)`, with `det = a · (b × c)`.
fn inverse_rows(a: Vec3, b: Vec3, c: Vec3) -> [Vec3; 3] {
    let det = a.cross(b).dot(c);
    // Cofactors of h (columns a,b,c). Row k of h⁻¹ = (cofactor of column k) / det.
    // For h with columns a,b,c, the inverse rows are:
    //   row0 = (b × c) / det
    //   row1 = (c × a) / det
    //   row2 = (a × b) / det
    let inv_det = 1.0 / det;
    [
        b.cross(c) * inv_det,
        c.cross(a) * inv_det,
        a.cross(b) * inv_det,
    ]
}

/// The many-particle system.
///
/// `positions`/`velocities`/`forces` are flat `3N` buffers (SoA). `masses`,
/// `species`, `charges`, `elements` are length `N`. `species` indexes per-atom
/// type parameters; `charges` (M4) defaults to 0.
#[derive(Clone)]
pub struct System {
    pub positions: Vec<f64>,
    pub velocities: Vec<f64>,
    pub forces: Vec<f64>,
    pub masses: Vec<f64>,
    pub species: Vec<u32>,
    pub charges: Vec<f64>,
    /// Element symbols for I/O (cosmetic; not used by the physics).
    pub elements: Vec<String>,
    pub box_size: SimBox,
    constraint_count: usize,
    conserved_com_dof: usize,
}

impl System {
    pub fn new(box_size: SimBox) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            species: Vec::new(),
            charges: Vec::new(),
            elements: Vec::new(),
            box_size,
            constraint_count: 0,
            conserved_com_dof: 0,
        }
    }

    pub fn n(&self) -> usize {
        self.masses.len()
    }

    /// Reserve `n` atoms and pre-zero velocities/forces.
    pub fn resize(&mut self, n: usize) {
        // Resizing invalidates physical manifold metadata.
        self.constraint_count = 0;
        self.conserved_com_dof = 0;
        self.positions.resize(3 * n, 0.0);
        self.velocities.resize(3 * n, 0.0);
        self.forces.resize(3 * n, 0.0);
        self.masses.resize(n, 1.0);
        self.species.resize(n, 0);
        self.charges.resize(n, 0.0);
        self.elements.resize(n, "X".to_string());
    }

    // --- SoA accessors (Vec3 bridge) ---

    #[inline]
    pub fn pos(&self, i: usize) -> Vec3 {
        Vec3::new(
            self.positions[3 * i],
            self.positions[3 * i + 1],
            self.positions[3 * i + 2],
        )
    }

    #[inline]
    pub fn set_pos(&mut self, i: usize, v: Vec3) {
        self.positions[3 * i] = v.0[0];
        self.positions[3 * i + 1] = v.0[1];
        self.positions[3 * i + 2] = v.0[2];
    }

    #[inline]
    pub fn add_pos(&mut self, i: usize, v: Vec3) {
        self.positions[3 * i] += v.0[0];
        self.positions[3 * i + 1] += v.0[1];
        self.positions[3 * i + 2] += v.0[2];
    }

    #[inline]
    pub fn vel(&self, i: usize) -> Vec3 {
        Vec3::new(
            self.velocities[3 * i],
            self.velocities[3 * i + 1],
            self.velocities[3 * i + 2],
        )
    }

    #[inline]
    pub fn set_vel(&mut self, i: usize, v: Vec3) {
        self.velocities[3 * i] = v.0[0];
        self.velocities[3 * i + 1] = v.0[1];
        self.velocities[3 * i + 2] = v.0[2];
    }

    #[inline]
    pub fn add_vel(&mut self, i: usize, v: Vec3) {
        self.velocities[3 * i] += v.0[0];
        self.velocities[3 * i + 1] += v.0[1];
        self.velocities[3 * i + 2] += v.0[2];
    }

    #[inline]
    pub fn force(&self, i: usize) -> Vec3 {
        Vec3::new(
            self.forces[3 * i],
            self.forces[3 * i + 1],
            self.forces[3 * i + 2],
        )
    }

    #[inline]
    pub fn set_force(&mut self, i: usize, v: Vec3) {
        self.forces[3 * i] = v.0[0];
        self.forces[3 * i + 1] = v.0[1];
        self.forces[3 * i + 2] = v.0[2];
    }

    /// Add `v` to the force on atom `i` (the accumulate primitive; used by the
    /// half-list kernel and M3/M4 composite models).
    #[inline]
    pub fn add_force(&mut self, i: usize, v: Vec3) {
        self.forces[3 * i] += v.0[0];
        self.forces[3 * i + 1] += v.0[1];
        self.forces[3 * i + 2] += v.0[2];
    }

    #[inline]
    pub fn mass(&self, i: usize) -> f64 {
        self.masses[i]
    }

    #[inline]
    pub fn species(&self, i: usize) -> u32 {
        self.species[i]
    }

    #[inline]
    pub fn charge(&self, i: usize) -> f64 {
        self.charges[i]
    }

    // --- raw SoA slices (for hot kernels / GPU upload) ---

    #[inline]
    pub fn positions_raw(&self) -> &[f64] {
        &self.positions
    }

    #[inline]
    pub fn velocities_raw(&self) -> &[f64] {
        &self.velocities
    }

    #[inline]
    pub fn forces_raw(&self) -> &[f64] {
        &self.forces
    }

    #[inline]
    pub fn forces_raw_mut(&mut self) -> &mut [f64] {
        &mut self.forces
    }

    /// Kinetic energy ½ Σ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.n() {
            let m = self.masses[i];
            let v = self.vel(i);
            ke += 0.5 * m * v.norm2();
        }
        ke
    }

    /// Center-of-mass linear momentum (should be ~0 for a properly initialized
    /// NVE run; conserved up to floating-point drift).
    pub fn total_momentum(&self) -> Vec3 {
        let mut p = Vec3::zero();
        for i in 0..self.n() {
            p += self.vel(i) * self.masses[i];
        }
        p
    }

    /// Independent kinetic degrees of freedom. Manual systems default to 3N.
    pub fn dof(&self) -> usize {
        (3 * self.n())
            .saturating_sub(self.constraint_count)
            .saturating_sub(self.conserved_com_dof)
    }
    /// Independent holonomic constraints (water contributes three).
    pub fn constraint_count(&self) -> usize {
        self.constraint_count
    }
    pub fn set_constraint_count(&mut self, count: usize) -> Result<()> {
        if count > 3 * self.n() || self.conserved_com_dof > 3 * self.n() - count {
            return Err(ModenaError::Physics(
                "constraint count exceeds available DOF".into(),
            ));
        }
        self.constraint_count = count;
        Ok(())
    }
    /// COM components fixed by initialized conserved momentum. Stochastic
    /// Langevin COM belongs to the ensemble and must use zero here.
    pub fn conserved_com_dof(&self) -> usize {
        self.conserved_com_dof
    }
    pub fn set_conserved_com_dof(&mut self, count: usize) -> Result<()> {
        if count > 3 || count > 3 * self.n() - self.constraint_count {
            return Err(ModenaError::Physics(
                "conserved COM DOF must be 0..=3 and fit the system".into(),
            ));
        }
        self.conserved_com_dof = count;
        Ok(())
    }

    /// Instantaneous temperature from KE: T = 2 KE / (dof * kB), kB = 1 in
    /// reduced units.
    pub fn temperature(&self) -> f64 {
        let dof = self.dof() as f64;
        if dof == 0.0 {
            0.0
        } else {
            2.0 * self.kinetic_energy() / dof
        }
    }

    /// Zero all forces (called before each force evaluation).
    pub fn zero_forces(&mut self) {
        for f in self.forces.iter_mut() {
            *f = 0.0;
        }
    }

    /// Wrap all positions into the primary image.
    pub fn wrap_positions(&mut self) {
        let b = self.box_size;
        for i in 0..self.n() {
            self.set_pos(i, b.wrap(self.pos(i)));
        }
    }
}
