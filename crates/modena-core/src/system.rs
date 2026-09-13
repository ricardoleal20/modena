//! The particle system: positions, velocities, forces, masses, species, box.
//!
//! M7 layout: the per-component position/velocity/force arrays are stored as
//! flat `Vec<f64>` of length `3N` (structure-of-arrays) for cache-line and
//! SIMD friendliness in the hot force kernel, and so downstream milestones
//! (bonded topology, Ewald FFT grid, GPU upload) operate on contiguous buffers.
//! [`Vec3`] remains the per-atom math helper; accessors [`pos`](Self::pos)/
//! [`set_pos`](Self::set_pos)/[`vel`](Self::vel)/... bridge the two views.
//! Masses, species, charges, and element labels stay per-atom (`N`).

use crate::vec3::Vec3;

/// Simulation box.
///
/// Orthorhombic in M1–M6 (stored as the per-axis edge `size`). M7 introduces
/// the [`min_image`](Self::min_image) seam: all minimum-image dispatch in the
/// engine routes through `SimBox::min_image`, so M4 can swap in a triclinic
/// cell (matrix + inverse + fractional wrap) behind the same method without
/// touching callers. The orthorhombic path is the fast path used until then.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SimBox {
    pub size: Vec3,
}

impl SimBox {
    pub fn cubic(edge: f64) -> Self {
        Self {
            size: Vec3::splat(edge),
        }
    }

    pub fn new(size: Vec3) -> Self {
        Self { size }
    }

    pub fn from_array(a: [f64; 3]) -> Self {
        Self {
            size: Vec3::new(a[0], a[1], a[2]),
        }
    }

    pub fn volume(&self) -> f64 {
        self.size.0[0] * self.size.0[1] * self.size.0[2]
    }

    /// Minimum-image of a displacement vector. The single dispatch point for
    /// periodic wrapping in the engine — M4's triclinic implementation replaces
    /// the body here; callers (force kernel, neighbor search) are unchanged.
    #[inline]
    pub fn min_image(&self, d: Vec3) -> Vec3 {
        // Orthorhombic fast path: dx -= L * round(dx / L) per axis.
        let mut r = d.0;
        let b = self.size.0;
        for k in 0..3 {
            r[k] -= b[k] * (r[k] / b[k]).round();
        }
        Vec3(r)
    }

    /// Wrap a position into the primary periodic image.
    pub fn wrap(&self, r: Vec3) -> Vec3 {
        let mut out = r.0;
        let b = self.size.0;
        for k in 0..3 {
            out[k] -= b[k] * (out[k] / b[k]).floor();
        }
        Vec3(out)
    }

    /// Scale one box axis by `factor` (used by M5 barostats; orthorhombic).
    pub fn scale_axis(&mut self, axis: usize, factor: f64) {
        self.size.0[axis] *= factor;
    }

    /// Isotropic scale of all three axes (used by M5 barostats).
    pub fn scale_uniform(&mut self, factor: f64) {
        for k in 0..3 {
            self.size.0[k] *= factor;
        }
    }
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
        }
    }

    pub fn n(&self) -> usize {
        self.masses.len()
    }

    /// Reserve `n` atoms and pre-zero velocities/forces.
    pub fn resize(&mut self, n: usize) {
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

    /// Number of positional degrees of freedom (3N). Thermostats subtract the
    /// COM constraints separately when needed.
    pub fn dof(&self) -> usize {
        3 * self.n()
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
