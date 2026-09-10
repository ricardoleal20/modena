//! The particle system: positions, velocities, forces, masses, species, box.

use crate::vec3::Vec3;

/// Orthorhombic simulation box (Phase 1). Triclinic is on the roadmap (M4).
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

    /// Wrap a position into the primary periodic image.
    pub fn wrap(&self, r: Vec3) -> Vec3 {
        let mut out = r.0;
        for (o_k, b_k) in out.iter_mut().zip(self.size.0.iter()) {
            *o_k -= b_k * (*o_k / b_k).floor();
        }
        Vec3(out)
    }
}

/// The many-particle system.
///
/// Arrays are length `n`. `species` indexes per-atom type parameters (e.g. LJ
/// `sigma`/`epsilon` cross-interaction tables). For Phase 1 single-component
/// systems all species are `0`.
#[derive(Clone)]
pub struct System {
    pub positions: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
    pub forces: Vec<Vec3>,
    pub masses: Vec<f64>,
    pub species: Vec<u32>,
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
            elements: Vec::new(),
            box_size,
        }
    }

    pub fn n(&self) -> usize {
        self.positions.len()
    }

    /// Reserve `n` atoms and pre-zero velocities/forces.
    pub fn resize(&mut self, n: usize) {
        self.positions.resize(n, Vec3::zero());
        self.velocities.resize(n, Vec3::zero());
        self.forces.resize(n, Vec3::zero());
        self.masses.resize(n, 1.0);
        self.species.resize(n, 0);
        self.elements.resize(n, "X".to_string());
    }

    /// Kinetic energy ½ Σ m v².
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.n() {
            let m = self.masses[i];
            ke += 0.5 * m * self.velocities[i].norm2();
        }
        ke
    }

    /// Center-of-mass linear momentum (should be ~0 for a properly initialized
    /// NVE run; conserved up to floating-point drift).
    pub fn total_momentum(&self) -> Vec3 {
        let mut p = Vec3::zero();
        for i in 0..self.n() {
            p += self.velocities[i] * self.masses[i];
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
            *f = Vec3::zero();
        }
    }

    /// Wrap all positions into the primary image.
    pub fn wrap_positions(&mut self) {
        let b = self.box_size;
        for p in self.positions.iter_mut() {
            *p = b.wrap(*p);
        }
    }
}
