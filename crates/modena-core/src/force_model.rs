//! Force models — the central extensibility seam.
//!
//! A `ForceModel` computes all forces on all atoms (and the potential energy
//! and scalar virial) for a given `System`. This is the trait an ML surrogate
//! potential (e.g. a NequIP/MACE-style learned model, or a coarse-grained
//! tabulated force) will implement to drop into the engine without touching the
//! integrator, thermostat, or I/O. Phase 1 ships one concrete model:
//! `PairForceModel`, which evaluates a `PairPotential` over a neighbor list.
//!
//! The trait is `&self` so stateless models are trivial; models needing scratch
//! space own it internally (`&mut self` is not required because scratch can be
//! thread-local or rebuilt per call).

use rayon::prelude::*;

use crate::error::Result;
use crate::neighbor::NeighborList;
use crate::potential::PairPotential;
use crate::system::System;
use crate::vec3::Vec3;

/// Symmetric pressure virial (scalar for Phase 1). `W = Σ_ij r_ij · F_ij`.
/// Pressure (reduced) = `(2 KE + W) / (3 V)`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Virial {
    pub w: f64,
}

/// The force/energy provider. Implement this to plug in a new potential —
/// analytical pair, many-body, bonded, or a learned ML surrogate.
pub trait ForceModel: Send + Sync {
    /// Compute forces, potential energy, and virial. The implementation reads
    /// `system.positions` and **overwrites** `system.forces` (every entry).
    /// Returns the total potential energy.
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64>;
}

/// Pairwise-additive force model: sums a `PairPotential` over all neighbor
/// pairs. Parallelized over atoms with rayon.
///
/// Uses a *full* neighbor list (each pair appears in both atoms' lists), so the
/// parallel per-atom pass writes only `forces[i]` — no races, no reaction-force
/// scatter. Each pair is evaluated twice, so energy and virial are halved. This
/// trades ~2× work for simplicity and correctness; a half-list + skin refactor
/// is on the performance roadmap (M7).
pub struct PairForceModel {
    pub potential: Box<dyn PairPotential>,
}

impl PairForceModel {
    pub fn new(potential: impl PairPotential + 'static) -> Self {
        Self {
            potential: Box::new(potential),
        }
    }
}

impl ForceModel for PairForceModel {
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        neighbors.build(system);

        let cutoff2 = self.potential.cutoff2();
        let n = system.n();
        let box_size = system.box_size.size;
        // Shared, read-only views for the parallel kernel.
        let positions = &system.positions;
        let species = &system.species;
        let nb = &*neighbors;
        let potential = &*self.potential;

        let results: Vec<(Vec3, f64, f64)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut f_i = Vec3::zero();
                let mut pe_i = 0.0;
                let mut w_i = 0.0;
                let ri = positions[i];
                let ti = species[i];
                for &j in nb.neighbors_of(i) {
                    let rj = positions[j];
                    let d = (ri - rj).min_image(box_size);
                    let r2 = d.norm2();
                    if r2 >= cutoff2 || r2 < 1e-14 {
                        continue;
                    }
                    let tj = species[j];
                    let (u, f_over_r) = potential.energy_force(r2, ti, tj);
                    let fv = d * f_over_r;
                    f_i += fv;
                    pe_i += u;
                    w_i += d.dot(fv);
                }
                (f_i, pe_i, w_i)
            })
            .collect();

        let mut pe = 0.0;
        let mut w = 0.0;
        let forces = &mut system.forces;
        for (i, (f_i, pe_i, w_i)) in results.into_iter().enumerate() {
            forces[i] = f_i;
            pe += pe_i;
            w += w_i;
        }
        // Each pair was counted from both ends.
        pe *= 0.5;
        w *= 0.5;

        virial.w = w;
        Ok(pe)
    }
}
