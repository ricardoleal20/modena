//! Force models — the central extensibility seam.
//!
//! A `ForceModel` computes all forces on all atoms (and the potential energy
//! and scalar virial) for a given `System`. This is the trait an ML surrogate
//! potential (e.g. a NequIP/MACE-style learned model, or a coarse-grained
//! tabulated force) will implement to drop into the engine without touching the
//! integrator, thermostat, or I/O. M7 ships one concrete model:
//! `PairForceModel`, which evaluates a `PairPotential` over a half neighbor list.
//!
//! The trait is `&self` so stateless models are trivial; models needing scratch
//! space own it internally (`&mut self` is not required because scratch can be
//! thread-local or rebuilt per call).
//!
//! M7 contract: a single `PairForceModel` **overwrites** `system.forces` (it
//! owns all forces). M3 introduces `CompositeForceModel`, at which point each
//! model switches to **accumulating** (`add_force`) and the integrator zeroes
//! forces + virial once before the composite pass.

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
    /// `system.positions` and **overwrites** `system.forces` (every entry) when
    /// it is the sole model; under the M3 composite contract it **accumulates**
    /// into `system.forces` instead. Returns the total potential energy.
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64>;
}

/// Pairwise-additive force model: sums a `PairPotential` over all neighbor
/// pairs. Parallelized over the half neighbor list with rayon.
///
/// Uses a *half* list (each pair `(i, j)` stored once) and applies Newton's
/// third law — `f_i += fv`, `f_j -= fv` — so each pair is evaluated exactly
/// once (half the pair work of the Phase-1 full list). Each rayon thread folds
/// into its own `3N` force accumulator, then the per-thread accumulators are
/// summed; this keeps the pair loop race-free without atomics.
pub struct PairForceModel {
    pub potential: Box<dyn PairPotential>,
}

impl PairForceModel {
    pub fn new(potential: impl PairPotential + 'static) -> Self {
        Self {
            potential: Box::new(potential),
        }
    }

    /// Construct from an already-boxed potential (used by the config dispatch,
    /// which builds the potential dynamically from `ModelSpec`).
    pub fn from_boxed(potential: Box<dyn PairPotential>) -> Self {
        Self { potential }
    }
}

impl ForceModel for PairForceModel {
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        // Skin-aware rebuild: cheap no-op when no atom has crossed the skin.
        neighbors.build(system);

        let cutoff2 = self.potential.cutoff2();
        let n = system.n();
        let box_size = system.box_size;
        // Read-only shared views for the parallel kernel.
        let positions = &system.positions;
        let species = &system.species;
        let pairs = neighbors.pairs();
        let potential = &*self.potential;

        if n == 0 {
            virial.w = 0.0;
            return Ok(0.0);
        }

        // Per-thread `(forces, pe, virial)` accumulators, combined by reduce.
        let (forces, pe, w) = pairs
            .par_iter()
            .fold(
                || (vec![0.0f64; 3 * n], 0.0f64, 0.0f64),
                |mut acc, &(i, j)| {
                    let ri =
                        Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
                    let rj =
                        Vec3::new(positions[3 * j], positions[3 * j + 1], positions[3 * j + 2]);
                    let d = box_size.min_image(ri - rj);
                    let r2 = d.norm2();
                    if r2 >= cutoff2 || r2 < 1e-14 {
                        return acc;
                    }
                    let ti = species[i];
                    let tj = species[j];
                    let (u, f_over_r) = potential.energy_force(r2, ti, tj);
                    let fv = d * f_over_r;
                    // Newton's third law: equal & opposite.
                    for k in 0..3 {
                        acc.0[3 * i + k] += fv.0[k];
                        acc.0[3 * j + k] -= fv.0[k];
                    }
                    acc.1 += u;
                    acc.2 += d.dot(fv);
                    acc
                },
            )
            .reduce(
                || (vec![0.0f64; 3 * n], 0.0f64, 0.0f64),
                |mut a, b| {
                    for k in 0..3 * n {
                        a.0[k] += b.0[k];
                    }
                    a.1 += b.1;
                    a.2 += b.2;
                    a
                },
            );

        // Single model → overwrite forces.
        system.forces.copy_from_slice(&forces);
        virial.w = w;
        Ok(pe)
    }
}
