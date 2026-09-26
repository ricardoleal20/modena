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
use crate::potential::{LennardJones, PairPotential};
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
    /// `system.positions` and **accumulates** into `system.forces` and `virial`
    /// (it does *not* zero them). The caller — the integrator, or
    /// [`crate::composite::CompositeForceModel`] — zeroes forces + virial once
    /// before the force pass. Returns this model's potential-energy contribution.
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
        let pairs = neighbors.pairs();

        if n == 0 {
            virial.w = 0.0;
            return Ok(0.0);
        }

        // Monomorphize the hot path for the shipped LJ potential (removes the
        // per-pair vtable call, lets LLVM inline the pair math — bit-identical,
        // same `energy_force` arithmetic). Fall back to the dyn path otherwise.
        let serial = pairs.len() < SERIAL_PAIR_THRESHOLD || rayon::current_num_threads() <= 1;
        let (pe, w) = if let Some(lj) = self
            .potential
            .as_any()
            .and_then(|a| a.downcast_ref::<LennardJones>())
        {
            pair_pass(
                lj,
                pairs,
                n,
                &system.positions,
                &system.species,
                box_size,
                cutoff2,
                serial,
                &mut system.forces,
            )
        } else {
            pair_pass(
                &*self.potential,
                pairs,
                n,
                &system.positions,
                &system.species,
                box_size,
                cutoff2,
                serial,
                &mut system.forces,
            )
        };

        // M3 accumulate contract: `pair_pass` accumulated directly into
        // `system.forces` (the caller zeroed once before the pass) and returned
        // the PE + virial contributions; compose via virial.
        virial.w += w;
        Ok(pe)
    }
}

/// Above this half-pair count the parallel path is faster than a serial pass
/// on the default rayon pool; below it, rayon's per-split allocations and
/// work-stealing churn make the *default* pool slower than one thread
/// (measured: 11 threads were 31% slower than 1 thread at N=864, and the
/// default pool burned ~35 s of system time per 5000 steps at N=2048). The
/// serial path is also what bit-identity at 1 thread is measured against.
const SERIAL_PAIR_THRESHOLD: usize = 100_000;

/// Evaluate the half pair list: accumulate forces directly into `forces_out`
/// (serial; bit-identical to the former single-worker fold) or via bounded
/// per-chunk scratch (parallel, `with_min_len` limits allocations to ~thread
/// count per step). Returns `(PE, virial)`.
///
/// Generic over the potential so the `LennardJones` fast path monomorphizes
/// `energy_force` (inline + constant folding); the dyn path is unchanged.
#[allow(clippy::too_many_arguments)]
fn pair_pass<P: PairPotential + ?Sized>(
    potential: &P,
    pairs: &[(usize, usize)],
    n: usize,
    positions: &[f64],
    species: &[u32],
    box_size: crate::system::SimBox,
    cutoff2: f64,
    serial: bool,
    forces_out: &mut [f64],
) -> (f64, f64) {
    if serial {
        // Accumulate straight into the caller-zeroed forces array: the
        // addition sequence per pair is identical to the former single-worker
        // rayon fold (0.0 + … per component), so 1-thread trajectories stay
        // byte-identical — the G1 gate.
        let mut pe = 0.0f64;
        let mut w = 0.0f64;
        for &(i, j) in pairs {
            let ri = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
            let rj = Vec3::new(positions[3 * j], positions[3 * j + 1], positions[3 * j + 2]);
            let d = box_size.min_image(ri - rj);
            let r2 = d.norm2();
            if r2 >= cutoff2 || r2 < 1e-14 {
                continue;
            }
            let (u, f_over_r) = potential.energy_force(r2, species[i], species[j]);
            let fv = d * f_over_r;
            for k in 0..3 {
                forces_out[3 * i + k] += fv.0[k];
                forces_out[3 * j + k] -= fv.0[k];
            }
            pe += u;
            w += d.dot(fv);
        }
        (pe, w)
    } else {
        // Parallel with bounded splits: each thread folds into its own 3N
        // accumulator; `with_min_len` caps the split count ~ thread count, so
        // the per-split `vec![0.0; 3n]` allocations no longer dominate (they
        // were ~35 s of system time per 5000 steps at N=2048, default pool).
        let chunk = (pairs.len() / rayon::current_num_threads()).max(1);
        let (forces, pe, w) = pairs
            .par_iter()
            .with_min_len(chunk)
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
                    let (u, f_over_r) = potential.energy_force(r2, species[i], species[j]);
                    let fv = d * f_over_r;
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
        // M3 accumulate: add the per-thread-reduced forces into the caller's
        // array (which the caller zeroed once).
        for (dst, &src) in forces_out.iter_mut().zip(forces.iter()).take(3 * n) {
            *dst += src;
        }
        (pe, w)
    }
}
