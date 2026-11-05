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
//! Every provider follows the accumulation contract: the driver zeroes force
//! and virial once before a full pass, and pair/composite/backend models add
//! their contributions. Empty systems do not erase previous contributions.

use rayon::prelude::*;
use std::sync::Mutex;

use crate::error::Result;
use crate::neighbor::NeighborList;
use crate::potential::{LennardJones, PairPotential};
use crate::system::System;
use crate::vec3::Vec3;

/// Scalar isotropic strain virial. For real-space pair interactions,
/// `W = Σ_ij r_ij · F_ij`; reciprocal/background terms include explicit cell
/// dependence. Constraint impulse stress is not included.
/// Pressure (reduced) = `(2 KE + W) / (3 V)`.
#[derive(Debug, Clone, Copy, Default)]
pub struct Virial {
    pub w: f64,
}

/// The force/energy provider. Implement this to plug in a new potential —
/// analytical pair, many-body, bonded, or a learned ML surrogate.
pub trait ForceModel: Send + Sync {
    /// The provider constructs its own neighbors and ignores the compatibility list.
    fn uses_provider_neighbors(&self) -> bool {
        false
    }
    /// Whether this model preserves total momentum for these atoms under a
    /// common continuous translation. Drivers may remove COM modes only when
    /// both the model and the controller preserve an initialized zero momentum.
    /// Providers opt in explicitly; an unknown model receives the full space.
    fn conserves_total_momentum(&self, _system: &System) -> bool {
        false
    }
    /// Maximum real-space neighbor radius required by this model. Drivers
    /// must reject an absent radius when they need to construct a list,
    /// unless `uses_provider_neighbors()` explicitly opts into an unused
    /// compatibility list. Such a list implies no physical provider cutoff.
    fn neighbor_cutoff(&self) -> Option<f64> {
        None
    }
    /// Compute forces, potential energy, and virial. The implementation reads
    /// `system.positions` and **accumulates** into `system.forces` and `virial`
    /// (it does *not* zero them). The caller — the integrator or driver —
    /// zeroes forces + virial once
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
/// once. Parallel passes reuse at most one `3N` accumulator per useful worker.
/// Concurrent callers use serial arithmetic when the workspace is leased;
/// no Rayon worker blocks on another caller's fork/join work.
pub struct PairForceModel {
    pub potential: Box<dyn PairPotential>,
    scratch: Mutex<PairScratch>,
}

#[derive(Default)]
struct PairScratch {
    partitions: Vec<PairPartition>,
}

#[derive(Default)]
struct PairPartition {
    forces: Vec<f64>,
    pe: f64,
    w: f64,
}

impl PairForceModel {
    pub fn new(potential: impl PairPotential + 'static) -> Self {
        Self {
            potential: Box::new(potential),
            scratch: Mutex::default(),
        }
    }

    /// Construct from an already-boxed potential (used by the config dispatch,
    /// which builds the potential dynamically from `ModelSpec`).
    pub fn from_boxed(potential: Box<dyn PairPotential>) -> Self {
        Self {
            potential,
            scratch: Mutex::default(),
        }
    }
}

impl ForceModel for PairForceModel {
    fn conserves_total_momentum(&self, _system: &System) -> bool {
        true
    }
    fn neighbor_cutoff(&self) -> Option<f64> {
        Some(self.potential.cutoff())
    }
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        // Skin-aware rebuild: cheap no-op when no atom has crossed the skin.
        neighbors.validate_pair_scales(system.n())?;
        neighbors.build(system);
        self.compute_with_built_neighbors(system, virial, neighbors)
    }
}

impl PairForceModel {
    /// Shared pair pass. Caller has built the list for this exact system;
    /// snapshots must consume their bounded list even with zero Verlet skin.
    pub(crate) fn compute_with_built_neighbors(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &NeighborList,
    ) -> Result<f64> {
        let cutoff2 = self.potential.cutoff2();
        let n = system.n();
        let box_size = system.box_size;
        let pairs = neighbors.pairs();

        if n == 0 {
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
                neighbors.scales(),
                n,
                &system.positions,
                &system.species,
                box_size,
                cutoff2,
                serial,
                &mut system.forces,
                &self.scratch,
            )
        } else {
            pair_pass(
                &*self.potential,
                pairs,
                neighbors.scales(),
                n,
                &system.positions,
                &system.species,
                box_size,
                cutoff2,
                serial,
                &mut system.forces,
                &self.scratch,
            )
        };

        // M3 accumulate contract: `pair_pass` accumulated directly into
        // `system.forces` (the caller zeroed once before the pass) and returned
        // the PE + virial contributions; compose via virial.
        virial.w += w;
        Ok(pe)
    }
}

impl PairForceModel {
    /// Deterministic force evaluation: identical numerics regardless of
    /// thread count or machine, because the pair sum is accumulated in a fixed
    /// serial order (the parallel path's chunked reduction is not
    /// order-deterministic at >1 thread — OpenMM solves the same problem with
    /// `CpuDeterministicForces`). Slightly slower at large N; for the paper,
    /// run-to-run identifiability is worth it when it matters.
    pub fn compute_deterministic(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        neighbors.validate_pair_scales(system.n())?;
        neighbors.build(system);
        let cutoff2 = self.potential.cutoff2();
        let n = system.n();
        if n == 0 {
            return Ok(0.0);
        }
        let (pe, w) = if let Some(lj) = self
            .potential
            .as_any()
            .and_then(|a| a.downcast_ref::<LennardJones>())
        {
            pair_pass(
                lj,
                neighbors.pairs(),
                neighbors.scales(),
                n,
                &system.positions,
                &system.species,
                system.box_size,
                cutoff2,
                true, // serial
                &mut system.forces,
                &self.scratch,
            )
        } else {
            pair_pass(
                &*self.potential,
                neighbors.pairs(),
                neighbors.scales(),
                n,
                &system.positions,
                &system.species,
                system.box_size,
                cutoff2,
                true, // serial
                &mut system.forces,
                &self.scratch,
            )
        };
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
/// model-owned per-chunk scratch (parallel). Returns `(PE, virial)`.
///
/// Generic over the potential so the `LennardJones` fast path monomorphizes
/// `energy_force` (inline + constant folding); the dyn path is unchanged.
#[allow(clippy::too_many_arguments)]
fn pair_pass<P: PairPotential + ?Sized>(
    potential: &P,
    pairs: &[(usize, usize)],
    scales: &crate::neighbor::PairScaleSet,
    n: usize,
    positions: &[f64],
    species: &[u32],
    box_size: crate::system::SimBox,
    cutoff2: f64,
    serial: bool,
    forces_out: &mut [f64],
    scratch: &Mutex<PairScratch>,
) -> (f64, f64) {
    // A nonblocking lease is essential: a caller may reenter this shared
    // model while the owner waits on work in the same Rayon pool. Poisoned
    // scratch also falls back to the safe, unchanged serial arithmetic.
    let workspace = if serial {
        None
    } else {
        scratch.try_lock().ok()
    };
    match workspace {
        None => {
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
                let scale = scales.scale(i, j).0;
                if scale == 0.0 {
                    continue;
                }
                let (u, f_over_r) = potential.energy_force(r2, species[i], species[j]);
                let (u, f_over_r) = (u * scale, f_over_r * scale);
                let fv = d * f_over_r;
                for k in 0..3 {
                    forces_out[3 * i + k] += fv.0[k];
                    forces_out[3 * j + k] -= fv.0[k];
                }
                pe += u;
                w += d.dot(fv);
            }
            (pe, w)
        }
        Some(mut workspace) => {
            let useful = rayon::current_num_threads().min(pairs.len()).max(1);
            let chunk = pairs.len().div_ceil(useful).max(1);
            let count = pairs.len().div_ceil(chunk);
            // Dropping excess partitions bounds retained buffers after pool
            // shrinkage; each remaining Vec retains its atom-count high water.
            workspace
                .partitions
                .resize_with(count, PairPartition::default);
            for partition in &mut workspace.partitions {
                partition.forces.resize(3 * n, 0.0);
            }
            workspace
                .partitions
                .par_iter_mut()
                .zip(pairs.par_chunks(chunk))
                .for_each(|(acc, pairs)| {
                    acc.forces.fill(0.0);
                    acc.pe = 0.0;
                    acc.w = 0.0;
                    for &(i, j) in pairs {
                        let ri =
                            Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
                        let rj =
                            Vec3::new(positions[3 * j], positions[3 * j + 1], positions[3 * j + 2]);
                        let d = box_size.min_image(ri - rj);
                        let r2 = d.norm2();
                        if r2 >= cutoff2 || r2 < 1e-14 {
                            continue;
                        }
                        let scale = scales.scale(i, j).0;
                        if scale == 0.0 {
                            continue;
                        }
                        let (u, f_over_r) = potential.energy_force(r2, species[i], species[j]);
                        let (u, f_over_r) = (u * scale, f_over_r * scale);
                        let fv = d * f_over_r;
                        for k in 0..3 {
                            acc.forces[3 * i + k] += fv.0[k];
                            acc.forces[3 * j + k] -= fv.0[k];
                        }
                        acc.pe += u;
                        acc.w += d.dot(fv);
                    }
                });
            let (mut pe, mut w) = (0.0, 0.0);
            for partition in &workspace.partitions {
                for (dst, &src) in forces_out.iter_mut().zip(&partition.forces) {
                    *dst += src;
                }
                pe += partition.pe;
                w += partition.w;
            }
            (pe, w)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prebuilt_bounded_zero_skin_list_reaches_the_parallel_pair_pass() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
        let mut s = crate::lattice::fcc(1.6796, [10, 5, 25], "Ar");
        let mut list = NeighborList::with_skin(2.5, 0.0);
        list.rebuild_bounded(&s, 200_000).unwrap();
        let count = list.n_pairs();
        assert!(count >= SERIAL_PAIR_THRESHOLD);
        let mut reference = s.clone();
        let mut v = Virial::default();
        let mut rv = Virial::default();
        // This immutable-list seam cannot perform an unrestricted rebuild.
        let u = pool.install(|| {
            model
                .compute_with_built_neighbors(&mut s, &mut v, &list)
                .unwrap()
        });
        let ru = model
            .compute_deterministic(&mut reference, &mut rv, &mut list)
            .unwrap();
        assert_eq!(list.n_pairs(), count);
        assert!((u - ru).abs() <= 1e-8 + 1e-12 * ru.abs());
        assert!((v.w - rv.w).abs() <= 1e-8 + 1e-12 * rv.w.abs());
        for (a, b) in s.forces.iter().zip(&reference.forces) {
            assert!((a - b).abs() <= 1e-10 + 1e-12 * b.abs());
        }
    }
}
