//! Spatial domain decomposition + ghost-atom exchange (M9).
//!
//! Partitions an orthorhombic periodic box into a 3-D grid of
//! `nr[0] × nr[1] × nr[2]` subdomains — one per parallel rank. Each rank *owns*
//! the atoms whose centers fall in its subdomain; atoms in the 26 neighbor
//! subdomains are imported as read-only *ghosts* so the local pair kernel sees
//! a complete neighborhood (a pair with true minimum-image distance `< cutoff`
//! has both partners present as owned-or-ghost on the owning rank, provided each
//! subdomain edge `≥ cutoff` — enforced by [`RankGrid::balanced`]).
//!
//! Pair work is partitioned by ownership of the lower-global-id atom: rank `r`
//! evaluates pair `(i, j)` (with `gid(i) < gid(j)`) iff atom `i` is owned by
//! `r`. The partner `j` may be owned by `r` (a *local* pair — Newton's third
//! law writes both `f_i` and `f_j`) or by another rank (a *ghost* on `r` — its
//! force is written to the local ghost copy and reduced back to its owner by
//! [`Communicator::reverse_ghost_reduce`]). Each pair has exactly one evaluator
//! and each atom exactly one force accumulator, so the decomposed
//! force/energy/virial summed over ranks equals the serial half-list
//! evaluation up to floating-point summation order (the tests assert ~1e-9
//! agreement) — the correctness anchor in `tests/decomp.rs`.
//!
//! The [`Communicator`] trait abstracts the transport. The default
//! [`LocalCommunicator`] simulates `P` ranks in one process (forces written
//! directly to a shared `System`, so the reverse reduce is a no-op); the `mpi`
//! feature supplies an `MpiCommunicator` backed by the `mpi` crate for real
//! distributed runs. Atom migration (re-ownership after a body crosses a
//! subdomain face) is ownership reassignment recomputed from current positions.

use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::potential::PairPotential;
use crate::system::{SimBox, System};
use crate::vec3::Vec3;

/// A 3-D rank grid `nr[0] × nr[1] × nr[2]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RankGrid {
    pub nr: [usize; 3],
}

impl RankGrid {
    /// The trivial single-rank grid.
    pub fn single() -> Self {
        Self { nr: [1, 1, 1] }
    }

    /// Number of ranks (subdomains).
    pub fn n_ranks(&self) -> usize {
        self.nr[0] * self.nr[1] * self.nr[2]
    }

    /// Linear rank id of the grid cell `(ix, iy, iz)`.
    pub fn rank_of(&self, cell: [usize; 3]) -> usize {
        cell[0] * self.nr[1] * self.nr[2] + cell[1] * self.nr[2] + cell[2]
    }

    /// Grid cell `(ix, iy, iz)` of a linear rank id.
    pub fn cell_of(&self, rank: usize) -> [usize; 3] {
        let [_, ny, nz] = self.nr;
        let ix = rank / (ny * nz);
        let iy = (rank / nz) % ny;
        let iz = rank % nz;
        [ix, iy, iz]
    }

    /// Choose a grid with as many ranks as possible up to `target`, subject to
    /// every subdomain edge being `>= min_edge` (so a single shell of neighbor
    /// subdomains is a sufficient ghost halo). Falls back to the single-rank grid
    /// if the box is too small to split along any axis.
    pub fn balanced(box_size: Vec3, min_edge: f64, target: usize) -> Self {
        let max_per_axis = [
            ((box_size.0[0] / min_edge).floor() as usize).max(1),
            ((box_size.0[1] / min_edge).floor() as usize).max(1),
            ((box_size.0[2] / min_edge).floor() as usize).max(1),
        ];
        // Greedily split the axis that is currently longest (relative to its
        // max) until we reach `target` ranks or no axis can split further.
        let mut nr = [1usize; 3];
        while nr[0] * nr[1] * nr[2] < target {
            // Pick the axis with the most headroom (current splits below max and
            // largest ratio box/splits).
            let mut best = None;
            let mut best_score = -1.0;
            for k in 0..3 {
                if nr[k] < max_per_axis[k] {
                    let edge = box_size.0[k] / (nr[k] + 1) as f64;
                    let score = edge; // prefer splitting the axis that keeps edges large
                    if score > best_score {
                        best_score = score;
                        best = Some(k);
                    }
                }
            }
            match best {
                Some(k) => nr[k] += 1,
                None => break, // no axis can split further
            }
        }
        Self { nr }
    }
}

/// Ownership assignment of atoms to ranks for an orthorhombic box.
///
/// Built from a [`SimBox`] (orthorhombic — triclinic decomposition is left for
/// a later roadmap item; the barostats are orthorhombic-only too), a pair
/// `cutoff`, and a [`RankGrid`]. The grid is constrained so each subdomain edge
/// `>= cutoff`; if the supplied grid violates this the constructor errors.
pub struct DomainDecomposition {
    grid: RankGrid,
    /// Per-axis subdomain edge length.
    sub: [f64; 3],
    box_size: SimBox,
    cutoff: f64,
}

impl DomainDecomposition {
    /// Construct a decomposition. Errors if the box is triclinic or any
    /// subdomain edge is below `cutoff` (which would make a one-shell ghost
    /// halo insufficient).
    pub fn new(box_size: SimBox, cutoff: f64, grid: RankGrid) -> Result<Self> {
        if box_size.is_triclinic() {
            return Err(ModenaError::Config(
                "domain decomposition requires an orthorhombic box".into(),
            ));
        }
        if !(cutoff > 0.0) {
            return Err(ModenaError::Config("cutoff must be > 0".into()));
        }
        let sub = [
            box_size.size.0[0] / grid.nr[0] as f64,
            box_size.size.0[1] / grid.nr[1] as f64,
            box_size.size.0[2] / grid.nr[2] as f64,
        ];
        for (k, &edge) in sub.iter().enumerate() {
            if edge + 1e-12 < cutoff {
                return Err(ModenaError::Config(format!(
                    "subdomain edge {} ({:.4}) < cutoff ({:.4}); use fewer ranks",
                    k, edge, cutoff
                )));
            }
        }
        Ok(Self {
            grid,
            sub,
            box_size,
            cutoff,
        })
    }

    pub fn grid(&self) -> RankGrid {
        self.grid
    }

    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Which rank owns the (sub)cell containing position `r`. Positions are
    /// wrapped into the primary image first.
    pub fn owner_of(&self, r: Vec3) -> usize {
        let w = self.box_size.wrap(r);
        let ix = (w.0[0] / self.sub[0]).floor() as usize;
        let iy = (w.0[1] / self.sub[1]).floor() as usize;
        let iz = (w.0[2] / self.sub[2]).floor() as usize;
        self.grid.rank_of([
            ix.min(self.grid.nr[0] - 1),
            iy.min(self.grid.nr[1] - 1),
            iz.min(self.grid.nr[2] - 1),
        ])
    }

    /// Owned atom indices for `rank`, in ascending global-id order.
    pub fn owned(&self, system: &System, rank: usize) -> Vec<usize> {
        let mut v: Vec<usize> = (0..system.n())
            .filter(|&i| self.owner_of(system.pos(i)) == rank)
            .collect();
        v.sort_unstable();
        v
    }

    /// Ghost atom indices for `rank`: atoms owned by any of the 26 neighbor
    /// subdomains (with periodic wrap), *excluding self*. These are the imported
    /// halo atoms. For a single-rank grid there are no neighbors, so this is
    /// empty (the rank already owns every atom).
    pub fn ghosts(&self, system: &System, rank: usize) -> Vec<usize> {
        let [ix, iy, iz] = self.grid.cell_of(rank);
        let [nx, ny, nz] = self.grid.nr;
        let neighbor_ranks: Vec<usize> = {
            let mut set = Vec::new();
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        if dx == 0 && dy == 0 && dz == 0 {
                            continue;
                        }
                        let jx = ((ix as i32 + dx + nx as i32) % nx as i32) as usize;
                        let jy = ((iy as i32 + dy + ny as i32) % ny as i32) as usize;
                        let jz = ((iz as i32 + dz + nz as i32) % nz as i32) as usize;
                        let r = self.grid.rank_of([jx, jy, jz]);
                        if r != rank {
                            set.push(r);
                        }
                    }
                }
            }
            set.sort_unstable();
            set.dedup();
            set
        };
        let mut g: Vec<usize> = (0..system.n())
            .filter(|&i| neighbor_ranks.contains(&self.owner_of(system.pos(i))))
            .collect();
        g.sort_unstable();
        g
    }

    /// Count atoms that changed owning rank after positions moved (migration
    /// diagnostic). In the in-process model ownership is recomputed from
    /// positions each call, so "migration" is implicit; this counts how many
    /// atoms would need to be sent/received in a distributed run.
    pub fn migration_count(&self, prev_owners: &[usize], system: &System) -> usize {
        prev_owners
            .iter()
            .enumerate()
            .filter(|&(i, &old)| old != self.owner_of(system.pos(i)))
            .count()
    }
}

/// Transport abstraction for a parallel decomposition. The default
/// [`LocalCommunicator`] simulates ranks in one process; the `mpi` feature
/// provides a real MPI backend.
///
/// All methods take the *linear rank* and operate on flattened `f64` buffers.
/// The contract is deliberately minimal: ghost import (forward), reverse ghost
/// force reduction, and scalar all-reduce (energy/virial).
pub trait Communicator: Send + Sync {
    /// Number of ranks.
    fn size(&self) -> usize;
    /// This rank's id.
    fn rank(&self) -> usize;

    /// Reverse-ghost reduce: each rank has accumulated forces on its *ghost*
    /// atoms (the partner side of cross-rank pairs). Sum those ghost forces back
    /// to their owners. For the in-process communicator this is a no-op (forces
    /// were written directly to the shared `System`); for MPI it is a halo
    /// exchange of force contributions.
    ///
    /// `forces` is the full `3N` global force buffer (in-process) or the local
    /// owned+ghost buffer (MPI); `ghost_local` / `ghost_global` identify the
    /// ghost slots and their owning global indices.
    fn reverse_ghost_reduce(
        &self,
        forces: &mut [f64],
        ghost_local: &[usize],
        ghost_global: &[usize],
    );

    /// Sum a scalar across all ranks (energy, virial).
    fn allreduce_sum(&self, value: f64) -> f64;
}

/// In-process communicator: all ranks share one `System`, so the reverse ghost
/// reduce is a no-op (forces are written directly to global indices) and the
/// all-reduce is identity. This is the default used by [`DomainForceModel`] and
/// by the correctness tests; it lets the decomposition run and be verified
/// without any parallel runtime.
#[derive(Default)]
pub struct LocalCommunicator {
    size: usize,
    rank: usize,
}

impl LocalCommunicator {
    /// A communicator for a single rank (the serial-equivalent case).
    pub fn single() -> Self {
        Self { size: 1, rank: 0 }
    }

    /// A communicator reporting `size` ranks, this being `rank`.
    pub fn new(size: usize, rank: usize) -> Self {
        assert!(rank < size, "rank must be < size");
        Self { size, rank }
    }
}

impl Communicator for LocalCommunicator {
    fn size(&self) -> usize {
        self.size
    }
    fn rank(&self) -> usize {
        self.rank
    }
    fn reverse_ghost_reduce(&self, _forces: &mut [f64], _gl: &[usize], _gg: &[usize]) {
        // No-op: in the in-process model forces are written directly to the
        // shared global force buffer, so ghost contributions are already in
        // place. A distributed backend would halo-exchange here.
    }
    fn allreduce_sum(&self, value: f64) -> f64 {
        value
    }
}

/// Evaluate the pair force/energy/virial for a single `rank` of a
/// [`DomainDecomposition`], using the in-process shared `System`. This is the
/// per-rank kernel underlying [`DomainForceModel`]; summed over all ranks it
/// equals the serial half-list evaluation up to floating-point summation order.
///
/// The rank builds a local view (owned atoms first, then ghosts), a local cell
/// list over that view (using the *global* box for minimum-image), and evaluates
/// every pair whose lower-global-id atom it owns — applying Newton's third law
/// to both atoms (the partner is either a co-owned local atom or a ghost whose
/// owner will read the accumulated force via the reverse reduce).
fn rank_pair_compute(
    decomp: &DomainDecomposition,
    comm: &dyn Communicator,
    potential: &dyn PairPotential,
    system: &mut System,
    virial: &mut Virial,
    rank: usize,
) -> f64 {
    let owned = decomp.owned(system, rank);
    let ghosts = decomp.ghosts(system, rank);
    let local_n = owned.len() + ghosts.len();
    if owned.is_empty() {
        return 0.0;
    }

    // Build a local subsystem: owned atoms (keeping global id as local index
    // basis via the map), then ghosts. The subsystem carries the GLOBAL box so
    // minimum-image wrapping is correct across the full periodic cell.
    let mut sub = System::new(system.box_size);
    sub.resize(local_n);
    let mut g2l = vec![usize::MAX; system.n()];
    for (local, &global) in owned.iter().chain(ghosts.iter()).enumerate() {
        g2l[global] = local;
        sub.set_pos(local, system.pos(global));
        sub.masses[local] = system.masses[global];
        sub.species[local] = system.species[global];
        sub.charges[local] = system.charges[global];
    }
    // g (local -> global) for the owned block + ghost block.
    let g: Vec<usize> = owned.iter().chain(ghosts.iter()).copied().collect();
    let n_owned = owned.len();

    let mut nl = NeighborList::new(decomp.cutoff());
    nl.build(&sub);
    let cutoff2 = potential.cutoff2();
    let box_size = sub.box_size;

    // Build symmetric adjacency from the half list so each owned atom sees all
    // its neighbors (owned and ghost).
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); local_n];
    for &(a, b) in nl.pairs() {
        adj[a].push(b);
        adj[b].push(a);
    }

    let mut pe = 0.0f64;
    let mut w = 0.0f64;
    let positions = &sub.positions;
    let species = &sub.species;
    let n_owned_local = n_owned;

    // Each owned local atom a evaluates pairs (a, b) with gid(a) < gid(b).
    for a in 0..n_owned_local {
        let ga = g[a];
        let ri = Vec3::new(positions[3 * a], positions[3 * a + 1], positions[3 * a + 2]);
        for &b in &adj[a] {
            let gb = g[b];
            if ga >= gb {
                continue; // pair owned by the rank owning the lower-gid atom
            }
            let rj = Vec3::new(positions[3 * b], positions[3 * b + 1], positions[3 * b + 2]);
            let d = box_size.min_image(ri - rj);
            let r2 = d.norm2();
            if r2 >= cutoff2 || r2 < 1e-14 {
                continue;
            }
            let (u, f_over_r) = potential.energy_force(r2, species[a], species[b]);
            let fv = d * f_over_r;
            // Newton's third law on both global atoms (in-process: direct writes).
            system.add_force(ga, fv);
            system.add_force(gb, -fv);
            pe += u;
            w += d.dot(fv);
        }
    }

    virial.w += w;
    // In-process: reverse reduce is a no-op (forces already global).
    comm.reverse_ghost_reduce(&mut system.forces, &[], &[]);
    pe
}

/// A pair force model evaluated through spatial domain decomposition.
///
/// Wraps a [`PairPotential`] and a [`RankGrid`]; `compute` loops over every rank
/// in the grid (using a [`LocalCommunicator`]), gathering owned + ghost atoms,
/// building a per-rank cell list, and evaluating the rank's pairs. The result
/// equals [`crate::force_model::PairForceModel`] up to floating-point summation
/// order for any grid (including multi-rank grids) — the decomposition only
/// partitions the work, verified in `tests/decomp.rs`.
///
/// This is the in-process reference: it demonstrates and tests the
/// decomposition algorithm (ownership, ghost halo, per-rank cell lists,
/// migration) without a parallel runtime. Production distributed runs use the
/// `mpi` feature's `MpiCommunicator` with the same per-rank kernel; production
/// throughput uses the `gpu` feature's `GpuPairForceModel`.
pub struct DomainForceModel {
    potential: Box<dyn PairPotential>,
    decomp: DomainDecomposition,
    comm: Box<dyn Communicator>,
}

impl DomainForceModel {
    pub fn new(
        potential: impl PairPotential + 'static,
        decomp: DomainDecomposition,
        comm: Box<dyn Communicator>,
    ) -> Self {
        Self {
            potential: Box::new(potential),
            decomp,
            comm,
        }
    }

    /// Convenience: single-rank decomposition (equivalent to the serial pair
    /// model, but routed through the decomposition kernel).
    pub fn single_rank(
        potential: impl PairPotential + 'static,
        box_size: SimBox,
        cutoff: f64,
    ) -> Self {
        let decomp = DomainDecomposition::new(box_size, cutoff, RankGrid::single())
            .expect("single-rank decomposition is always valid for an orthorhombic box");
        Self::new(potential, decomp, Box::new(LocalCommunicator::single()))
    }
}

impl ForceModel for DomainForceModel {
    fn neighbor_cutoff(&self) -> Option<f64> {
        Some(self.potential.cutoff())
    }
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        if !neighbors.pair_scales().is_empty() {
            return Err(ModenaError::Config(
                "domain force backend does not support special pair scales or exclusions".into(),
            ));
        }
        if system.box_size != self.decomp.box_size {
            return Err(ModenaError::Config(
                "domain force backend requires its original box".into(),
            ));
        }
        let mut pe = 0.0f64;
        for rank in 0..self.comm.size() {
            // The in-process kernel writes directly to the shared system; each
            // rank's contribution adds in. (For a real MPI backend this loop
            // collapses to a single `self.comm.rank()`.)
            pe += rank_pair_compute(
                &self.decomp,
                &*self.comm,
                &*self.potential,
                system,
                virial,
                rank,
            );
        }
        Ok(pe)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potential::LennardJones;
    use crate::system::SimBox;

    fn two_atom_system() -> System {
        let mut sys = System::new(SimBox::cubic(10.0));
        sys.resize(2);
        sys.set_pos(0, Vec3::new(0.0, 0.0, 0.0));
        sys.set_pos(1, Vec3::new(1.5, 0.0, 0.0));
        sys
    }

    #[test]
    fn rank_grid_balanced_respects_min_edge() {
        let grid = RankGrid::balanced(Vec3::new(10.0, 10.0, 10.0), 2.5, 8);
        // 10 / 2.5 = 4 per axis max → up to 64 ranks, but target 8.
        assert!(grid.n_ranks() <= 8);
        for k in 0..3 {
            let edge = 10.0 / grid.nr[k] as f64;
            assert!(edge + 1e-12 >= 2.5, "edge {edge} < min_edge");
        }
    }

    #[test]
    fn rank_grid_balanced_small_box_is_single() {
        let grid = RankGrid::balanced(Vec3::new(2.0, 2.0, 2.0), 2.5, 8);
        assert_eq!(grid.n_ranks(), 1);
    }

    #[test]
    fn owner_of_partitions_all_atoms_no_overlap() {
        let sys = two_atom_system();
        let decomp =
            DomainDecomposition::new(SimBox::cubic(10.0), 2.5, RankGrid { nr: [2, 1, 1] }).unwrap();
        let mut covered = 0;
        for r in 0..decomp.grid().n_ranks() {
            covered += decomp.owned(&sys, r).len();
        }
        assert_eq!(covered, sys.n());
        // Atom 0 at x=0 → rank 0; atom 1 at x=1.5 → rank 0 (subdomain edge 5).
        assert_eq!(decomp.owner_of(sys.pos(0)), 0);
        assert_eq!(decomp.owner_of(sys.pos(1)), 0);
    }

    #[test]
    fn owner_of_wraps_periodically() {
        let decomp =
            DomainDecomposition::new(SimBox::cubic(10.0), 2.5, RankGrid { nr: [2, 1, 1] }).unwrap();
        // Subdomain 0 = x in [0,5), subdomain 1 = x in [5,10). A negative x
        // wraps into [0,10) before binning.
        assert_eq!(decomp.owner_of(Vec3::new(4.9, 0.0, 0.0)), 0);
        assert_eq!(decomp.owner_of(Vec3::new(5.0, 0.0, 0.0)), 1);
        assert_eq!(decomp.owner_of(Vec3::new(9.9, 0.0, 0.0)), 1);
        // -0.1 wraps to 9.9 → rank 1.
        assert_eq!(decomp.owner_of(Vec3::new(-0.1, 0.0, 0.0)), 1);
    }

    #[test]
    fn rejects_subdomain_smaller_than_cutoff() {
        let err = DomainDecomposition::new(SimBox::cubic(4.0), 2.5, RankGrid { nr: [2, 1, 1] });
        assert!(err.is_err(), "edge 2.0 < cutoff 2.5 must error");
    }

    #[test]
    fn rejects_triclinic() {
        let tri = SimBox::triclinic([
            Vec3::new(10.0, 0.0, 0.0),
            Vec3::new(0.5, 10.0, 0.0),
            Vec3::new(0.0, 0.5, 10.0),
        ]);
        let err = DomainDecomposition::new(tri, 2.5, RankGrid::single());
        assert!(err.is_err());
    }

    #[test]
    fn domain_force_model_matches_serial_single_rank() {
        // The single-rank DomainForceModel must reproduce the serial pair model
        // exactly (it's the same kernel, just partitioned into one rank).
        let lj = LennardJones::new(1.0, 1.0, 2.5).unwrap();
        let mut sys = two_atom_system();

        let model = DomainForceModel::single_rank(lj.clone(), SimBox::cubic(10.0), 2.5);
        let mut virial = Virial::default();
        let mut nl = NeighborList::new(2.5);
        let pe = model.compute(&mut sys, &mut virial, &mut nl).unwrap();

        // Serial reference: direct LJ over the one pair.
        let r = 1.5_f64;
        let (u_ref, _) = lj.energy_force(r * r, 0, 0);
        assert!((pe - u_ref).abs() < 1e-12, "pe {pe} vs {u_ref}");
        let f0 = sys.force(0).0[0];
        let f1 = sys.force(1).0[0];
        assert!((f0 + f1).abs() < 1e-12, "Newton's 3rd law: f0 + f1 = 0");
        // Distinct forces (nonzero, equal & opposite).
        assert!(f0.abs() > 1e-6);
    }
}
