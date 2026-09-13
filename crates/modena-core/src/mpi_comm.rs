//! MPI backend for spatial domain decomposition (M9, `mpi` feature).
//!
//! Distributed-memory counterpart of [`crate::decomp`]. Each MPI rank owns a
//! local subset of atoms (those falling in its subdomain) plus a read-only ghost
//! halo imported from the 26 neighbor ranks. The driver performs the real
//! LAMMPS-style halo exchange:
//!
//! 1. **Forward exchange** — each rank sends its owned atoms (position, species,
//!    global id, packed into one `f64` buffer) to every neighbor rank; received
//!    atoms become local ghosts.
//! 2. **Local force eval** — a neighbor list over owned + ghosts (using the
//!    global box for minimum-image); each owned atom `i` evaluates pairs
//!    `(i, j)` with `gid(i) < gid(j)`. Force on the owned `i` is applied
//!    directly; force on a ghost `j` is accumulated into a ghost-force buffer.
//! 3. **Reverse exchange** — ghost-force buffers are sent back to the ghosts'
//!    owners, which sum them into their owned forces (Newton's third law across
//!    rank boundaries).
//! 4. **All-reduce** — potential energy and virial are summed across ranks.
//!
//! The result is bit-identical to the serial half-list evaluation: each pair is
//! evaluated once (by the rank owning the lower-gid atom) and each atom's force
//! is accumulated exactly once. Verified by `tests/mpi_decomp.rs` under
//! `mpiexec -n N`.
//!
//! Deadlock-free point-to-point exchange: with each neighbor, the higher-ranked
//! process sends first and the lower-ranked receives first (canonical ordering),
//! so no pair of blocking `send`/`receive_vec` calls ever wait on each other.

#![cfg(feature = "mpi")]

use crate::decomp::{DomainDecomposition, RankGrid};
use crate::neighbor::NeighborList;
use crate::potential::PairPotential;
use crate::system::System;
use crate::vec3::Vec3;
use mpi::collective::SystemOperation;
use mpi::traits::*;

/// Per-rank owned-atom data: a local [`System`] (carrying the *global* box so
/// minimum-image is correct) plus the global id of each owned atom.
pub struct OwnedAtoms {
    pub system: System,
    pub gids: Vec<u64>,
}

/// The 26 distinct neighbor ranks of `rank` in the grid (periodic, excluding
/// self), in ascending order.
fn neighbor_ranks(grid: RankGrid, rank: usize) -> Vec<usize> {
    let [ix, iy, iz] = grid.cell_of(rank);
    let [nx, ny, nz] = grid.nr;
    let mut set = Vec::new();
    for dx in -1i32..=1 {
        for dy in -1i32..=1 {
            for dz in -1i32..=1 {
                if dx == 0 && dy == 0 && dz == 0 {
                    continue;
                }
                let jx = ((ix as i32 + dx + nx as i32) % nx as i32) as usize;
                let jy = ((iy as i32 + dy + ny as i32) % ny as i32) as usize;
                let jz = ((iz as i32 + dz + nz as i32) % nz as i32) as usize;
                let r = grid.rank_of([jx, jy, jz]);
                if r != rank {
                    set.push(r);
                }
            }
        }
    }
    set.sort_unstable();
    set.dedup();
    set
}

/// Exchange a typed buffer with each neighbor (one value per neighbor). For each
/// neighbor `nb`: if `nb > rank`, send then receive; else receive then send.
/// Returns the received buffers in the same order as `neighbors`.
fn exchange_with_neighbors<C: Communicator, T: Equivalence + Clone>(
    world: &C,
    rank: usize,
    neighbors: &[usize],
    send_bufs: &[Vec<T>],
) -> Vec<Vec<T>> {
    let mut received: Vec<Option<Vec<T>>> = (0..neighbors.len()).map(|_| None).collect();
    for (idx, &nb) in neighbors.iter().enumerate() {
        let dest = world.process_at_rank(nb as mpi::Rank);
        if nb > rank {
            dest.send(&send_bufs[idx][..]);
            let (msg, _) = world.process_at_rank(nb as mpi::Rank).receive_vec::<T>();
            received[idx] = Some(msg);
        } else {
            let (msg, _) = world.process_at_rank(nb as mpi::Rank).receive_vec::<T>();
            // Received first; now send.
            dest.send(&send_bufs[idx][..]);
            received[idx] = Some(msg);
        }
    }
    received
        .into_iter()
        .map(|o| o.unwrap_or_default())
        .collect()
}

/// Pack owned atoms as `[x, y, z, species, gid]` per atom (5 `f64`s). `species`
/// and `gid` are small integers exactly representable as `f64` (< 2⁵³).
fn pack_owned(owned: &OwnedAtoms) -> Vec<f64> {
    let n = owned.system.n();
    let mut buf = Vec::with_capacity(5 * n);
    for i in 0..n {
        buf.push(owned.system.positions[3 * i]);
        buf.push(owned.system.positions[3 * i + 1]);
        buf.push(owned.system.positions[3 * i + 2]);
        buf.push(owned.system.species[i] as f64);
        buf.push(owned.gids[i] as f64);
    }
    buf
}

/// Run the decomposed pair force/energy/virial across all MPI ranks.
///
/// `owned` is this rank's owned atoms (with the global box). Returns the owned
/// forces (3 per owned atom, in `owned.gids` order), the global potential
/// energy, and the global scalar virial — all identical to the serial
/// half-list evaluation.
pub fn mpi_decomposed_pair_compute<C: Communicator, P: PairPotential>(
    world: &C,
    decomp: &DomainDecomposition,
    potential: &P,
    owned: &OwnedAtoms,
) -> (Vec<f64>, f64, f64) {
    let rank = world.rank() as usize;
    let grid = decomp.grid();
    let neighbors = neighbor_ranks(grid, rank);
    let n_owned = owned.system.n();

    // --- Forward exchange: send owned atoms to each neighbor, receive ghosts. ---
    let owned_pack = pack_owned(owned);
    let send_bufs: Vec<Vec<f64>> = neighbors.iter().map(|_| owned_pack.clone()).collect();
    let ghost_packs = exchange_with_neighbors::<C, f64>(world, rank, &neighbors, &send_bufs);

    // Unpack ghosts. Each ghost batch is stored with its source neighbor index
    // so the reverse exchange can route forces back.
    struct GhostBatch {
        gids: Vec<u64>,
        positions: Vec<f64>,
        species: Vec<u32>,
        forces: Vec<f64>, // 3 per ghost, accumulated for the reverse exchange
    }
    let mut ghosts: Vec<GhostBatch> = Vec::new();
    for pack in &ghost_packs {
        let ng = pack.len() / 5;
        let mut gids = Vec::with_capacity(ng);
        let mut positions = Vec::with_capacity(3 * ng);
        let mut species = Vec::with_capacity(ng);
        for k in 0..ng {
            positions.push(pack[5 * k]);
            positions.push(pack[5 * k + 1]);
            positions.push(pack[5 * k + 2]);
            species.push(pack[5 * k + 3] as u32);
            gids.push(pack[5 * k + 4] as u64);
        }
        ghosts.push(GhostBatch {
            gids,
            positions,
            species,
            forces: vec![0.0; 3 * ng],
        });
    }

    // --- Build the local view: owned atoms first, then ghosts. ---
    let n_ghost: usize = ghosts.iter().map(|g| g.gids.len()).sum();
    let local_n = n_owned + n_ghost;
    let mut sub = System::new(owned.system.box_size);
    sub.resize(local_n);
    let mut local_gids = Vec::with_capacity(local_n);
    for i in 0..n_owned {
        sub.set_pos(i, owned.system.pos(i));
        sub.masses[i] = owned.system.masses[i];
        sub.species[i] = owned.system.species[i];
        sub.charges[i] = owned.system.charges[i];
        local_gids.push(owned.gids[i]);
    }
    let mut ghost_offset = n_owned;
    // local ghost index -> (batch index, k-within-batch) for the reverse write.
    let mut ghost_index: Vec<(usize, usize)> = Vec::with_capacity(n_ghost);
    for (bidx, gb) in ghosts.iter().enumerate() {
        for (k, &gid) in gb.gids.iter().enumerate() {
            let local = ghost_offset;
            sub.set_pos(
                local,
                Vec3::new(
                    gb.positions[3 * k],
                    gb.positions[3 * k + 1],
                    gb.positions[3 * k + 2],
                ),
            );
            sub.species[local] = gb.species[k];
            local_gids.push(gid);
            ghost_index.push((bidx, k));
            ghost_offset += 1;
        }
    }

    // --- Local neighbor list over owned + ghosts (global box for min-image). ---
    let mut nl = NeighborList::new(decomp.cutoff());
    nl.build(&sub);
    let cutoff2 = potential.cutoff2();
    let box_size = sub.box_size;
    let positions = &sub.positions;
    let species = &sub.species;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); local_n];
    for &(a, b) in nl.pairs() {
        adj[a].push(b);
        adj[b].push(a);
    }

    let mut owned_forces = vec![0.0f64; 3 * n_owned];
    let mut pe = 0.0f64;
    let mut w = 0.0f64;
    for a in 0..n_owned {
        let ga = local_gids[a];
        let ri = Vec3::new(positions[3 * a], positions[3 * a + 1], positions[3 * a + 2]);
        for &b in &adj[a] {
            let gb = local_gids[b];
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
            for k in 0..3 {
                owned_forces[3 * a + k] += fv.0[k];
            }
            if b < n_owned {
                // Local owned pair: Newton on the local owned atom b.
                for k in 0..3 {
                    owned_forces[3 * b + k] -= fv.0[k];
                }
            } else {
                // Ghost b: accumulate -f into the ghost's force buffer.
                let (gbidx, gk) = ghost_index[b - n_owned];
                let base = 3 * gk;
                for k in 0..3 {
                    ghosts[gbidx].forces[base + k] -= fv.0[k];
                }
            }
            pe += u;
            w += d.dot(fv);
        }
    }

    // --- Reverse exchange: send ghost forces to owners, receive our owned
    //     atoms' force contributions from each neighbor. ---
    let rev_send: Vec<Vec<f64>> = ghosts.iter().map(|g| g.forces.clone()).collect();
    let rev_recv = exchange_with_neighbors::<C, f64>(world, rank, &neighbors, &rev_send);
    // Each neighbor's received buffer has 3 * n_owned entries (the forces for the
    // owned atoms we sent them as ghosts), in the order we packed them.
    for frecv in &rev_recv {
        for k in 0..frecv.len().min(owned_forces.len()) {
            owned_forces[k] += frecv[k];
        }
    }

    // --- All-reduce energy and virial. ---
    let mut pe_global = 0.0f64;
    let mut w_global = 0.0f64;
    world.all_reduce_into(&pe, &mut pe_global, SystemOperation::sum());
    world.all_reduce_into(&w, &mut w_global, SystemOperation::sum());

    (owned_forces, pe_global, w_global)
}
