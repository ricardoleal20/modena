//! MPI domain-decomposition integration test (M9, `mpi` feature).
//!
//! Run with:
//! ```sh
//! mpiexec -n 2 cargo test --features mpi --test mpi_decomp -- --nocapture
//! ```
//!
//! Each rank owns the atoms in its subdomain, exchanges a ghost halo with its
//! neighbors, evaluates its pairs, and reverse-reduces ghost forces. Rank 0
//! gathers every owned force and compares the assembled global force field +
//! energy + virial to the serial half-list reference — they must match exactly.

#![cfg(feature = "mpi")]

use modena_core::{
    decomp::{DomainDecomposition, RankGrid},
    force_model::{ForceModel, PairForceModel, Virial},
    lattice,
    mpi_comm::{mpi_decomposed_pair_compute, OwnedAtoms},
    neighbor::NeighborList,
    potential::LennardJones,
    system::System,
};
use mpi::traits::*;

/// Build the full system identically on every rank.
fn full_system() -> System {
    lattice::fcc(1.5, [5, 5, 5], "Ar") // box 7.5, 500 atoms
}

fn lj() -> LennardJones {
    LennardJones::new(1.0, 1.0, 2.5).unwrap()
}

#[test]
fn mpi_decomposition_matches_serial() {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let size = world.size() as usize;
    let rank = world.rank() as usize;
    // Graceful skip when run as a plain single-process `cargo test` (e.g. CI's
    // `--all-features` build job). The real 2-rank run is driven by a dedicated
    // `mpiexec -n 2` job; see the test module docs.
    if size != 2 {
        eprintln!("skipping MPI decomp test (size={size}); run under `mpiexec -n 2`");
        return;
    }

    let sys = full_system();
    // 2×1×1 grid: subdomain edge 7.5/2 = 3.75 ≥ cutoff 2.5.
    let grid = RankGrid { nr: [2, 1, 1] };
    let decomp = DomainDecomposition::new(sys.box_size, 2.5, grid).unwrap();

    // This rank's owned atoms (global box + global ids = original indices).
    let owned_idx = decomp.owned(&sys, rank);
    let mut local = System::new(sys.box_size);
    local.resize(owned_idx.len());
    let mut gids = Vec::with_capacity(owned_idx.len());
    for (dst, &src) in owned_idx.iter().enumerate() {
        local.set_pos(dst, sys.pos(src));
        local.masses[dst] = sys.masses[src];
        local.species[dst] = sys.species[src];
        local.charges[dst] = sys.charges[src];
        gids.push(src as u64);
    }
    let owned = OwnedAtoms {
        system: local,
        gids,
    };

    let (owned_forces, pe, w) = mpi_decomposed_pair_compute(&world, &decomp, &lj(), &owned);

    // Pack [gid, fx, fy, fz] per owned atom for the gather.
    let n_owned = owned_idx.len();
    let mut pack = Vec::with_capacity(4 * n_owned);
    for i in 0..n_owned {
        pack.push(owned_idx[i] as f64);
        pack.push(owned_forces[3 * i]);
        pack.push(owned_forces[3 * i + 1]);
        pack.push(owned_forces[3 * i + 2]);
    }

    // Gather to rank 0: ranks >0 send, rank 0 receives in order.
    let mut assembled: Vec<f64> = vec![0.0; 3 * sys.n()];
    if rank == 0 {
        // Rank 0's own contribution.
        for i in 0..n_owned {
            let gid = owned_idx[i] as usize;
            assembled[3 * gid] += owned_forces[3 * i];
            assembled[3 * gid + 1] += owned_forces[3 * i + 1];
            assembled[3 * gid + 2] += owned_forces[3 * i + 2];
        }
        for r in 1..size {
            let (recv, _): (Vec<f64>, _) = world.process_at_rank(r as mpi::Rank).receive_vec();
            for chunk in recv.chunks_exact(4) {
                let gid = chunk[0] as usize;
                assembled[3 * gid] += chunk[1];
                assembled[3 * gid + 1] += chunk[2];
                assembled[3 * gid + 2] += chunk[3];
            }
        }
    } else {
        world.process_at_rank(0).send(&pack[..]);
    }

    if rank == 0 {
        // Serial reference over the full system.
        let mut s = sys.clone();
        s.zero_forces();
        let mut vir = Virial::default();
        let mut nl = NeighborList::new(2.5);
        let model = PairForceModel::new(lj());
        let pe_ref = model.compute(&mut s, &mut vir, &mut nl).unwrap();

        // The decomposition is exact in real arithmetic: each pair is evaluated
        // once and each force accumulated once. The only difference vs the
        // serial rayon half-list is floating-point *summation order* (per-rank
        // accumulation + MPI all-reduce), so a relative tolerance with a floor
        // is the principled check — it rejects any real O(1) ownership /
        // double-counting bug while tolerating rounding on magnitudes ~1e4.
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(b.abs()).max(1.0);
        assert!(close(pe, pe_ref), "MPI energy {pe} vs serial {pe_ref}");
        assert!(close(w, vir.w), "MPI virial {w} vs serial {}", vir.w);
        for i in 0..sys.n() {
            for k in 0..3 {
                let mpi_f = assembled[3 * i + k];
                let ser_f = s.forces[3 * i + k];
                assert!(
                    close(mpi_f, ser_f),
                    "force atom {i} axis {k}: MPI {mpi_f} vs serial {ser_f}"
                );
            }
        }
    }

    // Keep `pack` alive until the send completes (blocking send already did).
    let _ = pack;
}
