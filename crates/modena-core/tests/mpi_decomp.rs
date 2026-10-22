//! MPI domain-decomposition integration test (M9, `mpi` feature).
//!
//! Run with:
//! ```sh
//! mpiexec -n 2 cargo test --features mpi --test mpi_decomp -- --nocapture
//! ```
//!
//! A single test exercises both correctness anchors of the distributed pair
//! driver, in sequence:
//!
//! 1. **Single-eval match** — rank 0 gathers every owned force and compares the
//!    assembled global force field + energy + virial to the serial half-list
//!    reference (relative tolerance for FP summation-order differences).
//! 2. **NVE conservation** — a velocity-Verlet trajectory driven by the MPI
//!    force driver across two ranks, with ghost re-exchange every step; total
//!    energy (all-reduced across ranks) must be conserved.
//!
//! Both checks live in one `#[test]` deliberately: MPI tests run as one binary
//! per rank under `mpiexec`, and a single test guarantees every rank executes
//! the same collective calls in the same order (parallel tests / differing
//! order across ranks deadlocks or segfaults under `MPI_THREAD_SINGLE`).

#![cfg(feature = "mpi")]

use modena_core::{
    decomp::{DomainDecomposition, RankGrid},
    force_model::{ForceModel, PairForceModel, Virial},
    lattice,
    mpi_comm::{mpi_decomposed_pair_compute, OwnedAtoms},
    neighbor::NeighborList,
    potential::LennardJones,
    system::System,
    vec3::Vec3,
};
use mpi::collective::SystemOperation;
use mpi::traits::*;

// Initialize MPI for this single test. Because there is only one `#[test]` in
// this binary (see module docs for why), `mpi::initialize` is called exactly
// once per process — no `OnceLock` needed. The `Universe` is held locally and
// dropped at the end of the test, **after a barrier**, so both ranks reach
// `MPI_Finalize` together (no cross-rank teardown race, which otherwise makes
// OpenMPI abort with a non-zero exit even when every assertion passed).

fn lj() -> LennardJones {
    LennardJones::new(1.0, 1.0, 2.5).unwrap()
}

/// Build this rank's owned-atom view from the full system + its owned indices.
fn build_owned(sys: &System, owned_idx: &[usize]) -> OwnedAtoms {
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
    OwnedAtoms {
        system: local,
        gids,
    }
}

#[test]
fn mpi_decomposition_matches_serial_and_conserves_nve() {
    let universe = mpi::initialize().expect("MPI initialization failed");
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

    let sys = lattice::fcc(1.5, [5, 5, 5], "Ar"); // 500 atoms, box 7.5
                                                  // 2×1×1 grid: subdomain edge 7.5/2 = 3.75 ≥ cutoff 2.5.
    let grid = RankGrid { nr: [2, 1, 1] };
    let decomp = DomainDecomposition::new(sys.box_size, 2.5, grid).unwrap();
    let owned_idx = decomp.owned(&sys, rank);

    // ───────────────────────────────────────────────────────────────────────
    // Part 1: single force/energy/virial evaluation == serial reference.
    // ───────────────────────────────────────────────────────────────────────
    let owned = build_owned(&sys, &owned_idx);
    let (owned_forces, pe, w) = mpi_decomposed_pair_compute(&world, &decomp, &lj(), &owned);

    // Pack [gid, fx, fy, fz] per owned atom; gather to rank 0.
    let n_owned = owned_idx.len();
    let mut pack = Vec::with_capacity(4 * n_owned);
    for i in 0..n_owned {
        pack.push(owned_idx[i] as f64);
        pack.push(owned_forces[3 * i]);
        pack.push(owned_forces[3 * i + 1]);
        pack.push(owned_forces[3 * i + 2]);
    }
    let mut assembled: Vec<f64> = vec![0.0; 3 * sys.n()];
    if rank == 0 {
        for i in 0..n_owned {
            let gid = owned_idx[i] as usize;
            assembled[3 * gid] += owned_forces[3 * i];
            assembled[3 * gid + 1] += owned_forces[3 * i + 1];
            assembled[3 * gid + 2] += owned_forces[3 * i + 2];
        }
        for r in 1..size {
            let (recv, _): (Vec<f64>, _) = world.process_at_rank(r as mpi::Rank).receive_vec();
            for chunk in recv.as_chunks::<4>().0 {
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
    // Keep `pack` alive until the blocking send completes.
    drop(pack);

    // ───────────────────────────────────────────────────────────────────────
    // Part 2: decomposed NVE trajectory — energy conserved across ranks.
    // ───────────────────────────────────────────────────────────────────────
    // Every rank holds the full system (positions + velocities) so it can
    // integrate its own atoms and seed velocities identically; only the force
    // evaluation is decomposed. Migration-free: with dt=0.004 and 300 steps the
    // max displacement (~0.36) is far below the subdomain edge (3.75), so no
    // atom crosses a boundary. Cross-rank migration is covered by the in-process
    // `ownership_migration_after_displacement` test.
    let mut sys = sys;
    let mut p_com = Vec3::zero();
    let mut total_m = 0.0;
    for i in 0..sys.n() {
        let v = Vec3::new(
            (((i as f64 * 12.9898).sin() * 43758.5453).fract() - 0.5) * 0.3,
            (((i as f64 * 78.233).sin() * 12543.123).fract() - 0.5) * 0.3,
            (((i as f64 * 39.345).sin() * 9981.777).fract() - 0.5) * 0.3,
        );
        sys.set_vel(i, v);
        p_com += v * sys.mass(i);
        total_m += sys.mass(i);
    }
    for i in 0..sys.n() {
        let v = sys.vel(i) - p_com * (1.0 / total_m);
        sys.set_vel(i, v);
    }

    // Initial force eval → initial total energy.
    let owned = build_owned(&sys, &owned_idx);
    let (owned_forces, pe, _) = mpi_decomposed_pair_compute(&world, &decomp, &lj(), &owned);
    sys.zero_forces();
    for (i, &gid) in owned_idx.iter().enumerate() {
        let f = Vec3::new(
            owned_forces[3 * i],
            owned_forces[3 * i + 1],
            owned_forces[3 * i + 2],
        );
        sys.add_force(gid, f);
    }
    // KE: each rank sums its owned KE, then all-reduce.
    let ke_local: f64 = owned_idx
        .iter()
        .map(|&i| 0.5 * sys.mass(i) * sys.vel(i).norm2())
        .sum();
    let mut ke_global = 0.0f64;
    world.all_reduce_into(&ke_local, &mut ke_global, SystemOperation::sum());
    let e0 = ke_global + pe;

    let dt = 0.004;
    let steps = 300;
    let mut e_max = e0;
    let mut e_min = e0;
    for _ in 0..steps {
        // Half-kick + drift on OWNED atoms only.
        for &i in &owned_idx {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
            let x = sys.pos(i) + v * dt;
            sys.set_pos(i, x);
        }
        sys.wrap_positions();

        // Re-exchange ghosts with updated positions and recompute forces.
        let owned = build_owned(&sys, &owned_idx);
        let (owned_forces, pe, _) = mpi_decomposed_pair_compute(&world, &decomp, &lj(), &owned);
        sys.zero_forces();
        for (i, &gid) in owned_idx.iter().enumerate() {
            let f = Vec3::new(
                owned_forces[3 * i],
                owned_forces[3 * i + 1],
                owned_forces[3 * i + 2],
            );
            sys.add_force(gid, f);
        }

        // Half-kick on owned atoms.
        for &i in &owned_idx {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
        }

        // Total energy across ranks.
        let ke_local: f64 = owned_idx
            .iter()
            .map(|&i| 0.5 * sys.mass(i) * sys.vel(i).norm2())
            .sum();
        let mut ke_global = 0.0f64;
        world.all_reduce_into(&ke_local, &mut ke_global, SystemOperation::sum());
        let e = ke_global + pe;
        e_max = e_max.max(e);
        e_min = e_min.min(e);
    }

    let drift = (e_max - e_min) / e0.abs().max(1e-12);
    // f64 + MPI all-reduce: drift is FP summation-order noise, well below 1e-3.
    assert!(
        drift < 1e-3,
        "MPI decomposed NVE drift too large: {drift:.3e}"
    );

    // Synchronize both ranks before the universe drops, so each reaches
    // `MPI_Finalize` together — no teardown race (which OpenMPI otherwise
    // reports as an abnormal exit despite all assertions passing).
    world.barrier();
    drop(world);
    drop(universe);
}
