//! MPI strong/weak scaling benchmark (M9, `mpi` feature).
//!
//! Measures wall-clock time of `mpi_decomposed_pair_compute` over a fixed
//! number of steps on an FCC Lennard-Jones system (rho ~ 0.8442, sigma = eps =
//! 1, cutoff = 2.5). Each rank builds the full system, extracts its owned
//! atoms, and calls the distributed force kernel `steps` times. Ghost exchange
//! happens every call (positions are static, so this benchmarks steady-state
//! communication + compute cost).
//!
//! Usage: `mpi_scaling <n_atoms_target> <steps>`
//!
//! ```sh
//! # Strong scaling (fixed problem size):
//! mpiexec --oversubscribe --allow-run-as-root -n 4 \
//!   target/release/examples/mpi_scaling 10000 100
//!
//! # Weak scaling (scale atoms with ranks):
//! mpiexec --oversubscribe --allow-run-as-root -n 2 \
//!   target/release/examples/mpi_scaling 2500 100
//! ```
//!
//! Rank 0 prints a single CSV line:
//! `atoms,ranks,steps,wall_s,steps_per_s,us_per_step`

#[cfg(not(feature = "mpi"))]
fn main() {
    eprintln!(
        "this example requires the `mpi` feature: cargo run --features mpi --example mpi_scaling"
    );
}

#[cfg(feature = "mpi")]
use modena_core::{
    decomp::{DomainDecomposition, RankGrid},
    lattice,
    mpi_comm::{mpi_decomposed_pair_compute, OwnedAtoms},
    potential::LennardJones,
    system::System,
};

#[cfg(feature = "mpi")]
use mpi::traits::*;

/// Build this rank's owned-atom view from the full system + its owned indices.
/// Same pattern as `tests/mpi_decomp.rs`.
#[cfg(feature = "mpi")]
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

/// Choose roughly cubic FCC replicas to reach ~`target` atoms (4 per cell).
#[cfg(feature = "mpi")]
fn fcc_replicas_for_target(target: usize) -> [u32; 3] {
    let cells = (target as f64 / 4.0).cbrt().round() as u32;
    let cells = cells.max(1);
    [cells, cells, cells]
}

#[cfg(feature = "mpi")]
fn main() {
    // ── MPI init (same pattern as tests/mpi_decomp.rs) ──────────────────────
    let universe = mpi::initialize().expect("MPI initialization failed");
    let world = universe.world();
    let size = world.size() as usize;
    let rank = world.rank() as usize;

    // ── Parse args ──────────────────────────────────────────────────────────
    let args: Vec<String> = std::env::args().collect();
    let n_target: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(10_000);
    let steps: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(100);

    // ── Build FCC LJ system (rho ~ 0.8442 → a = (4/rho)^(1/3) ~ 1.6796) ─────
    let cutoff = 2.5;
    let a = 1.6796;
    let replicas = fcc_replicas_for_target(n_target);
    let sys = lattice::fcc(a, replicas, "Ar");
    let n_atoms = sys.n();

    // ── Choose rank grid (subdomain edge >= cutoff) ─────────────────────────
    let grid = RankGrid::balanced(sys.box_size.size, cutoff, size);
    if rank == 0 {
        eprintln!(
            "system: {} atoms, box {:?}, grid {:?} ({} ranks requested, {} in grid)",
            n_atoms,
            sys.box_size.size.0,
            grid.nr,
            size,
            grid.n_ranks()
        );
    }

    if grid.n_ranks() != size {
        if rank == 0 {
            eprintln!(
                "WARNING: requested {} ranks but grid supports only {}; skipping",
                size,
                grid.n_ranks()
            );
        }
        world.barrier();
        drop(world);
        drop(universe);
        return;
    }

    let decomp = DomainDecomposition::new(sys.box_size, cutoff, grid).unwrap();
    let owned_idx = decomp.owned(&sys, rank);
    let owned = build_owned(&sys, &owned_idx);
    let lj = LennardJones::new(1.0, 1.0, cutoff).unwrap();

    // ── Warmup (untimed; first call may include one-time setup) ─────────────
    let _ = mpi_decomposed_pair_compute(&world, &decomp, &lj, &owned);

    // ── Timed loop ──────────────────────────────────────────────────────────
    world.barrier();
    let t0 = std::time::Instant::now();
    for _ in 0..steps {
        let _ = mpi_decomposed_pair_compute(&world, &decomp, &lj, &owned);
    }
    world.barrier();
    let wall = t0.elapsed().as_secs_f64();

    if rank == 0 {
        let steps_per_s = steps as f64 / wall;
        let us_per_step = wall / steps as f64 * 1.0e6;
        println!(
            "{},{},{},{:.6},{:.2},{:.2}",
            n_atoms, size, steps, wall, steps_per_s, us_per_step
        );
    }

    // ── Shutdown (barrier + explicit drop to avoid Finalize race) ───────────
    world.barrier();
    drop(world);
    drop(universe);
}
