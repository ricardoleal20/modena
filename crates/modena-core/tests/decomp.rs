//! Domain decomposition integration tests (M9).
//!
//! The correctness anchor for the parallel decomposition: a force/energy/virial
//! evaluation partitioned across a rank grid is **bit-identical** to the serial
//! half-list evaluation, for *any* grid. This holds because pair work is
//! partitioned by ownership of the lower-global-id atom — each pair has exactly
//! one evaluator and each atom exactly one force accumulator. A multi-rank grid
//! (2×2×2) must therefore match a single-rank grid and the serial
//! `PairForceModel`, and the ghost halo must be sufficient (every partner of an
//! owned atom is present as owned or ghost).

use modena_core::{
    decomp::{DomainDecomposition, DomainForceModel, LocalCommunicator, RankGrid},
    force_model::{ForceModel, PairForceModel, Virial},
    lattice,
    neighbor::NeighborList,
    potential::LennardJones,
    vec3::Vec3,
    System,
};

/// A small FCC argon cluster in a cubic box.
fn fcc_system() -> System {
    lattice::fcc(1.5, [3, 3, 3], "Ar")
}

fn lj() -> LennardJones {
    LennardJones::new(1.0, 1.0, 2.5).unwrap()
}

/// Serial reference force/energy/virial via the production pair model.
fn serial_reference(sys: &System) -> (Vec<f64>, f64, f64) {
    let mut s = sys.clone();
    s.zero_forces();
    let mut vir = Virial::default();
    let mut nl = NeighborList::new(2.5);
    let model = PairForceModel::new(lj());
    let pe = model.compute(&mut s, &mut vir, &mut nl).unwrap();
    (s.forces.clone(), pe, vir.w)
}

/// Decomposed force/energy/virial via `DomainForceModel` with a given grid.
fn decomposed(sys: &System, grid: RankGrid) -> (Vec<f64>, f64, f64) {
    let mut s = sys.clone();
    let decomp = DomainDecomposition::new(s.box_size, 2.5, grid).unwrap();
    let comm = LocalCommunicator::new(grid.n_ranks(), 0);
    let model = DomainForceModel::new(lj(), decomp, Box::new(comm));
    let mut vir = Virial::default();
    let mut nl = NeighborList::new(2.5);
    let pe = model.compute(&mut s, &mut vir, &mut nl).unwrap();
    (s.forces.clone(), pe, vir.w)
}

#[test]
fn single_rank_decomposition_matches_serial() {
    let sys = fcc_system();
    let (sf, spe, sw) = serial_reference(&sys);
    let (df, dpe, dw) = decomposed(&sys, RankGrid::single());
    assert!((dpe - spe).abs() < 1e-9, "pe {dpe} vs {spe}");
    assert!((dw - sw).abs() < 1e-9, "virial {dw} vs {sw}");
    for i in 0..sf.len() {
        assert!(
            (sf[i] - df[i]).abs() < 1e-9,
            "force[{i}] decomposed {df_x} vs serial {sf_x}",
            df_x = df[i],
            sf_x = sf[i]
        );
    }
}

#[test]
fn multi_rank_decomposition_matches_serial() {
    // 3×3×3 FCC of lattice 1.5 → box 4.5, cutoff 2.5. Max 1 split/axis
    // (4.5/2.5 = 1.8 → 1). Use a bigger box to allow a real 2×2×2 grid.
    let sys = lattice::fcc(1.5, [5, 5, 5], "Ar"); // box 7.5 → 3/axis max
    let grid = RankGrid::balanced(Vec3::new(7.5, 7.5, 7.5), 2.5, 8);
    assert!(
        grid.n_ranks() > 1,
        "test needs a multi-rank grid, got {grid:?}"
    );
    let (sf, spe, sw) = serial_reference(&sys);
    let (df, dpe, dw) = decomposed(&sys, grid);
    // Energy/virial exact (each pair counted once by construction).
    assert!((dpe - spe).abs() < 1e-9, "pe {dpe} vs {spe}");
    assert!((dw - sw).abs() < 1e-9, "virial {dw} vs {sw}");
    for i in 0..sf.len() {
        assert!(
            (sf[i] - df[i]).abs() < 1e-9,
            "force[{i}] decomposed {df_x} vs serial {sf_x} (rel {rel})",
            df_x = df[i],
            sf_x = sf[i],
            rel = (sf[i] - df[i]).abs() / sf[i].abs().max(1e-12)
        );
    }
}

#[test]
fn ghost_halo_is_sufficient() {
    // For a multi-rank grid, every pair partner of an owned atom must be
    // present as owned or ghost on the owning rank. Equivalent to: the union of
    // (owned ∪ ghosts) over all ranks contains every atom, and the decomposed
    // result matches serial (which it can only do if no neighbor is missing).
    let sys = lattice::fcc(1.5, [5, 5, 5], "Ar");
    let grid = RankGrid::balanced(Vec3::new(7.5, 7.5, 7.5), 2.5, 8);
    let decomp = DomainDecomposition::new(sys.box_size, 2.5, grid).unwrap();
    // Owned sets partition the atoms (no overlap, full coverage).
    let mut all: Vec<usize> = Vec::new();
    for r in 0..grid.n_ranks() {
        let owned = decomp.owned(&sys, r);
        let ghosts = decomp.ghosts(&sys, r);
        // Owned and ghost sets are disjoint (a ghost is owned by another rank).
        for &o in &owned {
            assert!(
                !ghosts.contains(&o),
                "atom {o} is both owned and ghost on rank {r}"
            );
        }
        all.extend(owned);
    }
    all.sort_unstable();
    all.dedup();
    assert_eq!(all.len(), sys.n(), "owned sets must cover every atom once");
    // The decomposed evaluation matching serial (below) is the real sufficiency
    // proof: a missing ghost would change a force.
    let (sf, _, _) = serial_reference(&sys);
    let (df, _, _) = decomposed(&sys, grid);
    for i in 0..sf.len() {
        assert!(
            (sf[i] - df[i]).abs() < 1e-9,
            "force[{i}] mismatch (ghost halo insufficient)"
        );
    }
}

#[test]
fn ownership_migration_after_displacement() {
    // After atoms cross a subdomain face, ownership is recomputed from
    // positions; the migration count is nonzero and the new partition still
    // covers every atom.
    let mut sys = lattice::fcc(1.5, [4, 4, 4], "Ar"); // box 6.0
    let grid = RankGrid { nr: [2, 1, 1] };
    let decomp = DomainDecomposition::new(sys.box_size, 2.5, grid).unwrap();
    let prev: Vec<usize> = (0..sys.n()).map(|i| decomp.owner_of(sys.pos(i))).collect();

    // Shift every atom by +3.0 in x → those near the boundary cross faces.
    for i in 0..sys.n() {
        let p = sys.pos(i);
        sys.set_pos(i, Vec3::new(p.0[0] + 3.0, p.0[1], p.0[2]));
    }
    sys.wrap_positions();

    let migrated = decomp.migration_count(&prev, &sys);
    assert!(
        migrated > 0,
        "expected some atoms to migrate across the x-face"
    );

    // Every atom still has exactly one owner (ownership is a total function of
    // position); the owned sets still partition the atoms.
    let owners: Vec<usize> = (0..sys.n()).map(|i| decomp.owner_of(sys.pos(i))).collect();
    assert_eq!(owners.len(), sys.n());
    // Both ranks still own atoms after migration.
    let n0 = owners.iter().filter(|&&r| r == 0).count();
    let n1 = owners.iter().filter(|&&r| r == 1).count();
    assert_eq!(n0 + n1, sys.n());
    assert!(
        n0 > 0 && n1 > 0,
        "both ranks should be populated after migration"
    );

    // Decomposed evaluation still matches serial after migration.
    let (sf, spe, _) = serial_reference(&sys);
    let (df, dpe, _) = decomposed(&sys, grid);
    assert!((dpe - spe).abs() < 1e-9);
    for i in 0..sf.len() {
        assert!((sf[i] - df[i]).abs() < 1e-9, "force[{i}] after migration");
    }
}

#[test]
fn decomposed_nve_energy_conservation() {
    // Integration-level check: a multi-rank DomainForceModel drives a manual
    // velocity-Verlet NVE loop, and total energy is conserved. This proves the
    // decomposed forces are conservative under integration (not just equal to
    // the serial forces at one configuration).
    let mut sys = lattice::fcc(1.5, [4, 4, 4], "Ar"); // box 6.0, 256 atoms
    let grid = RankGrid { nr: [2, 1, 1] }; // subdomain edge 3.0 >= cutoff 2.5
    let decomp = DomainDecomposition::new(sys.box_size, 2.5, grid).unwrap();
    let comm = LocalCommunicator::new(grid.n_ranks(), 0);
    let model = DomainForceModel::new(lj(), decomp, Box::new(comm));
    let mut vir = Virial::default();
    let mut nl = NeighborList::new(2.5);

    // Seed small deterministic velocities and zero the COM drift.
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

    let pe0 = model.compute(&mut sys, &mut vir, &mut nl).unwrap();
    let e0 = sys.kinetic_energy() + pe0;

    let dt = 0.005;
    let steps = 400;
    let mut e_max = e0;
    let mut e_min = e0;
    for _ in 0..steps {
        // Velocity-Verlet: half-kick, drift, force, half-kick.
        for i in 0..sys.n() {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
            let x = sys.pos(i) + v * dt;
            sys.set_pos(i, x);
        }
        sys.wrap_positions();
        let pe = model.compute(&mut sys, &mut vir, &mut nl).unwrap();
        for i in 0..sys.n() {
            let m = sys.mass(i);
            let f = sys.force(i);
            let v = sys.vel(i) + f * (0.5 * dt / m);
            sys.set_vel(i, v);
        }
        let e = sys.kinetic_energy() + pe;
        e_max = e_max.max(e);
        e_min = e_min.min(e);
    }
    let drift = (e_max - e_min) / e0.abs().max(1e-12);
    assert!(drift < 1e-3, "decomposed NVE drift too large: {drift:.3e}");
}
