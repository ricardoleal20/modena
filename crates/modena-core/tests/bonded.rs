//! Conservation contract for bonded interactions (Feature B, M3).
//!
//! For every bonded term — bond, angle, proper dihedral, improper — the reported
//! force must equal −∇U by central finite difference, and the per-term forces
//! must sum to zero (Newton's third law / translation invariance). The dihedral
//! gradient is the error-prone piece (a hand-derived exact analytic gradient of
//! the dihedral angle); this test is what makes it defensible.
//!
//! We also run a small all-bonded molecule NVE to confirm the composite
//! (pair + bonded) accumulate contract conserves energy end-to-end.

use modena_core::{
    Bond, BondedForceModel, BondedTopology, CompositeForceModel, Constraint, Constraints, Dihedral,
    DihedralTerm, ForceModel, Integrator, LennardJones, NeighborList, PairForceModel, SimBox,
    System, Vec3, VelocityVerlet, Virial,
};
use proptest::prelude::*;

/// A 5-atom "molecule": a bent chain 0-1-2-3-4 with a bond, an angle, and a
/// proper dihedral spanning 0-1-2-3, plus an improper on 1-2-3-4. Positions are
/// non-degenerate (no three collinear, no two coincident).
fn molecule() -> (System, BondedTopology) {
    let mut sys = System::new(SimBox::cubic(100.0));
    sys.resize(5);
    let pts = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.1, 0.2, 0.0),
        Vec3::new(2.0, 0.9, 0.3),
        Vec3::new(2.6, 1.8, -0.2),
        Vec3::new(3.5, 1.5, 0.5),
    ];
    for (i, p) in pts.iter().enumerate() {
        sys.set_pos(i, *p);
    }
    let mut topo = BondedTopology::new();
    topo.bonds.push(Bond {
        i: 0,
        j: 1,
        r0: 1.0,
        k: 200.0,
    });
    topo.bonds.push(Bond {
        i: 1,
        j: 2,
        r0: 1.2,
        k: 180.0,
    });
    topo.bonds.push(Bond {
        i: 2,
        j: 3,
        r0: 1.1,
        k: 220.0,
    });
    topo.bonds.push(Bond {
        i: 3,
        j: 4,
        r0: 1.0,
        k: 160.0,
    });
    topo.angles.push(modena_core::Angle {
        i: 0,
        j: 1,
        k: 2,
        theta0: 1.9,
        ktheta: 50.0,
    });
    topo.angles.push(modena_core::Angle {
        i: 1,
        j: 2,
        k: 3,
        theta0: 1.7,
        ktheta: 60.0,
    });
    topo.angles.push(modena_core::Angle {
        i: 2,
        j: 3,
        k: 4,
        theta0: 2.0,
        ktheta: 45.0,
    });
    topo.dihedrals.push(Dihedral {
        i: 0,
        j: 1,
        k: 2,
        l: 3,
        terms: vec![
            DihedralTerm {
                n: 1,
                k: 2.0,
                delta: 0.0,
            },
            DihedralTerm {
                n: 2,
                k: 1.5,
                delta: std::f64::consts::PI,
            },
        ],
    });
    topo.dihedrals.push(Dihedral {
        i: 1,
        j: 2,
        k: 3,
        l: 4,
        terms: vec![DihedralTerm {
            n: 2,
            k: 3.0,
            delta: 0.0,
        }],
    });
    (sys, topo)
}

fn bonded_pe(sys: &System, topo: &BondedTopology) -> f64 {
    let model = BondedForceModel::new(topo.clone());
    let mut nl = NeighborList::new(2.5);
    let mut v = Virial::default();
    let mut s = sys.clone();
    s.forces.fill(0.0);
    model.compute(&mut s, &mut v, &mut nl).unwrap()
}

fn bonded_forces(sys: &System, topo: &BondedTopology) -> Vec<Vec3> {
    let model = BondedForceModel::new(topo.clone());
    let mut nl = NeighborList::new(2.5);
    let mut v = Virial::default();
    let mut s = sys.clone();
    s.forces.fill(0.0);
    model.compute(&mut s, &mut v, &mut nl).unwrap();
    (0..s.n()).map(|i| s.force(i)).collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    /// Per-atom bonded force == −∇U, by central finite difference on total PE.
    #[test]
    fn bonded_force_is_negative_gradient_of_energy(
        jx in any::<f64>(), jy in any::<f64>(), jz in any::<f64>(),
    ) {
        let (mut sys, topo) = molecule();
        // Jitter every atom by a small random amount so the configuration is
        // generic (avoids measure-zero coincidences like θ = θ0 or φ at an extremum).
        let j = 0.08;
        for i in 0..sys.n() {
            let seed = (i as f64) * 0.37;
            sys.add_pos(i, Vec3::new(
                (jx * seed).rem_euclid(1.0) * j,
                (jy * seed).rem_euclid(1.0) * j,
                (jz * seed).rem_euclid(1.0) * j,
            ));
        }

        let forces = bonded_forces(&sys, &topo);
        let h = 1e-6_f64;
        for (i, force) in forces.iter().enumerate() {
            for k in 0..3 {
                let mut sp = sys.clone();
                let mut sm = sys.clone();
                sp.positions[3 * i + k] += h;
                sm.positions[3 * i + k] -= h;
                let pe_p = bonded_pe(&sp, &topo);
                let pe_m = bonded_pe(&sm, &topo);
                let fd_force = -(pe_p - pe_m) / (2.0 * h);
                let err = (force.0[k] - fd_force).abs();
                let scale = 1.0 + force.0[k].abs();
                prop_assert!(
                    err < 1e-3 * scale,
                    "atom {} axis {}: force={}, fd={}, err={}",
                    i, k, force.0[k], fd_force, err
                );
            }
        }
    }
}

#[test]
fn bonded_forces_sum_to_zero() {
    let (sys, topo) = molecule();
    let forces = bonded_forces(&sys, &topo);
    let net: Vec3 = forces.iter().copied().fold(Vec3::zero(), |a, b| a + b);
    assert!(
        net.norm() < 1e-9,
        "bonded forces must sum to zero (translation invariance), got {net:?}"
    );
}

/// A 4-atom molecule at a bonded equilibrium whose dihedral is **non-planar**
/// (φ = 90°). This is the key to a stable NVE test: the dihedral gradient
/// ∂φ/∂r ∝ 1/sinφ diverges at the coplanar geometries (φ = 0 or π), so a
/// planar test molecule (trans/cis) sits on a singularity and energy appears to
/// drift. By using a periodic term with n = 2, δ = 0 — whose minimum is at
/// φ = π/2 (gauche, 90°, far from coplanar) — the molecule rests at a smooth
/// minimum and velocity-Verlet conserves energy.
///
/// Geometry: bond 1–2 along z; atom 0 along x; atom 3 along y. All bonds 1.0,
/// all angles 90°, dihedral 0–1–2–3 = 90° = π/2.
fn equilibrium_molecule() -> (System, BondedTopology) {
    // Centered at (50, 50, 50) so the molecule is far from every periodic box
    // boundary — small vibrations never wrap (wrapping a 4-atom molecule in a
    // 100³ box breaks the cell-list neighbor search for its own LJ pairs).
    let c = Vec3::new(50.0, 50.0, 50.0);
    let mut sys = System::new(SimBox::cubic(100.0));
    sys.resize(4);
    let pts = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 1.0),
    ];
    for (i, p) in pts.iter().enumerate() {
        sys.set_pos(i, *p + c);
    }
    let mut topo = BondedTopology::new();
    for (a, b) in [(0, 1), (1, 2), (2, 3)] {
        topo.bonds.push(Bond {
            i: a,
            j: b,
            r0: 1.0,
            k: 50.0,
        });
    }
    for (a, b, c) in [(0, 1, 2), (1, 2, 3)] {
        topo.angles.push(modena_core::Angle {
            i: a,
            j: b,
            k: c,
            theta0: std::f64::consts::FRAC_PI_2, // 90°
            ktheta: 20.0,
        });
    }
    // n = 2, δ = 0 → U = k(1 + cos 2φ); minimum at φ = π/2 (the start).
    topo.dihedrals.push(Dihedral {
        i: 0,
        j: 1,
        k: 2,
        l: 3,
        terms: vec![DihedralTerm {
            n: 2,
            k: 2.0,
            delta: 0.0,
        }],
    });
    (sys, topo)
}

#[test]
fn composite_decomposes_into_pair_plus_bonded() {
    // The accumulate contract: the composite's forces, PE, and virial must equal
    // the sum of the pair-only and bonded-only contributions — no double-count,
    // no drop.
    let (sys, topo) = equilibrium_molecule();

    // Pair-only.
    let mut s_pair = sys.clone();
    s_pair.forces.fill(0.0);
    let mut nl = NeighborList::new(4.0);
    let mut vp = Virial::default();
    let pe_pair = PairForceModel::new(LennardJones::new(0.02, 1.0, 4.0).unwrap())
        .compute(&mut s_pair, &mut vp, &mut nl)
        .unwrap();

    // Bonded-only.
    let mut s_bond = sys.clone();
    s_bond.forces.fill(0.0);
    let mut nl2 = NeighborList::new(4.0);
    let mut vb = Virial::default();
    let pe_bond = BondedForceModel::new(topo.clone())
        .compute(&mut s_bond, &mut vb, &mut nl2)
        .unwrap();

    // Composite.
    let mut s_comp = sys.clone();
    s_comp.forces.fill(0.0);
    let mut nl3 = NeighborList::new(4.0);
    let mut vc = Virial::default();
    let pe_comp = CompositeForceModel::new(PairForceModel::new(
        LennardJones::new(0.02, 1.0, 4.0).unwrap(),
    ))
    .with_bonded(BondedForceModel::new(topo.clone()))
    .compute(&mut s_comp, &mut vc, &mut nl3)
    .unwrap();

    assert!(
        (pe_comp - (pe_pair + pe_bond)).abs() < 1e-9,
        "composite PE {pe_comp} != pair {pe_pair} + bonded {pe_bond}"
    );
    assert!(
        (vc.w - (vp.w + vb.w)).abs() < 1e-9,
        "composite virial {} != pair {} + bonded {}",
        vc.w,
        vp.w,
        vb.w
    );
    for i in 0..sys.n() {
        let sum = s_pair.force(i) + s_bond.force(i);
        let diff = (sum - s_comp.force(i)).norm();
        assert!(
            diff < 1e-9,
            "atom {i}: composite force != pair+bonded (diff {diff:.3e})"
        );
    }
}

#[test]
fn composite_pair_plus_bonded_conserves_energy() {
    // The molecule starts at rest (zero velocity ⇒ zero COM and angular
    // momentum) with a small position perturbation, so it vibrates in place and
    // never translates across the periodic boundary (which would break the
    // cell-list neighbor search for a tiny molecule in a huge box). Velocity-
    // Verlet then conserves energy for the pair+bonded composite.
    let (mut sys, topo) = equilibrium_molecule();
    // Excite a vibration: displace atom 0 and release from rest.
    sys.add_pos(0, Vec3::new(0.02, 0.0, 0.0));
    // Zero velocities.
    for i in 0..sys.n() {
        sys.set_vel(i, Vec3::zero());
    }

    let lj = LennardJones::new(0.02, 1.0, 4.0).unwrap();
    let pair = PairForceModel::new(lj);
    let bonded = BondedForceModel::new(topo);
    let model = CompositeForceModel::new(pair).with_bonded(bonded);

    let mut nl = NeighborList::new(4.0);
    let mut integ = VelocityVerlet::new();
    let mut v = Virial::default();
    sys.forces.fill(0.0);
    let _ = model.compute(&mut sys, &mut v, &mut nl);

    let dt = 0.0005;
    let mut e0 = None;
    let mut e_last = 0.0;
    let mut e_max: f64 = 0.0;
    for _ in 0..8000 {
        let e = integ.step(&mut sys, &model, &mut nl, dt).unwrap();
        if e0.is_none() {
            e0 = Some(e.total);
        }
        e_last = e.total;
        e_max = e_max.max(e.total);
    }
    let e0 = e0.unwrap();
    let drift = ((e_last - e0).abs().max((e_max - e0).abs())) / e0.abs().max(1e-12_f64);
    assert!(
        drift < 1e-3,
        "composite NVE drift {drift:.3e} exceeds 1e-3 (E0={e0:.6}, E_last={e_last:.6}, E_max={e_max:.6})"
    );
}

// ---------------------------------------------------------------------------
// SHAKE/RATTLE constraints (M3).
// ---------------------------------------------------------------------------

/// A free rigid dumbbell: two atoms tied by a single distance constraint, no
/// forces. Given an initial angular velocity (perpendicular to the bond), it
/// must rotate as a rigid rotor — bond length fixed, kinetic energy conserved.
/// This is the cleanest SHAKE+RATTLE check: nothing else is happening.
#[test]
fn shake_rattle_rigid_dumbbell_conserves_length_and_energy() {
    let mut sys = System::new(SimBox::cubic(100.0));
    sys.resize(2);
    // Centered, bond along x.
    let c = Vec3::new(50.0, 50.0, 50.0);
    sys.set_pos(0, c + Vec3::new(-0.5, 0.0, 0.0));
    sys.set_pos(1, c + Vec3::new(0.5, 0.0, 0.0));
    // Angular velocity about z: tangential (perpendicular to the bond).
    sys.set_vel(0, Vec3::new(0.0, 0.5, 0.0));
    sys.set_vel(1, Vec3::new(0.0, -0.5, 0.0));

    let cons = Constraints {
        bonds: vec![Constraint { i: 0, j: 1, d: 1.0 }],
    };
    // No forces: an empty bonded model contributes zero PE and zero force.
    let model = BondedForceModel::new(BondedTopology::new());
    let mut nl = NeighborList::new(2.5);
    let mut integ = VelocityVerlet::with_constraints(cons.clone());
    let mut v = Virial::default();
    sys.forces.fill(0.0);
    let _ = model.compute(&mut sys, &mut v, &mut nl);

    let dt = 0.001;
    let ke0 = sys.kinetic_energy();
    let mut max_len_err = 0.0_f64;
    let mut ke_drift = 0.0_f64;
    for _ in 0..5000 {
        let e = integ.step(&mut sys, &model, &mut nl, dt).unwrap();
        max_len_err = max_len_err.max(cons.max_residual(&sys));
        ke_drift = ke_drift.max((e.kinetic - ke0).abs());
    }
    assert!(
        max_len_err < 1e-8,
        "SHAKE bond-length residual {max_len_err:.3e} too large"
    );
    assert!(
        ke_drift / ke0 < 1e-4,
        "rigid-rotor KE drift {ke_drift:.3e} too large (ke0={ke0:.6})"
    );
}

/// A constrained molecule: the equilibrium 4-atom chain with its three bonds
/// held rigid by SHAKE/RATTLE, but the angle and dihedral terms still active
/// (so there is real bonded dynamics on the constraint manifold). Verifies bond
/// lengths stay fixed *and* total energy is conserved — the symplectic RATTLE
/// property.
#[test]
fn constrained_molecule_holds_bonds_and_conserves_energy() {
    let (mut sys, mut topo) = equilibrium_molecule();
    // Drop the harmonic bond terms — the constrained bonds are rigid, so a
    // harmonic term on the same pair would be redundant (and fight the
    // constraint). Keep the angles + dihedral as the only forces.
    topo.bonds.clear();
    let cons = Constraints {
        bonds: vec![
            Constraint { i: 0, j: 1, d: 1.0 },
            Constraint { i: 1, j: 2, d: 1.0 },
            Constraint { i: 2, j: 3, d: 1.0 },
        ],
    };
    // Excite angle/dihedral motion: a tangential kick on atom 0 (perpendicular
    // to bond 0-1, which is along x), zero COM velocity.
    sys.set_vel(0, Vec3::new(0.0, 0.3, 0.0));
    sys.set_vel(1, Vec3::new(0.0, -0.1, 0.0));
    sys.set_vel(2, Vec3::new(0.0, -0.1, 0.0));
    sys.set_vel(3, Vec3::new(0.0, -0.1, 0.0));

    let lj = LennardJones::new(0.02, 1.0, 4.0).unwrap();
    let pair = PairForceModel::new(lj);
    let bonded = BondedForceModel::new(topo);
    let model = CompositeForceModel::new(pair).with_bonded(bonded);

    let mut nl = NeighborList::new(4.0);
    let mut integ = VelocityVerlet::with_constraints(cons.clone());
    let mut v = Virial::default();
    sys.forces.fill(0.0);
    let _ = model.compute(&mut sys, &mut v, &mut nl);

    let dt = 0.0005;
    let mut e0 = None;
    let mut e_last = 0.0;
    let mut e_max = 0.0_f64;
    let mut max_len_err = 0.0_f64;
    for _ in 0..8000 {
        let e = integ.step(&mut sys, &model, &mut nl, dt).unwrap();
        if e0.is_none() {
            e0 = Some(e.total);
        }
        e_last = e.total;
        e_max = e_max.max(e.total);
        max_len_err = max_len_err.max(cons.max_residual(&sys));
    }
    let e0 = e0.unwrap();
    let drift = ((e_last - e0).abs().max((e_max - e0).abs())) / e0.abs().max(1e-12_f64);
    assert!(
        max_len_err < 1e-8,
        "constrained bond-length residual {max_len_err:.3e} too large"
    );
    assert!(
        drift < 1e-3,
        "constrained NVE drift {drift:.3e} exceeds 1e-3 (E0={e0:.6}, E_last={e_last:.6}, E_max={e_max:.6})"
    );
}
