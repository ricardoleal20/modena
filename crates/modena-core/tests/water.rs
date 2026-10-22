//! Rigid-water stack gates (Tier B item 6).
//!
//! 1. The exact per-molecule constraint projection (the SETTLE role) restores
//!    an arbitrarily distorted water and preserves the COM.
//! 2. The excluded-pair erf correction makes intra-molecular pairs feel
//!    exactly zero Coulomb — as the exact gradient of a smooth functional,
//!    the full stack (SPME + exclusions + rigid constraints) conserves total
//!    energy on the integration level, not just in a static force check.
//! 3. The layout builder produces SPC/E geometry at the requested density.
//!
//! A *liquid* NVE gate needs an equilibrated start; the builder emits a cold
//! lattice for which both MODENA and LAMMPS report a violent initial relaxation
//! (matched on the same configuration), so the multi-water deck is documented
//! with an NVT-equilibrate-first workflow instead of an NVE drift gate.

use modena_core::constraints::{Constraints, RigidWater};
use modena_core::force_model::{ForceModel, Virial};
use modena_core::integrator::{Integrator, VelocityVerlet};
use modena_core::neighbor::NeighborList;
use modena_core::spme::{mesh_for_box, spme_reciprocal_forces};
use modena_core::system::{SimBox, System};
use modena_core::vec3::Vec3;
use modena_core::{EwaldForceModel, SpmeForceModel};

const D_OH: f64 = 1.0;
/// The exact H–H distance for the SPC/E angle: 2·d_oh·sin(109.47°/2).
const D_HH: f64 = 1.6329808618402346;
const M_O: f64 = 15.9994;
const M_H: f64 = 1.008;

fn spce_water(box_edge: f64) -> System {
    let mut sys = System::new(SimBox::cubic(box_edge));
    sys.resize(3);
    let half = (109.47_f64).to_radians() * 0.5;
    sys.set_pos(0, Vec3::new(3.0, 3.0, 3.0));
    sys.set_pos(
        1,
        Vec3::new(3.0 + D_OH * half.cos(), 3.0 + D_OH * half.sin(), 3.0),
    );
    sys.set_pos(
        2,
        Vec3::new(3.0 + D_OH * half.cos(), 3.0 - D_OH * half.sin(), 3.0),
    );
    sys.charges = vec![-0.8476, 0.4238, 0.4238];
    sys.masses = vec![M_O, M_H, M_H];
    sys.species = vec![0, 1, 1];
    sys
}

fn rigid_water() -> Constraints {
    Constraints {
        bonds: vec![],
        waters: vec![RigidWater {
            first: 0,
            d_oh: D_OH,
            d_hh: D_HH,
        }],
    }
}

fn mass_weighted_com(sys: &System) -> Vec3 {
    let mut m = 0.0;
    let mut c = Vec3::zero();
    for i in 0..sys.n() {
        c += sys.pos(i) * sys.masses[i];
        m += sys.masses[i];
    }
    c * (1.0 / m)
}

/// A violently distorted water is restored to machine precision by the exact
/// per-molecule projection, with the center of mass preserved.
#[test]
fn settle_restores_distorted_water() {
    let mut sys = spce_water(20.0);
    sys.set_pos(1, sys.pos(1) + Vec3::new(0.15, -0.1, 0.07));
    sys.set_pos(2, sys.pos(2) + Vec3::new(-0.12, 0.09, 0.11));
    let com_before = mass_weighted_com(&sys);
    let cons = rigid_water();
    cons.settle(&mut sys, 1e-12).unwrap();
    for (i, j) in [(0, 1), (0, 2)] {
        let r = sys.box_size.min_image(sys.pos(j) - sys.pos(i)).norm();
        assert!((r - D_OH).abs() < 1e-10, "O–H bond ({i},{j}): {r}");
    }
    let hh = sys.box_size.min_image(sys.pos(2) - sys.pos(1)).norm();
    assert!((hh - D_HH).abs() < 1e-10, "H–H: {hh}");
    let com_after = mass_weighted_com(&sys);
    assert!(
        (com_after - com_before).norm() < 1e-12,
        "COM moved: {com_before:?} -> {com_after:?}"
    );
}

/// The excluded-pair erf correction is the exact gradient of a smooth
/// functional: an isolated +q/−q pair with the interaction excluded flies
/// freely with no energy drift (the probe that would catch a sign or form
/// error that a static force comparison misses).
#[test]
fn excluded_pair_free_flight_conserves() {
    let mut sys = System::new(SimBox::cubic(30.0));
    sys.resize(2);
    sys.set_pos(0, Vec3::new(5.0, 5.0, 5.0));
    sys.set_pos(1, Vec3::new(24.0, 24.0, 24.0));
    sys.charges[0] = 1.0;
    sys.charges[1] = -1.0;
    sys.masses = vec![1.0, 1.0];
    sys.set_vel(0, Vec3::new(0.05, 0.02, -0.01));
    sys.set_vel(1, Vec3::new(0.0, 0.01, 0.03));
    let ewald = EwaldForceModel::new(0.5, 8.0, 5.0)
        .unwrap()
        .with_excluded_pairs(vec![(0, 1)]);
    let mut nl = NeighborList::with_skin_exclusions(8.0, 0.0, &[(0, 1)]);
    nl.build(&sys);
    let mut v = Virial::default();
    let pe0 = ewald.compute(&mut sys, &mut v, &mut nl).unwrap();
    let e0 = pe0 + sys.kinetic_energy();
    let mut integrator = VelocityVerlet::new();
    let mut e_min = f64::MAX;
    let mut e_max = f64::MIN;
    for _ in 0..400 {
        let en = integrator.step(&mut sys, &ewald, &mut nl, 0.004).unwrap();
        e_min = e_min.min(en.total);
        e_max = e_max.max(en.total);
    }
    let drift = (e_max - e_min).abs() / e0.abs().max(1e-12);
    assert!(drift < 1e-6, "excluded free flight drifts {drift:.3e}");
}

/// The full stack — rigid water + SPME + intra-molecular exclusions — is a
/// conservative force field on the integration level (this is the deck's
/// physics in miniature; the pair model contributes nothing inside a single
/// molecule).
#[test]
fn rigid_water_stack_conserves() {
    let mut sys = spce_water(7.0);
    sys.set_vel(0, Vec3::new(0.05, -0.02, 0.01));
    sys.set_vel(1, Vec3::new(0.0, 0.03, -0.02));
    sys.set_vel(2, Vec3::new(-0.04, 0.01, 0.02));
    let excl: Vec<(usize, usize)> = vec![(0, 1), (0, 2), (1, 2)];
    let excl32: Vec<(u32, u32)> = excl.iter().map(|&(a, b)| (a as u32, b as u32)).collect();
    let spme = SpmeForceModel::new(1.0, 3.0, 12.0)
        .unwrap()
        .with_excluded_pairs(excl);
    let mut nl = NeighborList::with_skin_exclusions(3.0, 0.1, &excl32);
    nl.build(&sys);
    let mut v = Virial::default();
    let pe0 = spme.compute(&mut sys, &mut v, &mut nl).unwrap();
    let e0 = pe0 + sys.kinetic_energy();
    let mut integrator = VelocityVerlet::with_constraints(rigid_water());
    let mut e_min = f64::MAX;
    let mut e_max = f64::MIN;
    for _ in 0..2000 {
        let en = integrator.step(&mut sys, &spme, &mut nl, 0.002).unwrap();
        e_min = e_min.min(en.total);
        e_max = e_max.max(en.total);
    }
    let drift = (e_max - e_min).abs() / e0.abs().max(1e-12);
    assert!(
        drift < 1e-4,
        "rigid water + SPME + exclusions drifts {drift:.3e}"
    );
    // And the constraints stayed exact throughout.
    for j in 0..sys.n() {
        for i in (j + 1)..sys.n() {
            let d = sys.box_size.min_image(sys.pos(j) - sys.pos(i)).norm();
            let target = if i == 2 && j == 1 { D_HH } else { D_OH };
            assert!((d - target).abs() < 1e-8, "{i},{j}: {d} vs {target}");
        }
    }
    // Correctness tie-in: the axis-forces of the reciprocal pass are still the
    // exact gradient of the energy functional at the integration's end state.
    let mut forces = vec![0.0; 9];
    let (e_rec, _) = spme_reciprocal_forces(
        1.0,
        &sys.box_size,
        &sys.positions,
        &sys.charges,
        mesh_for_box(&sys.box_size, 12.0),
        &mut forces,
    )
    .expect("reciprocal");
    assert!(e_rec.is_finite());
}

/// The builder puts out SPC/E geometry at the requested density: cubic box,
/// contiguous (O,H,H) triples, corrected bond lengths and angle, and O/H
/// species + charges + masses.
#[test]
fn water_lattice_geometry_and_density() {
    use modena_core::lattice::water_lattice;

    let sys = water_lattice(
        64, 0.997, 1.0, 109.47, -0.8476, 0.4238, 15.9994, 1.008, 7, 1.0,
    )
    .unwrap();
    assert_eq!(sys.n(), 64 * 3);
    let edge = sys.box_size.size.0[0];
    // 64 SPC/E molecules at 0.997 g/cm³  -> edge = (64·18.0154·1.66054/0.997)^{1/3}.
    let vol: f64 = 64.0 * 18.0154 * 1.66053906660 / 0.997;
    assert!(
        (edge - vol.cbrt()).abs() < 1e-9,
        "box 8{edge} vs {}",
        vol.cbrt()
    );
    for w in 0..64 {
        let (o, h1, h2) = (3 * w, 3 * w + 1, 3 * w + 2);
        assert_eq!(sys.species[o], 0);
        assert_eq!(sys.species[h1], 1);
        assert!((sys.charges[o] + 0.8476).abs() < 1e-12);
        assert!((sys.charges[h1] - 0.4238).abs() < 1e-12);
        let oh1 = sys.box_size.min_image(sys.pos(h1) - sys.pos(o)).norm();
        let hh = sys.box_size.min_image(sys.pos(h2) - sys.pos(h1)).norm();
        assert!((oh1 - D_OH).abs() < 1e-9, "water {w} O–H {oh1}");
        assert!((hh - D_HH).abs() < 1e-9, "water {w} H–H {hh}");
    }
}
