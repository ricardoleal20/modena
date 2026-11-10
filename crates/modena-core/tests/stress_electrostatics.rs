#[allow(dead_code)]
#[path = "support/stress.rs"]
mod support;
use modena_core::neighbor::PairScale;
use modena_core::{EwaldForceModel, ForceModel, NeighborList, SimBox, System, Vec3, Virial};
use support::{affine, complete, near};

const EPS: [f64; 4] = [1e-4, 3e-5, 1e-5, 3e-6];
const RC: f64 = 3.5;
const KC: f64 = 8.;
fn fixture(triclinic: bool, neutral: bool, inside: bool, outside: bool) -> System {
    let cell = if triclinic {
        SimBox::triclinic([
            Vec3::new(8., 0., 0.),
            Vec3::new(0.8, 9., 0.),
            Vec3::new(0.3, 0.1, 10.),
        ])
    } else {
        SimBox::from_array([8., 9., 10.])
    };
    let mut s = System::new(cell);
    let positions: &[[f64; 3]] = if outside {
        &[[0.1, 0.1, 0.1], [0.59, 0.4, 0.3]]
    } else if inside {
        &[[0.13, 0.27, 0.39], [0.31, 0.34, 0.48], [0.32, 0.51, 0.58]]
    } else {
        &[[0.13, 0.27, 0.39], [0.61, 0.44, 0.73], [0.32, 0.81, 0.58]]
    };
    s.resize(positions.len());
    let h = cell.basis();
    for (i, p) in positions.iter().enumerate() {
        s.set_pos(i, h[0] * p[0] + h[1] * p[1] + h[2] * p[2]);
    }
    s.charges = if outside {
        vec![1., -0.7]
    } else if neutral {
        vec![1., -0.7, -0.3]
    } else {
        vec![1., -0.7, 0.2]
    };
    s
}
fn model(scaled: bool, outside: bool) -> EwaldForceModel {
    let m = EwaldForceModel::new(0.8, RC, KC).unwrap();
    if scaled {
        m.with_pair_scales(if outside {
            vec![PairScale {
                i: 0,
                j: 1,
                lj: 1.,
                coulomb: 0.,
            }]
        } else {
            vec![
                PairScale {
                    i: 0,
                    j: 1,
                    lj: 1.,
                    coulomb: 0.4,
                },
                PairScale {
                    i: 1,
                    j: 2,
                    lj: 1.,
                    coulomb: 0.,
                },
            ]
        })
        .unwrap()
    } else {
        m
    }
}
// Enumerate a deliberately generous fixed integer cube, independently of the
// production adaptive bounds; exact equality proves no cutoff-mode crossings.
fn modes(s: &System) -> Vec<[i32; 3]> {
    let inv = s.box_size.inv_basis_rows();
    for h in s.box_size.basis() {
        assert!(h.norm() * KC / (2. * std::f64::consts::PI) < 16.);
    }
    let mut modes = Vec::new();
    for i in -16..=16 {
        for j in -16..=16 {
            for k in -16..=16 {
                if i > 0 || (i == 0 && j > 0) || (i == 0 && j == 0 && k > 0) {
                    let v = (inv[0] * f64::from(i) + inv[1] * f64::from(j) + inv[2] * f64::from(k))
                        * (2. * std::f64::consts::PI);
                    if v.norm2() < KC * KC {
                        modes.push([i, j, k]);
                    }
                }
            }
        }
    }
    modes
}
fn branches(s: &System) -> Vec<(bool, Vec3)> {
    let mut out = Vec::new();
    for i in 0..s.n() {
        for j in i + 1..s.n() {
            let raw = s.pos(i) - s.pos(j);
            let d = s.box_size.min_image(raw);
            let shift = raw - d;
            let inv = s.box_size.inv_basis_rows();
            out.push((
                d.norm() < RC,
                Vec3::new(
                    inv[0].dot(shift).round(),
                    inv[1].dot(shift).round(),
                    inv[2].dot(shift).round(),
                ),
            ));
        }
    }
    out
}
#[test]
fn ewald_all_six_affine_derivatives_fixed_modes_and_smooth_branches() {
    rejected_half_cell_outside_fixture_has_image_tie();
    let fixtures = [
        (false, true, false, false, false),
        (false, false, false, false, true),
        (false, false, true, false, true),
        (false, false, false, true, true),
        (true, true, false, false, false),
        (true, false, false, false, true),
        (true, false, true, false, true),
        (true, false, false, true, true),
    ];
    for &(tri, neutral, inside, outside, _) in &fixtures {
        let source = fixture(tri, neutral, inside, outside);
        let expected_modes = modes(&source);
        let expected_branches = branches(&source);
        // Prove all branches/mode sets for all four steps BEFORE derivatives.
        for c in 0..6 {
            for eps in EPS {
                for sign in [-1., 1.] {
                    let s = affine(&source, c, sign * eps);
                    assert_eq!(modes(&s), expected_modes, "reciprocal membership");
                    assert_eq!(branches(&s), expected_branches, "image/cutoff branch");
                }
            }
        }
        println!("branch-proof tri={tri} neutral={neutral} inside={inside} outside={outside} modes={} branches={expected_branches:?}", expected_modes.len());
    }
    for (tri, neutral, inside, outside, scaled) in fixtures {
        let source = fixture(tri, neutral, inside, outside);
        let m = model(scaled, outside);
        let mut s = source.clone();
        s.forces.fill(0.25);
        let mut nl = NeighborList::new(RC);
        let e = m.evaluate(&mut s, &mut nl).unwrap();
        let tensor = complete(&e.tensor);
        near(tensor.trace(), e.scalar_virial, 5e-12);
        let mut warm = source.clone();
        warm.forces.fill(0.25);
        let warm_e = m.evaluate(&mut warm, &mut nl).unwrap();
        assert_eq!(e, warm_e);
        assert_eq!(s.forces, warm.forces);
        let mut legacy = source.clone();
        legacy.forces.fill(0.25);
        let mut w = Virial { w: 7. };
        let energy = m.compute(&mut legacy, &mut w, &mut nl).unwrap();
        near(energy, e.energy, 5e-12);
        near(w.w - 7., e.scalar_virial, 5e-12);
        assert_eq!(legacy.forces, s.forces);
        near(m.energy(&source, &mut nl).unwrap(), e.energy, 5e-12); // empty correction force slice
        for (c, &expected) in tensor.0.iter().enumerate() {
            let mut values = Vec::new();
            for eps in EPS {
                let up = m
                    .energy(&affine(&source, c, eps), &mut NeighborList::new(RC))
                    .unwrap();
                let down = m
                    .energy(&affine(&source, c, -eps), &mut NeighborList::new(RC))
                    .unwrap();
                let fd = -(up - down) / (2. * eps);
                values.push(fd);
                println!("tri={tri} neutral={neutral} inside={inside} outside={outside} scaled={scaled} component={c} eps={eps:.8e} analytic={expected:.16e} fd={fd:.16e} error={:.8e}",(fd-expected).abs());
            }
            let tol = 2e-7 * (1. + expected.abs());
            assert!(
                values.windows(2).any(|v| (v[0] - expected).abs() <= tol
                    && (v[1] - expected).abs() <= tol
                    && (v[0] - v[1]).abs() <= 2. * tol),
                "component {c} {values:?} expected {expected}"
            );
        }
    }
}
#[test]
fn net_charge_background_is_isotropic_and_self_has_zero_cell_derivative() {
    let mut s = System::new(SimBox::from_array([8., 9., 10.]));
    s.resize(1);
    s.charges[0] = 0.5;
    // No real pairs or reciprocal modes: only self and background remain.
    let m = EwaldForceModel::new(0.8, RC, 0.01).unwrap();
    let e = m.evaluate(&mut s, &mut NeighborList::new(RC)).unwrap();
    let t = complete(&e.tensor);
    let bg = -std::f64::consts::PI * 0.5_f64.powi(2) / (2. * 720. * 0.8_f64.powi(2));
    for (a, b) in t.0.into_iter().zip([bg, bg, bg, 0., 0., 0.]) {
        near(a, b, 2e-14);
    }
    near(e.scalar_virial, 3. * bg, 2e-14);
    let self_energy = -0.8 / std::f64::consts::PI.sqrt() * 0.25;
    near(e.energy, self_energy + bg, 2e-14);
}

/// Rejected spec fixture: a half-cell separation is an image tie and cannot
/// certify a smooth Hamiltonian derivative. Controller ruling changes .6 to .59.
fn rejected_half_cell_outside_fixture_has_image_tie() {
    let mut s = fixture(false, false, false, true);
    s.set_pos(1, Vec3::new(4.8, 3.6, 3.));
    let d = s.pos(0) - s.pos(1);
    near(d.norm2(), (d + Vec3::new(8., 0., 0.)).norm2(), 2e-14);
    assert_ne!(
        branches(&affine(&s, 3, 1e-4)),
        branches(&affine(&s, 3, -1e-4))
    );
}
