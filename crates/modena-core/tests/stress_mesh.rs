//! Independent mesh-Hamiltonian derivatives: no production oracle helpers.
#[path = "../examples/stress_validation/mod.rs"]
mod stress_validation;
use modena_core::neighbor::PairScale;
use modena_core::spme::{
    mesh_for_box, spme_reciprocal_energy, spme_reciprocal_forces, spme_reciprocal_tensor,
};
use modena_core::stress::{SymTensor6, TensorAvailability};
use modena_core::{ForceModel, NeighborList, SimBox, SpmeForceModel, System, Vec3, Virial};
use std::f64::consts::PI;
use stress_validation::mesh_oracle::{
    energy_strains, grid_energy, grid_potential, mesh_energy, strain_cell, uniform_mesh_energy,
};

const CELL: [[f64; 3]; 3] = [[8., 0., 0.], [0., 9., 0.], [0., 0., 10.]];
const FRACTIONAL: [[f64; 3]; 3] = [[0.13, 0.27, 0.39], [0.61, 0.44, 0.73], [0.32, 0.81, 0.58]];
const CHARGES: [f64; 3] = [1., -0.7, 0.2];
const MESHES: [[usize; 3]; 4] = [[5, 5, 5], [5, 9, 5], [4, 6, 8], [6, 6, 6]];
const EPS: [f64; 4] = [1e-4, 3e-5, 1e-5, 3e-6];
fn near(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t * (1. + b.abs()),
        "{a:.16e} != {b:.16e} tolerance {t}"
    );
}
fn positions(s: &[[f64; 3]], h: [[f64; 3]; 3]) -> Vec<f64> {
    s.iter()
        .flat_map(|s| (0..3).map(move |a| (0..3).map(|b| h[b][a] * s[b]).sum::<f64>()))
        .collect()
}
fn adjacent(rows: &[[f64; 6]; 4], expected: [f64; 6], tol: f64) {
    for (c, e) in expected.into_iter().enumerate() {
        let t = tol * (1. + e.abs());
        assert!(
            rows.windows(2).any(|r| (r[0][c] - e).abs() <= t
                && (r[1][c] - e).abs() <= t
                && (r[0][c] - r[1][c]).abs() <= 2. * t),
            "component {c} expected {e} rows {rows:?}"
        );
    }
}
#[test]
fn h0_independent_energy_coordinate_forces_and_scalar_trace() {
    for mesh in MESHES {
        let p = positions(&FRACTIONAL, CELL);
        let b = SimBox::from_array([8., 9., 10.]);
        let mut forces = vec![0.; 9];
        let (u, w) = spme_reciprocal_forces(0.8, &b, &p, &CHARGES, mesh, &mut forces).unwrap();
        let oracle = mesh_energy(&CHARGES, &FRACTIONAL, CELL, mesh, 0.8).unwrap();
        for lambda in [0.97, 1.03] {
            let scaled = CELL.map(|v| v.map(|x| x * lambda));
            let real = mesh_energy(&CHARGES, &FRACTIONAL, scaled, mesh, 0.8).unwrap();
            let complex = uniform_mesh_energy(
                &CHARGES,
                &FRACTIONAL,
                CELL,
                mesh,
                0.8,
                num_complex::Complex::new(lambda, 0.),
            )
            .unwrap();
            near(real, complex.re, 1e-11);
            near(complex.im, 0., 1e-11);
            println!("uniform extension mesh={mesh:?} lambda={lambda} real_energy={real:.16e} complex_energy={complex:?}");
        }
        let complex = uniform_mesh_energy(
            &CHARGES,
            &FRACTIONAL,
            CELL,
            mesh,
            0.8,
            num_complex::Complex::new(1., 1e-20),
        )
        .unwrap();
        near(u, complex.re, 1e-11);
        near(w, -complex.im / 1e-20, 1e-11);
        println!(
            "H0 independent complex energy mesh={mesh:?} W={w:.16e} derivative={:.16e}",
            -complex.im / 1e-20
        );
        near(u, oracle, 1e-11);
        near(
            u,
            spme_reciprocal_energy(0.8, 8., &b, &p, &CHARGES, mesh).unwrap(),
            1e-11,
        );
        for (i, &f) in forces.iter().enumerate() {
            let mut up = FRACTIONAL;
            let mut down = FRACTIONAL;
            up[i / 3][i % 3] += 1e-5 / CELL[i % 3][i % 3];
            down[i / 3][i % 3] -= 1e-5 / CELL[i % 3][i % 3];
            let fd = -(mesh_energy(&CHARGES, &up, CELL, mesh, 0.8).unwrap()
                - mesh_energy(&CHARGES, &down, CELL, mesh, 0.8).unwrap())
                / 2e-5;
            near(f, fd, 2e-7);
            println!("H0 force mesh={mesh:?} coordinate={i} actual={f:.16e} fd={fd:.16e}");
        }
        // Independent isotropic energy derivative certifies legacy scalar W first.
        let mut traces = [0.; 4];
        for (row, eps) in EPS.into_iter().enumerate() {
            let up = CELL.map(|v| v.map(|x| x * (1. + eps)));
            let down = CELL.map(|v| v.map(|x| x * (1. - eps)));
            traces[row] = -(mesh_energy(&CHARGES, &FRACTIONAL, up, mesh, 0.8).unwrap()
                - mesh_energy(&CHARGES, &FRACTIONAL, down, mesh, 0.8).unwrap())
                / (2. * eps);
        }
        assert!(traces
            .windows(2)
            .any(|v| (v[0] - w).abs() <= 2e-6 * (1. + w.abs())
                && (v[1] - w).abs() <= 2e-6 * (1. + w.abs())
                && (v[0] - v[1]).abs() <= 4e-6 * (1. + w.abs())));
        println!("H0 mesh={mesh:?} U={u:.16e} independent_U={oracle:.16e} W={w:.16e} four_trace_fd={traces:?}");
    }
}
#[test]
fn nyquist_extension_cancels_only_affected_shear() {
    let mesh = [4, 6, 8];
    for double in [false, true] {
        let grid: Vec<f64> = (0..4)
            .flat_map(|x| {
                (0..6).flat_map(move |y| {
                    (0..8).map(move |z| {
                        (PI * x as f64
                            + if double {
                                PI * y as f64
                            } else {
                                2. * PI * y as f64 / 6.
                            }
                            + 2. * PI * z as f64 / 8.)
                            .cos()
                    })
                })
            })
            .collect();
        for eps in EPS {
            let d = energy_strains(CELL, eps, |h| grid_energy(&grid, h, mesh, 0.8)).unwrap();
            near(d[3], 0., 2e-6);
            near(d[4], 0., 2e-6);
            if double {
                near(d[5], 0., 2e-6);
            } else {
                assert!(d[5].abs() > 1e-4);
            }
            println!("mixed Nyquist double={double} eps={eps:.8e} energy_derivatives={d:?}");
        }
        let h = [[8., 0.37, -0.16], [0.8, 9., 0.23], [0.3, 0.1, 10.]];
        let potential = grid_potential(&grid, h, mesh, 0.8).unwrap();
        let real = potential.iter().map(|z| z[0] * z[0]).sum::<f64>().sqrt();
        let imag = potential.iter().map(|z| z[1] * z[1]).sum::<f64>().sqrt();
        assert!(imag <= 1e-12 * real.max(1.));
        println!("complex potential double={double} arbitrary_H={h:?} real_norm={real:.16e} imag_norm={imag:.16e}");
    }
}
#[test]
fn oracle_rejects_invalid_or_unbounded_inputs() {
    assert!(mesh_energy(&[1.], &[], CELL, [5; 3], 0.8).is_err());
    assert!(mesh_energy(&[f64::NAN], &[[0.; 3]], CELL, [5; 3], 0.8).is_err());
    assert!(grid_energy(&[0.; 125], [[0.; 3]; 3], [5; 3], 0.8).is_err());
    assert!(grid_energy(&[0.; 125], CELL, [5; 3], 0.).is_err());
    assert!(grid_energy(&[0.; 125], CELL, [6; 3], 0.8).is_err());
    assert!(mesh_energy(&[], &[], CELL, [8; 3], 0.8).is_err());
}

// Full energy controls use no production electrostatics helpers. The real
// branch and selected lattice image are frozen by a geometric preflight.
#[derive(Clone)]
struct Fixture {
    name: &'static str,
    s: Vec<[f64; 3]>,
    q: Vec<f64>,
    scales: Vec<(usize, usize, f64)>,
}
fn fixtures() -> Vec<Fixture> {
    let main = Fixture {
        name: "charged",
        s: FRACTIONAL.to_vec(),
        q: CHARGES.to_vec(),
        scales: vec![],
    };
    let mut scaled = main.clone();
    scaled.name = "charged-scaled";
    scaled.scales = vec![(0, 1, 0.4), (1, 2, 0.)];
    let mut inside = scaled.clone();
    inside.name = "inside-scaled";
    inside.s = vec![[0.13, 0.27, 0.39], [0.31, 0.34, 0.48], [0.32, 0.51, 0.58]];
    let outside = Fixture {
        name: "outside-excluded",
        s: vec![[0.1, 0.1, 0.1], [0.59, 0.4, 0.3]],
        q: vec![1., -0.7],
        scales: vec![(0, 1, 0.)],
    };
    vec![main, scaled, inside, outside]
}
fn nearest(d: [f64; 3], h: [[f64; 3]; 3]) -> ([i32; 3], f64, f64) {
    let mut best = ([0; 3], f64::INFINITY, f64::INFINITY);
    for x in -2..=2 {
        for y in -2..=2 {
            for z in -2..=2 {
                let shift = [x, y, z];
                let ds = std::array::from_fn(|a| d[a] - f64::from(shift[a]));
                let r2 = positions(&[ds], h).iter().map(|v| v * v).sum::<f64>();
                if r2 < best.1 {
                    best = (shift, r2, best.1);
                } else if r2 < best.2 {
                    best.2 = r2;
                }
            }
        }
    }
    best
}
fn branches(f: &Fixture, h: [[f64; 3]; 3]) -> Vec<([i32; 3], bool)> {
    let mut out = Vec::new();
    for i in 0..f.s.len() {
        for j in i + 1..f.s.len() {
            let d = std::array::from_fn(|a| f.s[i][a] - f.s[j][a]);
            let (shift, r2, second) = nearest(d, h);
            assert!(second - r2 > 64. * f64::EPSILON * 100.);
            assert!((r2 - 3.5 * 3.5).abs() > 1e-4);
            out.push((shift, r2 < 3.5 * 3.5));
        }
    }
    out
}
fn geometry_preflight() {
    // Rejected original outside fixture retains its exact-half-cell tie.
    let (_, first, second) = nearest([-0.5, -0.3, -0.2], CELL);
    near(first, second, 2e-14);
    for f in fixtures() {
        let expected = branches(&f, CELL);
        for s in &f.s {
            for x in s {
                let t = (4. * x).fract();
                assert!(t > 1e-4 && t < 1. - 1e-4);
            }
        }
        for c in 0..6 {
            for eps in EPS {
                for sign in [-1., 1.] {
                    assert_eq!(branches(&f, strain_cell(CELL, c, sign * eps)), expected);
                }
            }
        }
        println!(
            "geometry-proof fixture={} all48 images/cutoffs={expected:?}",
            f.name
        );
    }
    for mesh in MESHES {
        for s in FRACTIONAL {
            for a in 0..3 {
                let t = (s[a] * mesh[a] as f64).fract();
                assert!(t > 1e-4 && t < 1. - 1e-4);
            }
        }
    }
}
fn other_energy(f: &Fixture, h: [[f64; 3]; 3]) -> f64 {
    let mut u = 0.;
    for i in 0..f.s.len() {
        for j in i + 1..f.s.len() {
            let d = std::array::from_fn(|a| f.s[i][a] - f.s[j][a]);
            let (_, r2, _) = nearest(d, h);
            let r = r2.sqrt();
            let scale = f
                .scales
                .iter()
                .find(|(a, b, _)| *a == i && *b == j)
                .map_or(1., |p| p.2);
            let qq = f.q[i] * f.q[j];
            if r < 3.5 {
                u += qq * (scale * libm::erfc(0.8 * r) - (1. - scale) * libm::erf(0.8 * r)) / r;
            } else {
                u -= qq * (1. - scale) / r;
            }
        }
    }
    let volume = stress_validation::mesh_oracle::inverse(h).unwrap().0;
    u -= 0.8 / PI.sqrt() * f.q.iter().map(|x| x * x).sum::<f64>();
    u -= PI * f.q.iter().sum::<f64>().powi(2) / (2. * volume * 0.8 * 0.8);
    u
}
fn full_energy(f: &Fixture, h: [[f64; 3]; 3], mesh: [usize; 3]) -> Result<f64, String> {
    Ok(mesh_energy(&f.q, &f.s, h, mesh, 0.8)? + other_energy(f, h))
}
fn reference_modes(cutoff: f64) -> Vec<[i32; 3]> {
    // Full +/- enumeration, independent of the production half-space cache.
    let mut modes = Vec::new();
    for x in -20..=20 {
        for y in -20..=20 {
            for z in -20..=20 {
                let n = [x, y, z];
                let k2 = (0..3)
                    .map(|a| (2. * PI * f64::from(n[a]) / CELL[a][a]).powi(2))
                    .sum::<f64>();
                if k2 > 0. && k2 < cutoff * cutoff {
                    modes.push(n);
                }
            }
        }
    }
    modes
}
fn reference_energy(f: &Fixture, h: [[f64; 3]; 3], modes: &[[i32; 3]]) -> Result<f64, String> {
    let (volume, inv) = stress_validation::mesh_oracle::inverse(h)?;
    let mut u = other_energy(f, h);
    for n in modes {
        let k2 = (0..3)
            .map(|a| (2. * PI * (0..3).map(|b| inv[b][a] * f64::from(n[b])).sum::<f64>()).powi(2))
            .sum::<f64>();
        let mut re = 0.;
        let mut im = 0.;
        for (s, q) in f.s.iter().zip(&f.q) {
            let phase = 2. * PI * (0..3).map(|a| s[a] * f64::from(n[a])).sum::<f64>();
            let (sin, cos) = phase.sin_cos();
            re += q * cos;
            im += q * sin;
        }
        u += 2. * PI / volume * (-k2 / (4. * 0.8 * 0.8)).exp() / k2 * (re * re + im * im);
    }
    Ok(u)
}
fn energy_forces<F: Fn(&Fixture) -> f64>(f: &Fixture, energy: F) -> Vec<f64> {
    (0..f.s.len() * 3)
        .map(|i| {
            let mut up = f.clone();
            let mut down = f.clone();
            up.s[i / 3][i % 3] += 1e-5 / CELL[i % 3][i % 3];
            down.s[i / 3][i % 3] -= 1e-5 / CELL[i % 3][i % 3];
            -(energy(&up) - energy(&down)) / 2e-5
        })
        .collect()
}

#[test]
fn fixed_fixture_geometry_precedes_full_derivatives() {
    geometry_preflight();
}

#[test]
fn arbitrary_cell_complex_potential_has_energy_and_conjugacy() {
    let h = [[8., 0.37, -0.16], [0.8, 9., 0.23], [0.3, 0.1, 10.]];
    for mesh in MESHES {
        let grid: Vec<f64> = (0..mesh.iter().product())
            .map(|i| ((i as f64 + 1.) * 0.371).sin() + 0.3 * ((i as f64 + 1.) * 0.683).cos())
            .collect();
        let p = grid_potential(&grid, h, mesh, 0.8).unwrap();
        let energy = grid_energy(&grid, h, mesh, 0.8).unwrap();
        near(
            energy,
            grid.iter().zip(&p).map(|(q, z)| q * z[0]).sum(),
            1e-12,
        );
        let real = p.iter().map(|z| z[0] * z[0]).sum::<f64>().sqrt();
        let imag = p.iter().map(|z| z[1] * z[1]).sum::<f64>().sqrt();
        assert!(imag <= 1e-12 * real.max(1.));
        // The complete complex convolution is the gradient dU/d(grid)/2.
        for i in [0, grid.len() / 2, grid.len() - 1] {
            let mut up = grid.clone();
            let mut down = grid.clone();
            up[i] += 1e-5;
            down[i] -= 1e-5;
            let fd = (grid_energy(&up, h, mesh, 0.8).unwrap()
                - grid_energy(&down, h, mesh, 0.8).unwrap())
                / 2e-5;
            near(fd, 2. * p[i][0], 2e-7);
        }
        println!("arbitrary-cell mesh={mesh:?} energy={energy:.16e} real_norm={real:.16e} imaginary_norm={imag:.16e}");
    }
}

// First-order scalar energy arithmetic, independent of production force/stress.
#[derive(Clone, Copy)]
struct Dual {
    u: f64,
    d: f64,
}
impl Dual {
    fn constant(u: f64) -> Self {
        Self { u, d: 0. }
    }
    fn plus(self, b: Self) -> Self {
        Self {
            u: self.u + b.u,
            d: self.d + b.d,
        }
    }
    fn scale(self, b: f64) -> Self {
        Self {
            u: self.u * b,
            d: self.d * b,
        }
    }
    fn times(self, b: Self) -> Self {
        Self {
            u: self.u * b.u,
            d: self.d * b.u + self.u * b.d,
        }
    }
    fn inverse(self) -> Self {
        Self {
            u: 1. / self.u,
            d: -self.d / (self.u * self.u),
        }
    }
    fn erf(self) -> Self {
        Self {
            u: libm::erf(self.u),
            d: 2. / PI.sqrt() * (-self.u * self.u).exp() * self.d,
        }
    }
    fn erfc(self) -> Self {
        Self {
            u: libm::erfc(self.u),
            d: -2. / PI.sqrt() * (-self.u * self.u).exp() * self.d,
        }
    }
}
fn other_uniform(f: &Fixture, lambda: f64) -> Dual {
    let mut u = Dual::constant(0.);
    for i in 0..f.s.len() {
        for j in i + 1..f.s.len() {
            let d = std::array::from_fn(|a| f.s[i][a] - f.s[j][a]);
            let (_, r2, _) = nearest(d, CELL);
            let r0 = r2.sqrt();
            let r = Dual {
                u: r0 * lambda,
                d: r0,
            };
            let scale = f
                .scales
                .iter()
                .find(|(a, b, _)| *a == i && *b == j)
                .map_or(1., |p| p.2);
            let screened = if r0 < 3.5 {
                r.scale(0.8)
                    .erfc()
                    .scale(scale)
                    .plus(r.scale(0.8).erf().scale(-(1. - scale)))
            } else {
                Dual::constant(-(1. - scale))
            };
            u = u.plus(screened.times(r.inverse()).scale(f.q[i] * f.q[j]));
        }
    }
    let l = Dual { u: lambda, d: 1. };
    let volume = l.times(l).times(l).scale(720.);
    u = u.plus(
        volume
            .inverse()
            .scale(-PI * f.q.iter().sum::<f64>().powi(2) / (2. * 0.8 * 0.8)),
    );
    u.plus(Dual::constant(
        -0.8 / PI.sqrt() * f.q.iter().map(|x| x * x).sum::<f64>(),
    ))
}
fn tensor(t: TensorAvailability) -> SymTensor6 {
    let TensorAvailability::Complete(t) = t else {
        panic!("missing SPME tensor")
    };
    t
}
fn system(f: &Fixture) -> System {
    let mut s = System::new(SimBox::from_array([8., 9., 10.]));
    s.resize(f.q.len());
    s.positions = positions(&f.s, CELL);
    s.charges = f.q.clone();
    s
}
fn model(f: &Fixture, kcut: f64) -> SpmeForceModel {
    SpmeForceModel::new(0.8, 3.5, kcut)
        .unwrap()
        .with_pair_scales(
            f.scales
                .iter()
                .map(|&(i, j, coulomb)| PairScale {
                    i,
                    j,
                    coulomb,
                    lj: 1.,
                })
                .collect(),
        )
        .unwrap()
}
#[test]
fn reciprocal_six_strains_integer_wraps_and_legacy_arithmetic() {
    for mesh in MESHES {
        let source = system(&fixtures()[0]);
        let mut f = vec![0.125; 9];
        let (u, s) = spme_reciprocal_tensor(
            0.8,
            &source.box_size,
            &source.positions,
            &source.charges,
            &mut f,
            mesh,
        )
        .unwrap();
        let mut old_f = vec![0.125; 9];
        let (old_u, w) = spme_reciprocal_forces(
            0.8,
            &source.box_size,
            &source.positions,
            &source.charges,
            mesh,
            &mut old_f,
        )
        .unwrap();
        assert_eq!(u, old_u);
        assert_eq!(f, old_f);
        near(s.trace(), w, 5e-12);
        let mut rows = [[0.; 6]; 4];
        for (row, eps) in EPS.into_iter().enumerate() {
            rows[row] = energy_strains(CELL, eps, |h| {
                mesh_energy(&CHARGES, &FRACTIONAL, h, mesh, 0.8)
            })
            .unwrap();
            println!(
                "reciprocal mesh={mesh:?} eps={eps:.8e} analytic={:?} fd={:?}",
                s.0, rows[row]
            );
        }
        adjacent(&rows, s.0, 2e-6);
        let shifts = [[2., -1., 3.], [-2., 1., 0.], [1., 2., -2.]];
        let moved: Vec<[f64; 3]> = FRACTIONAL
            .iter()
            .zip(shifts)
            .map(|(s, n)| std::array::from_fn(|a| s[a] + n[a]))
            .collect();
        let mut moved_f = vec![0.125; 9];
        let (moved_u, moved_s) = spme_reciprocal_tensor(
            0.8,
            &source.box_size,
            &positions(&moved, CELL),
            &source.charges,
            &mut moved_f,
            mesh,
        )
        .unwrap();
        near(u, moved_u, 1e-10);
        for (a, b) in s.0.into_iter().zip(moved_s.0) {
            near(a, b, 1e-10);
        }
        for (a, b) in f.into_iter().zip(moved_f) {
            near(a, b, 1e-10);
        }
    }
}
#[test]
fn full_model_h0_energy_forces_scalar_and_all_six_strains() {
    geometry_preflight();
    for f in fixtures() {
        let source = system(&f);
        let m = model(&f, 1.2);
        let mesh = mesh_for_box(&source.box_size, 1.2);
        assert_eq!(mesh, [4; 3]);
        let mut s = source.clone();
        let e = m.evaluate(&mut s, &mut NeighborList::new(3.5)).unwrap();
        let t = tensor(e.tensor.clone());
        let mut legacy = source.clone();
        let mut w = Virial { w: 7. };
        let old = m
            .compute(&mut legacy, &mut w, &mut NeighborList::new(3.5))
            .unwrap();
        assert_eq!(old, e.energy);
        near(w.w - 7., e.scalar_virial, 5e-12);
        near(t.trace(), e.scalar_virial, 5e-12);
        assert_eq!(legacy.forces, s.forces);
        let oracle = full_energy(&f, CELL, mesh).unwrap();
        near(oracle, e.energy, 1e-11);
        for l in [0.97, 1., 1.03] {
            let h = CELL.map(|v| v.map(|x| x * l));
            near(other_uniform(&f, l).u, other_energy(&f, h), 1e-11);
        }
        let rec = uniform_mesh_energy(
            &f.q,
            &f.s,
            CELL,
            mesh,
            0.8,
            num_complex::Complex::new(1., 1e-20),
        )
        .unwrap();
        let other = other_uniform(&f, 1.);
        near(e.energy, rec.re + other.u, 1e-11);
        near(e.scalar_virial, -rec.im / 1e-20 - other.d, 1e-11);
        println!("full H0 fixture={} U={:.16e} independent_U={oracle:.16e} W={:.16e} independent_W={:.16e}",f.name,e.energy,e.scalar_virial,-rec.im/1e-20-other.d);
        let fd_forces = energy_forces(&f, |f| full_energy(f, CELL, mesh).unwrap());
        for (&a, &b) in s.forces.iter().zip(&fd_forces) {
            near(a, b, 2e-7);
        }
        println!(
            "full H0 forces fixture={} actual={:?} independent_fd={fd_forces:?}",
            f.name, s.forces
        );
        let mut rows = [[0.; 6]; 4];
        for (row, eps) in EPS.into_iter().enumerate() {
            rows[row] = energy_strains(CELL, eps, |h| full_energy(&f, h, mesh)).unwrap();
            println!(
                "full six-strain fixture={} eps={eps:.8e} analytic={:?} fd={:?}",
                f.name, t.0, rows[row]
            );
        }
        adjacent(&rows, t.0, 2e-6);
    }
}
#[test]
fn reciprocal_rejects_triclinic_and_known_empty_model_is_complete() {
    let tri = SimBox::triclinic([
        Vec3::new(8., 0., 0.),
        Vec3::new(0.8, 9., 0.),
        Vec3::new(0.3, 0.1, 10.),
    ]);
    assert!(spme_reciprocal_tensor(0.8, &tri, &[], &[], &mut [], [4; 3]).is_err());
    let mut s = System::new(SimBox::from_array([8., 9., 10.]));
    let e = SpmeForceModel::new(0.8, 3.5, 1.2)
        .unwrap()
        .evaluate(&mut s, &mut NeighborList::new(3.5))
        .unwrap();
    assert_eq!(e, modena_core::stress::ModelEvaluation::zero());
}
fn reference_values(f: &Fixture, cutoff: f64) -> Vec<f64> {
    let modes = reference_modes(cutoff);
    let u = reference_energy(f, CELL, &modes).unwrap();
    let forces = energy_forces(f, |f| reference_energy(f, CELL, &modes).unwrap());
    let strain = energy_strains(CELL, 1e-5, |h| reference_energy(f, h, &modes)).unwrap();
    println!("independent Ewald fixture={} cutoff={cutoff} mode_count={} U={u:.16e} forces={forces:?} tensor={strain:?}",f.name,modes.len());
    [vec![u], forces, strain.to_vec()].concat()
}
fn mesh_values(f: &Fixture, n: usize) -> Vec<f64> {
    let source = system(f);
    let mesh = [n; 3];
    let mut forces = vec![0.; source.forces.len()];
    let (ur, sr) = spme_reciprocal_tensor(
        0.8,
        &source.box_size,
        &source.positions,
        &source.charges,
        &mut forces,
        mesh,
    )
    .unwrap();
    // To exercise actual full-model real/correction/self/background, subtract
    // the same small reciprocal evaluation from an actual complete model pass.
    // Large explicit grids never use the bounded direct mesh oracle.
    let mut full = source.clone();
    let e = model(f, 1.2)
        .evaluate(&mut full, &mut NeighborList::new(3.5))
        .unwrap();
    let mut small_f = vec![0.; forces.len()];
    let (us, ss) = spme_reciprocal_tensor(
        0.8,
        &source.box_size,
        &source.positions,
        &source.charges,
        &mut small_f,
        [4; 3],
    )
    .unwrap();
    let u = e.energy - us + ur;
    let t = tensor(e.tensor).plus(ss.scaled(-1.)).plus(sr);
    for (i, x) in forces.iter_mut().enumerate() {
        *x += full.forces[i] - small_f[i];
    }
    [vec![u], forces, t.0.to_vec()].concat()
}
#[test]
fn charged_scaled_translation_converges_on_fixed_ladder_to_independent_ewald() {
    geometry_preflight();
    let mut statuses = Vec::new();
    for f in fixtures() {
        let mut moved = f.clone();
        for s in &mut moved.s {
            for a in 0..3 {
                s[a] += [0.31, -0.27, 0.19][a] / CELL[a][a];
            }
        }
        let mut references = Vec::new();
        for current in [&f, &moved] {
            let rows = [8., 10., 12.].map(|cutoff| reference_values(current, cutoff));
            for adjacent in rows.windows(2) {
                for (a, b) in adjacent[0].iter().zip(&adjacent[1]) {
                    near(*a, *b, 2e-10);
                }
            }
            references.push(rows[2].clone());
        }
        for (a, b) in references[0].iter().zip(&references[1]) {
            near(*a, *b, 2e-10);
        }
        let mut errors = Vec::new();
        let mut translation = Vec::new();
        for n in [8, 16, 32] {
            let values = [mesh_values(&f, n), mesh_values(&moved, n)];
            let error = values
                .iter()
                .zip(&references)
                .flat_map(|(a, b)| a.iter().zip(b).map(|(x, y)| (x - y).abs()))
                .fold(0., f64::max);
            let movement = values[0]
                .iter()
                .zip(&values[1])
                .map(|(x, y)| (x - y).abs())
                .fold(0., f64::max);
            println!("mesh convergence fixture={} mesh={n} max_EFS_error={error:.16e} max_translation_change={movement:.16e} original_EFS={:?} translated_EFS={:?}",f.name,values[0],values[1]);
            errors.push(error);
            translation.push(movement);
        }
        let status =
            classify_mesh_ladder(errors.try_into().unwrap(), translation.try_into().unwrap());
        println!("scientific_accuracy fixture={} status={status:?} tolerance=2e-4; production pressure remains Diagnostic", f.name);
        // The controller retains the measured outside force failure. A green
        // harness certifies classification, never blanket scientific acceptance.
        let expected = if f.name == "outside-excluded" {
            ScientificStatus::Failed
        } else {
            ScientificStatus::Passed
        };
        assert_eq!(
            status, expected,
            "new outcome requires a reviewed protocol epoch"
        );
        statuses.push(status);
    }
    let aggregate = scientific_aggregate(&statuses);
    println!("scientific_accuracy aggregate={aggregate:?}; outside32 force failure retained");
    assert_eq!(aggregate, ScientificStatus::Failed);
}

#[test]
#[allow(clippy::excessive_precision)] // Preserve the original printed observations verbatim.
fn retained_four_case_ladder_preserves_scientific_failure() {
    // Actual fixed-protocol observations from the first 8/9 run, retained at
    // target/b07-evidence/task4-spme-tensor-20261110/green-stress-mesh.log.
    // This regression tests status propagation without retrying that campaign.
    let observations = [
        (
            [
                0.033318724963932844,
                0.002618902157218845,
                0.00016453094863359483,
            ],
            [
                0.025203214839077545,
                0.0032484181719485122,
                0.00013877802985862858,
            ],
        ),
        (
            [
                0.03331872500629632,
                0.0026189021517154885,
                0.00016453093757912163,
            ],
            [
                0.02520321483907755,
                0.003248418171948514,
                0.00013877802985862685,
            ],
        ),
        (
            [
                0.03305855096522901,
                0.0028620501045207303,
                0.0001752401058369088,
            ],
            [
                0.024780657380824495,
                0.0029317841015066744,
                0.00016415162260114552,
            ],
        ),
        (
            [
                0.038419837852351105,
                0.0027527128078139347,
                0.00024207810466960128,
            ],
            [
                0.021198804997027053,
                0.0023688585537265158,
                0.0004829819977220807,
            ],
        ),
    ];
    let statuses =
        observations.map(|(errors, translation)| classify_mesh_ladder(errors, translation));
    assert_eq!(
        statuses,
        [
            ScientificStatus::Passed,
            ScientificStatus::Passed,
            ScientificStatus::Passed,
            ScientificStatus::Failed
        ]
    );
    assert_eq!(scientific_aggregate(&statuses), ScientificStatus::Failed);
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum ScientificStatus {
    Passed,
    Failed,
}
fn classify_mesh_ladder(errors: [f64; 3], translation: [f64; 3]) -> ScientificStatus {
    let finite = errors
        .iter()
        .chain(&translation)
        .all(|x| x.is_finite() && *x >= 0.);
    let decreases = errors
        .windows(2)
        .chain(translation.windows(2))
        .all(|v| v[1] < v[0]);
    if finite && decreases && errors[2] <= 2e-4 && translation[2] <= 2e-4 {
        ScientificStatus::Passed
    } else {
        ScientificStatus::Failed
    }
}
fn scientific_aggregate(statuses: &[ScientificStatus]) -> ScientificStatus {
    if !statuses.is_empty() && statuses.iter().all(|s| *s == ScientificStatus::Passed) {
        ScientificStatus::Passed
    } else {
        ScientificStatus::Failed
    }
}
