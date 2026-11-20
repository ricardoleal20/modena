use modena_core::potential::PairPotential;
use modena_core::stress::{SymTensor6, TensorAvailability};
use modena_core::system::{SimBox, System};
use modena_core::vec3::Vec3;

pub struct Harmonic;
impl PairPotential for Harmonic {
    fn cutoff(&self) -> f64 {
        8.0
    }
    fn cutoff2(&self) -> f64 {
        64.0
    }
    fn energy_force(&self, r2: f64, _: u32, _: u32) -> (f64, f64) {
        (0.5 * r2, -1.0)
    }
}
pub fn atoms(p: &[[f64; 3]], m: &[f64]) -> System {
    let mut s = System::new(SimBox::cubic(20.0));
    s.resize(p.len());
    s.positions = p.iter().flatten().copied().collect();
    s.masses = m.to_vec();
    s
}
pub fn complete(x: &TensorAvailability) -> SymTensor6 {
    match x {
        TensorAvailability::Complete(t) => *t,
        _ => panic!("tensor unavailable"),
    }
}
pub fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1.0 + b.abs()), "{a} != {b}");
}

/// Independently apply each normal or symmetric half-shear generator to
/// coordinates AND cell vectors; no tensor production helpers are used.
pub fn affine(source: &System, component: usize, eta: f64) -> System {
    let deform = |v: Vec3| {
        let mut out = v.0;
        if component < 3 {
            out[component] += eta * v.0[component];
        } else {
            let (a, b) = [(0, 1), (0, 2), (1, 2)][component - 3];
            out[a] += 0.5 * eta * v.0[b];
            out[b] += 0.5 * eta * v.0[a];
        }
        Vec3(out)
    };
    let mut s = source.clone();
    s.box_size = SimBox::triclinic(source.box_size.basis().map(deform));
    for i in 0..s.n() {
        s.set_pos(i, deform(source.pos(i)));
    }
    s
}

/// Frozen oblique free rotor from the B07 specification.
pub fn rotor(m1: f64, m2: f64) -> (System, modena_core::constraints::Constraints) {
    use modena_core::constraints::{Constraint, Constraints};
    let r = Vec3::new(1., 2., 2.) * (1.3 / 3.);
    let u = Vec3::new(2., -1., 0.) * (1.3 * 0.7 / 5_f64.sqrt());
    let center = Vec3::new(5., 6., 7.);
    let velocity = Vec3::new(0.2, -0.1, 0.3);
    let mut s = atoms(
        &[
            (center - r * (m2 / (m1 + m2))).0,
            (center + r * (m1 / (m1 + m2))).0,
        ],
        &[m1, m2],
    );
    s.set_vel(0, velocity - u * (m2 / (m1 + m2)));
    s.set_vel(1, velocity + u * (m1 / (m1 + m2)));
    (
        s,
        Constraints {
            bonds: vec![Constraint { i: 0, j: 1, d: 1.3 }],
            waters: vec![],
        },
    )
}

/// Fixed Rodrigues rotation used by the B07 water fixtures.
pub fn water_rotation(v: Vec3) -> Vec3 {
    let axis = Vec3::new(1., 2., 3.) * (1. / 14_f64.sqrt());
    v * 0.41_f64.cos()
        + axis.cross(v) * 0.41_f64.sin()
        + axis * (axis.dot(v) * (1. - 0.41_f64.cos()))
}

/// Oblique rigid triangle with the specification's principal-axis rotation.
pub fn water(masses: [f64; 3]) -> (System, modena_core::constraints::Constraints) {
    use modena_core::constraints::{Constraints, RigidWater};
    let half = 109.47_f64.to_radians() / 2.;
    let p = [
        Vec3::zero(),
        Vec3::new(half.cos(), half.sin(), 0.),
        Vec3::new(half.cos(), -half.sin(), 0.),
    ];
    let center = (p[0] * masses[0] + p[1] * masses[1] + p[2] * masses[2])
        * (1. / masses.iter().sum::<f64>());
    let omega = water_rotation(Vec3::new(0., 0., 0.7));
    let mut s = atoms(&[[0.; 3]; 3], &masses);
    for (i, x) in p.into_iter().enumerate() {
        let rho = water_rotation(x - center);
        s.set_pos(i, Vec3::new(5., 6., 7.) + rho);
        s.set_vel(i, Vec3::new(0.2, -0.1, 0.3) + omega.cross(rho));
    }
    (
        s,
        Constraints {
            bonds: vec![],
            waters: vec![RigidWater {
                first: 0,
                d_oh: 1.,
                d_hh: 2. * half.sin(),
            }],
        },
    )
}
