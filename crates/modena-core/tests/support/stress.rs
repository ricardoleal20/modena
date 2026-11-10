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
