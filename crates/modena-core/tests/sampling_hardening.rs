//! Independent analytical sampling and periodic-geometry regressions from the audit.
use modena_core::thermostat::{Thermostat, VelocityRescaling};
use modena_core::{
    Constraint, Constraints, ForceModel, Integrator, LangevinBAOAB, NeighborList, Result, SimBox,
    System, Vec3, VelocityVerlet, Virial,
};

struct ZeroForce;
impl ForceModel for ZeroForce {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.0)
    }
}

fn gas(n: usize) -> System {
    let mut s = System::new(SimBox::cubic(1000.0));
    s.resize(n);
    s.velocities.fill(1.0);
    s
}

#[test]
fn csvr_one_step_has_canonical_conditional_moments() {
    // Exact CIR transition: E[K'|K]=cK+(1-c)fT/2,
    // Var[K'|K]=(1-c)^2 f T^2/2 + 2c(1-c)KT.
    let mut s = gas(4);
    let f = 12.0;
    let k = 2.0;
    let t = 1.4;
    let c = (-1.0_f64).exp();
    let expected_mean = c * k + (1.0 - c) * f * t / 2.0;
    let expected_var = (1.0 - c).powi(2) * f * t * t / 2.0 + 2.0 * c * (1.0 - c) * k * t;
    let mut thermostat = VelocityRescaling::new(t, 1.0, 42);
    let (mut sum, mut sum2) = (0.0, 0.0);
    for _ in 0..80_000 {
        s.velocities.fill((2.0 * k / f).sqrt());
        thermostat.apply(&mut s, 1.0);
        let kp = s.kinetic_energy();
        sum += kp;
        sum2 += kp * kp;
    }
    let mean = sum / 80_000.0;
    let var = sum2 / 80_000.0 - mean * mean;
    assert!(
        (mean - expected_mean).abs() < 0.045,
        "mean={mean}, canonical={expected_mean}"
    );
    assert!(
        (var - expected_var).abs() < 0.14,
        "variance={var}, canonical={expected_var}"
    );
}

#[test]
fn csvr_stationary_kinetic_energy_has_canonical_mean_and_variance() {
    for seed in [7, 42, 91] {
        let mut s = gas(4);
        let mut thermostat = VelocityRescaling::new(1.0, 1.0, seed);
        let (mut sum, mut sum2) = (0.0, 0.0);
        for step in 0..41_000 {
            thermostat.apply(&mut s, 1.0);
            if step >= 1000 {
                let k = s.kinetic_energy();
                sum += k;
                sum2 += k * k;
            }
        }
        let mean = sum / 40_000.0;
        let variance = sum2 / 40_000.0 - mean * mean;
        // K~Gamma(f/2,T), hence mean=6 and variance=6.
        assert!(
            (mean - 6.0).abs() < 0.085,
            "seed={seed}, mean K={mean}, canonical=6"
        );
        assert!(
            (variance - 6.0).abs() < 0.25,
            "seed={seed}, variance K={variance}, canonical=6"
        );
    }
}

#[test]
fn csvr_decorrelated_signed_factor_can_reverse_velocity_direction() {
    let mut s = gas(1);
    let mut thermostat = VelocityRescaling::new(1.0, 1.0, 13);
    let mut negative = 0;
    for _ in 0..10_000 {
        s.velocities.fill(1.0);
        thermostat.apply(&mut s, 100.0);
        negative += usize::from(s.velocities[0] < 0.0);
    }
    assert!(
        (4500..5500).contains(&negative),
        "negative signed draws={negative}, expected half"
    );
}

#[test]
fn baoab_unrelated_dimer_preserves_free_atom_maxwell_variance() {
    for seed in [7, 42, 91] {
        for constrained in [false, true] {
            let mut s = gas(3);
            s.velocities.fill(0.0);
            s.masses[2] = 2.0;
            s.set_pos(0, Vec3::new(100.0, 100.0, 100.0));
            s.set_pos(1, Vec3::new(101.0, 100.0, 100.0));
            s.set_pos(2, Vec3::new(200.0, 200.0, 200.0));
            let mut cons = Constraints::new();
            cons.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
            let mut integrator = if constrained {
                LangevinBAOAB::with_constraints(100.0, 1.0, seed, cons)
            } else {
                LangevinBAOAB::new(100.0, 1.0, seed)
            };
            let mut neighbors = NeighborList::new(2.5);
            let mut sum = Vec3::zero();
            let mut sum2 = Vec3::zero();
            for step in 0..26_000 {
                integrator
                    .step(&mut s, &ZeroForce, &mut neighbors, 0.01)
                    .unwrap();
                if step >= 1000 {
                    let v = s.vel(2);
                    sum += v;
                    sum2 += Vec3::new(v.0[0] * v.0[0], v.0[1] * v.0[1], v.0[2] * v.0[2]);
                }
            }
            for k in 0..3 {
                let mean = sum.0[k] / 25_000.0;
                let variance = sum2.0[k] / 25_000.0 - mean * mean;
                assert!((variance-0.5).abs()<0.035,
                    "seed={seed}, constrained={constrained}, axis={k}, variance={variance}, T/m=0.5");
            }
        }
    }
}

#[test]
fn failed_constraint_projection_is_an_integration_error() {
    let mut s = gas(2);
    s.positions.fill(0.0); // coincident reference: no possible SHAKE direction
    let mut cons = Constraints::new();
    cons.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
    let mut integrator = VelocityVerlet::with_constraints(cons);
    assert!(
        integrator
            .step(&mut s, &ZeroForce, &mut NeighborList::new(2.5), 0.01)
            .is_err(),
        "singular constraint projection silently succeeded"
    );
}

fn skew_cell() -> SimBox {
    SimBox::triclinic([
        Vec3::new(10.0, 0.0, 0.0),
        Vec3::new(4.9, 10.0, 0.0),
        Vec3::new(0.0, 0.0, 10.0),
    ])
}

#[test]
fn triclinic_image_is_nearest_euclidean_lattice_image() {
    let image = skew_cell().min_image(Vec3::new(6.37, 3.0, 0.0));
    assert!(
        (image - Vec3::new(-3.63, 3.0, 0.0)).norm() < 1e-12,
        "nearest image={image:?}"
    );
}

#[test]
fn triclinic_search_covers_translations_beyond_adjacent_fractional_cells() {
    let b = SimBox::triclinic([
        Vec3::new(10.0, 0.0, 0.0),
        Vec3::new(49.0, 10.0, 0.0),
        Vec3::new(0.0, 0.0, 10.0),
    ]);
    let d = Vec3::new(19.6, 3.0, 0.0); // fractional (.49,.3,0), nearest shift is 2*a.
    assert!((b.min_image(d) - Vec3::new(-0.4, 3.0, 0.0)).norm() < 1e-12);
}

#[test]
fn triclinic_scaling_updates_basis_inverse_and_volume() {
    let mut b = skew_cell();
    b.scale_uniform(2.0);
    assert!(
        (b.volume() - 8000.0).abs() < 1e-9,
        "scaled volume={}",
        b.volume()
    );
    assert!((b.wrap(Vec3::new(25.0, 2.0, 0.0)) - Vec3::new(5.0, 2.0, 0.0)).norm() < 1e-12);
    b.scale_axis(1, 0.5);
    assert!((b.volume() - 4000.0).abs() < 1e-9);
}

/// Assessment of the pressure limitation, not a claim of validated molecular
/// stress. A centered free rotor needs W_constraint=-2K_rot, because its
/// centripetal internal forces cancel the rotational kinetic pressure.
#[test]
fn rotating_rigid_dimer_exposes_missing_constraint_pressure() {
    let mut s = System::new(SimBox::cubic(100.0));
    s.resize(2);
    let center = Vec3::new(50.0, 50.0, 50.0);
    s.set_pos(0, center + Vec3::new(-0.5, 0.0, 0.0));
    s.set_pos(1, center + Vec3::new(0.5, 0.0, 0.0));
    s.set_vel(0, Vec3::new(0.0, 0.5, 0.0));
    s.set_vel(1, Vec3::new(0.0, -0.5, 0.0));
    let mut cons = Constraints::new();
    cons.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
    let mut integrator = VelocityVerlet::with_constraints(cons);
    let mut neighbors = NeighborList::new(2.5);
    let k = s.kinetic_energy();
    let dt = 0.001;
    let before = s.clone();
    let e = integrator
        .step(&mut s, &ZeroForce, &mut neighbors, dt)
        .unwrap();
    let mut impulse_virial = 0.0;
    for i in 0..s.n() {
        let midpoint = (before.pos(i) + s.pos(i)) * 0.5 - center;
        impulse_virial += midpoint.dot(s.vel(i) - before.vel(i)) * s.masses[i] / dt;
    }
    assert!(
        (impulse_virial + 2.0 * k).abs() < 1e-7,
        "rotor impulse W={impulse_virial}, expected={}",
        -2.0 * k
    );
    assert_eq!(
        e.virial, 0.0,
        "reported virial currently includes model forces only"
    );
    assert!(
        e.pressure > 0.0,
        "unreported rotational constraint pressure must remain an explicit limitation"
    );
    eprintln!("rotor K={k}, independent impulse W={impulse_virial}, reported W={}, reported P={}, complete rotor P=0",e.virial,e.pressure);
}
