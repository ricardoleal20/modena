//! Bounded numerical validation. Independent expectations are recorded beside observations.
use modena_core::thermostat::{Thermostat, VelocityRescaling};
use modena_core::{
    Constraint, Constraints, EwaldForceModel, ForceModel, Integrator, LangevinBAOAB, LennardJones,
    NeighborList, PairForceModel, PairPotential, Result, SimBox, SpmeForceModel, System, Vec3,
    VelocityVerlet, Virial,
};
use serde_json::{json, Value};

struct ZeroForce;
impl ForceModel for ZeroForce {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.0)
    }
}
fn observation(
    rows: &mut Vec<Value>,
    name: String,
    observed: f64,
    expected: f64,
    tolerance: f64,
    basis: &str,
) {
    rows.push(json!({"name":name,"observed":observed,"expected":expected,
        "absolute_tolerance":tolerance,"independent_basis":basis}));
}
fn gas(n: usize) -> System {
    let mut s = System::new(SimBox::cubic(1000.0));
    s.resize(n);
    s.velocities.fill(1.0);
    s
}
fn eval(m: &dyn ForceModel, s: &System) -> (f64, f64, Vec<f64>) {
    let mut s = s.clone();
    s.zero_forces();
    let mut v = Virial::default();
    let u = m
        .compute(&mut s, &mut v, &mut NeighborList::new(4.2))
        .unwrap();
    (u, v.w, s.forces)
}
fn strain(s: &System, factor: f64) -> System {
    let mut s = s.clone();
    s.box_size.scale_uniform(factor);
    for x in &mut s.positions {
        *x *= factor;
    }
    s
}
fn main() {
    let extended = std::env::args().any(|x| x == "--extended");
    let samples = if extended { 400_000 } else { 40_000 };
    let mut rows = Vec::new();
    for seed in [7, 42, 91] {
        let mut s = gas(4);
        let mut thermostat = VelocityRescaling::new(1.0, 1.0, seed);
        let (mut sum, mut sum2) = (0.0, 0.0);
        for step in 0..samples + 1000 {
            thermostat.apply(&mut s, 1.0);
            if step >= 1000 {
                let k = s.kinetic_energy();
                sum += k;
                sum2 += k * k;
            }
        }
        let mean = sum / samples as f64;
        let variance = sum2 / samples as f64 - mean * mean;
        observation(
            &mut rows,
            format!("csvr_seed{seed}_mean"),
            mean,
            6.0,
            0.085,
            "K~Gamma(f/2,T); f=12,T=1",
        );
        observation(
            &mut rows,
            format!("csvr_seed{seed}_variance"),
            variance,
            6.0,
            0.25,
            "Var K=f*T^2/2; fixed acceptance band, correlated samples",
        );
        for constrained in [false, true] {
            let mut s = gas(3);
            s.velocities.fill(0.0);
            s.masses[2] = 2.0;
            s.set_pos(0, Vec3::splat(100.0));
            s.set_pos(1, Vec3::new(101.0, 100.0, 100.0));
            s.set_pos(2, Vec3::splat(200.0));
            let mut c = Constraints::new();
            c.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
            let mut integrator = if constrained {
                LangevinBAOAB::with_constraints(100.0, 1.0, seed, c)
            } else {
                LangevinBAOAB::new(100.0, 1.0, seed)
            };
            let mut nl = NeighborList::new(2.5);
            let (mut sum, mut sum2) = (Vec3::zero(), Vec3::zero());
            for step in 0..samples + 1000 {
                integrator.step(&mut s, &ZeroForce, &mut nl, 0.01).unwrap();
                if step >= 1000 {
                    let v = s.vel(2);
                    sum += v;
                    sum2 += Vec3::new(v.0[0] * v.0[0], v.0[1] * v.0[1], v.0[2] * v.0[2]);
                }
            }
            for axis in 0..3 {
                let mean = sum.0[axis] / samples as f64;
                let var = sum2.0[axis] / samples as f64 - mean * mean;
                observation(
                    &mut rows,
                    format!("baoab_seed{seed}_constrained{constrained}_axis{axis}_mean"),
                    mean,
                    0.0,
                    0.025,
                    "Independent free atom Maxwell mean=0; constraints on separate dimer",
                );
                observation(
                    &mut rows,
                    format!("baoab_seed{seed}_constrained{constrained}_axis{axis}_variance"),
                    var,
                    0.5,
                    0.035,
                    "Maxwell Var(v)=T/m=1/2; fixed band, correlated samples",
                );
            }
        }
    }
    // Nonlinear rigid triangle: three rotational modes, Krot~Gamma(3/2,T).
    for seed in [7, 42, 91] {
        let mut s = System::new(SimBox::cubic(100.0));
        s.resize(3);
        for i in 0..3 {
            let a = 2.0 * std::f64::consts::PI * i as f64 / 3.0;
            s.set_pos(i, Vec3::splat(50.0) + Vec3::new(a.cos(), a.sin(), 0.0));
        }
        let c = Constraints {
            bonds: vec![],
            waters: vec![modena_core::constraints::RigidWater {
                first: 0,
                d_oh: 3.0_f64.sqrt(),
                d_hh: 3.0_f64.sqrt(),
            }],
        };
        let mut integrator = LangevinBAOAB::with_constraints(2.0, 1.0, seed, c.clone());
        let mut nl = NeighborList::new(2.5);
        let (mut sum, mut sum2) = (0.0, 0.0);
        let mut residual = 0.0_f64;
        for step in 0..samples + 1000 {
            integrator.step(&mut s, &ZeroForce, &mut nl, 0.01).unwrap();
            if step >= 1000 {
                let com = s.total_momentum() * (1.0 / 3.0);
                let k = (0..3).map(|i| 0.5 * (s.vel(i) - com).norm2()).sum::<f64>();
                sum += k;
                sum2 += k * k;
                residual = residual.max(c.max_residual(&s));
            }
        }
        let mean = sum / samples as f64;
        let var = sum2 / samples as f64 - mean * mean;
        let corr = (-0.04_f64).exp();
        let inflation = (1.0 + corr) / (1.0 - corr);
        observation(
            &mut rows,
            format!("rigid_rotation_seed{seed}_mean"),
            mean,
            1.5,
            5.0 * (1.5 * inflation / samples as f64).sqrt() + 0.01,
            "Three quadratic rotational modes: Gamma(3/2,1); conservative OU correlation band",
        );
        observation(
            &mut rows,
            format!("rigid_rotation_seed{seed}_variance"),
            var,
            1.5,
            5.0 * (13.5 * inflation / samples as f64).sqrt() + 0.01,
            "Gamma variance=3/2; centered-square variance=13.5, conservative OU correlation band",
        );
        observation(
            &mut rows,
            format!("rigid_rotation_seed{seed}_constraint_residual"),
            residual,
            0.0,
            1e-10,
            "Three fixed analytical pair distances sqrt(3)",
        );
    }
    let lj = LennardJones::new(1.0, 1.0, 4.2).unwrap();
    let h = 1e-6_f64;
    for r in [1.1_f64, 1.4, 2.1] {
        let (_, f) = lj.energy_force(r * r, 0, 0);
        let fd = -(lj.energy_force((r + h).powi(2), 0, 0).0
            - lj.energy_force((r - h).powi(2), 0, 0).0)
            / (2.0 * h);
        observation(
            &mut rows,
            format!("lj_force_r{r}"),
            f * r,
            fd,
            1e-7,
            "Independent central energy derivative; h=1e-6",
        );
    }
    let mut pair = System::new(SimBox::cubic(20.0));
    pair.resize(2);
    pair.set_pos(0, Vec3::new(9.4, 10.0, 10.0));
    pair.set_pos(1, Vec3::new(10.6, 10.0, 10.0));
    pair.set_vel(0, Vec3::new(0.0, 0.1, 0.0));
    pair.set_vel(1, Vec3::new(0.0, -0.1, 0.0));
    let model = PairForceModel::new(lj);
    let mut bands = Vec::new();
    for dt in [0.004_f64, 0.002, 0.001] {
        let mut s = pair.clone();
        let mut v = Virial::default();
        let mut nl = NeighborList::with_skin(4.2, 0.3);
        let e0 = model.compute(&mut s, &mut v, &mut nl).unwrap() + s.kinetic_energy();
        let mut vv = VelocityVerlet::new();
        let mut band = 0.0_f64;
        for _ in 0..(2.0 / dt).round() as usize {
            let e = vv.step(&mut s, &model, &mut nl, dt).unwrap();
            band = band.max((e.total - e0).abs());
        }
        observation(
            &mut rows,
            format!("nve_energy_band_dt{dt}"),
            band,
            0.0,
            0.001,
            "Two-body smooth LJ Hamiltonian; bounded energy error",
        );
        bands.push(band);
    }
    for i in 0..2 {
        observation(
            &mut rows,
            format!("nve_dt_halving_band_ratio{i}"),
            bands[i] / bands[i + 1],
            4.0,
            0.3,
            "Velocity Verlet global energy error O(dt^2), same physical duration=2",
        );
    }
    let mut s = System::new(SimBox::cubic(9.0));
    s.resize(3);
    s.charges = vec![1.0, -0.7, -0.3];
    for (i, p) in [
        Vec3::new(1.1, 2.3, 3.4),
        Vec3::new(2.7, 3.2, 3.8),
        Vec3::new(6.2, 5.8, 7.0),
    ]
    .into_iter()
    .enumerate()
    {
        s.set_pos(i, p);
    }
    let direct = EwaldForceModel::new(0.8, 4.2, 6.4).unwrap();
    let (reference, _, reference_f) = eval(&direct, &s);
    for kcut in [4.8, 6.4] {
        let m = SpmeForceModel::new(0.8, 4.2, kcut).unwrap();
        let (u, w, f) = eval(&m, &s);
        let fd = -(eval(&m, &strain(&s, 1.0 + h)).0 - eval(&m, &strain(&s, 1.0 - h)).0) / (2.0 * h);
        observation(
            &mut rows,
            format!("spme_kcut{kcut}_energy"),
            u,
            reference,
            0.005,
            "Direct reciprocal Ewald sum, same kappa=.8 real rc=4.2; neutral 3-charge fixture",
        );
        let max_force = f
            .iter()
            .zip(&reference_f)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        observation(
            &mut rows,
            format!("spme_kcut{kcut}_force_max_error"),
            max_force,
            0.0,
            0.005,
            "Direct Ewald independent reciprocal summation",
        );
        observation(
            &mut rows,
            format!("spme_kcut{kcut}_strain_virial"),
            w,
            fd,
            1e-7,
            "Central homogeneous strain -dU/dln(lambda), full box+positions scaled; h=1e-6",
        );
    }
    let (_, w, _) = eval(&direct, &s);
    let fd = -(eval(&direct, &strain(&s, 1.0 + h)).0 - eval(&direct, &strain(&s, 1.0 - h)).0)
        / (2.0 * h);
    observation(
        &mut rows,
        "direct_ewald_strain_virial".into(),
        w,
        fd,
        1e-7,
        "Central full-cell strain derivative",
    );
    let mut system = gas(4);
    let mut integrator = LangevinBAOAB::new(1.23, 0.7, 42);
    let mut nl = NeighborList::new(2.5);
    for _ in 0..13 {
        integrator
            .step(&mut system, &ZeroForce, &mut nl, 0.01)
            .unwrap();
    }
    let state: Value = serde_json::from_str(
        &serde_json::to_string(&integrator.checkpoint_state().unwrap()).unwrap(),
    )
    .unwrap();
    let mut split = system.clone();
    let mut restored = LangevinBAOAB::new(1.23, 0.7, 42);
    restored.restore_state(&state).unwrap();
    let mut split_nl = NeighborList::new(2.5);
    for _ in 0..17 {
        integrator
            .step(&mut system, &ZeroForce, &mut nl, 0.01)
            .unwrap();
        restored
            .step(&mut split, &ZeroForce, &mut split_nl, 0.01)
            .unwrap();
    }
    let seam = system
        .positions
        .iter()
        .zip(&split.positions)
        .chain(system.velocities.iter().zip(&split.velocities))
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    observation(
        &mut rows,
        "langevin_json_restart_seam_max_error".into(),
        seam,
        0.0,
        0.0,
        "Unsplit deterministic trajectory with identical stochastic stream; 13+17 steps",
    );
    println!(
        "{}",
        json!({"schema_version":1,"campaign":if extended {"extended"}else{"smoke"},"samples_per_seed":samples,"seeds":[7,42,91],"electrostatics":{"box":9.0,"kappa":0.8,"real_cutoff":4.2,"kcuts":[4.8,6.4],"meshes":[modena_core::spme::mesh_for_box(&s.box_size,4.8),modena_core::spme::mesh_for_box(&s.box_size,6.4)],"precision":"f64","reference":"direct Ewald kcut6.4"},"sampling_note":"Correlated moments with fixed regression tolerances; no claim of formal uncertainty or converged liquid observables","observations":rows})
    );
}
