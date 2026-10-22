//! Read-only audit probes using the public MODENA API. No engine changes.
use modena_core::{
    lattice, Constraint, Constraints, ForceModel, InputConfig, Integrator, LangevinBAOAB,
    LennardJones, NeighborList, PairForceModel, Result, SimBox, Simulation, System, Vec3, Virial,
};
use modena_core::thermostat::{Thermostat, VelocityRescaling};
use std::path::Path;

struct ZeroForce;
impl ForceModel for ZeroForce {
    fn compute(&self, _: &mut System, _: &mut Virial, _: &mut NeighborList) -> Result<f64> {
        Ok(0.0)
    }
}

fn bussi() {
    let mut s = lattice::fcc(3.0, [1, 1, 1], "X");
    s.velocities.fill(1.0); // KE=Nf/2: exactly target T=1 on the engine's DOF.
    let mut t = VelocityRescaling::new(1.0, 1.0, 42);
    let mut sum = 0.0;
    let mut sum2 = 0.0;
    let samples = 200_000;
    for i in 0..samples + 1000 {
        t.apply(&mut s, 1.0);
        if i >= 1000 {
            let temp = s.temperature();
            sum += temp;
            sum2 += temp * temp;
        }
    }
    let nf = s.dof() as f64;
    println!("bussi Nf={nf:.0} samples={samples} target_T=1 measured_mean_T={:.8} measured_T_variance={:.8} predicted_bug_mean_T={:.8}",
        sum/samples as f64, sum2/samples as f64-(sum/samples as f64).powi(2), (nf+1.0)/(2.0*nf));
}

fn baoab(constrained: bool) {
    let mut s = System::new(SimBox::new(Vec3::new(1000.0, 1000.0, 1000.0)));
    s.resize(3);
    s.set_pos(0, Vec3::new(100.0, 100.0, 100.0));
    s.set_pos(1, Vec3::new(101.0, 100.0, 100.0));
    s.set_pos(2, Vec3::new(200.0, 200.0, 200.0)); // free atom: never constrained.
    let mut cons = Constraints::new();
    cons.bonds.push(Constraint { i: 0, j: 1, d: 1.0 });
    let mut integrator = if constrained {
        LangevinBAOAB::with_constraints(100.0, 1.0, 42, cons)
    } else {
        LangevinBAOAB::new(100.0, 1.0, 42)
    };
    let mut nb = NeighborList::new(2.5);
    let samples = 100_000;
    let mut sum = 0.0;
    let mut mean = Vec3::zero();
    for i in 0..samples + 1000 {
        integrator.step(&mut s, &ZeroForce, &mut nb, 0.01).unwrap();
        if i >= 1000 {
            let v = s.vel(2);
            sum += v.dot(v) / 3.0;
            mean += v;
        }
    }
    let c = (-1.0_f64).exp();
    println!("baoab constrained={constrained} free_atom_samples={samples} target_v_component_variance=1 measured_mean_v2_per_component={:.8} predicted_mean_v2={:.8} mean_velocity={:?}",
        sum/samples as f64, if constrained {(1.0+c)/(3.0+c)} else {1.0}, mean*(1.0/samples as f64));
}

fn barostat(enabled: bool) {
    let base = r#"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = [2, 2, 2]
element = "Ar"
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.005
[run]
steps = 1
initial_temperature = 1.0
seed = 42
"#;
    let text = if enabled { format!("{base}\n[barostat]\ntype=\"berendsen\"\ntarget_pressure=0.0\ntau=0.1\ncompressibility=1.0\n") } else { base.into() };
    let config = InputConfig::from_str(&text).unwrap();
    let mut sim = Simulation::from_config(&config, Path::new(".")).unwrap();
    let out = sim.run().unwrap();
    let cached = sim.system().forces.clone();
    let mut s = sim.system().clone();
    s.forces.fill(0.0);
    let model = PairForceModel::new(LennardJones::new(1.0, 1.0, 2.5).unwrap());
    let mut virial = Virial::default();
    let pe = model.compute(&mut s, &mut virial, &mut NeighborList::new(2.5)).unwrap();
    let diff2: f64 = cached.iter().zip(&s.forces).map(|(a,b)|(a-b).powi(2)).sum();
    let norm2: f64 = s.forces.iter().map(|f|f*f).sum();
    println!("barostat enabled={enabled} reported_PE={:.10} recomputed_PE={pe:.10} abs_PE_error={:.10} relative_force_l2_error={:.8}",
        out.final_energies.potential, (out.final_energies.potential-pe).abs(), (diff2/norm2.max(1e-30)).sqrt());
}

fn main() {
    bussi();
    baoab(false);
    baoab(true);
    barostat(false);
    barostat(true);
}
