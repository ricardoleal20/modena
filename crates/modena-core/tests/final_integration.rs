//! Final review regressions: mesh COM, bounded alloy populations, and restart imports.
use modena_core::{
    random_alloy, AlloySpecies, EwaldForceModel, ForceModel, InputConfig, NeighborList, SimBox,
    Simulation, SpmeForceModel, System, Vec3, Virial,
};
use std::path::{Path, PathBuf};

static FIXTURE_ID: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let d = std::env::temp_dir().join(format!(
            "modena_final_{}_{}_{}",
            std::process::id(),
            FIXTURE_ID.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&d).unwrap();
        std::fs::write(d.join("charges.xyz"), "3\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nX 1.1 2.3 3.4 0 1.0 1.0\nX 2.7 3.2 3.8 0 1.0 -0.7\nX 6.2 5.8 7.0 0 1.0 -0.3\n").unwrap();
        Self(d)
    }
    fn sim(&self, text: &str) -> Simulation {
        Simulation::from_config_without_outputs(&InputConfig::from_str(text).unwrap(), &self.0)
            .unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn config(method: &str, controller: &str) -> String {
    format!(
        r#"
[structure]
type="xyz"
file="charges.xyz"
box_size=[9.0,9.0,9.0]
[model]
type="lennard_jones"
sigma=1.0
epsilon=0.0
cutoff=4.2
[integrator]
{controller}
dt=0.001
[electrostatics]
method="{method}"
kappa=0.8
real_cutoff=4.2
kspace_cutoff=4.8
[run]
steps=10
initial_temperature=1.0
seed=51
log_every=1
"#
    )
}
const VV: &str = "type=\"velocity_verlet\"";
const NHC: &str = "type=\"nose_hoover_chain\"\ntarget_temperature=1.0\ntau=0.1\nchain_len=3";
const LAN: &str = "type=\"langevin_baoab\"\ngamma=0.7\ntarget_temperature=1.0\nseed=771";
fn eval(model: &dyn ForceModel, sys: &System) -> (f64, Vec3) {
    let mut s = sys.clone();
    s.zero_forces();
    let u = model
        .compute(
            &mut s,
            &mut Virial::default(),
            &mut NeighborList::with_skin(4.2, 0.3),
        )
        .unwrap();
    let f = (0..s.n()).fold(Vec3::zero(), |f, i| f + s.force(i));
    (u, f)
}
#[test]
fn charged_mesh_net_force_is_common_translation_energy_derivative() {
    let fixture = Fixture::new();
    let s = fixture.sim(&config("spme", VV)).system().clone();
    for (k, want_x) in [(4.8, -0.002196221738), (6.4, 0.000145507097)] {
        let model = SpmeForceModel::new(0.8, 4.2, k).unwrap();
        let (_, f) = eval(&model, &s);
        assert!((f.0[0] - want_x).abs() < 1e-11, "net force {f:?}");
        let h = 1e-5;
        let mut hi = s.clone();
        let mut lo = s.clone();
        for i in 0..s.n() {
            hi.set_pos(i, s.pos(i) + Vec3::new(h, 0.0, 0.0));
            lo.set_pos(i, s.pos(i) - Vec3::new(h, 0.0, 0.0));
        }
        let derivative = -(eval(&model, &hi).0 - eval(&model, &lo).0) / (2.0 * h);
        assert!(
            (derivative - f.0[0]).abs() < 2e-10,
            "net force {} vs -dU/dx {derivative}",
            f.0[0]
        );
    }
    assert!(
        eval(&EwaldForceModel::new(0.8, 4.2, 6.4).unwrap(), &s)
            .1
            .norm()
            < 1e-14
    );
}
#[test]
fn charged_mesh_controllers_initialize_and_keep_full_kinetic_space() {
    let fixture = Fixture::new();
    let csvr = format!(
        "{}\n[thermostat]\ntype=\"velocity_rescaling\"\ntarget_temperature=1.0\ntau=0.1\nseed=231",
        config("spme", VV)
    );
    let mc=format!("{csvr}\n[barostat]\ntype=\"monte_carlo\"\ntarget_temperature=1.0\ntarget_pressure=0.1\ninterval=3\nmax_log_volume=0.02\nseed=93");
    for c in [
        config("spme", VV),
        csvr,
        config("spme", NHC),
        config("spme", LAN),
        mc,
    ] {
        let mut s = fixture.sim(&c);
        assert_eq!(s.metadata().conserved_com_dof, 0, "{c}");
        assert_eq!(s.metadata().kinetic_dof, 9);
        assert!((s.system().kinetic_energy() - 4.5).abs() < 1e-12);
        assert!(
            s.system().total_momentum().norm() > 1e-6,
            "full space must draw COM modes"
        );
        s.run().unwrap();
        assert_eq!(s.metadata().kinetic_dof, 9);
    }
}
#[test]
fn constrained_charged_mesh_has_eight_modes_and_valid_tangent_velocities() {
    let fixture = Fixture::new();
    // Distance between recorded atoms 0 and 1: sqrt(1.6²+0.9²+0.4²).
    for controller in [VV, NHC, LAN] {
        let c = format!(
            "{}\n[[constraints.bonds]]\natoms=[0,1]\nlength=1.8788294228055935",
            config("spme", controller)
        );
        let mut s = fixture.sim(&c);
        assert_eq!(s.metadata().constraint_count, 1);
        assert_eq!(s.metadata().conserved_com_dof, 0);
        assert_eq!(s.metadata().kinetic_dof, 8);
        assert!((s.system().kinetic_energy() - 4.0).abs() < 1e-12);
        s.run().unwrap();
        let sys = s.system();
        let d = sys.box_size.min_image(sys.pos(0) - sys.pos(1));
        assert!((d.norm() - 1.8788294228055935).abs() < 1e-9);
        assert!(d.dot(sys.vel(0) - sys.vel(1)).abs() < 1e-9);
    }
}
#[test]
fn momentum_conserving_and_zero_charge_controls_keep_com_reduction() {
    let fixture = Fixture::new();
    for method in ["ewald", "spme"] {
        if method == "spme" {
            let path = fixture.0.join("charges.xyz");
            std::fs::write(
                path,
                "3\nplain neutral atoms\nX 1.1 2.3 3.4\nX 2.7 3.2 3.8\nX 6.2 5.8 7.0\n",
            )
            .unwrap();
        }
        let mut s = fixture.sim(&config(method, VV));
        assert_eq!(s.metadata().conserved_com_dof, 3);
        assert_eq!(s.metadata().kinetic_dof, 6);
        assert!(s.system().total_momentum().norm() < 1e-14);
        s.run().unwrap();
        assert!(s.system().total_momentum().norm() < 1e-13);
        let langevin = fixture.sim(&config(method, LAN));
        assert_eq!(langevin.metadata().conserved_com_dof, 0);
    }
}
#[test]
fn public_velocity_initializer_respects_declared_com_manifold() {
    for (n, com) in [(1, 0), (3, 0), (3, 3)] {
        let mut s = System::new(SimBox::cubic(9.0));
        s.resize(n);
        s.set_conserved_com_dof(com).unwrap();
        modena_core::runner::init_velocities(&mut s, 1.0, 51);
        assert!(
            (s.kinetic_energy()
                - if n == 1 {
                    1.5
                } else if com == 0 {
                    4.5
                } else {
                    3.0
                })
            .abs()
                < 1e-12
        );
        assert_eq!(s.conserved_com_dof(), com);
        if com == 0 {
            assert!(s.total_momentum().norm() > 1e-6);
        } else {
            assert!(s.total_momentum().norm() < 1e-14);
        }
    }
}
#[test]
fn charged_mesh_split_continuation_preserves_full_space() {
    let fixture = Fixture::new();
    let csvr = format!(
        "{}\n[thermostat]\ntype=\"velocity_rescaling\"\ntarget_temperature=1.0\ntau=0.1\nseed=231",
        config("spme", VV)
    );
    for c in [
        config("spme", VV),
        csvr,
        config("spme", NHC),
        config("spme", LAN),
    ] {
        let mut whole = fixture.sim(&c);
        whole.run().unwrap();
        let mut first = fixture.sim(&c.replace("steps=10", "steps=4"));
        first.run().unwrap();
        first.write_checkpoint(&fixture.0.join("cp.chk")).unwrap();
        let mut second = fixture.sim(&c.replace("steps=10", "steps=6\nrestore_file=\"cp.chk\""));
        assert_eq!(second.metadata().kinetic_dof, 9);
        second.run().unwrap();
        assert_eq!(whole.system().positions, second.system().positions);
        assert_eq!(whole.system().velocities, second.system().velocities);
        assert_eq!(whole.system().forces, second.system().forces);
    }
}
#[test]
fn unsupported_mesh_com_checkpoint_is_rejected_and_state_only_reassesses_model() {
    let fixture = Fixture::new();
    let c = config("spme", VV);
    let s = fixture.sim(&c);
    s.write_checkpoint(&fixture.0.join("cp.chk")).unwrap();
    let path = fixture.0.join("cp.chk");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut cp: serde_json::Value = serde_json::from_str(text.split_once('\n').unwrap().1).unwrap();
    // Recreate old unsupported COM3 metadata, with exactly zero imported momentum.
    cp["continuation"]["conserved_com_dof"] = serde_json::json!(3);
    cp["velocities"] = serde_json::json!([1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    std::fs::write(&path, format!("MODENA_RESTART 2\n{cp}\n")).unwrap();
    let resume = c.replace("steps=10", "steps=0\nrestore_file=\"cp.chk\"");
    assert!(
        Simulation::from_config_without_outputs(
            &InputConfig::from_str(&resume).unwrap(),
            &fixture.0
        )
        .is_err(),
        "COM3 mesh checkpoint must fail closed"
    );
    let imported = fixture.sim(&resume.replace(
        "restore_file=\"cp.chk\"",
        "restore_file=\"cp.chk\"\nrestart_mode=\"state_only\"",
    ));
    assert_eq!(imported.metadata().conserved_com_dof, 0);
    assert_eq!(imported.metadata().kinetic_dof, 9);
    // New direct Ewald protocol can conserve this zero-momentum import.
    let direct = fixture.sim(
        &resume
            .replace("method=\"spme\"", "method=\"ewald\"")
            .replace(
                "restore_file=\"cp.chk\"",
                "restore_file=\"cp.chk\"\nrestart_mode=\"state_only\"",
            ),
    );
    assert_eq!(direct.metadata().conserved_com_dof, 3);
    let constrained = fixture.sim(&format!(
        "{}\n[[constraints.bonds]]\natoms=[0,1]\nlength=1.8788294228055935",
        resume.replace(
            "restore_file=\"cp.chk\"",
            "restore_file=\"cp.chk\"\nrestart_mode=\"state_only\""
        )
    ));
    assert_eq!(constrained.metadata().conserved_com_dof, 0);
    assert_eq!(constrained.metadata().kinetic_dof, 8);
}
fn mixture(fractions: &[f64]) -> Vec<AlloySpecies> {
    fractions
        .iter()
        .enumerate()
        .map(|(i, &c)| AlloySpecies::neutral(&format!("S{i}"), i as f64 + 1.0, c))
        .collect()
}
#[test]
fn tiny_alloy_populations_are_bounded_deterministic_and_exact() {
    for (fractions, want) in [
        (
            vec![0.13, 0.13, 0.13, 0.13, 0.13, 0.35],
            vec![1, 1, 1, 0, 0, 1],
        ),
        (vec![0.0, 0.5, 0.5, 0.0], vec![0, 2, 2, 0]),
        (
            vec![0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.16],
            vec![1, 1, 1, 0, 0, 0, 1],
        ),
    ] {
        let species = mixture(&fractions);
        let s = random_alloy(3.0, [1, 1, 1], &species, 42).unwrap();
        let repeated = random_alloy(3.0, [1, 1, 1], &species, 42).unwrap();
        assert_eq!(s.n(), 4);
        assert_eq!(s.species, repeated.species);
        let mut counts = vec![0; fractions.len()];
        for i in 0..s.n() {
            let sp = s.species[i] as usize;
            counts[sp] += 1;
            assert_eq!(s.masses[i], sp as f64 + 1.0);
            assert_eq!(s.elements[i], format!("S{sp}"));
        }
        assert_eq!(counts, want);
    }
}
#[test]
fn exact_review_alloy_constructs_and_runs_common_config() {
    let raw = include_str!("fixtures/final_alloy.toml");
    let c = InputConfig::from_str(raw).unwrap();
    let mut s = Simulation::from_config_without_outputs(&c, Path::new(".")).unwrap();
    assert_eq!(s.system().n(), 4);
    assert!(s.system().species.iter().all(|&sp| sp < 6));
    assert_eq!(s.run().unwrap().steps, 2);
}
#[test]
fn alloy_invalid_input_and_capacity_return_errors() {
    for bad in [vec![f64::NAN], vec![-0.1, 1.1], vec![0.4, 0.4]] {
        assert!(random_alloy(3.0, [1, 1, 1], &mixture(&bad), 42).is_err());
    }
    for replicas in [[0, 1, 1], [u32::MAX, u32::MAX, u32::MAX]] {
        assert!(random_alloy(3.0, replicas, &mixture(&[1.0]), 42).is_err());
    }
}

#[test]
fn charged_mesh_csvr_and_nhc_consume_the_full_nine_mode_target() {
    use modena_core::{Integrator, NoseHooverChain, Thermostat, VelocityRescaling, VelocityVerlet};
    let fixture = Fixture::new();
    for nhc in [false, true] {
        let mut c = config("spme", if nhc { NHC } else { VV }).replace("steps=10", "steps=1");
        if !nhc {
            c.push_str("\n[thermostat]\ntype=\"velocity_rescaling\"\ntarget_temperature=1.0\ntau=0.1\nseed=231");
        }
        let mut runner = fixture.sim(&c);
        // The independent kinetic dimension is 3 * 3 atoms = 9. Feed that
        // literal dimension to the real controller outside runner setup.
        let mut reference = runner.system().clone();
        reference.set_conserved_com_dof(0).unwrap();
        let model = SpmeForceModel::new(0.8, 4.2, 4.8).unwrap();
        let mut neighbors = NeighborList::with_skin(4.2, 0.3);
        if nhc {
            NoseHooverChain::new(1.0, 0.1, 3)
                .step(&mut reference, &model, &mut neighbors, 0.001)
                .unwrap();
        } else {
            VelocityVerlet::new()
                .step(&mut reference, &model, &mut neighbors, 0.001)
                .unwrap();
            VelocityRescaling::new(1.0, 0.1, 231).apply(&mut reference, 0.001);
        }
        runner.run().unwrap();
        assert_eq!(
            runner.system().velocities,
            reference.velocities,
            "controller target is not nine modes; nhc={nhc}"
        );
    }
}
#[test]
fn charged_rigid_water_uses_six_modes_per_molecule_without_com_removal() {
    let c = r#"units="real"
[structure]
type="water"
n_waters=2
density=0.05
[model]
type="lennard_jones"
cutoff=2.5
[[model.types]]
sigma=1.0
epsilon=0.1
[[model.types]]
sigma=1.0
epsilon=0.1
[integrator]
type="velocity_verlet"
dt=0.00000001
[electrostatics]
method="spme"
kappa=0.8
real_cutoff=4.2
kspace_cutoff=4.8
[run]
steps=2
initial_temperature=1.0
seed=51
"#;
    let mut s =
        Simulation::from_config_without_outputs(&InputConfig::from_str(c).unwrap(), Path::new("."))
            .unwrap();
    assert_eq!(s.metadata().constraint_count, 6);
    assert_eq!(s.metadata().kinetic_dof, 12);
    assert_eq!(s.metadata().conserved_com_dof, 0);
    assert!((s.system().kinetic_energy() - 6.0 * 0.001987204259 / 0.1).abs() < 1e-12);
    s.run().unwrap();
}
