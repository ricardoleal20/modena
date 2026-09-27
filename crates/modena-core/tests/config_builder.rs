//! The Rust "modena-language" (`Simulation::builder()`, M-RL-0): fluent config
//! must produce the same validated `InputConfig` as a hand-written TOML file.

use modena_core::{ConservationVerdict, Simulation};

#[test]
fn builder_runs_nve_and_conserves() {
    let mut sim = Simulation::builder()
        .fcc(1.6796, [2, 2, 2], "Ar")
        .lj(1.0, 1.0, 2.5)
        .velocity_verlet(0.005)
        .steps(400)
        .temperature(0.8)
        .seed(42)
        .build()
        .expect("builder config is valid");
    let out = sim.run().expect("run completes");
    assert_ne!(
        out.conservation.verdict,
        ConservationVerdict::Drifting,
        "NVE conserved"
    );
    assert!(out.final_energies.potential.is_finite());
}

#[test]
fn builder_to_config_matches_toml_semantics() {
    let a = Simulation::builder()
        .fcc(1.6796, [2, 2, 2], "Ar")
        .lj(1.0, 1.0, 2.5)
        .velocity_verlet(0.005)
        .steps(100)
        .to_config()
        .expect("valid");
    let toml = r#"
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
steps = 100
"#;
    let b = modena_core::InputConfig::from_str(toml).expect("parses");
    assert_eq!(a.run.steps, b.run.steps);
    assert_eq!(a.run.neighbor_skin, b.run.neighbor_skin);
    match (a.structure, b.structure) {
        (
            modena_core::StructureSpec::Fcc {
                lattice_constant: la,
                ..
            },
            modena_core::StructureSpec::Fcc {
                lattice_constant: lb,
                ..
            },
        ) => {
            assert!((la - lb).abs() < 1e-12);
        }
        _ => panic!("both must be fcc"),
    }
}

#[test]
fn builder_requires_core_sections() {
    // Structure-only is invalid.
    let err = Simulation::builder()
        .fcc(1.6796, [1, 1, 1], "Ar")
        .to_config()
        .unwrap_err()
        .to_string();
    assert!(err.contains("missing model"), "got: {err}");
    // Model + structure but no integrator is invalid too.
    let err = Simulation::builder()
        .fcc(1.6796, [1, 1, 1], "Ar")
        .lj(1.0, 1.0, 2.5)
        .to_config()
        .unwrap_err()
        .to_string();
    assert!(err.contains("missing integrator"), "got: {err}");
}
