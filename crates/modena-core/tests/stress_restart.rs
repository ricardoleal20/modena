use modena_core::io::restart::{load, load_with_stress, save, save_with_stress, StressCheckpoint};
use modena_core::stress::{MissingTensorReason, TensorAvailability, ValidationStatus};
use modena_core::RunEvent;
use modena_core::{InputConfig, Simulation};
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;
use std::sync::atomic::{AtomicUsize, Ordering};

fn unique_test_path(label: &str) -> PathBuf {
    static SEQ: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "modena_{label}_{}_{}",
        std::process::id(),
        SEQ.fetch_add(1, Ordering::Relaxed)
    ))
}

fn rotor(dir: &Path, steps: usize, restore: Option<&str>) -> Simulation {
    Simulation::from_config_without_outputs(&rotor_config(dir, steps, restore), dir).unwrap()
}

fn rotor_config(dir: &Path, steps: usize, restore: Option<&str>) -> InputConfig {
    std::fs::write(dir.join("rotor.xyz"), "2\nrotor\nX 1 1 1\nX 2.2 1 1\n").unwrap();
    let restore = restore
        .map(|path| format!("restore_file = {path:?}"))
        .unwrap_or_default();
    let input = format!(
        r#"
[structure]
type = "xyz"
file = "rotor.xyz"
box_size = [8.0,8.0,8.0]
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "langevin_baoab"
gamma = 0.7
target_temperature = 1.0
dt = 0.003
seed = 71
[[constraints.bonds]]
atoms = [0,1]
length = 1.2
[run]
steps = {steps}
initial_temperature = 1.0
seed = 51
log_every = 1
{restore}
"#
    );
    InputConfig::from_str(&input).unwrap()
}

fn make_stress_checkpoint_after_rotor_step() -> StressCheckpoint {
    let dir = tempfile::tempdir().unwrap();
    let mut sim = rotor(dir.path(), 1, None);
    sim.run().unwrap();
    let file = dir.path().join("rotor.chk");
    sim.write_checkpoint(&file).unwrap();
    load_with_stress(file).unwrap()
}

#[test]
fn legacy_v2_reader_rejects_v3_explicitly() {
    let cp = make_stress_checkpoint_after_rotor_step();
    let path = unique_test_path("b07-v3");
    save_with_stress(&path, &cp).unwrap();
    let restored = load_with_stress(&path).unwrap();
    assert_eq!(restored.mechanical.positions, cp.mechanical.positions);
    assert_eq!(restored.stress, cp.stress);
    assert!(load(&path)
        .unwrap_err()
        .to_string()
        .contains("load_with_stress"));
    save(&path, &cp.mechanical).unwrap();
    assert!(load_with_stress(&path).unwrap().stress.is_none());
    std::fs::remove_file(path).unwrap();
}

fn saved(sim: &Simulation, file: &Path) -> Vec<u8> {
    sim.write_checkpoint(file).unwrap();
    std::fs::read(file).unwrap()
}

fn legacy(tensor: &TensorAvailability) -> bool {
    matches!(tensor, TensorAvailability::Missing(why) if why.iter().any(|m| m.reason == MissingTensorReason::LegacyRestart))
}

#[test]
fn constrained_baoab_seventeen_plus_twenty_three_exact_continuation() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let mut direct = rotor(dir, 23, None);
    let stop = AtomicBool::new(false);
    let phase = direct
        .run_streaming(Some(&stop), |e| {
            if matches!(e, RunEvent::Pressure { step: 17, .. }) {
                stop.store(true, Ordering::Release);
            }
        })
        .unwrap();
    assert_eq!(phase.steps, 17);
    let seam = saved(&direct, &dir.join("step17.chk"));
    let cp = load_with_stress(dir.join("step17.chk")).unwrap();
    let stress = cp.stress.as_ref().unwrap();
    assert_eq!(stress.last_stages.len(), 7);
    assert!(stress.last_stages.iter().any(|r| !r.impulses.is_empty()));
    assert!(stress.last_bath.is_some());
    assert_ne!(
        cp.mechanical.continuation.as_ref().unwrap().forces,
        vec![0.; 6]
    );
    let mut resumed = rotor(dir, 23, Some("step17.chk"));
    assert_eq!(saved(&resumed, &dir.join("loaded17.chk")), seam);
    save(dir.join("v2.chk"), &cp.mechanical).unwrap();
    let mut old = rotor(dir, 23, Some("v2.chk"));
    assert!(legacy(&old.pressure_observation().unwrap().model));
    assert!(old.interval_pressure().is_none());
    let mut pressure = Vec::new();
    direct
        .run_streaming(None, |e| {
            if let RunEvent::Pressure {
                endpoint, interval, ..
            } = e
            {
                pressure.push((endpoint.clone(), interval.clone()));
            }
        })
        .unwrap();
    let evidence = std::env::var_os("MODENA_RESTART_TEST_EVIDENCE");
    let retain = |name: &str, value: &serde_json::Value| {
        if let Some(dir) = &evidence {
            std::fs::write(
                Path::new(dir).join(name),
                serde_json::to_vec(value).unwrap(),
            )
            .unwrap();
        }
    };
    retain("direct-observations.json", &serde_json::json!(pressure));
    for (sim, from_v2, out) in [
        (&mut resumed, false, "loaded40.chk"),
        (&mut old, true, "legacy40.chk"),
    ] {
        let mut rows = Vec::new();
        sim.run_streaming(None, |e| {
            if let RunEvent::Pressure {
                step,
                endpoint,
                interval,
                ..
            } = e
            {
                if evidence.is_some() {
                    rows.push(serde_json::json!([endpoint, interval]));
                }
                let (expected, history) = &pressure[*step - 17];
                if !from_v2 || *step > 17 {
                    assert_eq!(endpoint, expected);
                }
                match (from_v2, *step) {
                    (true, 17) => assert!(interval.is_none()),
                    (true, 18) => assert!(legacy(&interval.as_ref().unwrap().total)),
                    _ => {
                        assert_eq!(interval, history);
                        assert!(matches!(
                            interval.as_ref().unwrap().total,
                            TensorAvailability::Complete(_)
                        ));
                    }
                }
            }
        })
        .unwrap();
        retain(
            if from_v2 {
                "legacy-observations.json"
            } else {
                "loaded-observations.json"
            },
            &serde_json::json!(rows),
        );
        assert_eq!(sim.current_step(), 40);
        assert_eq!(sim.system().positions, direct.system().positions);
        assert_eq!(sim.system().velocities, direct.system().velocities);
        assert_eq!(sim.system().forces, direct.system().forces);
        // Compare all bytes, including signs of zero, ChaCha cursor, thermostat,
        // barostat, physical signature, model cache, endpoint, stages and OU bath.
        assert_eq!(
            saved(sim, &dir.join(out)),
            saved(&direct, &dir.join("direct40.chk"))
        );
    }
    let legacy = std::fs::read(dir.join("v2.chk")).unwrap();
    let mechanical: serde_json::Value =
        serde_json::from_slice(legacy.splitn(2, |b| *b == b'\n').nth(1).unwrap()).unwrap();
    assert_eq!(mechanical["continuation"].as_object().unwrap().len(), 14);
    println!("17+23 steps, 24 observations per continuation, exact serialized seam and final equality; v2 {} bytes, v3 {} bytes, step 17 RNG word position {}",legacy.len(),seam.len(),mechanical["continuation"]["integrator"]["rng"]["word_pos"]);
    if let Some(destination) = std::env::var_os("MODENA_RESTART_TEST_EVIDENCE") {
        for name in [
            "step17.chk",
            "loaded17.chk",
            "v2.chk",
            "direct40.chk",
            "loaded40.chk",
            "legacy40.chk",
        ] {
            std::fs::copy(dir.join(name), Path::new(&destination).join(name)).unwrap();
        }
    }
}

#[test]
fn corrupted_or_partial_stress_is_rejected_and_failed_writes_preserve_bytes() {
    use serde_json::json;
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("state.chk");
    let cp = make_stress_checkpoint_after_rotor_step();
    save_with_stress(&path, &cp).unwrap();
    let original = std::fs::read(&path).unwrap();
    let mut bad = cp.clone();
    bad.stress
        .as_mut()
        .unwrap()
        .last_bath
        .as_mut()
        .unwrap()
        .reaction[0]
        .0[1] = f64::NAN;
    assert!(save_with_stress(&path, &bad).is_err());
    assert_eq!(std::fs::read(&path).unwrap(), original);
    let value: serde_json::Value =
        serde_json::from_slice(original.splitn(2, |b| *b == b'\n').nth(1).unwrap()).unwrap();
    for (key, replacement) in [
        ("/stress/estimator_version", json!(2)),
        ("/stress/step", json!(8)),
        ("/stress/dt", json!(0)),
        ("/stress/endpoint/volume", json!(1)),
        ("/stress/model/tensor/value", json!([1, 2])),
        ("/stress/model/scalar_virial", json!(23)),
        ("/stress/endpoint/scalar", json!(11)),
        ("/stress/interval/scalar", json!(11)),
        ("/stress/interval/dt", json!(5)),
        ("/stress/last_bath", json!(null)),
        ("/stress/last_bath/reaction", json!([])),
        ("/stress/last_stages/0/stage/kind", json!("vv_final")),
        ("/stress/last_stages/0/tensor_impulse/0", json!(93)),
        ("/stress/last_stages", json!([])),
        ("/stress/interval", json!(null)),
        ("/mechanical/continuation/forces", json!([])),
        (
            "/stress/physical_signature/pressure_estimator/endpoint_timing",
            json!("before"),
        ),
        ("/stress/physical_signature/model", json!({"pin":"changed"})),
        (
            "/stress/physical_signature/barostat",
            json!({"type":"unknown"}),
        ),
        ("/stress/solver/tolerance", json!(1e-12)),
    ] {
        let mut bad = value.clone();
        *bad.pointer_mut(key).unwrap() = replacement;
        std::fs::write(&path, format!("MODENA_RESTART 3\n{bad}\n")).unwrap();
        assert!(load_with_stress(&path).is_err(), "accepted malformed {key}");
    }
    for node in [
        "",
        "/mechanical",
        "/mechanical/continuation",
        "/stress",
        "/stress/model",
        "/stress/endpoint",
        "/stress/last_stages/0",
        "/stress/last_bath",
    ] {
        let mut bad = value.clone();
        bad.pointer_mut(node).unwrap()["unknown"] = json!(true);
        std::fs::write(&path, format!("MODENA_RESTART 3\n{bad}\n")).unwrap();
        assert!(load_with_stress(&path).is_err(), "unknown field at {node}");
    }
    let mut claimed = value;
    for key in ["/stress/endpoint/validation", "/stress/interval/validation"] {
        *claimed.pointer_mut(key).unwrap() =
            json!({"kind":"validated","value":{"pub_record":"unverified"}});
    }
    std::fs::write(&path, format!("MODENA_RESTART 3\n{claimed}\n")).unwrap();
    let s = load_with_stress(&path).unwrap().stress.unwrap();
    assert_eq!(s.endpoint.validation, ValidationStatus::Diagnostic);
    assert_eq!(s.interval.unwrap().validation, ValidationStatus::Diagnostic);
}

#[test]
fn runner_rejects_changed_solver_and_controller_parameters() {
    use serde_json::json;
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("changed.chk");
    let cp = make_stress_checkpoint_after_rotor_step();
    let config = rotor_config(dir.path(), 1, Some("changed.chk"));
    let mut solver = cp.clone();
    let stress = solver.stress.as_mut().unwrap();
    stress.solver.tolerance = 1e-12;
    stress.physical_signature["pressure_estimator"]["solver"]["tolerance"] = json!(1e-12);
    save_with_stress(&path, &solver).unwrap();
    assert!(Simulation::from_config_without_outputs(&config, dir.path())
        .err()
        .unwrap()
        .to_string()
        .contains("solver"));
    save_with_stress(&path, &cp).unwrap();
    let mut changed = config.clone();
    if let modena_core::io::input::ModelSpec::LennardJones { cutoff, .. } = &mut changed.model {
        *cutoff = 1.9;
    } else {
        unreachable!()
    }
    assert!(Simulation::from_config_without_outputs(&changed, dir.path()).is_err());
    let mut bad = cp;
    bad.mechanical.continuation.as_mut().unwrap().integrator["gamma"] = json!(0.11);
    save_with_stress(&path, &bad).unwrap();
    assert!(Simulation::from_config_without_outputs(&config, dir.path()).is_err());
}

#[test]
fn post_controller_endpoint_and_pre_move_interval_survive_restore() {
    use serde_json::json;
    let dir = tempfile::tempdir().unwrap();
    let mut config = rotor_config(dir.path(), 1, None);
    config.constraints = None;
    config.integrator =
        serde_json::from_value(json!({"type":"velocity_verlet","dt":0.003})).unwrap();
    config.thermostat = serde_json::from_value(
        json!({"type":"velocity_rescaling","target_temperature":1.0,"tau":0.1,"seed":29}),
    )
    .unwrap();
    config.barostat = serde_json::from_value(
        json!({"type":"berendsen","target_pressure":0.25,"tau":2.0,"compressibility":0.1}),
    )
    .unwrap();
    let mut direct = Simulation::from_config_without_outputs(&config, dir.path()).unwrap();
    let old_volume = direct.system().box_size.volume();
    direct.run().unwrap();
    let interval = direct.interval_pressure().unwrap();
    assert_eq!(interval.volume, old_volume);
    assert!(
        matches!(&interval.total, TensorAvailability::Missing(why) if why.iter().any(|r|r.reason==MissingTensorReason::PostStepTransform))
    );
    assert_ne!(direct.pressure_observation().unwrap().volume, old_volume);
    let seam = saved(&direct, &dir.path().join("moved.chk"));
    config.run.restore_file = Some("moved.chk".into());
    let mut restored = Simulation::from_config_without_outputs(&config, dir.path()).unwrap();
    assert_eq!(saved(&restored, &dir.path().join("restored.chk")), seam);
    direct.run().unwrap();
    restored.run().unwrap();
    assert_eq!(
        saved(&restored, &dir.path().join("restored.chk")),
        saved(&direct, &dir.path().join("moved.chk"))
    );
}
