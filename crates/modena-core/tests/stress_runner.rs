use modena_core::barostat::{Barostat, BerendsenBarostat};
use modena_core::force_model::Virial;
#[allow(dead_code)]
#[path = "support/stress.rs"]
mod support;

#[test]
fn anisotropic_scalar_call_fails_before_mutation() {
    let mut system = support::atoms(&[[1., 2., 3.]], &[1.]);
    let before = system.clone();
    let mut controller = BerendsenBarostat::new(1., 1., 0.1).anisotropic();
    let error = controller
        .apply(&mut system, &Virial { w: 9. }, 0.01)
        .unwrap_err();
    assert!(error.to_string().contains("tensor pressure"));
    assert_eq!(system.positions, before.positions);
    assert_eq!(system.box_size, before.box_size);
}

use modena_core::barostat::{NoneBarostat, ParrinelloRahman};
use modena_core::stress::{
    endpoint_pressure, MissingTensor, MissingTensorReason, ModelEvaluation, PressureObservation,
    SymTensor6, TensorAvailability, ValidationStatus,
};
use modena_core::{InputConfig, RunEvent, Simulation};
use std::fs;

const REAL: &str = r#"
units = "real"
[structure]
type = "fcc"
lattice_constant = 4.0
replicas = [2,2,2]
mass = 10.0
[model]
type = "lennard_jones"
sigma = 2.0
epsilon = 3.0
cutoff = 3.9
[[constraints.bonds]]
atoms = [0,1]
length = 2.8
[integrator]
type = "velocity_verlet"
dt = 0.0003
[run]
steps = 3
initial_temperature = 120.0
log_every = 2
pressure_jsonl = "sub/pressure.jsonl"
log_file = "legacy.csv"
"#;

#[test]
fn experimental_diagonal_controls_use_the_three_endpoint_components_and_reject_invalid_states_atomically(
) {
    let source = support::atoms(&[[1., 2., 3.], [2., 3., 4.]], &[1., 1.]);
    let v = source.box_size.volume();
    let p = endpoint_pressure(
        &source,
        &ModelEvaluation {
            energy: 0.,
            scalar_virial: 6. * v,
            tensor: TensorAvailability::Complete(SymTensor6([
                v,
                2. * v,
                3. * v,
                0.4 * v,
                -0.2 * v,
                0.1 * v,
            ])),
        },
        None,
    )
    .unwrap();
    let dt: f64 = 0.001;
    for (kind, mut controller) in [
        (
            "berendsen",
            Box::new(BerendsenBarostat::new(0., 2., 0.2).anisotropic()) as Box<dyn Barostat>,
        ),
        (
            "piston",
            Box::new(ParrinelloRahman::new(0., 2., 1.).anisotropic()) as Box<dyn Barostat>,
        ),
    ] {
        let mut system = source.clone();
        let state = controller.checkpoint_state().unwrap();
        let mut missing = p.clone();
        missing.total = TensorAvailability::Missing(vec![MissingTensor {
            component: "model".into(),
            reason: MissingTensorReason::LegacyForceModel,
        }]);
        let mut wrong_volume = p.clone();
        wrong_volume.volume *= 2.;
        let mut nonfinite = p.clone();
        nonfinite.total = TensorAvailability::Complete(SymTensor6([0., 0., 0., 0., 0., f64::NAN]));
        for bad in [missing, wrong_volume, nonfinite] {
            assert!(controller
                .apply_with_pressure(&mut system, &bad, dt)
                .is_err());
        }
        system.set_constraint_count(1).unwrap();
        assert!(controller.apply_with_pressure(&mut system, &p, dt).is_err());
        system.set_constraint_count(0).unwrap();
        assert!(controller
            .apply(&mut system, &Virial { w: 9. }, dt)
            .is_err());
        assert_eq!(system.positions, source.positions);
        assert_eq!(system.velocities, source.velocities);
        assert_eq!(system.box_size, source.box_size);
        assert_eq!(state, controller.checkpoint_state().unwrap());
        controller.apply_with_pressure(&mut system, &p, dt).unwrap();
        for (k, diagonal) in [1., 2., 3.].into_iter().enumerate() {
            let mu = if kind == "berendsen" {
                (1. + 0.2 * dt / 2. * diagonal).cbrt()
            } else {
                (0.5 * dt * dt * 3. * v * diagonal / (source.dof() as f64 * 4.)).exp()
            };
            support::near(
                system.box_size.size.0[k],
                source.box_size.size.0[k] * mu,
                1e-14,
            );
            support::near(system.pos(1).0[k], source.pos(1).0[k] * mu, 1e-14);
        }
        println!(
            "{kind}: initial diagonal = [1,2,3], final box = {:?}",
            system.box_size.size
        );
    }
    assert!(NoneBarostat
        .apply_with_pressure(&mut source.clone(), &p, dt)
        .is_err());
}

#[test]
fn real_jsonl_and_events_use_display_units_and_the_legacy_csv_cadence() {
    let dir = tempfile::tempdir().unwrap();
    let cfg = InputConfig::from_str(REAL).unwrap();
    let factors = cfg.unit_factors.unwrap();
    let mut sim = Simulation::from_config(&cfg, dir.path()).unwrap();
    let initial = sim.pressure_observation().unwrap().clone();
    assert!(sim.interval_pressure().is_none());
    let mut samples = Vec::new();
    sim.run_streaming(None, |event| {
        if let RunEvent::Pressure {
            step,
            time,
            endpoint,
            interval,
        } = event
        {
            samples.push((*step, *time, endpoint.clone(), interval.clone()));
        }
    })
    .unwrap();
    let rows = fs::read_to_string(dir.path().join("sub/pressure.jsonl")).unwrap();
    let json: Vec<serde_json::Value> = rows
        .lines()
        .map(|s| serde_json::from_str(s).unwrap())
        .collect();
    assert_eq!(json.len(), 2);
    assert_eq!(samples.len(), 2);
    for ((step, time, e, i), row) in samples.iter().zip(&json) {
        assert_eq!(row["step"], *step);
        assert_eq!(row["time"], *time);
        assert_eq!(row["endpoint"], serde_json::to_value(e).unwrap());
        assert_eq!(row["interval"], serde_json::to_value(i).unwrap());
        assert_eq!(e.validation, ValidationStatus::Diagnostic);
    }
    assert_eq!((samples[0].0, samples[0].1), (0, 0.));
    assert_eq!(samples[1].0, 2);
    support::near(samples[1].1, 0.0006, 1e-15);
    assert!(json[0]["interval"].is_null());
    fn check_scaled(x: SymTensor6, y: SymTensor6, scale: f64) {
        for k in 0..6 {
            support::near(x.0[k] * scale, y.0[k], 1e-13);
        }
    }
    fn check_endpoint(x: &PressureObservation, y: &PressureObservation, e: f64, p: f64) {
        check_scaled(x.kinetic, y.kinetic, e);
        check_scaled(support::complete(&x.model), support::complete(&y.model), e);
        check_scaled(x.constraint, y.constraint, e);
        check_scaled(support::complete(&x.total), support::complete(&y.total), p);
        support::near(y.scalar.unwrap(), x.scalar.unwrap() * p, 1e-13);
    }
    check_endpoint(&initial, &samples[0].2, factors.epsilon, factors.pressure);
    support::near(samples[0].2.volume, 512., 1e-13);
    // End precisely on the second logged sample to compare the internal getters.
    let mut two = cfg.clone();
    two.run.steps = 2;
    let mut sim2 = Simulation::from_config_without_outputs(&two, dir.path()).unwrap();
    let legacy = sim2.run().unwrap().final_energies;
    check_endpoint(
        sim2.pressure_observation().unwrap(),
        &samples[1].2,
        factors.epsilon,
        factors.pressure,
    );
    let x = sim2.interval_pressure().unwrap();
    let y = samples[1].3.as_ref().unwrap();
    assert!(x.constraint_average.0.iter().any(|c| c.abs() > 1e-5));
    check_scaled(x.kinetic_average, y.kinetic_average, factors.epsilon);
    check_scaled(x.constraint_average, y.constraint_average, factors.epsilon);
    check_scaled(
        support::complete(&x.model_average),
        support::complete(&y.model_average),
        factors.epsilon,
    );
    check_scaled(
        support::complete(&x.total),
        support::complete(&y.total),
        factors.pressure,
    );
    support::near(y.volume, x.volume * factors.sigma.powi(3), 1e-13);
    support::near(
        y.scalar.unwrap(),
        x.scalar.unwrap() * factors.pressure,
        1e-13,
    );
    support::near(y.dt, 0.0003, 1e-15);
    let csv = fs::read_to_string(dir.path().join("legacy.csv")).unwrap();
    assert_eq!(
        csv.lines().next().unwrap(),
        "step,time,kinetic,potential,total,temperature,pressure,virial"
    );
    support::near(
        csv.lines()
            .nth(2)
            .unwrap()
            .split(',')
            .nth(6)
            .unwrap()
            .parse()
            .unwrap(),
        legacy.pressure,
        1e-9,
    );
    println!(
        "JSONL {} bytes / {} rows; CSV {} bytes; step 2 endpoint total {:?}; interval total {:?}",
        rows.len(),
        json.len(),
        csv.len(),
        samples[1].2.total,
        y.total
    );
}

#[test]
fn thermostat_recomputes_rotor_endpoint_and_leaves_pre_controller_interval_and_forces() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("rotor.xyz"),"2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nX 5 5 5 0 1 0\nX 6 7 7 0 2 0\n").unwrap();
    let config = r#"
[structure]
type = "xyz"
file = "rotor.xyz"
box_size = [20.0,20.0,20.0]
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.0
[[constraints.bonds]]
atoms = [0,1]
length = 3.0
[integrator]
type = "velocity_verlet"
dt = 0.001
[run]
steps = 1
initial_temperature = 1.0
seed = 42
"#;
    let mut reference =
        Simulation::from_config(&InputConfig::from_str(config).unwrap(), dir.path()).unwrap();
    reference.run().unwrap();
    let target = 4. * reference.system().temperature();
    let hot = format!("{config}\n[thermostat]\ntype = \"berendsen\"\ntau = 0.001\ntarget_temperature = {target:.17}\n");
    let mut scaled =
        Simulation::from_config(&InputConfig::from_str(&hot).unwrap(), dir.path()).unwrap();
    let e = scaled.run().unwrap().final_energies;
    let a = reference.pressure_observation().unwrap();
    let b = scaled.pressure_observation().unwrap();
    assert_eq!(reference.system().forces, scaled.system().forces);
    assert_eq!(reference.system().positions, scaled.system().positions);
    assert_eq!(reference.interval_pressure(), scaled.interval_pressure());
    for k in 0..6 {
        support::near(a.kinetic.0[k] * 4., b.kinetic.0[k], 1e-12);
        support::near(a.constraint.0[k] * 4., b.constraint.0[k], 1e-12);
    }
    assert!(a.constraint.trace().abs() > 1.);
    support::near(e.pressure, 2. * e.kinetic / (3. * b.volume), 1e-13);
    assert!((b.scalar.unwrap() - e.pressure).abs() > 1e-5);
    println!(
        "velocity x2: K {:?} -> {:?}, reaction {:?} -> {:?}, legacy {}, endpoint {}",
        a.kinetic,
        b.kinetic,
        a.constraint,
        b.constraint,
        e.pressure,
        b.scalar.unwrap()
    );
}

#[test]
fn pressure_output_path_failure_and_output_disabled_mode() {
    let cfg = InputConfig::from_str(REAL).unwrap();
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("sub"), "regular file").unwrap();
    assert!(Simulation::from_config(&cfg, dir.path()).is_err());
    Simulation::from_config_without_outputs(&cfg, dir.path())
        .unwrap()
        .run()
        .unwrap();
    assert!(!dir.path().join("legacy.csv").exists());
}

#[test]
fn pressure_output_collision_precedes_other_writers() {
    for (other, alias) in [
        ("log_file", "same.out"),
        ("trajectory_file", "nested/../same.out"),
        ("restart_file", "./same.out"),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let mut cfg = InputConfig::from_str(REAL).unwrap();
        cfg.run.log_file = Some("untouched.csv".into());
        cfg.run.trajectory_file = Some("untouched.xyz".into());
        cfg.run.pressure_jsonl = Some(alias.into());
        match other {
            "log_file" => cfg.run.log_file = Some("same.out".into()),
            "trajectory_file" => cfg.run.trajectory_file = Some("same.out".into()),
            _ => cfg.run.restart_file = Some("same.out".into()),
        }
        Simulation::from_config_without_outputs(&cfg, dir.path()).unwrap();
        let err = Simulation::from_config(&cfg, dir.path())
            .err()
            .expect("aliased output accepted");
        assert!(
            matches!(err,modena_core::ModenaError::Config(ref reason) if reason.contains("pressure_jsonl") && reason.contains(other)),
            "{err}"
        );
        assert!(!dir.path().join("untouched.csv").exists());
        assert!(!dir.path().join("untouched.xyz").exists());
        if let Ok(m) = fs::metadata(dir.path().join("same.out")) {
            assert_eq!(m.len(), 0);
        }
    }
}

#[test]
fn pressure_output_collision_preserves_restore_checkpoint() {
    let dir = tempfile::tempdir().unwrap();
    let mut cfg = InputConfig::from_str(REAL).unwrap();
    let checkpoint = dir.path().join("good.chk");
    Simulation::from_config_without_outputs(&cfg, dir.path())
        .unwrap()
        .write_checkpoint(&checkpoint)
        .unwrap();
    let intact = fs::read(&checkpoint).unwrap();
    cfg.run.restore_file = Some("good.chk".into());
    cfg.run.pressure_jsonl = Some(checkpoint.to_str().unwrap().into());
    Simulation::from_config_without_outputs(&cfg, dir.path()).unwrap();
    let rejected = Simulation::from_config(&cfg, dir.path()).err();
    assert_eq!(
        fs::read(&checkpoint).unwrap().len(),
        intact.len(),
        "good checkpoint truncated"
    );
    assert_eq!(fs::read(&checkpoint).unwrap(), intact);
    assert!(
        matches!(rejected,Some(modena_core::ModenaError::Config(reason)) if reason.contains("restore_file"))
    );
    assert!(!dir.path().join("legacy.csv").exists());
}

#[cfg(unix)]
#[test]
fn pressure_output_collision_rejects_filesystem_aliases() {
    use std::os::unix::fs::symlink;
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("keep"), b"original").unwrap();
    fs::hard_link(dir.path().join("keep"), dir.path().join("hard")).unwrap();
    symlink("keep", dir.path().join("soft")).unwrap();
    symlink(".", dir.path().join("aliasdir")).unwrap();
    symlink("uncreated", dir.path().join("dangling")).unwrap();
    for (pressure, log) in [
        ("hard", "keep"),
        ("soft", "keep"),
        ("aliasdir/sub/../new", "new"),
        ("dangling", "uncreated"),
    ] {
        let mut cfg = InputConfig::from_str(REAL).unwrap();
        cfg.run.pressure_jsonl = Some(pressure.into());
        cfg.run.log_file = Some(log.into());
        let err = Simulation::from_config(&cfg, dir.path())
            .err()
            .expect("aliased output accepted");
        assert!(
            matches!(err,modena_core::ModenaError::Config(ref reason) if reason.contains("log_file")),
            "{err}"
        );
        assert_eq!(fs::read(dir.path().join("keep")).unwrap(), b"original");
        for fresh in ["new", "uncreated"] {
            if let Ok(m) = fs::metadata(dir.path().join(fresh)) {
                assert_eq!(m.len(), 0);
            }
        }
        assert!(fs::symlink_metadata(dir.path().join("dangling"))
            .unwrap()
            .is_symlink());
    }
}

#[test]
fn pressure_output_collision_detects_new_case_alias_on_insensitive_volume() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("VolumeCase"), b"volume probe").unwrap();
    if !dir.path().join("volumecase").exists() {
        println!("case-sensitive test volume");
        return;
    }
    let mut cfg = InputConfig::from_str(REAL).unwrap();
    cfg.run.pressure_jsonl = Some("Sample.jsonl".into());
    cfg.run.log_file = Some("sample.jsonl".into());
    cfg.run.trajectory_file = Some("untouched.xyz".into());
    let err = Simulation::from_config(&cfg, dir.path())
        .err()
        .expect("new case alias accepted");
    assert!(
        matches!(err,modena_core::ModenaError::Config(ref reason) if reason.contains("log_file")),
        "{err}"
    );
    assert_eq!(
        fs::metadata(dir.path().join("Sample.jsonl"))
            .expect("public probe must not be unlinked")
            .len(),
        0
    );
    assert!(!dir.path().join("untouched.xyz").exists());
    println!("case-insensitive test volume; new case alias rejected, public probe retained empty");
}

#[cfg(unix)]
fn assert_invalid_pressure_spelling(pressure: &str) {
    let dir = tempfile::tempdir().unwrap();
    let mut cfg = InputConfig::from_str(REAL).unwrap();
    cfg.run.pressure_jsonl = Some(pressure.into());
    let saved = dir.path().join("saved");
    let marker = dir.path().join("marker");
    fs::write(&saved, b"preserve these saved bytes").unwrap();
    fs::write(&marker, b"regular intermediate file").unwrap();
    assert_eq!(
        fs::metadata(dir.path().join(pressure)).unwrap_err().kind(),
        std::io::ErrorKind::NotADirectory
    );
    let result = Simulation::from_config(&cfg, dir.path());
    assert_eq!(
        fs::read(&saved).unwrap(),
        b"preserve these saved bytes",
        "truncated through {pressure}"
    );
    assert_eq!(fs::read(&marker).unwrap(), b"regular intermediate file");
    assert!(
        matches!(result,Err(modena_core::ModenaError::Io(e)) if e.kind()==std::io::ErrorKind::NotADirectory)
    );
    assert!(!dir.path().join("legacy.csv").exists());
}
#[cfg(unix)]
#[test]
fn pressure_output_preserves_intermediate_regular_file_semantics() {
    assert_invalid_pressure_spelling("marker/../saved");
}
#[cfg(unix)]
#[test]
fn pressure_output_preserves_trailing_directory_separator() {
    assert_invalid_pressure_spelling("saved/");
}
#[cfg(unix)]
#[test]
fn pressure_output_preserves_trailing_directory_dot() {
    assert_invalid_pressure_spelling("saved/.");
}

#[test]
fn pressure_output_collision_after_creating_counterpart_parent() {
    let dir = tempfile::tempdir().unwrap();
    let mut cfg = InputConfig::from_str(REAL).unwrap();
    cfg.run.pressure_jsonl = Some("saved".into());
    cfg.run.log_file = Some("new-dir/../saved".into());
    fs::write(dir.path().join("saved"), b"old data").unwrap();
    let rejected = Simulation::from_config(&cfg, dir.path()).err();
    assert!(
        matches!(rejected,Some(modena_core::ModenaError::Config(reason)) if reason.contains("log_file"))
    );
    assert_eq!(fs::read(dir.path().join("saved")).unwrap(), b"old data");
}
