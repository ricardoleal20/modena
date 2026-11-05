#![cfg(unix)]

use modena_core::{
    describe_ase_emt, AseEmtForceModel, AseEmtSpec, ForceModel, NeighborList, SimBox, System,
    UnitFactors, Virial,
};
use std::os::unix::fs::PermissionsExt;

fn spec(python: &str) -> AseEmtSpec {
    AseEmtSpec {
        python: python.into(),
        ase_version: "3.29.0".into(),
        model_sha256: "a".repeat(64),
        atomic_numbers: vec![29],
        composition: "total".into(),
        pbc: [true; 3],
        asap_cutoff: false,
        timeout_ms: 5000,
    }
}
fn fixture(evaluate: &str) -> (tempfile::TempDir, String) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("python");
    let script = format!(
        r##"#!/usr/bin/env python3
import sys,json,time
request=json.load(sys.stdin)
identity={{"protocol":"modena-ase-emt-v1","provider":"ase_emt","model":"EMT","task":"energy_forces_stress","composition":"total","units":{{"positions":"angstrom","energy":"eV","forces":"eV/angstrom","stress":"eV/angstrom^3"}},"capabilities":["energy","forces","stress","full_pbc","provider_neighbors","conservative","total_momentum"],"ase_version":"3.29.0","numpy_version":"2.5.4","model_sha256":"{}","asap_cutoff":False}}
if request['action']=='info':
 print(json.dumps({{"protocol":"modena-ase-emt-v1","action":"info","identity":identity}}))
else:
 {}
"##,
        "a".repeat(64),
        evaluate
    );
    std::fs::write(&path, script).unwrap();
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    (dir, path.to_string_lossy().into_owned())
}
fn system() -> System {
    let mut s = System::new(SimBox::cubic(4.0));
    s.resize(1);
    s.atomic_numbers = vec![29];
    s.forces = vec![2.0, 3.0, 4.0];
    s
}
const VALID: &str = "print(json.dumps({'protocol':'modena-ase-emt-v1','action':'evaluate','identity':identity,'evaluation':{'energy_ev':2.0,'forces_ev_per_angstrom':[[1.0,2.0,3.0]],'stress_ev_per_angstrom3':[[0.1,0.0,0.0],[0.0,0.2,0.0],[0.0,0.0,0.3]]}}))";
#[test]
fn rejects_unsupported_specs_before_launch() {
    for (field, value) in [
        ("ase_version", ""),
        ("model_sha256", "bad"),
        ("composition", "residual"),
    ] {
        let mut v = serde_json::to_value(spec("/does/not/exist")).unwrap();
        v[field] = value.into();
        let error = AseEmtForceModel::new(
            serde_json::from_value(v).unwrap(),
            UnitFactors::new(1.0, 1.0, 1.0),
        )
        .err()
        .unwrap()
        .to_string();
        assert!(error.contains(field), "{error}");
    }
    for z in [vec![], vec![0], vec![1], vec![6], vec![255]] {
        let mut s = spec("/missing");
        s.atomic_numbers = z;
        assert!(AseEmtForceModel::new(s, UnitFactors::new(1.0, 1.0, 1.0))
            .err()
            .unwrap()
            .to_string()
            .contains("atomic_numbers"));
    }
    let mut s = spec("/missing");
    s.pbc[1] = false;
    assert!(AseEmtForceModel::new(s, UnitFactors::new(1.0, 1.0, 1.0))
        .err()
        .unwrap()
        .to_string()
        .contains("pbc"));
    assert!(AseEmtForceModel::new(spec("/missing"), UnitFactors::new(0.0, 1.0, 1.0)).is_err());
}
#[test]
fn discovery_and_pins_are_checked() {
    let (_d, p) = fixture(VALID);
    let id = describe_ase_emt(&p, 5000).unwrap();
    assert_eq!(id.ase_version, "3.29.0");
    let mut s = spec(&p);
    s.model_sha256 = "b".repeat(64);
    assert!(AseEmtForceModel::new(s, UnitFactors::new(1.0, 1.0, 1.0)).is_err());
    let mut s = spec(&p);
    s.ase_version = "0.0".into();
    assert!(AseEmtForceModel::new(s, UnitFactors::new(1.0, 1.0, 1.0)).is_err());
}
#[test]
fn valid_evaluation_converts_and_accumulates_without_neighbor_cutoff() {
    let (_d, p) = fixture(VALID);
    let model = AseEmtForceModel::new(spec(&p), UnitFactors::new(2.0, 3.0, 1.0)).unwrap();
    let mut s = system();
    let e = model.evaluate(&s).unwrap();
    assert_eq!(e.energy_ev, 2.0);
    let mut w = Virial { w: 7.0 };
    let energy = model
        .compute(&mut s, &mut w, &mut NeighborList::new(1.0))
        .unwrap();
    let c = 1.602176634e-19 * 6.02214076e23 / 4184.0 / 3.0;
    assert!((energy - 2.0 * c).abs() < 1e-12);
    assert!((s.forces[0] - (2.0 + 2.0 * c)).abs() < 1e-12);
    assert!((w.w - (7.0 - 512.0 * 0.6 * c)).abs() < 1e-10);
    assert!(model.uses_provider_neighbors());
    assert_eq!(model.neighbor_cutoff(), None);
}
#[test]
fn failed_results_leave_accumulators_unchanged() {
    for body in [
        "sys.exit(7)", "time.sleep(30)", "print('not json')",
        "print('x'*65537, file=sys.stderr)", "print('x'*17000000)",
        "print(json.dumps({'protocol':'modena-ase-emt-v1','action':'evaluate','identity':identity,'evaluation':{'energy_ev':1.0}}))",
        "identity changed",
    ] {
        let body = if body.starts_with("identity") {
            format!("identity['numpy_version']='changed'; {VALID}")
        } else { body.into() };
        let (_dir, python) = fixture(&body);
        let mut spec = spec(&python);
        // Discovery starts another interpreter before evaluation. Leave startup
        // margin under concurrent suites; the blocking actor is six times longer.
        spec.timeout_ms = 5000;
        let model = AseEmtForceModel::new(spec, UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
        let mut system = system();
        let before = system.forces.clone();
        let mut virial = Virial { w: 7.0 };
        let start = std::time::Instant::now();
        let error = model.compute(&mut system, &mut virial, &mut NeighborList::new(1.0)).unwrap_err();
        if body == "time.sleep(30)" {
            assert!(error.to_string().contains("timeout after 5000 ms"), "{error}");
            assert!(start.elapsed() < std::time::Duration::from_secs(15), "blocking child was not stopped promptly");
        }
        assert_eq!(system.forces, before);
        assert_eq!(virial.w, 7.0);
    }
}
#[test]
fn rejects_wrong_atomic_identity_and_charges() {
    let (_d, p) = fixture(VALID);
    let m = AseEmtForceModel::new(spec(&p), UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
    let mut s = system();
    s.atomic_numbers[0] = 28;
    assert!(m.evaluate(&s).is_err());
    s.atomic_numbers[0] = 29;
    s.charges[0] = 1.0;
    assert!(m.evaluate(&s).is_err());
}
#[test]
fn installed_reference_provider_is_optional_and_explicit() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP installed ASE: set MODENA_ASE_PYTHON");
        return;
    };
    let id = describe_ase_emt(&python, 30000).unwrap();
    let mut config = spec(&python);
    config.model_sha256 = id.model_sha256;
    config.ase_version = id.ase_version;
    config.atomic_numbers = vec![29; 4];
    config.timeout_ms = 30000;
    let m = AseEmtForceModel::new(config, UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
    let mut s = System::new(SimBox::cubic(3.6));
    s.resize(4);
    s.atomic_numbers = vec![29; 4];
    s.positions = vec![0.03, 0.0, 0.0, 0.0, 1.8, 1.8, 1.8, 0.0, 1.8, 1.8, 1.8, 0.0];
    let a = m.evaluate(&s).unwrap();
    assert!(a.energy_ev.is_finite());
    s.positions[0] += 0.1;
    let b = m.evaluate(&s).unwrap();
    assert!((a.energy_ev - b.energy_ev).abs() > 1e-8);
}

#[test]
fn malformed_properties_and_conversion_overflow_are_atomic() {
    // Each mutation reaches the real deserialize/validation/accumulation path.
    for body in [
        VALID.replace("[[1.0,2.0,3.0]]", "[]"),
        VALID.replace("[[1.0,2.0,3.0]]", "[[1.0,2.0]]"),
        VALID.replace("2.0,'forces", "float('nan'),'forces"),
        VALID.replace("[[1.0,2.0,3.0]]", "[[float('inf'),2.0,3.0]]"),
        VALID.replace("[[1.0,2.0,3.0]]", "[[1e308,2.0,3.0]]"),
        VALID.replace("[0.1,0.0,0.0]", "[0.1,0.1,0.0]"),
        VALID.replace("[0.1,0.0,0.0]", "[0.1,0.0]"),
        VALID.replace("'energy_ev':2.0,", ""),
        VALID.replace("'forces_ev_per_angstrom':[[1.0,2.0,3.0]],", ""),
        VALID.replace("'stress_ev_per_angstrom3'", "'missing_stress'"),
        format!("identity['units']['energy']='kcal/mol'; {VALID}"),
        format!("identity['capabilities']=[]; {VALID}"),
        VALID.replace("'action':'evaluate'", "'action':'info'"),
        VALID.replace("'energy_ev':2.0", "'extra':1,'energy_ev':2.0"),
    ] {
        let (_dir, python) = fixture(&body);
        let model = AseEmtForceModel::new(spec(&python), UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
        let mut s = system();
        let before = s.forces.clone();
        let mut w = Virial { w: 7.0 };
        assert!(
            model
                .compute(&mut s, &mut w, &mut NeighborList::new(1.0))
                .is_err(),
            "{body}"
        );
        assert_eq!(s.forces, before);
        assert_eq!(w.w, 7.0);
    }
}

#[test]
fn snapshot_request_contains_current_scaled_positions_cell_and_z() {
    let body = r#"assert request['atomic_numbers']==[29]
 assert request['positions_angstrom']==[[0.5,1.0,1.5]]
 assert request['cell_angstrom']==[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]
 assert request['pbc']==[True]*3
 "#
    .to_owned()
        + VALID;
    let (_dir, python) = fixture(&body);
    let m = AseEmtForceModel::new(spec(&python), UnitFactors::new(2.0, 1.0, 1.0)).unwrap();
    let mut s = system();
    s.positions = vec![0.25, 0.5, 0.75];
    assert!(m.evaluate(&s).is_ok());
    s.positions[0] += 1.0;
    assert!(
        m.evaluate(&s).is_err(),
        "a cached request/result hid changed coordinates"
    );
}

#[test]
fn real_provider_matches_direct_ase_and_derivatives() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP direct ASE parity: set MODENA_ASE_PYTHON");
        return;
    };
    let id = describe_ase_emt(&python, 30000).unwrap();
    let mut config = spec(&python);
    config.model_sha256 = id.model_sha256;
    config.ase_version = id.ase_version;
    config.atomic_numbers = vec![29; 4];
    config.timeout_ms = 30000;
    let model = AseEmtForceModel::new(config, UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
    for cell in [
        [[3.6, 0.0, 0.0], [0.0, 3.6, 0.0], [0.0, 0.0, 3.6]],
        [[3.6, 0.0, 0.0], [0.3, 3.6, 0.0], [0.2, 0.1, 3.6]],
    ] {
        let mut s = System::new(SimBox::triclinic(cell.map(modena_core::Vec3)));
        s.resize(4);
        s.atomic_numbers = vec![29; 4];
        s.positions = vec![0.03, 0.0, 0.0, 0.0, 1.8, 1.8, 1.8, 0.0, 1.8, 1.8, 1.8, 0.0];
        let reference_code = r#"import json,sys
from ase import Atoms
from ase.calculators.emt import EMT
cell=json.loads(sys.argv[1])
a=Atoms(numbers=[29]*4,positions=[[0.03,0,0],[0,1.8,1.8],[1.8,0,1.8],[1.8,1.8,0]],cell=cell,pbc=True)
a.calc=EMT(asap_cutoff=False)
print(json.dumps({'energy_ev':a.get_potential_energy(),'forces_ev_per_angstrom':a.get_forces().tolist(),'stress_ev_per_angstrom3':a.get_stress(voigt=False).tolist()}))"#;
        let out = std::process::Command::new(&python)
            .arg("-c")
            .arg(reference_code)
            .arg(serde_json::to_string(&cell).unwrap())
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let direct: modena_core::ProviderEvaluation = serde_json::from_slice(&out.stdout).unwrap();
        let got = model.evaluate(&s).unwrap();
        assert!((got.energy_ev - direct.energy_ev).abs() < 1e-10);
        for (a, b) in got
            .forces_ev_per_angstrom
            .iter()
            .flatten()
            .zip(direct.forces_ev_per_angstrom.iter().flatten())
        {
            assert!((a - b).abs() < 1e-10);
        }
        let direct_w = -s.box_size.volume()
            * (0..3)
                .map(|i| direct.stress_ev_per_angstrom3[i][i])
                .sum::<f64>();
        let mut accum = s.clone();
        accum.zero_forces();
        let mut w = Virial::default();
        let energy = model
            .compute(&mut accum, &mut w, &mut NeighborList::new(1.0))
            .unwrap();
        let conversion = 1.602176634e-19 * 6.02214076e23 / 4184.0;
        assert!((energy / conversion - direct.energy_ev).abs() < 1e-10);
        assert!((w.w / conversion - direct_w).abs() < 1e-9);
        let mut plus = s.clone();
        let mut minus = s.clone();
        plus.positions[0] += 1e-5;
        minus.positions[0] -= 1e-5;
        let fd = -(model.evaluate(&plus).unwrap().energy_ev
            - model.evaluate(&minus).unwrap().energy_ev)
            / 2e-5;
        assert!((fd - got.forces_ev_per_angstrom[0][0]).abs() < 2e-6);
        let mut plus = s.clone();
        let mut minus = s.clone();
        plus.box_size.scale_uniform(1.0 + 1e-6);
        minus.box_size.scale_uniform(1.0 - 1e-6);
        for x in &mut plus.positions {
            *x *= 1.0 + 1e-6;
        }
        for x in &mut minus.positions {
            *x *= 1.0 - 1e-6;
        }
        let fd_w = -(model.evaluate(&plus).unwrap().energy_ev
            - model.evaluate(&minus).unwrap().energy_ev)
            / 2e-6;
        assert!((fd_w - direct_w).abs() < 2e-6);
        let mut translated = s.clone();
        for p in translated.positions.chunks_exact_mut(3) {
            p[0] += 0.71;
            p[1] -= 0.37;
            p[2] += 1.12;
        }
        let translated = model.evaluate(&translated).unwrap();
        assert!((translated.energy_ev - got.energy_ev).abs() < 1e-10);
        for (a, b) in translated
            .forces_ev_per_angstrom
            .iter()
            .flatten()
            .zip(got.forces_ev_per_angstrom.iter().flatten())
        {
            assert!((a - b).abs() < 1e-10);
        }
        for k in 0..3 {
            assert!(
                got.forces_ev_per_angstrom
                    .iter()
                    .map(|f| f[k])
                    .sum::<f64>()
                    .abs()
                    < 1e-10
            );
        }
    }
}

#[test]
fn malformed_json_diagnostics_are_bounded() {
    let body="print(json.dumps({'protocol':'modena-ase-emt-v1','action':'evaluate','identity':identity,'x'*100000:1}))";
    let (_dir, python) = fixture(body);
    let model = AseEmtForceModel::new(spec(&python), UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
    let error = model.evaluate(&system()).unwrap_err().to_string();
    assert!(
        error.len() <= 8192,
        "unbounded diagnostic length {}",
        error.len()
    );
}

fn runner_toml(python: &str, hash: &str, version: &str) -> String {
    format!(
        r#"units = "real"
[structure]
type = "fcc"
lattice_constant = 3.6
replicas = [1, 1, 1]
element = "Cu"
mass = 63.546
[model]
type = "ase_emt"
python = {python:?}
ase_version = {version:?}
model_sha256 = {hash:?}
atomic_numbers = [29, 29, 29, 29]
composition = "total"
[integrator]
type = "velocity_verlet"
dt = 0.0001
[run]
steps = 2
initial_temperature = 20.0
log_every = 100
"#
    )
}

#[test]
fn external_toml_runs_total_provider_and_reports_real_stress_metadata() {
    use modena_core::{io::input::InputConfig, runner::Simulation};
    let (_dir, python) = fixture(RUNNER_VALID);
    let config = InputConfig::from_str(&runner_toml(&python, &"a".repeat(64), "3.29.0")).unwrap();
    let mut sim =
        Simulation::from_config_without_outputs(&config, std::path::Path::new(".")).unwrap();
    assert_eq!(sim.system().atomic_numbers, vec![29; 4]);
    let metadata = serde_json::to_value(sim.metadata()).unwrap();
    assert_eq!(metadata["pressure_status"], "unconstrained_provider_stress");
    assert_eq!(metadata["provider"]["provider"], "ase_emt");
    let summary = sim.run().unwrap();
    assert_eq!(summary.conservation.diagnostics.unwrap().samples, 3);
    assert!(
        (summary.final_energies.potential - 2.0 * modena_core::external::EV_TO_KCAL_MOL).abs()
            < 1e-10
    );
    let expected_w = -3.6_f64.powi(3) * 0.6 * modena_core::external::EV_TO_KCAL_MOL;
    assert!((summary.final_energies.virial - expected_w).abs() < 1e-10);
    let expected_pressure = (2.0 * summary.final_energies.kinetic + expected_w)
        / (3.0 * 3.6_f64.powi(3))
        * config.unit_factors.unwrap().pressure;
    assert!((summary.final_energies.pressure - expected_pressure).abs() < 1e-9);
}

#[test]
fn external_toml_rejects_unsupported_coupling_before_launch() {
    use modena_core::io::input::InputConfig;
    let valid = runner_toml("/must/not/launch", &"a".repeat(64), "3.29.0");
    for (text, expected) in [
        (
            valid.replace("units = \"real\"", "units = \"reduced\""),
            "real",
        ),
        (
            valid.replace("composition = \"total\"", "composition = \"residual\""),
            "composition",
        ),
        (
            valid.replace("mass = 63.546", "mass = 63.546\ncharge = 1.0"),
            "charges",
        ),
        (format!("{valid}\n[constraints]\nbonds = []"), "constraints"),
        (
            format!(
                "{valid}\n[electrostatics]\nkappa = 0.3\nreal_cutoff = 3.0\nkspace_cutoff = 2.0"
            ),
            "electrostatics",
        ),
        (valid.replace("mass = 63.546\n", ""), "mass"),
        (
            valid.replace("[29, 29, 29, 29]", "[0, 29, 29, 29]"),
            "atomic_numbers",
        ),
        (format!("{valid}\n[topology]\nbonds = []"), "topology"),
        (
            format!("{valid}\n[barostat]\ntype = \"berendsen\"\ntarget_pressure = 1.0"),
            "barostat",
        ),
        (
            valid.replace(
                "composition = \"total\"",
                "composition = \"total\"\npbc = [true, false, true]",
            ),
            "pbc",
        ),
    ] {
        let error = InputConfig::from_str(&text).unwrap_err().to_string();
        assert!(error.contains(expected), "expected {expected}: {error}");
    }
}

const RUNNER_VALID: &str = "n=len(request['atomic_numbers']); print(json.dumps({'protocol':'modena-ase-emt-v1','action':'evaluate','identity':identity,'evaluation':{'energy_ev':2.0,'forces_ev_per_angstrom':[[0.0,0.0,0.0]]*n,'stress_ev_per_angstrom3':[[0.1,0.0,0.0],[0.0,0.2,0.0],[0.0,0.0,0.3]]}}))";

fn build_runner(text: &str, dir: &std::path::Path) -> modena_core::Simulation {
    let config = modena_core::InputConfig::from_str(text).unwrap();
    modena_core::Simulation::from_config_without_outputs(&config, dir).unwrap()
}

#[test]
fn provider_evaluation_error_emits_failed_at_unlogged_integration_step() {
    let (dir, python) = fixture(RUNNER_VALID);
    let mut sim = build_runner(&runner_toml(&python, &"a".repeat(64), "3.29.0"), dir.path());
    // A working real subprocess during setup becomes an actual child failure.
    std::fs::write(&python, "#!/usr/bin/env python3\nimport sys\nsys.exit(7)\n").unwrap();
    let mut events = Vec::new();
    let error = sim
        .run_streaming(None, |event| events.push(event.clone()))
        .unwrap_err();
    assert!(error.to_string().contains("7"), "{error}");
    assert!(
        matches!(events.last(), Some(modena_core::RunEvent::Failed { step: 1, error }) if error.contains("7")),
        "{events:?}"
    );
    assert!(!events
        .iter()
        .any(|e| matches!(e, modena_core::RunEvent::Finished { .. })));
}

#[test]
fn external_checkpoint_checks_atomic_numbers_runtime_identity_and_transport_independence() {
    use modena_core::{io::restart, InputConfig, Simulation};
    let (dir, python) = fixture(RUNNER_VALID);
    let text = runner_toml(&python, &"a".repeat(64), "3.29.0");
    let mut sim = build_runner(&text, dir.path());
    sim.run().unwrap();
    let path = dir.path().join("state.chk");
    sim.write_checkpoint(&path).unwrap();
    let cp = restart::load(&path).unwrap();
    assert_eq!(
        cp.continuation.as_ref().unwrap().atomic_numbers,
        vec![29; 4]
    );
    let resume = format!("{text}\nrestore_file = \"state.chk\"\n");
    let restored = build_runner(&resume, dir.path());
    assert_eq!(restored.system().positions, sim.system().positions);
    assert_eq!(restored.system().velocities, sim.system().velocities);
    assert_eq!(restored.system().forces, sim.system().forces);
    assert_eq!(restored.current_step(), 2);
    // Executable location and timeout have no physical meaning.
    let alternate = dir.path().join("alternate-python");
    std::fs::copy(&python, &alternate).unwrap();
    let transport_only = resume
        .replace(&python, alternate.to_str().unwrap())
        .replace(
            "composition = \"total\"",
            "composition = \"total\"\ntimeout_ms = 10000",
        );
    assert_eq!(
        build_runner(&transport_only, dir.path()).system().forces,
        sim.system().forces
    );
    for changed in [
        resume.replace("[29, 29, 29, 29]", "[28, 29, 29, 29]"),
        resume.replace(&"a".repeat(64), &"b".repeat(64)),
        resume.replace("3.29.0", "3.28.0"),
        resume.replace(
            "composition = \"total\"",
            "composition = \"total\"\nasap_cutoff = true",
        ),
    ] {
        let result = InputConfig::from_str(&changed)
            .and_then(|c| Simulation::from_config_without_outputs(&c, dir.path()));
        let error = result
            .err()
            .expect("changed physical identity accepted")
            .to_string();
        assert!(
            error.contains("physical identity")
                || error.contains("pin mismatch")
                || error.contains("asap_cutoff"),
            "{error}"
        );
    }
    // A valid newly discovered runtime identity (same source/version pins,
    // different NumPy) must not restore old mechanics.
    let script = std::fs::read_to_string(&python).unwrap();
    std::fs::write(&python, script.replace("2.5.4", "2.5.5")).unwrap();
    let c = InputConfig::from_str(&resume).unwrap();
    let error = Simulation::from_config_without_outputs(&c, dir.path())
        .err()
        .unwrap()
        .to_string();
    assert!(error.contains("physical identity"), "{error}");
    std::fs::write(&python, script).unwrap();
    // An external checkpoint cannot use the legacy native unspecified-Z default.
    let mut cp = cp;
    cp.continuation.as_mut().unwrap().atomic_numbers.clear();
    restart::save(&path, &cp).unwrap();
    assert!(Simulation::from_config_without_outputs(&c, dir.path()).is_err());
}

const NATIVE_TOML: &str = r#"[structure]
type = "fcc"
lattice_constant = 1.6
replicas = [1,1,1]
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
"#;

#[test]
fn older_native_v2_checkpoint_keeps_signature_and_unspecified_z_compatibility() {
    use modena_core::io::restart;
    let dir = tempfile::tempdir().unwrap();
    let mut sim = build_runner(NATIVE_TOML, dir.path());
    sim.run().unwrap();
    let path = dir.path().join("state.chk");
    sim.write_checkpoint(&path).unwrap();
    let cp = restart::load(&path).unwrap();
    let state = cp.continuation.unwrap();
    // Hand-written legacy physical-v2 identity: no atomic-number/provider keys.
    let expected = serde_json::json!({"schema":"modena-physical-v2","model":{"type":"lennard_jones","sigma":1.0,"epsilon":1.0,"cutoff":2.5,"types":[],"overrides":[],"mixing":"lorentz_berthelot"},"integrator":{"type":"velocity_verlet","dt":0.005},"thermostat":null,"barostat":null,"electrostatics":null,"topology":null,"constraints":null,"water_geometry":null,"initialization_seed":42,"initial_temperature":1.0,"masses":[1.0,1.0,1.0,1.0],"species":[0,0,0,0],"charges":[0.0,0.0,0.0,0.0],"elements":["X","X","X","X"],"output_unit_factors":null});
    assert_eq!(state.physical_signature, expected);
    let text = std::fs::read_to_string(&path).unwrap();
    let (header, payload) = text.split_once('\n').unwrap();
    let mut payload: serde_json::Value = serde_json::from_str(payload).unwrap();
    payload["continuation"]
        .as_object_mut()
        .unwrap()
        .remove("atomic_numbers");
    std::fs::write(
        &path,
        format!("{header}\n{}\n", serde_json::to_string(&payload).unwrap()),
    )
    .unwrap();
    let restored = build_runner(
        &format!("{NATIVE_TOML}\nrestore_file = \"state.chk\"\n"),
        dir.path(),
    );
    assert_eq!(restored.system().positions, sim.system().positions);
    assert_eq!(restored.system().forces, sim.system().forces);
    assert!(restored.metadata().provider.is_none());
}

#[test]
fn checkpoint_rejects_nonempty_atomic_number_buffer_with_wrong_length() {
    let dir = tempfile::tempdir().unwrap();
    let sim = build_runner(NATIVE_TOML, dir.path());
    let path = dir.path().join("state.chk");
    sim.write_checkpoint(&path).unwrap();
    let mut cp = modena_core::io::restart::load(&path).unwrap();
    cp.continuation.as_mut().unwrap().atomic_numbers = vec![29];
    assert!(modena_core::io::restart::save(&path, &cp).is_err());
}

#[test]
fn real_ase_runner_split_continuation_and_selected_z_mutation() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP ASE runner/restart/selected-Z: set MODENA_ASE_PYTHON");
        return;
    };
    let id = describe_ase_emt(&python, 30000).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let text = runner_toml(&python, &id.model_sha256, &id.ase_version);
    let mut first = build_runner(&text, dir.path());
    let a = first.run().unwrap();
    assert_eq!(a.conservation.diagnostics.unwrap().samples, 3);
    assert_eq!(a.metadata.provider.as_ref(), Some(&id));
    first.write_checkpoint(&dir.path().join("ase.chk")).unwrap();
    let mut resumed = build_runner(&format!("{text}\nrestore_file = \"ase.chk\""), dir.path());
    let b = resumed.run().unwrap();
    let mut full = build_runner(&text.replace("steps = 2", "steps = 4"), dir.path());
    let c = full.run().unwrap();
    assert_eq!(resumed.current_step(), 4);
    assert_eq!(resumed.system().positions, full.system().positions);
    assert_eq!(resumed.system().velocities, full.system().velocities);
    assert_eq!(resumed.system().forces, full.system().forces);
    assert_eq!(b.final_energies.total, c.final_energies.total);
    assert_eq!(b.final_energies.virial, c.final_energies.virial);
    // Z is explicit model input, independent of the cosmetic Cu label/species.
    let mut changed = build_runner(
        &text.replace("[29, 29, 29, 29]", "[28, 29, 29, 29]"),
        dir.path(),
    );
    let changed_summary = changed.run().unwrap();
    assert!((changed_summary.conservation.e0 - a.conservation.e0).abs() > 1e-6);
    assert_eq!(changed.system().species, vec![0; 4]);
    assert_eq!(changed.system().elements, vec!["Cu"; 4]);
    let resumed_changed = format!(
        "{}\nrestore_file = \"ase.chk\"",
        text.replace("[29, 29, 29, 29]", "[28, 29, 29, 29]")
    );
    let c = modena_core::InputConfig::from_str(&resumed_changed).unwrap();
    assert!(modena_core::Simulation::from_config_without_outputs(&c, dir.path()).is_err());
    eprintln!(
        "EXECUTED ASE {} runner, split restart and selected-Z mutation",
        id.ase_version
    );
}

#[test]
fn external_runner_checks_z_count_before_process_launch_and_keeps_ff_species_separate() {
    let invalid = runner_toml("/must/not/launch", &"a".repeat(64), "3.29.0")
        .replace("[29, 29, 29, 29]", "[29]");
    let config = modena_core::InputConfig::from_str(&invalid).unwrap();
    let error =
        modena_core::Simulation::from_config_without_outputs(&config, std::path::Path::new("."))
            .err()
            .unwrap()
            .to_string();
    assert!(error.contains("atomic_numbers"), "{error}");
    let (dir, python) = fixture(RUNNER_VALID);
    let text = runner_toml(&python, &"a".repeat(64), "3.29.0").replace(
        "type = \"fcc\"\nlattice_constant = 3.6\nreplicas = [1, 1, 1]\nelement = \"Cu\"\nmass = 63.546",
        "type = \"random_alloy\"\nlattice_constant = 3.6\nreplicas = [1,1,1]\n[[structure.species]]\nelement = \"Cu\"\nmass = 63.546\nconcentration = 0.0\n[[structure.species]]\nelement = \"Ni\"\nmass = 58.6934\nconcentration = 1.0"
    );
    let sim = build_runner(&text, dir.path());
    assert_eq!(sim.system().species, vec![1; 4]);
    assert_eq!(sim.system().elements, vec!["Ni"; 4]);
    assert_eq!(sim.system().atomic_numbers, vec![29; 4]);
}

#[test]
fn legacy_real_fixed_cell_signature_migrates_only_exact_old_pressure_factor() {
    use modena_core::{io::restart, InputConfig, Simulation};
    for integrator in [
        "type = \"velocity_verlet\"\ndt = 0.0001",
        "type = \"langevin_baoab\"\ndt = 0.0001\ngamma = 1.0\ntarget_temperature = 10.0",
        "type = \"nose_hoover_chain\"\ndt = 0.0001\ntarget_temperature = 10.0",
    ] {
        let dir = tempfile::tempdir().unwrap();
        let real = format!("units = \"real\"\n{NATIVE_TOML}")
            .replace("lattice_constant = 1.6", "lattice_constant = 3.2")
            .replace("sigma = 1.0", "sigma = 2.0")
            .replace("epsilon = 1.0", "epsilon = 3.0")
            .replace("type = \"velocity_verlet\"\ndt = 0.005", integrator);
        let mut sim = build_runner(&real, dir.path());
        sim.run().unwrap();
        let path = dir.path().join("old.chk");
        sim.write_checkpoint(&path).unwrap();
        let mut old = restart::load(&path).unwrap();
        old.continuation.as_mut().unwrap().physical_signature["output_unit_factors"][5] =
            serde_json::json!(3.0 / (2.0 * 2.0 * 2.0) * 69477.8);
        restart::save(&path, &old).unwrap();
        let resume = format!("{real}\nrestore_file = \"old.chk\"");
        let mut resumed = build_runner(&resume, dir.path());
        resumed
            .write_checkpoint(&dir.path().join("new.chk"))
            .unwrap();
        let new = restart::load(dir.path().join("new.chk")).unwrap();
        assert_eq!(
            new.continuation.as_ref().unwrap().physical_signature["output_unit_factors"][5],
            serde_json::json!(25713.158612344094)
        );
        resumed.run().unwrap();
        sim.run().unwrap();
        assert_eq!(resumed.system().positions, sim.system().positions);
        assert_eq!(resumed.system().velocities, sim.system().velocities);
        for (index, value) in [(5, 26054.17500000001), (0, 2.01), (4, 500.0)] {
            let mut bad = old.clone();
            bad.continuation.as_mut().unwrap().physical_signature["output_unit_factors"][index] =
                serde_json::json!(value);
            restart::save(&path, &bad).unwrap();
            let config = InputConfig::from_str(&resume).unwrap();
            assert!(Simulation::from_config_without_outputs(&config, dir.path()).is_err());
        }
        let mut bad = old.clone();
        bad.continuation.as_mut().unwrap().physical_signature["initialization_seed"] =
            serde_json::json!(43);
        restart::save(&path, &bad).unwrap();
        assert!(Simulation::from_config_without_outputs(
            &InputConfig::from_str(&resume).unwrap(),
            dir.path()
        )
        .is_err());
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("old.chk");
    let coupled = format!("units = \"real\"\n{NATIVE_TOML}\n[barostat]\ntype = \"monte_carlo\"\ntarget_pressure = 1.0\ntarget_temperature = 10.0").replace("type = \"velocity_verlet\"\ndt = 0.005","type = \"langevin_baoab\"\ndt = 0.0001\ngamma = 1.0\ntarget_temperature = 10.0");
    let sim = build_runner(&coupled, dir.path());
    sim.write_checkpoint(&path).unwrap();
    let mut old = restart::load(&path).unwrap();
    old.continuation.as_mut().unwrap().physical_signature["output_unit_factors"][5] =
        serde_json::json!(69477.8);
    restart::save(&path, &old).unwrap();
    let mut config = InputConfig::from_str(&coupled).unwrap();
    config.run.restore_file = Some("old.chk".into());
    assert!(Simulation::from_config_without_outputs(&config, dir.path()).is_err());
    config.run.restart_mode = modena_core::io::input::RestartMode::StateOnly;
    assert!(Simulation::from_config_without_outputs(&config, dir.path()).is_ok());
}

#[test]
fn real_pressure_targets_and_inverse_pressure_compressibility_use_standard_atmospheres() {
    use modena_core::io::input::BarostatSpec;
    let text = format!("units = \"real\"\n{NATIVE_TOML}")
        .replace("sigma = 1.0", "sigma = 2.0")
        .replace("epsilon = 1.0", "epsilon = 3.0")
        .replace(
            "type = \"velocity_verlet\"\ndt = 0.005",
            "type = \"langevin_baoab\"\ndt = 0.0001\ngamma = 1.0\ntarget_temperature = 10.0",
        );
    let mc=modena_core::InputConfig::from_str(&format!("{text}\n[barostat]\ntype = \"monte_carlo\"\ntarget_pressure = 25713.158612344094\ntarget_temperature = 10.0")).unwrap();
    let Some(BarostatSpec::MonteCarlo {
        target_pressure, ..
    }) = mc.barostat
    else {
        panic!()
    };
    assert!((target_pressure - 1.0).abs() < 1e-14);
    let ber=modena_core::InputConfig::from_str(&format!("{text}\n[barostat]\ntype = \"berendsen\"\ntarget_pressure = 25713.158612344094\ncompressibility = 0.000038890593531357555")).unwrap();
    let Some(BarostatSpec::Berendsen {
        target_pressure,
        compressibility,
        ..
    }) = ber.barostat
    else {
        panic!()
    };
    assert!((target_pressure - 1.0).abs() < 1e-14);
    assert!((compressibility - 1.0).abs() < 1e-12, "{compressibility}");
}
