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
        timeout_ms: 2000,
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
    let id = describe_ase_emt(&p, 2000).unwrap();
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
    for body in ["sys.exit(7)","time.sleep(3)","print('not json')",
        "print('x'*65537, file=sys.stderr)", "print('x'*17000000)",
      "print(json.dumps({'protocol':'modena-ase-emt-v1','action':'evaluate','identity':identity,'evaluation':{'energy_ev':1.0}}))",
      "identity changed"] {
        let body=if body.starts_with("identity") {format!("identity['numpy_version']='changed'; {VALID}")} else {body.into()};
        let (_d,p)=fixture(&body);let mut spec=spec(&p); spec.timeout_ms=1000;let model=AseEmtForceModel::new(spec,UnitFactors::new(1.0,1.0,1.0)).unwrap();
        let mut s=system();let before=s.forces.clone();let mut w=Virial{w:7.0};assert!(model.compute(&mut s,&mut w,&mut NeighborList::new(1.0)).is_err(),"{body}");assert_eq!(s.forces,before);assert_eq!(w.w,7.0);
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
