#![cfg(unix)]

use modena_core::{
    AseEmtForceModel, AseEmtSpec, ForceModel, NeighborList, SimBox, System, UnitFactors, Virial,
};
use std::os::unix::fs::PermissionsExt;
use std::time::{Duration, Instant};

fn config(python: &str) -> AseEmtSpec {
    serde_json::from_value(serde_json::json!({"python":python,"ase_version":"3.29.0","model_sha256":"a".repeat(64),"atomic_numbers":[29],"composition":"total","transport":"persistent","timeout_ms":2000})).expect("persistent transport must parse")
}
fn system() -> System {
    let mut s = System::new(SimBox::cubic(4.0));
    s.resize(1);
    s.atomic_numbers = vec![29];
    s.forces = vec![2.0, 3.0, 4.0];
    s
}
fn fixture(body: &str) -> (tempfile::TempDir, String) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("python");
    let script = format!(
        r##"#!/usr/bin/env python3
import sys,json,time,os,struct
from pathlib import Path
root=Path(__file__).parent
with (root/'starts').open('a') as f:f.write(str(os.getpid())+'\n')
identity={{"protocol":"modena-ase-emt-v1","provider":"ase_emt","model":"EMT","task":"energy_forces_stress","composition":"total","units":{{"positions":"angstrom","energy":"eV","forces":"eV/angstrom","stress":"eV/angstrom^3"}},"capabilities":["energy","forces","stress","full_pbc","provider_neighbors","conservative","total_momentum"],"ase_version":"3.29.0","numpy_version":"2.5.4","model_sha256":"{}","asap_cutoff":False}}
def send(value):
 data=json.dumps(value).encode();sys.stdout.buffer.write(struct.pack('>I',len(data))+data);sys.stdout.buffer.flush()
while True:
 header=sys.stdin.buffer.read(4)
 if not header:break
 frame=json.loads(sys.stdin.buffer.read(struct.unpack('>I',header)[0]))
 request=frame['request']
 with (root/'nonces').open('a') as f:f.write(json.dumps(frame['nonce'])+'\n')
 with (root/'requests').open('a') as f:f.write(request['action']+'\n')
 while (root/('block-'+request['action'])).exists() and not (root/('release-'+request['action'])).exists():time.sleep(0.005)
 response={{'protocol':'modena-ase-emt-v1','action':request['action'],'identity':identity}}
 if request['action']=='evaluate':
  response['evaluation']={{'energy_ev':request['positions_angstrom'][0][0]+2.0,'forces_ev_per_angstrom':[[request['positions_angstrom'][0][0]+1.0,2.0,3.0]]*len(request['atomic_numbers']),'stress_ev_per_angstrom3':[[0.1,0,0],[0,0.2,0],[0,0,0.3]]}}
  {}
 send({{'session_protocol':'modena-ase-emt-session-v1','request_id':frame['request_id'],'nonce':frame['nonce'],'response':response}})
"##,
        "a".repeat(64),
        body
    );
    std::fs::write(&path, script).unwrap();
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    (dir, path.to_string_lossy().into_owned())
}
fn model(python: &str) -> AseEmtForceModel {
    AseEmtForceModel::new(config(python), UnitFactors::new(1.0, 1.0, 1.0)).unwrap()
}
fn assert_reaped(dir: &std::path::Path) {
    let pid = std::fs::read_to_string(dir.join("starts")).unwrap();
    assert_eq!(pid.lines().count(), 1, "unexpected process restart");
    let output = std::process::Command::new("kill")
        .args(["-0", pid.trim()])
        .output()
        .unwrap();
    assert!(!output.status.success(), "provider child survived cleanup");
}

fn append_after_reply(python: &str, code: &str) {
    let mut script = std::fs::read_to_string(python).unwrap();
    script.push_str(code);
    std::fs::write(python, script).unwrap();
}

#[test]
fn persistent_rejects_unsolicited_future_replies_before_their_request() {
    for after_id in [1, 2] {
        let (dir, python) = fixture("pass");
        append_after_reply(
            &python,
            &format!(
                r#" if frame['request_id']=={after_id}:
  fake={{'session_protocol':'modena-ase-emt-session-v1','request_id':frame['request_id']+1,'response':{{'protocol':'modena-ase-emt-v1','action':'evaluate','identity':identity,'evaluation':{{'energy_ev':-777.0,'forces_ev_per_angstrom':[[0.0]*3],'stress_ev_per_angstrom3':[[0.0]*3]*3}}}}}}
  if 'nonce' in frame:fake['nonce']=frame['nonce']
  send(fake)
  (root/'forged-sent').touch()
"#
            ),
        );
        let result = AseEmtForceModel::new(config(&python), UnitFactors::new(1.0, 1.0, 1.0));
        if let Ok(m) = result {
            if after_id == 2 && m.evaluate(&system()).is_err() {
                assert_reaped(dir.path());
                continue;
            }
            let start = Instant::now();
            while !dir.path().join("forged-sent").exists() {
                assert!(start.elapsed() < Duration::from_secs(5));
                std::thread::sleep(Duration::from_millis(5));
            }
            let mut s = system();
            let before = s.forces.clone();
            let mut w = Virial { w: 7.0 };
            assert!(
                m.compute(&mut s, &mut w, &mut NeighborList::new(1.0))
                    .is_err(),
                "accepted a reply fabricated before request {}",
                after_id + 1
            );
            assert_eq!(s.forces, before);
            assert_eq!(w.w, 7.0);
        }
        assert_reaped(dir.path());
    }
}

#[test]
fn persistent_reaps_idle_failures_without_another_request_or_drop() {
    for failure in [
        "sys.stderr.write('x'*70000);sys.stderr.flush()",
        "os.close(sys.stdout.fileno())",
        "sys.stdout.buffer.write(struct.pack('>I',17000000));sys.stdout.buffer.flush()",
        "send({'malformed':'envelope'})",
    ] {
        let (dir, python) = fixture("pass");
        append_after_reply(
            &python,
            &format!(" time.sleep(0.1)\n {failure}\n time.sleep(20)\n"),
        );
        let m = model(&python);
        let pid = std::fs::read_to_string(dir.path().join("starts")).unwrap();
        let start = Instant::now();
        loop {
            if !std::process::Command::new("kill")
                .args(["-0", pid.trim()])
                .output()
                .unwrap()
                .status
                .success()
            {
                break;
            }
            assert!(
                start.elapsed() < Duration::from_secs(3),
                "idle failure did not kill/reap: {failure}"
            );
            std::thread::sleep(Duration::from_millis(10));
        }
        assert!(m
            .evaluate(&system())
            .unwrap_err()
            .to_string()
            .contains("terminal"));
        assert_reaped(dir.path());
    }
}

#[test]
fn persistent_reuses_handshake_process_and_validates_before_sending() {
    let (dir, python) = fixture("pass");
    let m = model(&python);
    let mut s = system();
    let first = m.evaluate(&s).unwrap();
    s.positions[0] = 0.25;
    let changed = m.evaluate(&s).unwrap();
    assert_eq!(changed.energy_ev - first.energy_ev, 0.25);
    assert_eq!(
        changed.forces_ev_per_angstrom[0][0] - first.forces_ev_per_angstrom[0][0],
        0.25
    );
    s.charges[0] = 1.0;
    assert!(m.evaluate(&s).is_err());
    s.charges[0] = 0.0;
    assert!(
        m.evaluate(&s).is_ok(),
        "invalid unsent snapshots must not poison a session"
    );
    assert_eq!(
        std::fs::read_to_string(dir.path().join("requests")).unwrap(),
        "info\nevaluate\nevaluate\nevaluate\n"
    );
    let nonces: Vec<[u8; 32]> = std::fs::read_to_string(dir.path().join("nonces"))
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(nonces.len(), 4);
    for i in 0..nonces.len() {
        for j in 0..i {
            assert_ne!(nonces[i], nonces[j], "request nonce was reused");
        }
    }
    drop(m);
    assert_reaped(dir.path());
}

#[test]
fn persistent_protocol_and_provider_failures_are_atomic_terminal_and_reaped() {
    for body in [
        "frame['request_id']-=1", "identity['numpy_version']='changed'", "sys.exit(7)",
        "send({'session_protocol':'other','request_id':frame['request_id'],'response':response});time.sleep(20)",
        "response['protocol']='other'", "response['action']='info'", "del response['identity']",
        "response['error']='calculator failed';del response['evaluation']",
        "response['x'*100000]=1",
        "sys.stdout.buffer.write(b'\\x00\\x00');sys.stdout.buffer.flush();sys.exit(0)",
        "sys.stdout.buffer.write(struct.pack('>I',4)+b'oops');sys.stdout.buffer.flush();time.sleep(20)",
        "sys.stdout.buffer.write(struct.pack('>I',17000000));sys.stdout.buffer.flush();time.sleep(20)",
        "sys.stderr.write('x'*70000);sys.stderr.flush();time.sleep(20)",
        "response['evaluation']['forces_ev_per_angstrom']=[]", "response['evaluation']['energy_ev']=float('nan')",
        "response['evaluation']['stress_ev_per_angstrom3'][0][1]=0.2",
        "response['evaluation']['forces_ev_per_angstrom'][0][0]=1e308",
        "time.sleep(20)",
    ] {
        let (dir, python)=fixture(body);
        let m=model(&python);
        let mut s=system();
        let before=s.forces.clone();
        let mut w=Virial { w:7.0 };
        let start=Instant::now();
        let error=m.compute(&mut s,&mut w,&mut NeighborList::new(1.0)).unwrap_err();
        assert!(error.to_string().len()<8192,"unbounded diagnostic");
        assert!(start.elapsed()<Duration::from_secs(8),"hung session: {body}");
        assert_eq!(s.forces,before,"{body}");
        assert_eq!(w.w,7.0,"{body}");
        let calls=std::fs::read(dir.path().join("requests")).unwrap();
        let error=m.evaluate(&s).unwrap_err().to_string();
        assert!(error.contains("terminal"),"{body}: {error}");
        assert_eq!(std::fs::read(dir.path().join("requests")).unwrap(),calls);
        assert_reaped(dir.path());
    }
}

#[test]
fn persistent_deadline_includes_writes_to_a_nonreading_worker() {
    // Worker completes handshake, then never drains stdin. A large snapshot
    // exceeds pipe capacity so a blocking write would defeat the deadline.
    let (dir, python) = fixture("pass");
    let script = std::fs::read_to_string(&python).unwrap();
    // Insert the stall after the first complete reply.
    let script = script.replace(
        "'response':response})",
        "'response':response})\n time.sleep(20)",
    );
    std::fs::write(&python, script).unwrap();
    let mut cfg = config(&python);
    cfg.atomic_numbers = vec![29; 20000];
    let m = AseEmtForceModel::new(cfg, UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
    let mut s = System::new(SimBox::cubic(4.0));
    s.resize(20000);
    s.atomic_numbers = vec![29; 20000];
    let start = Instant::now();
    assert!(m.evaluate(&s).unwrap_err().to_string().contains("timeout"));
    assert!(start.elapsed() < Duration::from_secs(8));
    assert_reaped(dir.path());
}

#[test]
fn persistent_bad_handshake_and_pin_mismatch_reap_the_child() {
    for mutation in [
        "identity['units']['energy']='wrong'",
        "identity['ase_version']='wrong'",
        "identity['model_sha256']='b'*64",
        "response['action']='evaluate'",
        "frame['request_id']=0",
    ] {
        let (dir, python) = fixture("pass");
        let script = std::fs::read_to_string(&python).unwrap().replace(
            " if request['action']=='evaluate':",
            &format!(" {mutation}\n if request['action']=='evaluate':"),
        );
        std::fs::write(&python, script).unwrap();
        assert!(
            AseEmtForceModel::new(config(&python), UnitFactors::new(1.0, 1.0, 1.0)).is_err(),
            "{mutation}"
        );
        assert_reaped(dir.path());
    }
}

#[test]
fn persistent_accumulation_failure_is_atomic_and_terminal() {
    let (dir, python) = fixture("pass");
    let m = model(&python);
    let mut s = system();
    s.forces[1] = f64::NAN;
    let mut w = Virial { w: 7.0 };
    assert!(m
        .compute(&mut s, &mut w, &mut NeighborList::new(1.0))
        .is_err());
    assert!(s.forces[1].is_nan());
    assert_eq!(s.forces[0], 2.0);
    assert_eq!(s.forces[2], 4.0);
    assert_eq!(w.w, 7.0);
    assert!(m
        .evaluate(&system())
        .unwrap_err()
        .to_string()
        .contains("terminal"));
    assert_reaped(dir.path());
}

#[test]
fn persistent_post_evaluation_conversion_failure_is_terminal() {
    let (dir, python) = fixture("pass");
    let factors = UnitFactors {
        sigma: 1e-110,
        ..UnitFactors::new(1.0, 1.0, 1.0)
    };
    let m = AseEmtForceModel::new(config(&python), factors).unwrap();
    let mut s = system();
    s.box_size = SimBox::cubic(4e110);
    // The scaled cell is finite, whereas native-volume conversion overflows.
    assert!(m.evaluate(&s).is_ok());
    let mut w = Virial { w: 7.0 };
    let before = s.forces.clone();
    assert!(m
        .compute(&mut s, &mut w, &mut NeighborList::new(1.0))
        .is_err());
    assert_eq!(s.forces, before);
    assert_eq!(w.w, 7.0);
    assert!(m.evaluate(&s).unwrap_err().to_string().contains("terminal"));
    assert_reaped(dir.path());
}

#[test]
fn persistent_real_changed_snapshots_match_direct_ase() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP: set MODENA_ASE_PYTHON for direct ASE session reference");
        return;
    };
    let id = modena_core::describe_ase_emt(&python, 30000).unwrap();
    let mut cfg = config(&python);
    cfg.ase_version = id.ase_version;
    cfg.model_sha256 = id.model_sha256;
    cfg.atomic_numbers = vec![29; 4];
    cfg.timeout_ms = 30000;
    let m = AseEmtForceModel::new(cfg, UnitFactors::new(1.0, 1.0, 1.0)).unwrap();
    let mut last: Option<f64> = None;
    for (x, length, skew) in [(0.03, 3.6, 0.0), (0.17, 3.6, 0.0), (0.17, 3.7, 0.2)] {
        let cell = [[length, 0.0, 0.0], [skew, length, 0.0], [0.1, 0.3, length]];
        let positions = [
            [x, 0.0, 0.0],
            [0.0, 1.8, 1.8],
            [1.8, 0.0, 1.8],
            [1.8, 1.8, 0.0],
        ];
        let mut s = System::new(SimBox::triclinic(cell.map(modena_core::Vec3)));
        s.resize(4);
        s.atomic_numbers = vec![29; 4];
        s.positions = positions.into_iter().flatten().collect();
        let code = r#"import sys,json
from ase import Atoms
from ase.calculators.emt import EMT
a=Atoms(numbers=[29]*4,positions=json.loads(sys.argv[1]),cell=json.loads(sys.argv[2]),pbc=True)
a.calc=EMT(asap_cutoff=False)
print(json.dumps({'energy_ev':a.get_potential_energy(),'forces_ev_per_angstrom':a.get_forces().tolist(),'stress_ev_per_angstrom3':a.get_stress(voigt=False).tolist()}))"#;
        let out = std::process::Command::new(&python)
            .args([
                "-c",
                code,
                &serde_json::to_string(&positions).unwrap(),
                &serde_json::to_string(&cell).unwrap(),
            ])
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let expected: modena_core::ProviderEvaluation =
            serde_json::from_slice(&out.stdout).unwrap();
        let got = m.evaluate(&s).unwrap();
        assert!((got.energy_ev - expected.energy_ev).abs() < 1e-10);
        for (a, b) in got
            .forces_ev_per_angstrom
            .iter()
            .flatten()
            .chain(got.stress_ev_per_angstrom3.iter().flatten())
            .zip(
                expected
                    .forces_ev_per_angstrom
                    .iter()
                    .flatten()
                    .chain(expected.stress_ev_per_angstrom3.iter().flatten()),
            )
        {
            assert!((a - b).abs() < 1e-10);
        }
        if let Some(e) = last {
            assert!((got.energy_ev - e).abs() > 1e-8);
        }
        last = Some(got.energy_ev);
    }
}

fn runner_toml(python: &str) -> String {
    format!(
        r#"units="real"
[structure]
type="fcc"
lattice_constant=3.6
replicas=[1,1,1]
mass=63.546
[model]
type="ase_emt"
python={python:?}
ase_version="3.29.0"
model_sha256="{}"
atomic_numbers=[29,29,29,29]
composition="total"
transport="persistent"
timeout_ms=5000
[integrator]
type="velocity_verlet"
dt=0.0001
[run]
steps=2
initial_temperature=10.0
"#,
        "a".repeat(64)
    )
}

#[test]
fn persistent_checkpoint_transport_independence_and_safe_step_cancellation() {
    use std::sync::atomic::{AtomicBool, Ordering};
    let (dir, python) = fixture("pass");
    let text = runner_toml(&python);
    let build = |text: &str| {
        modena_core::Simulation::from_config_without_outputs(
            &modena_core::InputConfig::from_str(text).unwrap(),
            dir.path(),
        )
        .unwrap()
    };
    let mut first = build(&text);
    first.run().unwrap();
    first
        .write_checkpoint(&dir.path().join("state.chk"))
        .unwrap();
    let mut resumed = build(&format!("{text}\nrestore_file=\"state.chk\"\n"));
    resumed.run().unwrap();
    let mut full = build(&text.replace("steps=2", "steps=4"));
    full.run().unwrap();
    assert_eq!(resumed.system().positions, full.system().positions);
    assert_eq!(resumed.system().velocities, full.system().velocities);
    assert_eq!(resumed.system().forces, full.system().forces);
    // Change transport in the checkpoint signature's input configuration: the
    // native physical signature intentionally contains neither transport nor path.
    let cp = modena_core::io::restart::load(dir.path().join("state.chk")).unwrap();
    let signature = &cp.continuation.unwrap().physical_signature;
    assert!(signature["model"].get("transport").is_none());
    assert!(signature["model"].get("python").is_none());
    let stop = AtomicBool::new(false);
    let before = full.system().positions.clone();
    let summary = full
        .run_streaming(Some(&stop), |e| {
            if matches!(e, modena_core::RunEvent::Started { .. }) {
                stop.store(true, Ordering::Release);
            }
        })
        .unwrap();
    assert_eq!(summary.steps, 0);
    assert_eq!(full.system().positions, before);
    assert_eq!(
        std::fs::read_to_string(dir.path().join("starts"))
            .unwrap()
            .lines()
            .count(),
        3
    );
}

#[test]
fn persistent_setup_cancellation_stops_at_existing_safe_boundaries_and_reaps() {
    use std::sync::atomic::{AtomicBool, Ordering};
    for action in ["pre", "info", "evaluate"] {
        let (dir, python) = fixture("pass");
        let cfg = modena_core::InputConfig::from_str(&runner_toml(&python)).unwrap();
        let stop = AtomicBool::new(action == "pre");
        if action != "pre" {
            std::fs::write(dir.path().join(format!("block-{action}")), "").unwrap();
        }
        std::thread::scope(|scope| {
            let handle = scope
                .spawn(|| modena_core::Simulation::from_config_with_stop(&cfg, dir.path(), &stop));
            if action != "pre" {
                let count = if action == "info" { 1 } else { 2 };
                let start = Instant::now();
                while std::fs::read_to_string(dir.path().join("requests"))
                    .unwrap_or_default()
                    .lines()
                    .count()
                    < count
                {
                    assert!(start.elapsed() < Duration::from_secs(10));
                    std::thread::sleep(Duration::from_millis(5));
                }
                stop.store(true, Ordering::Release);
                std::fs::write(dir.path().join(format!("release-{action}")), "").unwrap();
            }
            assert!(matches!(
                handle.join().unwrap(),
                Err(modena_core::ModenaError::Cancelled)
            ));
        });
        let actions = std::fs::read_to_string(dir.path().join("requests")).unwrap_or_default();
        assert_eq!(
            actions,
            match action {
                "pre" => "",
                "info" => "info\n",
                _ => "info\nevaluate\n",
            }
        );
        if action == "pre" {
            assert!(!dir.path().join("starts").exists());
        } else {
            assert_reaped(dir.path());
        }
    }
}
