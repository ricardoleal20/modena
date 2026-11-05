use std::process::Command;

#[test]
fn provider_info_reports_transport_failure_without_non_json_stdout() {
    let output = Command::new(env!("CARGO_BIN_EXE_modena"))
        .args([
            "provider-info",
            "--python",
            "/nonexistent/modena-python",
            "--timeout-ms",
            "100",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    let error = String::from_utf8_lossy(&output.stderr);
    assert!(!error.contains("unrecognized subcommand"), "{error}");
    assert!(error.contains("provider"), "{error}");
}

#[test]
fn installed_provider_cli_discovery_and_run_metadata() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP: set MODENA_ASE_PYTHON for actual ASE CLI reference");
        return;
    };
    let output = Command::new(env!("CARGO_BIN_EXE_modena"))
        .args([
            "provider-info",
            "--python",
            &python,
            "--timeout-ms",
            "120000",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let identity: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(identity["protocol"], "modena-ase-emt-v1");
    assert_eq!(identity["provider"], "ase_emt");
    assert_eq!(identity["composition"], "total");
    assert_eq!(identity["asap_cutoff"], false);
    assert_eq!(identity["model_sha256"].as_str().unwrap().len(), 64);
    let d = std::env::temp_dir().join(format!("modena_cli_ase_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    let config = format!(
        r#"units="real"
[structure]
type="fcc"
lattice_constant=3.61
replicas=[1,1,1]
mass=63.546
[model]
type="ase_emt"
python={python:?}
ase_version={version}
model_sha256={hash}
atomic_numbers=[29,29,29,29]
composition="total"
timeout_ms=120000
transport="persistent"
[integrator]
type="velocity_verlet"
dt=0.0001
[run]
steps=2
initial_temperature=10.0
log_every=100
"#,
        version = identity["ase_version"],
        hash = identity["model_sha256"]
    );
    std::fs::write(d.join("in.toml"), config).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_modena"))
        .args(["run", "--diagnostics-json", "-i"])
        .arg(d.join("in.toml"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.contains("ase_emt"));
    assert!(text.contains(identity["ase_version"].as_str().unwrap()));
    let metadata = text
        .lines()
        .find_map(|line| line.trim().strip_prefix("run metadata = "))
        .unwrap();
    let metadata: serde_json::Value = serde_json::from_str(metadata).unwrap();
    assert_eq!(metadata["provider"], identity);
    assert_eq!(metadata["pressure_status"], "unconstrained_provider_stress");
    assert!(text
        .lines()
        .any(|line| line.trim().starts_with("energy diagnostics = ")));
    std::fs::remove_dir_all(d).unwrap();
}

#[cfg(unix)]
#[test]
fn cold_copper_cli_pressure_matches_direct_si_ase_from_other_cwd_with_venv_symlink() {
    let Ok(python) = std::env::var("MODENA_ASE_PYTHON") else {
        eprintln!("SKIP actual cold Cu SI pressure/relative CLI: set MODENA_ASE_PYTHON");
        return;
    };
    let info = Command::new(env!("CARGO_BIN_EXE_modena"))
        .args(["provider-info", "--python", &python])
        .output()
        .unwrap();
    assert!(info.status.success());
    let identity: serde_json::Value = serde_json::from_slice(&info.stdout).unwrap();
    let root =
        std::env::temp_dir().join(format!("modena_pressure_relative_{}", std::process::id()));
    let config_dir = root.join("config");
    let other = root.join("other");
    std::fs::create_dir_all(&config_dir).unwrap();
    std::fs::create_dir_all(&other).unwrap();
    let link = config_dir.join("venv");
    let _ = std::fs::remove_file(&link);
    std::os::unix::fs::symlink(
        std::path::Path::new(&python)
            .parent()
            .unwrap()
            .parent()
            .unwrap(),
        &link,
    )
    .unwrap();
    let config = format!(
        r#"units="real"
[structure]
type="fcc"
lattice_constant=3.6
replicas=[1,1,1]
mass=63.546
[model]
type="ase_emt"
python="./venv/bin/python"
ase_version={version}
model_sha256={hash}
atomic_numbers=[29,29,29,29]
composition="total"
[integrator]
type="velocity_verlet"
dt=0.0001
[run]
steps=0
initial_temperature=0.0
log_file="pressure.csv"
"#,
        version = identity["ase_version"],
        hash = identity["model_sha256"]
    );
    std::fs::write(config_dir.join("in.toml"), config).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_modena"))
        .current_dir(&other)
        .args(["run", "--diagnostics-json", "-i"])
        .arg(config_dir.join("in.toml"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let csv = std::fs::read_to_string(config_dir.join("pressure.csv")).unwrap();
    let mut rows = csv.lines();
    let header = rows.next().unwrap();
    let index = header
        .split(',')
        .position(|name| name == "pressure")
        .unwrap();
    let pressure: f64 = rows
        .next()
        .unwrap()
        .split(',')
        .nth(index)
        .unwrap()
        .parse()
        .unwrap();
    let direct=Command::new(&python).args(["-c",r#"from ase import Atoms
from ase.calculators.emt import EMT
a=Atoms(numbers=[29]*4,positions=[[0,0,0],[0,1.8,1.8],[1.8,0,1.8],[1.8,1.8,0]],cell=[3.6]*3,pbc=True)
a.calc=EMT(asap_cutoff=False)
print(-sum(a.get_stress()[:3])/3 * 1.602176634e-19 / 1e-30 / 101325)
"#]).output().unwrap();
    assert!(
        direct.status.success(),
        "{}",
        String::from_utf8_lossy(&direct.stderr)
    );
    let expected: f64 = String::from_utf8(direct.stdout)
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert!(
        (pressure - expected).abs() < 1e-7,
        "CLI {pressure} atm vs independent SI/ASE {expected} atm"
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let metadata = stdout
        .lines()
        .find_map(|line| line.trim().strip_prefix("run metadata = "))
        .unwrap();
    let metadata: serde_json::Value = serde_json::from_str(metadata).unwrap();
    assert_eq!(metadata["provider"], identity);
    eprintln!("EXECUTED cold Cu CLI/SI pressure: actual={pressure:.12} expected={expected:.12} atm; relative venv symlink preserved");
    std::fs::remove_dir_all(root).unwrap();
}

#[cfg(unix)]
#[test]
fn dotted_bare_interpreter_uses_child_path_and_never_input_directory_decoy() {
    use std::os::unix::fs::PermissionsExt;
    let root = std::env::temp_dir().join(format!("modena_dotted_path_{}", std::process::id()));
    let bin = root.join("bin");
    let input = root.join("input");
    let cwd = root.join("other-cwd");
    for dir in [&bin, &input, &cwd] {
        std::fs::create_dir_all(dir).unwrap();
    }
    let marker = root.join("path-provider.calls");
    let decoy_marker = root.join("decoy.calls");
    let provider = bin.join(".emt-python");
    let script = format!(
        r#"#!/usr/bin/env python3
import sys,json
with open({marker:?},'a') as f:f.write('path-provider\n')
r=json.load(sys.stdin)
id={{"protocol":"modena-ase-emt-v1","provider":"ase_emt","model":"EMT","task":"energy_forces_stress","composition":"total","units":{{"positions":"angstrom","energy":"eV","forces":"eV/angstrom","stress":"eV/angstrom^3"}},"capabilities":["energy","forces","stress","full_pbc","provider_neighbors","conservative","total_momentum"],"ase_version":"3.29.0","numpy_version":"2.5.4","model_sha256":"{hash}","asap_cutoff":False}}
out={{'protocol':'modena-ase-emt-v1','action':r['action'],'identity':id}}
if r['action']=='evaluate':out['evaluation']={{'energy_ev':0.0,'forces_ev_per_angstrom':[[0.0]*3]*4,'stress_ev_per_angstrom3':[[0.0]*3]*3}}
print(json.dumps(out))
"#,
        marker = marker.to_str().unwrap(),
        hash = "a".repeat(64)
    );
    std::fs::write(&provider, script).unwrap();
    std::fs::set_permissions(&provider, std::fs::Permissions::from_mode(0o755)).unwrap();
    let decoy = input.join(".emt-python");
    std::fs::write(&decoy,format!("#!/usr/bin/env python3\nwith open({:?},'a') as f:f.write('wrong executable')\nraise SystemExit(7)\n",decoy_marker.to_str().unwrap())).unwrap();
    std::fs::set_permissions(decoy, std::fs::Permissions::from_mode(0o755)).unwrap();
    let mut paths = vec![bin];
    paths.extend(std::env::split_paths(
        &std::env::var_os("PATH").unwrap_or_default(),
    ));
    let path = std::env::join_paths(paths).unwrap();
    let info = Command::new(env!("CARGO_BIN_EXE_modena"))
        .current_dir(&cwd)
        .env("PATH", &path)
        .args(["provider-info", "--python", ".emt-python"])
        .output()
        .unwrap();
    assert!(
        info.status.success(),
        "{}",
        String::from_utf8_lossy(&info.stderr)
    );
    let identity: serde_json::Value = serde_json::from_slice(&info.stdout).unwrap();
    let config = format!(
        r#"units="real"
[structure]
type="fcc"
lattice_constant=3.6
replicas=[1,1,1]
mass=63.546
[model]
type="ase_emt"
python=".emt-python"
ase_version="3.29.0"
model_sha256="{}"
atomic_numbers=[29,29,29,29]
composition="total"
[integrator]
type="velocity_verlet"
dt=0.0001
[run]
steps=0
initial_temperature=0.0
"#,
        "a".repeat(64)
    );
    std::fs::write(input.join("in.toml"), config).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_modena"))
        .current_dir(&cwd)
        .env("PATH", &path)
        .args(["run", "--diagnostics-json", "-i"])
        .arg(input.join("in.toml"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).unwrap();
    let metadata = stdout
        .lines()
        .find_map(|line| line.trim().strip_prefix("run metadata = "))
        .unwrap();
    let metadata: serde_json::Value = serde_json::from_str(metadata).unwrap();
    assert_eq!(metadata["provider"], identity);
    assert_eq!(std::fs::read_to_string(marker).unwrap().lines().count(), 3);
    assert!(
        !decoy_marker.exists(),
        "selected the same-named config-directory decoy"
    );
    std::fs::remove_dir_all(root).unwrap();
}
