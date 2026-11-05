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
