#[test]
fn check_preserves_configured_output_sentinels() {
    let d = std::env::temp_dir().join(format!(
        "modena_cli_check_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&d).unwrap();
    for file in ["traj.xyz", "log.csv", "restart.chk"] {
        std::fs::write(d.join(file), b"sentinel\n").unwrap();
    }
    let text = r#"
[structure]
type="fcc"
lattice_constant=2.4
replicas=[1,1,1]
[model]
type="lennard_jones"
sigma=1.0
epsilon=1.0
cutoff=2.5
[integrator]
type="velocity_verlet"
dt=0.003
[run]
steps=1
trajectory_file="traj.xyz"
log_file="log.csv"
restart_file="restart.chk"
"#;
    std::fs::write(d.join("in.toml"), text).unwrap();
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_modena"))
        .args(["check", "--input"])
        .arg(d.join("in.toml"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for file in ["traj.xyz", "log.csv", "restart.chk"] {
        assert_eq!(
            std::fs::read(d.join(file)).unwrap(),
            b"sentinel\n",
            "check modified {file}"
        );
    }
    let missing = text
        .replace("traj.xyz", "new.xyz")
        .replace("log.csv", "new.csv")
        .replace("restart.chk", "new.chk");
    std::fs::write(d.join("in.toml"), missing).unwrap();
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_modena"))
        .args(["check", "--input"])
        .arg(d.join("in.toml"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for file in ["new.xyz", "new.csv", "new.chk"] {
        assert!(!d.join(file).exists(), "check created {file}");
    }
    std::fs::remove_dir_all(d).unwrap();
}

#[test]
fn constrained_run_labels_pressure_as_diagnostic_and_nvt_energy_as_exchange() {
    let d = std::env::temp_dir().join(format!("modena_cli_diagnostic_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    std::fs::write(
        d.join("atoms.xyz"),
        "3\nfixture\nX 1 1 1\nX 2.2 1 1\nX 5 5 5\n",
    )
    .unwrap();
    let text = r#"
[structure]
type="xyz"
file="atoms.xyz"
box_size=[8.0,8.0,8.0]
[model]
type="lennard_jones"
sigma=1.0
epsilon=1.0
cutoff=2.5
[integrator]
type="langevin_baoab"
dt=0.001
gamma=1.0
target_temperature=1.0
[[constraints.bonds]]
atoms=[0,1]
length=1.2
[run]
steps=3
"#;
    std::fs::write(d.join("in.toml"), text).unwrap();
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_modena"))
        .args(["run", "--input"])
        .arg(d.join("in.toml"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(
        text.contains("diagnostic") && text.contains("constraint impulse stress"),
        "{text}"
    );
    assert!(text.contains("energy exchanged"), "{text}");
    assert!(!text.contains("conservation  drift"), "{text}");
    std::fs::remove_dir_all(d).unwrap();
}
