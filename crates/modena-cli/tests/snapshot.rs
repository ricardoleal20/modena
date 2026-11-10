use serde_json::{json, Value};
use std::io::Write;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};

fn request() -> Value {
    json!({"protocol":"modena-native-lj-v1", "units":{"length":"angstrom","energy":"eV","force":"eV/angstrom"},
        "atomic_numbers":[18,18], "positions_angstrom":[[0.,0.,0.],[4.,0.,0.]],
        "cell_angstrom":[[20.,0.,0.],[0.,20.,0.],[0.,0.,20.]], "pbc":[true,true,true],
        "sigma_angstrom":3.4,"epsilon_ev":0.0103,"cutoff_angstrom":8.5,"skin_angstrom":0.3})
}

fn invoke(input: &[u8]) -> Output {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!(
        "modena_snapshot_{}_{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_modena"))
        .current_dir(&dir)
        .arg("snapshot")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    // A bounded reader can close its pipe before a deliberately oversized input completes.
    let _ = child.stdin.take().unwrap().write_all(input);
    let output = child.wait_with_output().unwrap();
    assert_eq!(
        std::fs::read_dir(&dir).unwrap().count(),
        0,
        "snapshot wrote files"
    );
    std::fs::remove_dir(dir).unwrap();
    output
}

fn evaluate(input: &Value) -> Value {
    let out = invoke(&serde_json::to_vec(input).unwrap());
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).unwrap()
}

#[test]
fn snapshot_matches_independent_shifted_force_pair_and_finite_difference() {
    let mut req = request();
    let out = evaluate(&req);
    // Independent direct derivative of 4 epsilon ((sigma/r)^12-(sigma/r)^6).
    let u = |r: f64| 4.0 * 0.0103 * ((3.4 / r).powi(12) - (3.4 / r).powi(6));
    let f = |r: f64| 24.0 * 0.0103 / r * (2.0 * (3.4 / r).powi(12) - (3.4 / r).powi(6));
    let expected = u(4.0) - u(8.5) + (4.0 - 8.5) * f(8.5);
    assert!((out["energy_ev"].as_f64().unwrap() - expected).abs() < 1e-14);
    assert_eq!(out["protocol"], req["protocol"]);
    assert_eq!(out["units"], req["units"]);
    let force = out["forces_ev_per_angstrom"][0][0].as_f64().unwrap();
    assert!((force + f(4.0) - f(8.5)).abs() < 1e-14);
    assert_eq!(
        out["forces_ev_per_angstrom"][1][0].as_f64().unwrap(),
        -force
    );
    assert!(out.get("stress").is_none());
    let h = 1e-5;
    req["positions_angstrom"][0][0] = json!(h);
    let plus = evaluate(&req)["energy_ev"].as_f64().unwrap();
    req["positions_angstrom"][0][0] = json!(-h);
    let minus = evaluate(&req)["energy_ev"].as_f64().unwrap();
    assert!((force + (plus - minus) / (2. * h)).abs() < 1e-10);
    req["positions_angstrom"][0][0] = json!(20.0);
    assert_eq!(evaluate(&req), out);
    req["positions_angstrom"][1][0] = json!(8.5);
    assert_eq!(evaluate(&req)["energy_ev"], 0.0);
}

#[test]
fn snapshot_rejects_unsupported_and_malformed_requests_without_stdout() {
    let invalid = [
        ("protocol", json!("wrong")),
        ("units", json!({})),
        ("pbc", json!([true, false, true])),
        ("atomic_numbers", json!([18, 29])),
        ("atomic_numbers", json!([0, 0])),
        ("atomic_numbers", json!([119, 119])),
        ("atomic_numbers", json!([])),
        ("positions_angstrom", json!([[0., 0., 0.]])),
        ("positions_angstrom", json!([[0., 0.], [4., 0., 0.]])),
        ("positions_angstrom", json!([[0., 0., 0.], [20., 0., 0.]])),
        ("positions_angstrom", json!([[0., 0., 0.], [1e-8, 0., 0.]])),
        (
            "cell_angstrom",
            json!([[20., 1., 0.], [0., 20., 0.], [0., 0., 20.]]),
        ),
        (
            "cell_angstrom",
            json!([[0., 0., 0.], [0., 20., 0.], [0., 0., 20.]]),
        ),
        ("sigma_angstrom", json!(-1.)),
        ("epsilon_ev", json!(-1.)),
        ("cutoff_angstrom", json!(10.)),
        ("skin_angstrom", json!(-0.1)),
        ("extra", json!(true)),
        ("velocities", json!([])),
        ("sigma_angstrom", json!(1e200)),
    ];
    for (key, value) in invalid {
        let mut req = request();
        req[key] = value;
        let out = invoke(&serde_json::to_vec(&req).unwrap());
        assert!(!out.status.success(), "accepted {key}: {req}");
        assert!(out.stdout.is_empty());
        assert!(!String::from_utf8_lossy(&out.stderr).contains("unrecognized subcommand"));
    }
    for input in [
        b"{} garbage".as_slice(),
        b"NaN",
        b"{}",
        &vec![b' '; 8 * 1024 * 1024 + 1],
    ] {
        let out = invoke(input);
        assert!(!out.status.success());
        assert!(out.stdout.is_empty());
    }
}
