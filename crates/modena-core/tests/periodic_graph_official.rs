//! Explicit external gates: ignored compilation alone is never G1/G2 acceptance.
#[path = "support/periodic_oracle.rs"]
#[allow(dead_code)]
mod periodic_oracle;

use periodic_oracle::*;
use serde::Deserialize;
use serde_json::json;
use std::io::Write;
use std::process::{Command, Stdio};

#[derive(Debug, Deserialize)]
struct Response {
    fixture_id: String,
    edge_count: usize,
    versions: serde_json::Value,
    source_sha256: String,
    max_vector_error: f64,
    source_path: String,
    interpreter: String,
    platform: String,
    machine: String,
    request_sha256: String,
    fixture_sha256: String,
    cases_sha256: String,
    helper_sha256: String,
}
fn sha256(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|b| b.is_ascii_hexdigit())
}
#[test]
#[ignore = "requires explicitly selected pinned official ASE/NumPy interpreter"]
fn official_ase_images() {
    let python = std::env::var_os("MODENA_ASE_PYTHON")
        .expect("set MODENA_ASE_PYTHON to the pinned ASE=3.29.0 NumPy=2.5.4 interpreter");
    let all = fixtures();
    // Preflight the whole suite before invoking any official helper.
    for f in all.iter().filter(|f| f.official) {
        assert!(
            cutoff_margin(f) >= 1e-8,
            "{}: unresolved cutoff margin",
            f.id
        );
        assert_eq!(
            actual(f),
            fixed_cube(f),
            "{}: independent tuple preflight",
            f.id
        );
    }
    let helper = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/support/ase_graph_reference.py"
    );
    for f in all.iter().filter(|f| f.official) {
        let edges = actual(f);
        let request =
            serde_json::to_vec(&json!({"fixture": f, "edges": wire_edges(&edges)})).unwrap();
        let mut child = Command::new(&python)
            .arg(helper)
            .env("PYTHONDONTWRITEBYTECODE", "1")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap_or_else(|e| panic!("{}: ASE interpreter: {e}", f.id));
        child.stdin.take().unwrap().write_all(&request).unwrap();
        let output = child
            .wait_with_output()
            .expect("wait for official ASE helper");
        let stderr = String::from_utf8_lossy(&output.stderr[..output.stderr.len().min(8192)]);
        assert!(
            output.status.success(),
            "{}: ASE helper {}\nstderr (first 8192 bytes): {}",
            f.id,
            output.status,
            stderr
        );
        let response: Response = serde_json::from_slice(&output.stdout)
            .unwrap_or_else(|e| panic!("{}: malformed/missing ASE response: {e}\n{stderr}", f.id));
        assert_eq!(response.fixture_id, f.id);
        assert_eq!(response.edge_count, edges.len());
        assert_eq!(response.versions["ase"], "3.29.0");
        assert_eq!(response.versions["numpy"], "2.5.4");
        assert!(response.versions["python"]
            .as_str()
            .is_some_and(|s| !s.is_empty()));
        assert!(response.max_vector_error.is_finite() && response.max_vector_error >= 0.0);
        for digest in [
            &response.source_sha256,
            &response.request_sha256,
            &response.fixture_sha256,
            &response.cases_sha256,
            &response.helper_sha256,
        ] {
            assert!(sha256(digest));
        }
        for identity in [
            &response.source_path,
            &response.interpreter,
            &response.platform,
            &response.machine,
        ] {
            assert!(!identity.is_empty());
        }
        let mut manifest: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        manifest["fixed_cube_cutoff_margin"] = json!(cutoff_margin(f));
        println!("{}", serde_json::to_string(&manifest).unwrap());
    }
}

#[test]
#[ignore = "requires explicitly selected pinned official ASE/NumPy interpreter"]
fn official_ase_rejects_wrong_tuple_with_python_optimization() {
    let python = std::env::var_os("MODENA_ASE_PYTHON")
        .expect("set MODENA_ASE_PYTHON to the pinned ASE=3.29.0 NumPy=2.5.4 interpreter");
    let mut f = fixtures().into_iter().find(|f| f.id == "multi4").unwrap();
    f.id = "optimized_tuple_identity_control".into();
    f.positions[1] = f.positions[0];
    assert!(
        cutoff_margin(&f) >= 1e-8,
        "optimized tuple control cutoff margin"
    );
    let edges = actual(&f);
    assert_eq!(edges, [(0, 1, [0; 3]), (1, 0, [0; 3])]);
    let helper = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/support/ase_graph_reference.py"
    );
    let run = |edges: &[Tuple]| {
        let request =
            serde_json::to_vec(&json!({"fixture": f, "edges": wire_edges(edges)})).unwrap();
        let mut child = Command::new(&python)
            .arg(helper)
            .env("PYTHONDONTWRITEBYTECODE", "1")
            .env("PYTHONOPTIMIZE", "1")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("spawn optimized official ASE control");
        child.stdin.take().unwrap().write_all(&request).unwrap();
        child
            .wait_with_output()
            .expect("wait for optimized official ASE control")
    };
    let positive = run(&edges);
    assert!(
        positive.status.success(),
        "optimized positive control failed: {}",
        String::from_utf8_lossy(&positive.stderr[..positive.stderr.len().min(8192)])
    );
    let response: Response = serde_json::from_slice(&positive.stdout).unwrap();
    assert_eq!(response.edge_count, 2);
    assert_eq!(response.max_vector_error, 0.0);
    // Same count and exactly identical Cartesian vectors, but one tuple is
    // duplicated and the reversed index tuple is missing. Vector allclose alone
    // cannot detect this control, even at zero tolerance.
    let wrong = [edges[0], edges[0]];
    assert!(edges
        .iter()
        .chain(&wrong)
        .all(|&e| displacement(&f, e) == [0.0; 3]));
    let negative = run(&wrong);
    assert!(
        !negative.status.success(),
        "PYTHONOPTIMIZE=1 accepted wrong integer tuples with the same count and zero vectors"
    );
    let stderr = String::from_utf8_lossy(&negative.stderr[..negative.stderr.len().min(8192)]);
    assert!(stderr.contains("integer-image tuple mismatch"), "{stderr}");
    assert!(
        !stderr.contains("assert_allclose"),
        "tuple guard must fail before vector comparison: {stderr}"
    );
    assert!(
        negative.stdout.is_empty(),
        "tuple mismatch must not publish a success response"
    );
    println!("optimized tuple control: valid graph accepted; same-count/same-vector wrong tuples rejected: {stderr}");
}

fn mace_request() -> serde_json::Value {
    let all = fixtures();
    let mut cases = Vec::new();
    for id in ["mace3", "gauge", "deform", "displaced"] {
        let f = all.iter().find(|f| f.id == id).unwrap();
        let edges = actual(f);
        let margin = cutoff_margin(f);
        assert!(margin >= 1e-8, "{id}: model cutoff margin");
        assert_eq!(edges, fixed_cube(f), "{id}: independent tuple preflight");
        cases.push(json!({"fixture": f, "edges": wire_edges(&edges), "margin": margin}));
    }
    let base = all.iter().find(|f| f.id == "displaced").unwrap();
    let base_edges = actual(base);
    let mut differences = Vec::new();
    for atom in 0..base.positions.len() {
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut f = base.clone();
                f.positions[atom][axis] += sign * 1e-5;
                f.id = format!("fd_{atom}_{axis}_{sign}");
                // Every displaced snapshot owns a newly built native graph.
                let edges = actual(&f);
                assert_eq!(edges, base_edges, "FD neighbor set changed");
                assert_eq!(edges, fixed_cube(&f), "FD independent tuple preflight");
                let margin = cutoff_margin(&f);
                assert!(margin >= 1e-8, "FD cutoff margin");
                differences.push(json!({"fixture": f, "edges": wire_edges(&edges),
                    "margin": margin, "atom": atom, "axis": axis, "sign": sign}));
            }
        }
    }
    json!({"schema": 1, "cases": cases, "finite_differences": differences,
        "finite_difference_h": 1e-5})
}

fn mace_run(request: &serde_json::Value) -> std::process::Output {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .unwrap();
    let python = std::env::var_os("MODENA_MACE_PYTHON")
        .expect("set MODENA_MACE_PYTHON to the pinned isolated Python 3.12.8 interpreter");
    assert!(root.join(&python).is_file(), "resolve MACE interpreter");
    // Do not resolve the executable symlink: Python must retain its venv location.
    let python_path = root.join(python);
    let python = python_path
        .parent()
        .unwrap()
        .canonicalize()
        .unwrap()
        .join(python_path.file_name().unwrap());
    let helper = root.join("crates/modena-core/tests/support/mace_parity.py");
    // Resolve the existing parent so the missing-helper RED reaches Python.
    let helper = helper
        .parent()
        .unwrap()
        .canonicalize()
        .unwrap()
        .join(helper.file_name().unwrap());
    let mut child = Command::new(&python)
        .arg("-O")
        .arg(helper)
        .current_dir(&root)
        .env("PYTHONDONTWRITEBYTECODE", "1")
        .env("OMP_NUM_THREADS", "1")
        .env("MKL_NUM_THREADS", "1")
        .env("OPENBLAS_NUM_THREADS", "1")
        .env("WANDB_MODE", "disabled")
        .env("HF_HUB_OFFLINE", "1")
        .env("CUDA_VISIBLE_DEVICES", "")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn official initialized MACE helper");
    let write_result = child
        .stdin
        .take()
        .unwrap()
        .write_all(&serde_json::to_vec(request).unwrap());
    let output = child
        .wait_with_output()
        .expect("wait for official MACE helper");
    if output.status.success() {
        write_result.expect("write complete MACE request");
    }
    output
}

#[test]
#[ignore = "requires explicitly selected pinned official MACE Python environment"]
fn official_mace_initialized_graph_energy_forces() {
    let request = mace_request();
    let output = mace_run(&request);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "MACE helper {}\n{stderr}",
        output.status
    );
    let result: serde_json::Value = serde_json::from_slice(&output.stdout)
        .unwrap_or_else(|e| panic!("malformed MACE response: {e}\n{stderr}"));
    assert_eq!(result["status"], "passed");
    assert_eq!(result["cases"].as_array().unwrap().len(), 4);
    assert_eq!(result["finite_difference_components"], 9);
    for key in [
        "model_state_sha256",
        "request_sha256",
        "lock_sha256",
        "cases_sha256",
        "helper_sha256",
    ] {
        assert!(result[key].as_str().is_some_and(sha256), "missing {key}");
    }
    println!("{}", serde_json::to_string(&result).unwrap());

    // Optimized Python must reject malformed tuples before any model forward.
    for (kind, expected) in [
        ("tuple", "integer-image tuple mismatch"),
        ("fractional", "exact integer shift required"),
        ("schema", "request schema mismatch"),
    ] {
        let mut bad = request.clone();
        if kind == "tuple" {
            bad["cases"][0]["edges"][1] = bad["cases"][0]["edges"][0].clone();
        } else if kind == "fractional" {
            bad["cases"][0]["edges"][0]["edge_cell_shift"][0] = json!(0.5);
        } else {
            bad["schema"] = json!(2);
        }
        let negative = mace_run(&bad);
        assert!(!negative.status.success(), "optimized MACE accepted {kind}");
        let stderr = String::from_utf8_lossy(&negative.stderr);
        assert!(stderr.contains(expected), "{stderr}");
        assert!(
            negative.stdout.is_empty(),
            "negative control published success"
        );
        println!("optimized MACE {kind} guard: rejected before evaluation");
    }
}
