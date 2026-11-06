//! Explicit external gate: ignored compilation alone is never G1 acceptance.
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
