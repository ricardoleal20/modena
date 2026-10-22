import hashlib
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys

root = Path(__file__).parent
out = root / "root-verification"
labels = ["workspace", "fmt", "clippy-default", "clippy-gpu", "clippy-mpi", "release", "scientific-build", "bench-compile", "python", "tooling", "numerical", "gpu-adapter", "mpi-build", "mpi-2rank", "original-probe-build", "original-probes"]
labels += [f"example-{name}-{action}" for name in ("molecular-nvt", "continue-nvt", "ionic-spme", "mc-npt", "water-nvt") for action in ("check", "run")]
records = [json.loads((out / f"{name}.json").read_text()) for name in labels]
assert all(r["returncode"] == 0 and r["source_status"] == "" and r["source_commit"] == "e51ff1f47bc72a08296f984604add1f2649df46e" for r in records)
by_label = {r["label"]: r for r in records}
counts = by_label["workspace"]["rust_test_counts"]
rust_counts = {"passed": sum(c[0] for c in counts), "failed": sum(c[1] for c in counts), "ignored": sum(c[2] for c in counts)}
assert rust_counts == {"passed": 266, "failed": 0, "ignored": 1}
for name in ("python", "tooling"):
    text = (out / f"{name}.stderr").read_text()
    assert "Ran 10 tests" in text and text.rstrip().endswith("OK") and "skipped=" not in text
assert by_label["gpu-adapter"]["rust_test_counts"] == [[7, 0, 0]]
for name in ("gpu-adapter", "mpi-2rank"):
    text = (out / f"{name}.stdout").read_text() + (out / f"{name}.stderr").read_text()
    assert "skipping" not in text.lower()
assert by_label["mpi-2rank"]["rust_test_counts"] == [[1, 0, 0], [1, 0, 0]]
numeric = json.loads((out / "scientific/validation.json").read_text())
assert numeric["passed"] and len(numeric["observations"]) == 67 and all(row["passed"] for row in numeric["observations"])

text = (out / "original-probes.stdout").read_text()
rows = []
def row(name, observed, expected, tolerance, basis):
    passed = math.isfinite(observed) and abs(observed - expected) <= tolerance
    rows.append(dict(name=name, observed=observed, expected=expected, absolute_tolerance=tolerance, independent_basis=basis, passed=passed))
mean, variance = re.search(r"measured_mean_T=([0-9.e+-]+) measured_T_variance=([0-9.e+-]+)", text).groups()
row("csvr_mean_T", float(mean), 1, .01, "Gamma kinetic distribution, f=12")
row("csvr_variance_T", float(variance), 1/6, .01, "Gamma temperature variance=2/f")
for constrained, value in re.findall(r"baoab constrained=(false|true).*?measured_mean_v2_per_component=([0-9.e+-]+)", text):
    row("baoab_free_atom_" + constrained, float(value), 1, .02, "Maxwell T/m=1; free spectator with or without independent dimer")
for enabled, pe, force in re.findall(r"barostat enabled=(false|true).*?abs_PE_error=([0-9.e+-]+) relative_force_l2_error=([0-9.e+-]+)", text):
    row("barostat_pe_" + enabled, float(pe), 0, 1e-9, "independent fresh force evaluation, printed precision")
    row("barostat_force_" + enabled, float(force), 0, 1e-7, "independent fresh force evaluation, printed precision")
assert len(rows) == 8 and all(r["passed"] for r in rows)
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
probe = dict(passed=True, observations=rows, source_commit=records[0]["source_commit"], library_sha256=sha("target/release/libmodena_core.rlib"), probe_source_sha256=sha("docs/scientific-hardening/evidence/original-audit/probes.rs"), limitation="raw original bug labels remain historical; printed cache precision and correlated moments are not formal uncertainty")
(out / "original-probe-observations.json").write_text(json.dumps(probe, indent=2) + "\n")

archive = Path("docs/scientific-hardening/evidence/final-review-probe")
manifest = json.loads((archive / "manifest.json").read_text())
assert all(sha(archive / name) == value for name, value in manifest["files"].items())
warnings = {name: [line for line in (out / f"{name}.stderr").read_text().splitlines() if line.startswith("warning:")] for name in labels}
assert not any(warnings.values()), warnings
summary = dict(passed=True, verified_source_commit=records[0]["source_commit"], source_status="", rust=rust_counts, python=10, tooling=10, gpu={"passed":7,"adapter_available":True,"skipped":0}, mpi={"unique_tests":1,"ranks":2,"passed_per_rank":1,"skipped":0}, numerical_observations=67, original_probe_observations=8, example_check_run_actions=10, checks=records, platform=platform.platform(), python_version=sys.version, python_executable=sys.executable, rustc=subprocess.check_output(["rtk","proxy","rustc","--version"],text=True).strip(), mpi_version=subprocess.check_output(["rtk","proxy","mpiexec","--version"],text=True).strip(), cpu=subprocess.check_output(["rtk","proxy","sysctl","-n","machdep.cpu.brand_string"],text=True).strip(), binaries={name:sha(path) for name,path in {"modena":"target/release/modena", "modena-mcp":"target/release/modena-mcp", "core_default_release":"target/release/libmodena_core.rlib", "scientific_validation":"target/release/examples/scientific_validation"}.items()}, archive_fixture_hashes_verified=len(manifest["files"]))
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({key:summary[key] for key in ("passed","verified_source_commit","rust","python","tooling","gpu","mpi","numerical_observations","original_probe_observations","example_check_run_actions","archive_fixture_hashes_verified")},indent=2))
