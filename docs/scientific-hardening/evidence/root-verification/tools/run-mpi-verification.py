import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).parent
out = root / "root-verification"
rows = [json.loads(line) for line in (out / "mpi-build.stdout").read_text().splitlines() if line.startswith("{")]
artifacts = [r for r in rows if r.get("reason") == "compiler-artifact" and r.get("target", {}).get("name") == "mpi_decomp" and "test" in r["target"]["kind"] and "mpi" in r.get("features", []) and r.get("executable")]
assert len(artifacts) == 1, artifacts
artifact = artifacts[0]
binary = Path(artifact["executable"])
assert binary.is_file() and os.access(binary, os.X_OK)
record = {"cargo_artifact": artifact, "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(), "ranks": 2, "unique_tests": 1}
(out / "mpi-executable.json").write_text(json.dumps(record, indent=2) + "\n")
command = ["rtk", "proxy", sys.executable, str(root / "verify-command.py"), "mpi-2rank", "rtk", "proxy", "mpiexec", "--oversubscribe", "-n", "2", str(binary), "--nocapture"]
raise SystemExit(subprocess.run(command).returncode)
