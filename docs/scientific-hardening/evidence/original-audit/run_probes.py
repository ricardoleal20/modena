"""Build the default core and run read-only public-API scientific probes."""
import json
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
build = subprocess.run(
    ["cargo", "build", "-p", "modena-core", "--message-format=json"],
    cwd=ROOT, capture_output=True, text=True, check=True,
)
artifacts = [json.loads(line) for line in build.stdout.splitlines() if line.startswith("{")]
library = next(
    Path(filename)
    for message in artifacts
    if message.get("reason") == "compiler-artifact"
    and message.get("target", {}).get("name") == "modena_core"
    for filename in message["filenames"]
    if filename.endswith(".rlib")
)
with tempfile.TemporaryDirectory(prefix="modena-science-probes-") as temporary:
    binary = Path(temporary) / "probes"
    subprocess.run([
        "rustc", "--edition=2021", "-O", str(HERE / "probes.rs"),
        "--extern", "modena_core=" + str(library),
        "-L", "dependency=" + str(ROOT / "target/debug/deps"),
        "-o", str(binary),
    ], cwd=ROOT, check=True)
    result = subprocess.run([str(binary)], cwd=ROOT, text=True, capture_output=True, check=True)
    (HERE / "PROBE-RESULTS.txt").write_text(result.stdout)
    (HERE / "PROBE-BUILD.txt").write_text(
        "library=" + str(library) + "\n" + build.stderr +
        subprocess.check_output(["rustc", "--version"], text=True)
    )
    print(result.stdout)
