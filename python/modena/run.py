"""Transport A: drive the `modena` CLI from Python (zero dependencies).

`check()` validates a config; `run()` executes it and parses the thermo log
back into arrays plus the conservation verdict.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

from .config import Config


def _binary() -> str:
    import shutil

    return shutil.which("modena") or str(Path(__file__).resolve().parents[2] / "target" / "release" / "modena")


@dataclass
class RunResult:
    """Result of `run()` — thermo series + the conservation verdict."""

    steps_rows: list[dict]  # parsed energy CSV rows
    conservation: str       # from the CLI's final line (Healthy/Warning/Drifting)
    perf: Optional[str]     # e.g. "752.77 µs/step"
    toml: str               # the rendered config (reproducibility)


def check(cfg: Config, binary: Optional[str] = None) -> str:
    """Validate `cfg` against the engine (writes nothing). Returns STDOUT."""
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        f.write(cfg.to_toml())
        path = Path(f.name)
    try:
        r = subprocess.run(
            [binary or _binary(), "check", "-i", str(path)],
            capture_output=True, text=True, timeout=120,
        )
        if r.returncode != 0:
            raise RuntimeError(f"modena check failed (rc={r.returncode}):\n{r.stdout}\n{r.stderr}")
        return r.stdout
    finally:
        path.unlink(missing_ok=True)


def run(cfg: Config, binary: Optional[str] = None) -> RunResult:
    """Run the simulation, parse the log CSV + verdict, return everything."""
    with tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False) as f:
        f.write(cfg.to_toml())
        path = Path(f.name)
    log_path = None
    try:
        with tempfile.TemporaryDirectory() as td:
            log_path = Path(td) / "energy.csv"
            # Render with the log file inside the tmpdir (no duplicate keys).
            cfg_with_log = replace(cfg, run=replace(cfg.run, log_file=str(log_path)))
            toml_text = cfg_with_log.to_toml()
            path.write_text(toml_text)
            r = subprocess.run(
                [binary or _binary(), "run", "-i", str(path)],
                capture_output=True, text=True, timeout=3600,
            )
            out = r.stdout + r.stderr
            rows: list[dict] = []
            conservation = "unknown"
            perf = None
            if log_path.exists():
                with open(log_path) as lf:
                    rows = list(csv.DictReader(lf))
            import re
            for line in out.splitlines():
                if "conservation" in line:
                    m = re.search(r"\((Healthy|Warning|Drifting)\)", line)
                    conservation = m.group(1) if m else "unknown"
                if "µs/step" in line or "us/step" in line:
                    perf = line.split("perf")[-1].strip()
            if r.returncode != 0:
                raise RuntimeError(f"modena run failed (rc={r.returncode}):\n{out}")
            return RunResult(rows, conservation.strip(), perf, toml_text)
    finally:
        path.unlink(missing_ok=True)