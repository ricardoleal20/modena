"""Transport A: drive the `modena` CLI from Python (zero dependencies).

`check()` validates a config; `run()` executes it and parses the thermo log
back into arrays plus the conservation verdict.
"""

from __future__ import annotations

import copy
import csv
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

from .config import Config


def _binary() -> str:
    import shutil

    if found := shutil.which("modena"):
        return found
    root = Path(__file__).resolve().parents[2]
    targets = [root / "target"]
    if target_dir := os.environ.get("CARGO_TARGET_DIR"):
        targets.insert(0, Path(target_dir))
    filename = "modena.exe" if os.name == "nt" else "modena"
    for target in targets:
        for profile in ("release", "debug"):
            candidate = target / profile / filename
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate.resolve())
    raise FileNotFoundError(
        "Cannot find the MODENA CLI. Install modena on PATH, build it with "
        "`cargo build -p modena-cli --release`, or pass binary='/path/to/modena'."
    )


@dataclass
class RunResult:
    """Result of `run()` — thermo series + the conservation verdict."""

    steps_rows: list[dict]  # parsed energy CSV rows
    conservation: str       # from the CLI's final line (Healthy/Warning/Drifting)
    perf: Optional[str]     # e.g. "752.77 µs/step"
    toml: str               # the rendered config (reproducibility)


def _resolve_paths(cfg: Config, base_dir: Optional[Path] = None) -> Config:
    """Resolve sources and intentional outputs before relocating TOML to tempdir."""
    cfg = copy.deepcopy(cfg)
    base = Path(base_dir or Path.cwd()).resolve()
    for obj, names in ((cfg.structure, ("file",)), (cfg.model, ("file",)),
                       (cfg.run, ("restore_file", "restart_file", "trajectory_file", "log_file"))):
        for name in names:
            value = getattr(obj, name, None)
            if value is not None:
                setattr(obj, name, str((base / value).resolve()))
    return cfg


def check(cfg: Config, binary: Optional[str] = None, *, base_dir: Optional[Path] = None) -> str:
    """Validate `cfg` against the engine (writes nothing). Returns STDOUT."""
    cfg = _resolve_paths(cfg, base_dir)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "input.toml"
        path.write_text(cfg.to_toml(), encoding="utf-8")
        r = subprocess.run(
            [binary or _binary(), "check", "-i", str(path)],
            capture_output=True, text=True, timeout=120,
        )
        if r.returncode != 0:
            raise RuntimeError(f"modena check failed (rc={r.returncode}):\n{r.stdout}\n{r.stderr}")
        return r.stdout


def run(cfg: Config, binary: Optional[str] = None, *, base_dir: Optional[Path] = None) -> RunResult:
    """Run the simulation and return numeric thermo rows plus its verdict.

    Source, trajectory, and checkpoint paths use `base_dir` (default: cwd).
    The thermo CSV is temporary; a configured `run.log_file` is left intact.
    `RunResult.toml` records the resolved input actually passed to the CLI.
    """
    cfg = _resolve_paths(cfg, base_dir)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "input.toml"
        log_path = Path(td) / "energy.csv"
        # Render with the log file inside the tmpdir (no duplicate keys).
        cfg_with_log = replace(cfg, run=replace(cfg.run, log_file=str(log_path)))
        toml_text = cfg_with_log.to_toml()
        path.write_text(toml_text, encoding="utf-8")
        r = subprocess.run(
            [binary or _binary(), "run", "-i", str(path)],
            capture_output=True, text=True, timeout=3600,
        )
        out = r.stdout + r.stderr
        if r.returncode != 0:
            raise RuntimeError(f"modena run failed (rc={r.returncode}):\n{out}")
        rows: list[dict] = []
        conservation = "unknown"
        perf = None
        if log_path.exists():
            with log_path.open(encoding="utf-8", newline="") as lf:
                rows = [
                    {key: int(value) if key == "step" else float(value)
                     for key, value in row.items()}
                    for row in csv.DictReader(lf)
                ]
        for line in out.splitlines():
            if "conservation" in line:
                m = re.search(r"\((Healthy|Warning|Drifting)\)", line)
                conservation = m.group(1) if m else "unknown"
            if "µs/step" in line or "us/step" in line:
                perf = line.split("perf")[-1].strip()
        return RunResult(rows, conservation, perf, toml_text)
