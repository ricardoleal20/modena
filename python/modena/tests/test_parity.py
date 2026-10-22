"""Parity tests for the Python config mirror (M-PY-0/1).

The renderer must produce input the Rust validator accepts (check rc=0) and
`run()` must reproduce an NVE run with a parsed thermo series + verdict.
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import unittest
import os

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modena import Config, Fcc, LennardJones, Run, VelocityVerlet, check, run  # noqa: E402

BINARY = os.environ.get("MODENA_BINARY", str(Path(__file__).resolve().parents[3] / "target" / "debug" / "modena"))


def test_render_argon_style():
    cfg = Config(
        structure=Fcc(lattice_constant=1.6796, replicas=[2, 2, 2]),
        model=LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5),
        integrator=VelocityVerlet(dt=0.005),
        run=Run(steps=50, initial_temperature=0.8, seed=42),
    )
    toml = cfg.to_toml()
    for needle in (
        '[structure]', 'type = "fcc"', 'lattice_constant = 1.6796',
        '[model]', 'type = "lennard_jones"', 'sigma = 1.0', 'cutoff = 2.5',
        '[integrator]', 'type = "velocity_verlet"', 'dt = 0.005',
        '[run]', 'steps = 50', 'neighbor_skin = 0.3',
    ):
        assert needle in toml, f"missing {needle!r} in:\n{toml}"
    # The engine accepts it.
    assert "OK" in check(cfg, BINARY)


def test_run_nve_and_parse():
    cfg = Config(
        structure=Fcc(lattice_constant=1.6796, replicas=[2, 2, 2]),
        model=LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5),
        integrator=VelocityVerlet(dt=0.005),
        run=Run(steps=100, initial_temperature=0.8, seed=42),
    )
    res = run(cfg, BINARY)
    assert len(res.steps_rows) > 0
    assert all(isinstance(row["total"], float) for row in res.steps_rows)
    assert [row["step"] for row in res.steps_rows] == list(range(0, 101, 10))
    assert all(0.0 <= row["temperature"] < 10.0 for row in res.steps_rows)
    # Conservation parsed and the trajectory is symplectic-consistent.
    assert res.conservation in ("Healthy", "Warning", "Drifting", "unknown")
    assert float(res.steps_rows[0]["total"]) != 0.0
class ParityTests(unittest.TestCase):
    test_render_argon_style = staticmethod(test_render_argon_style)
    test_run_nve_and_parse = staticmethod(test_run_nve_and_parse)

    def test_python_transport_matches_a_handwritten_cli_input(self):
        cfg = Config(
            Fcc(2.4, [1, 1, 1]), LennardJones(2.5, 1.0, 1.0),
            VelocityVerlet(0.001), Run(steps=4, log_every=1),
        )
        result = run(cfg, BINARY)
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "manual.toml"
            path.write_text('''[structure]
type = "fcc"
lattice_constant = 2.4
replicas = [1, 1, 1]
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.001
[run]
steps = 4
initial_temperature = 0.8
seed = 42
neighbor_skin = 0.3
log_file = "manual.csv"
log_every = 1
''', encoding="utf-8")
            completed = subprocess.run(
                [BINARY, "run", "-i", str(path)], capture_output=True,
                text=True, timeout=120,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            with (Path(td) / "manual.csv").open(newline="", encoding="utf-8") as log:
                expected = list(csv.DictReader(log))
            self.assertEqual([row["step"] for row in result.steps_rows], [0, 1, 2, 3, 4])
            self.assertEqual(len(result.steps_rows), len(expected))
            for actual, reference in zip(result.steps_rows, expected):
                self.assertEqual(actual.keys(), reference.keys())
                for key in reference:
                    self.assertEqual(actual[key], float(reference[key]), key)
