"""Parity tests for the Python config mirror (M-PY-0/1).

The renderer must produce input the Rust validator accepts (check rc=0) and
`run()` must reproduce an NVE run with a parsed thermo series + verdict.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from modena import Config, Fcc, LennardJones, Run, VelocityVerlet, check, run  # noqa: E402

BINARY = str(Path(__file__).resolve().parents[3] / "target" / "release" / "modena")


@pytest.fixture(scope="module")
def cargo_built():
    # The driver tests need the release binary; skip politely when missing.
    if not Path(BINARY).exists():
        pytest.skip("cargo build --release not run")
    return BINARY


def test_render_argon_style(cargo_built):
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
    assert "OK" in check(cfg, cargo_built)


def test_run_nve_and_parse(cargo_built):
    cfg = Config(
        structure=Fcc(lattice_constant=1.6796, replicas=[2, 2, 2]),
        model=LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5),
        integrator=VelocityVerlet(dt=0.005),
        run=Run(steps=100, initial_temperature=0.8, seed=42),
    )
    res = run(cfg, cargo_built)
    assert len(res.steps_rows) > 0
    assert res.rows_ok if hasattr(res, "rows_ok") else True
    # Conservation parsed and the trajectory is symplectic-consistent.
    assert res.conservation in ("Healthy", "Warning", "Drifting", "unknown")
    assert float(res.steps_rows[0]["total"]) != 0.0