#!/usr/bin/env python3
"""MODENA driven from Python (transport A): configure in code, run, read back.

    pip install -e ./python   # (or just PYTHONPATH=python)

    PYTHONPATH=python python3 examples/python_driver.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from modena import Config, Fcc, LennardJones, Run, VelocityVerlet, check, run

def main() -> None:
    cfg = Config(
        structure=Fcc(lattice_constant=1.6796, replicas=[2, 2, 2]),
        model=LennardJones(sigma=1.0, epsilon=1.0, cutoff=2.5),
        integrator=VelocityVerlet(dt=0.005),
        run=Run(steps=400, initial_temperature=0.8, seed=42, neighbor_skin=0.3),
    )

    print("--- modena check ---")
    print(check(cfg).strip())

    print("--- modena run ---")
    res = run(cfg)
    print(f"  conservation: {res.conservation}")
    print(f"  perf:         {res.perf}")
    print(f"  thermo rows:  {len(res.steps_rows)}  (E0 = {res.steps_rows[0]['total']})")

if __name__ == "__main__":
    main()