"""Exploratory local LJ comparison; explicit steps, raw logs, no file I/O.

Run after `rtk cargo build --release --bin modena`, with no concurrent builds/tests:
    rtk proxy python3 docs/audits/2026-10-03-modena-vs-lammps/benchmark_lj.py
This is a local serial baseline, not an HPC or best-configuration comparison.
"""
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MODENA = ROOT / "target/release/modena"
LAMMPS = shutil.which("lmp_serial")
STEPS = 1000
REPEATS = 3
A = 1.6796
RHO = 4 / A**3

def decks(cells, steps):
    modena = f'''units = "reduced"
[structure]
type = "fcc"
lattice_constant = {A}
replicas = {list(cells)}
element = "Ar"
[model]
type = "lennard_jones"
sigma = 1.0
epsilon = 1.0
cutoff = 2.5
[integrator]
type = "velocity_verlet"
dt = 0.005
[run]
steps = {steps}
initial_temperature = 1.0
seed = 87287
neighbor_skin = 0.3
log_every = {steps}
trajectory_every = {steps}
'''
    lammps = f'''units lj
atom_style atomic
boundary p p p
lattice fcc {RHO:.17g}
region box block 0 {cells[0]} 0 {cells[1]} 0 {cells[2]}
create_box 1 box
create_atoms 1 box
mass 1 1.0
pair_style lj/smooth/linear 2.5
pair_coeff * * 1.0 1.0 2.5
compute audit_temp all temp
compute_modify audit_temp extra/dof 0
velocity all create 1.0 87287 dist gaussian mom yes rot no temp audit_temp
neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes
fix 1 all nve
timestep 0.005
thermo {steps}
thermo_style custom step temp pe ke etotal press
thermo_modify temp audit_temp norm no
run {steps}
'''
    return modena, lammps

def execute(engine, cells, steps, work):
    md, lp = decks(cells, steps)
    env = os.environ.copy()
    env.pop("RAYON_NUM_THREADS", None)
    env["OMP_NUM_THREADS"] = "1"
    if engine.startswith("modena"):
        if engine == "modena_1t":
            env["RAYON_NUM_THREADS"] = "1"
        path = work / "in.toml"
        path.write_text(md)
        cmd = [str(MODENA), "run", "-i", str(path)]
    else:
        path = work / "in.lammps"
        path.write_text(lp)
        cmd = [LAMMPS, "-log", "none", "-in", str(path)]
    start = time.perf_counter()
    p = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True, check=True)
    wall = time.perf_counter() - start
    if engine.startswith("modena"):
        assert f"finished {steps} steps" in p.stdout
    else:
        assert re.search(rf"Loop time of .* for {steps} steps with", p.stdout)
    return wall, p.stdout + p.stderr

assert MODENA.exists() and LAMMPS, "Build MODENA and provide lmp_serial."
rows = []
result = {
    "sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
    "platform": platform.platform(),
    "cpu_count": os.cpu_count(),
    "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
    "lammps": next(line.strip() for line in subprocess.check_output([LAMMPS, "-h"], text=True).splitlines() if "Large-scale" in line),
    "steps": STEPS,
    "repeats": REPEATS,
    "density": RHO,
    "potential": "LJ shifted-force, sigma=epsilon=mass=1, cutoff=2.5",
    "precision": "CPU f64",
    "timing": "whole process wall time; includes setup; no trajectory/log files",
    "caveats": ["shared desktop, no idle guarantee", "three repeats only", "independent RNG initial velocities", "LAMMPS serial; default MODENA uses multiple cores", "not statistical property validation"],
    "rows": rows,
}
with tempfile.TemporaryDirectory(prefix="modena-lj-audit-") as directory:
    work = Path(directory)
    for cells in [(6, 6, 6), (8, 8, 8), (10, 20, 25)]:
        n = 4 * cells[0] * cells[1] * cells[2]
        engines = ["lammps_serial", "modena_1t", "modena_default"]
        for engine in engines:
            execute(engine, cells, 200, work) # discarded warmup
        samples = {engine: [] for engine in engines}
        drifts = {engine: [] for engine in engines}
        loads = []
        for repeat in range(REPEATS):
            loads.append(os.getloadavg())
            order = engines[repeat:] + engines[:repeat]
            for engine in order:
                wall, raw = execute(engine, cells, STEPS, work)
                (HERE / f"LJ-{n}-{engine}-r{repeat}.txt").write_text(raw)
                samples[engine].append(wall)
                if engine.startswith("modena"):
                    drifts[engine].append(float(re.search(r"drift = (\S+)", raw).group(1)))
                else:
                    thermo = [line.split() for line in raw.splitlines()
                              if re.match(rf"^\s*(0|{STEPS})\s", line) and len(line.split()) == 6]
                    e0, ef = float(thermo[0][4]), float(thermo[-1][4])
                    drifts[engine].append((ef-e0)/abs(e0))
        for engine in engines:
            times = samples[engine]
            row = {"atoms": n, "engine": engine, "wall_s": times,
                   "median_ms_per_step": statistics.median(times)*1000/STEPS,
                   "relative_endpoint_energy_drift": drifts[engine],
                   "loadavg_before_repeats": loads}
            rows.append(row)
            print(n, engine, f"{row['median_ms_per_step']:.4f} ms/step", flush=True)
(HERE / "LJ-RESULTS.json").write_text(json.dumps(result, indent=2)+"\n")
