#!/usr/bin/env python3
"""Failure-sensitive whole-process LJ benchmark; Python >=3.11, stdlib only."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REPLICAS = {864: (6, 6, 6), 2048: (8, 8, 8), 5000: (5, 10, 25),
            10000: (10, 10, 25), 20000: (10, 20, 25)}


def capture(command):
    try:
        result = subprocess.run(command, text=True, capture_output=True, timeout=10)
        return result.stdout.strip() or result.stderr.strip()
    except (OSError, subprocess.TimeoutExpired):
        return None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decks(atoms, steps, seed):
    reps = REPLICAS[atoms]
    modena = f'''units = "reduced"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = {list(reps)}
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
seed = {seed}
neighbor_skin = 0.3
trajectory_every = {steps + 1}
log_every = {steps + 1}
'''
    lammps = f'''units lj
atom_style atomic
boundary p p p
lattice fcc {4 / 1.6796 ** 3:.17g}
region box block 0 {reps[0]} 0 {reps[1]} 0 {reps[2]}
create_box 1 box
create_atoms 1 box
mass 1 1.0
pair_style lj/smooth/linear 2.5
pair_coeff 1 1 1.0 1.0 2.5
velocity all create 1.0 {seed} dist gaussian mom yes rot no
neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes
fix 1 all nve
timestep 0.005
thermo {steps}
thermo_style custom step temp pe ke etotal press
thermo_modify norm no flush yes
run {steps}
'''
    return modena, lammps


def resources(text):
    if platform.system() == 'Darwin':
        user = re.search(r'([\d.]+)\s+user', text)
        system = re.search(r'([\d.]+)\s+sys', text)
        rss = re.search(r'(\d+)\s+maximum resident set size', text)
        return (float(user[1]) if user else None, float(system[1]) if system else None,
                int(rss[1]) if rss else None)
    user = re.search(r'User time \(seconds\): ([\d.]+)', text)
    system = re.search(r'System time \(seconds\): ([\d.]+)', text)
    rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', text)
    return (float(user[1]) if user else None, float(system[1]) if system else None,
            int(rss[1]) * 1024 if rss else None)


def drift(stdout, engine):
    if engine == 'modena':
        value = re.search(r'conservation\s+drift\s*=\s*([+\-\d.eE]+)', stdout)
        return float(value[1]) if value else None
    rows = re.findall(r'^\s*\d+\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*$', stdout, re.M)
    return (float(rows[-1][3]) - float(rows[0][3])) / abs(float(rows[0][3])) if len(rows) >= 2 and float(rows[0][3]) else None


def run_one(out, engine, binary, atoms, steps, seed, threads, repetition, warmup):
    label = f'{engine}-{atoms}-t{threads}-' + ('warmup' if warmup else f'r{repetition}')
    deck = out / (label + ('.toml' if engine == 'modena' else '.lmp'))
    deck.write_text(decks(atoms, steps, seed)[engine == 'lammps'])
    command = ([binary, 'run', '-i', str(deck)] if engine == 'modena' else
               [binary, '-in', str(deck), '-log', 'none'])
    env = os.environ.copy()
    env['RAYON_NUM_THREADS'] = str(threads)
    env['OMP_NUM_THREADS'] = '1'
    timer = '/usr/bin/time'
    timed = [timer, '-l' if platform.system() == 'Darwin' else '-v', *command] if Path(timer).exists() else command
    start = time.perf_counter()
    try:
        result = subprocess.run(timed, cwd=out, env=env, text=True, capture_output=True)
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
    except OSError as error:
        code, stdout, stderr = 127, '', str(error)
    wall = time.perf_counter() - start
    (out / (label + '.stdout')).write_text(stdout)
    (out / (label + '.stderr')).write_text(stderr)
    user, system, rss = resources(stderr)
    return dict(engine=engine, atoms=atoms, steps=steps, seed=seed, threads=threads,
                precision='f64', dof=3 * atoms - 3, conserved_com_dof=3, repetition=repetition,
                warmup=warmup, returncode=code, wall_seconds=wall, user_cpu_seconds=user,
                system_cpu_seconds=system, rss_bytes=rss, relative_energy_drift=drift(stdout, engine),
                input=deck.name, input_sha256=sha(deck), command=command,
                stdout=label + '.stdout', stderr=label + '.stderr')


def render(runs):
    print('engine        threads atoms  n   median_ms/step')
    keys = sorted({(r['engine'], r['threads'], r['atoms']) for r in runs if not r['warmup'] and r['returncode'] == 0})
    for engine, threads, atoms in keys:
        group = [r for r in runs if (r['engine'], r['threads'], r['atoms']) == (engine, threads, atoms) and not r['warmup'] and r['returncode'] == 0]
        median = statistics.median(r['wall_seconds'] * 1000 / r['steps'] for r in group)
        print(f'{engine:13s} {threads:7d} {atoms:5d} {len(group):2d} {median:17.6f}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modena', default=str(ROOT / 'target/release/modena'))
    parser.add_argument('--lammps', help='optional serial binary, same resources as --threads 1 only')
    parser.add_argument('--ladder', nargs='+', type=int, choices=REPLICAS, default=[864, 2048, 20000])
    parser.add_argument('--threads', nargs='+', type=int, default=[1])
    parser.add_argument('--steps', type=int, default=1000)
    parser.add_argument('--warmup-steps', type=int, default=200)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--seed', type=int, default=87287)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.steps, args.warmup_steps, args.runs, *args.threads) <= 0:
        parser.error('steps, warmup, runs and thread counts must be positive')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'benchmark.json').exists():
        parser.error('output already contains benchmark.json; use a fresh directory')
    binaries = {'modena': str(Path(args.modena).resolve())}
    if args.lammps:
        binaries['lammps'] = shutil.which(args.lammps) or str(Path(args.lammps).resolve())
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), git_commit=capture(['git', '-C', str(ROOT), 'rev-parse', 'HEAD']),
        git_status=capture(['git', '-C', str(ROOT), 'status', '--porcelain']),
        rustc=capture(['rustc', '--version']), python=sys.version, python_executable=sys.executable,
        platform=platform.platform(), machine=platform.machine(), logical_cpus=os.cpu_count(),
        hardware=capture(['sysctl', '-n', 'machdep.cpu.brand_string']) if platform.system() == 'Darwin' else capture(['lscpu']),
        binaries={k: {'path': v, 'sha256': sha(v) if Path(v).is_file() else None, 'version': capture([v, '--version'] if k == 'modena' else [v, '-h'])} for k, v in binaries.items()}, load_before=os.getloadavg() if hasattr(os, 'getloadavg') else None)
    artifact = dict(schema_version=1, provenance=provenance, scope='whole process: startup, initialization and integration; no trajectory or thermo files; stdout captured',
        warmup='separate fresh process, discarded; not a continued trajectory',
        hamiltonian=dict(potential='LJ shifted-force', sigma=1, epsilon=1, mass=1, cutoff=2.5,
            skin=0.3, lattice_constant=1.6796, dt=0.005, temperature=1, units='reduced',
            initialization='Gaussian, conserved COM removed, rescaled to 3N-3 DOF; engine RNG algorithms differ'),
        comparison='serial LAMMPS vs MODENA threads=1 only; other thread counts are unmatched', runs=[])
    for atoms in args.ladder:
        for engine, binary in binaries.items():
            for threads in (args.threads if engine == 'modena' else [1]):
                for repetition in range(args.runs + 1):
                    row = run_one(out, engine, binary, atoms, args.warmup_steps if repetition == 0 else args.steps,
                                  args.seed, threads, repetition, repetition == 0)
                    artifact['runs'].append(row)
                    (out / 'benchmark.json').write_text(json.dumps(artifact, indent=2) + '\n')
                    if row['returncode']:
                        print(f"FAILED {engine}: exit {row['returncode']}; retained logs in {out}", file=sys.stderr)
                        return 1
    render(artifact['runs'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
