#!/usr/bin/env python3
"""Matched ordinary-release full-CLI pair-scratch experiment (Python >=3.11).

Retains every AB/BA sample, inputs, logs, per-child CPU/RSS and source/binary
provenance. Static multithread gates live in pair_scratch_profile; chaotic
multithread trajectories are not compared componentwise after integration.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REPLICAS = {864: [6, 6, 6], 2048: [8, 8, 8], 5000: [10, 5, 25], 10000: [10, 10, 25]}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def capture(command):
    return subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip()


def deck(atoms, steps):
    return f'''units = "reduced"
[structure]
type = "fcc"
lattice_constant = 1.6796
replicas = {REPLICAS[atoms]}
element = "Ar"
mass = 1.0
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
log_every = 100
log_file = "thermo.csv"
trajectory_every = {steps + 1}
restart_every = {steps}
restart_file = "final.checkpoint"
'''


def finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def verify_run(directory, stdout, atoms, steps):
    if re.findall(r'MODENA: finished (\d+) steps', stdout) != [str(steps)]:
        raise ValueError('completed steps do not match requested steps')
    header, payload = (directory / 'final.checkpoint').read_text().split('\n', 1)
    def reject(value):
        raise ValueError(f'nonfinite JSON: {value}')
    state = json.loads(payload, parse_constant=reject)
    if header != 'MODENA_RESTART 2' or state['step'] != steps or state['atoms'] != atoms:
        raise ValueError('checkpoint version/steps/atoms mismatch')
    expected_edges = [1.6796 * r for r in REPLICAS[atoms]]
    edges = state.get('box_edges')
    if (not isinstance(edges, list) or len(edges) != 3 or state.get('box_vectors') is not None
            or any(not finite_number(x) or not math.isclose(x, y, rel_tol=0, abs_tol=1e-12)
                   for x, y in zip(edges, expected_edges))):
        raise ValueError('checkpoint FCC box mismatch')
    if not finite_number(state.get('time')) or not math.isclose(state['time'], steps*0.005, rel_tol=1e-12, abs_tol=1e-10):
        raise ValueError('checkpoint integration time mismatch')
    continuation = state.get('continuation')
    if not isinstance(continuation, dict):
        raise ValueError('missing or malformed checkpoint continuation')
    # Validate saved state, independently of its declared physical signature.
    # Current native FCC checkpoints write explicit unspecified atomic numbers
    # (zero), although the general Rust reader defaults older missing Z to [].
    identities = [
        ('masses', 1.0, finite_number),
        ('species', 0, lambda x: type(x) is int),
        ('atomic_numbers', 0, lambda x: type(x) is int),
        ('charges', 0.0, finite_number),
        ('elements', 'Ar', lambda x: type(x) is str),
    ]
    for key, expected_value, valid_type in identities:
        values = continuation.get(key)
        if (not isinstance(values, list) or len(values) != atoms
                or any(not valid_type(x) or x != expected_value for x in values)):
            raise ValueError(f'checkpoint continuation {key} mismatch')
    for key, expected_value in [('constraint_count', 0), ('conserved_com_dof', 3)]:
        if type(continuation.get(key)) is not int or continuation[key] != expected_value:
            raise ValueError(f'checkpoint continuation {key} mismatch')
    # Velocity Verlet and the absent thermostat/barostat are stateless here.
    for key in ['integrator', 'thermostat', 'barostat']:
        if key not in continuation or continuation[key] is not None:
            raise ValueError(f'checkpoint continuation {key} state mismatch')
    for values in [state['positions'], state['velocities'], continuation['forces']]:
        if len(values) != 3 * atoms or not all(finite_number(x) for x in values):
            raise ValueError('invalid positions/velocities/forces')
    if not all(finite_number(continuation[k]) for k in ['potential_energy', 'virial']):
        raise ValueError('invalid energy/virial')
    sig = continuation['physical_signature']
    expected = dict(schema='modena-physical-v2', initialization_seed=87287, initial_temperature=1.0,
                    masses=[1.0]*atoms, species=[0]*atoms, charges=[0.0]*atoms,
                    elements=['Ar']*atoms, integrator=dict(type='velocity_verlet', dt=0.005),
                    thermostat=None, barostat=None, electrostatics=None, topology=None,
                    constraints=None, water_geometry=None, output_unit_factors=None)
    if any(k not in sig or sig[k] != v for k, v in expected.items()):
        raise ValueError('checkpoint physical signature mismatch')
    model = sig['model']
    if any(model.get(k) != v for k, v in dict(type='lennard_jones', sigma=1.0, epsilon=1.0, cutoff=2.5).items()):
        raise ValueError('checkpoint LJ physics mismatch')
    if model.get('types') or model.get('overrides'):
        raise ValueError('unexpected typed LJ physics')
    with (directory / 'thermo.csv').open() as stream:
        thermo = list(csv.reader(stream))
    columns = ['step', 'time', 'kinetic', 'potential', 'total', 'temperature', 'pressure', 'virial']
    expected_steps = list(range(0, steps + 1, 100))
    if thermo[0] != columns or len(thermo) != len(expected_steps) + 1:
        raise ValueError('thermo header/cadence mismatch')
    for step, row in zip(expected_steps, thermo[1:]):
        if len(row) != 8 or int(row[0]) != step or not all(math.isfinite(float(x)) for x in row):
            raise ValueError('invalid thermo steps/values')
    # Hash actual files: decoded float equality considers +0.0 == -0.0.
    return dict(checkpoint_sha256=sha(directory / 'final.checkpoint'),
                thermo_sha256=sha(directory / 'thermo.csv'))


def run_one(out, binary, atoms, threads, steps, repetition, variant, warmup, rows):
    name = f'n{atoms}-t{threads}-' + ('warmup' if warmup else f'r{repetition}') + f'-{variant}'
    directory = out / name
    directory.mkdir()
    input_path = directory / 'input.toml'
    input_path.write_text(deck(atoms, steps))
    env_values = dict(RAYON_NUM_THREADS=str(threads), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    command = [str(binary), 'run', '-i', str(input_path)]
    row = dict(directory=name, atoms=atoms, threads=threads, steps=steps, repetition=repetition,
               variant=variant, warmup=warmup, command=command, environment=env_values,
               input_sha256=sha(input_path), load_before=list(os.getloadavg()))
    rows.append(row)
    start = time.perf_counter()
    with (directory / 'stdout').open('w') as stdout, (directory / 'stderr').open('w') as stderr:
        process = subprocess.Popen(command, cwd=directory, env={**os.environ, **env_values},
                                   stdout=stdout, stderr=stderr)
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    row.update(returncode=process.returncode, wall_seconds=time.perf_counter()-start,
               user_seconds=usage.ru_utime, system_seconds=usage.ru_stime,
               max_rss_bytes=usage.ru_maxrss * (1 if platform.system() == 'Darwin' else 1024),
               load_after=list(os.getloadavg()))
    if process.returncode:
        raise ValueError(f'{name}: command exited {process.returncode}')
    evidence = verify_run(directory, (directory / 'stdout').read_text(), atoms, steps)
    row.update(completed_steps=steps, checkpoint_sha256=sha(directory / 'final.checkpoint'),
               thermo_sha256=sha(directory / 'thermo.csv'), physics_verified=True)
    return evidence


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', required=True, type=Path)
    parser.add_argument('--after', required=True, type=Path)
    parser.add_argument('--sizes', type=int, nargs='+', choices=REPLICAS, default=list(REPLICAS))
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 4, 8])
    parser.add_argument('--steps', type=int, default=1000)
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)
    if min(args.steps, args.repetitions, *args.threads) <= 0 or len(set(args.sizes)) != len(args.sizes) or len(set(args.threads)) != len(args.threads):
        parser.error('positive steps/repetitions/threads and unique sizes/threads required')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    artifact = dict(schema_version=1, passed=False, runs=[], summary=[], command=sys.argv,
                    created_utc=datetime.now(timezone.utc).isoformat(),
                    scope='fresh full CLI: startup, setup, integration, identical output policy; warm OS caches possible',
                    numerical_policy='1T exact checkpoint/thermo; >1T static E/W atol1e-8 rtol1e-12, all F atol1e-10 rtol1e-12 in isolated profile; scientific gates separate')
    try:
        binaries = dict(before=args.before.resolve(strict=True), after=args.after.resolve(strict=True))
        source = {p: sha(ROOT / p) for p in capture(['git', '-C', str(ROOT), 'ls-files']).splitlines() if (ROOT / p).is_file()}
        artifact['provenance'] = dict(git_commit=capture(['git', '-C', str(ROOT), 'rev-parse', 'HEAD']),
            git_status=capture(['git', '-C', str(ROOT), 'status', '--porcelain']),
            source_files=source, source_manifest_sha256=hashlib.sha256(json.dumps(source, sort_keys=True).encode()).hexdigest(),
            binaries={k:dict(path=str(v), sha256=sha(v)) for k,v in binaries.items()},
            driver_sha256=sha(__file__), rustc=capture(['rustc', '-Vv']),
            rustflags=os.environ.get('RUSTFLAGS'), cargo_encoded_rustflags=os.environ.get('CARGO_ENCODED_RUSTFLAGS'),
            python=sys.version, platform=platform.platform(), logical_cpus=os.cpu_count(),
            hardware=capture(['sysctl', '-n', 'machdep.cpu.brand_string']) if platform.system() == 'Darwin' else platform.processor())
        for atoms in args.sizes:
            for threads in args.threads:
                for repetition in range(args.repetitions + 1):
                    warmup = repetition == 0
                    order = ['before', 'after'] if repetition % 2 else ['after', 'before']
                    reference = None
                    for variant in order:
                        evidence = run_one(out, binaries[variant], atoms, threads,
                            min(50, args.steps) if warmup else args.steps, repetition, variant, warmup, artifact['runs'])
                        if reference is not None and threads == 1 and evidence != reference:
                            raise ValueError('one-thread exact checkpoint/thermo parity failed')
                        reference = evidence
                samples = [r for r in artifact['runs'] if r['atoms']==atoms and r['threads']==threads and not r['warmup']]
                medians = {v:statistics.median(r['wall_seconds'] for r in samples if r['variant']==v) for v in binaries}
                ratios = [next(r['wall_seconds'] for r in samples if r['repetition']==i and r['variant']=='before') /
                          next(r['wall_seconds'] for r in samples if r['repetition']==i and r['variant']=='after') for i in range(1,args.repetitions+1)]
                artifact['summary'].append(dict(atoms=atoms, threads=threads, measured_pairs=args.repetitions,
                    median_seconds=medians, before_over_after=medians['before']/medians['after'], paired_ratios=ratios))
        artifact['passed'] = True
    except (OSError, ValueError, KeyError, TypeError, IndexError, subprocess.SubprocessError) as error:
        artifact['error'] = str(error)
        print(f'FAILED: {error}; evidence: {out}', file=sys.stderr)
    finally:
        artifact['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (out / 'pair_scratch.json').write_text(json.dumps(artifact, indent=2, allow_nan=False)+'\n')
    if artifact['passed']:
        print(json.dumps(artifact['summary'], indent=2))
    return 0 if artifact['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
