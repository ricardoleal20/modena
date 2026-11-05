#!/usr/bin/env python3
"""Reproduce pinned EMT full-CLI transport timings; Python >=3.11, stdlib only.

Every run starts a fresh CLI/session. OS caches may be warm. No isolated
steady-evaluation or solver-speed claim follows from this experiment.
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
THREADS = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', RAYON_NUM_THREADS='1')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def strict_json(text):
    def reject(value):
        raise ValueError(f'nonfinite JSON: {value}')
    return json.loads(text, parse_constant=reject)


def capture(command, env=None):
    result = subprocess.run(command, env=env, text=True, capture_output=True, timeout=30)
    if result.returncode or not result.stdout.strip():
        raise ValueError(f'metadata command failed: {command!r}: {result.stderr}')
    return result.stdout.strip()


def validate_identity(identity):
    if not isinstance(identity, dict):
        raise ValueError('provider metadata must be an object')
    expected = dict(protocol='modena-ase-emt-v1', provider='ase_emt', model='EMT',
                    task='energy_forces_stress', composition='total', asap_cutoff=False)
    if any(identity.get(k) != v for k, v in expected.items()):
        raise ValueError('invalid provider identity')
    if identity.get('units') != dict(positions='angstrom', energy='eV',
                                     forces='eV/angstrom', stress='eV/angstrom^3'):
        raise ValueError('invalid provider units')
    if not re.fullmatch('[0-9a-f]{64}', identity.get('model_sha256', '')):
        raise ValueError('invalid provider model hash')
    for key in ['ase_version', 'numpy_version']:
        if not isinstance(identity.get(key), str) or not identity[key]:
            raise ValueError(f'missing provider {key}')
    if not {'energy', 'forces', 'stress', 'full_pbc'}.issubset(identity.get('capabilities', [])):
        raise ValueError('missing provider capabilities')


def deck(identity, python, transport, replica, steps):
    q = json.dumps
    return f'''units = "real"
[structure]
type = "fcc"
lattice_constant = 3.61
replicas = [{replica}, {replica}, {replica}]
element = "Cu"
mass = 63.546
[model]
type = "ase_emt"
python = {q(python)}
ase_version = {q(identity['ase_version'])}
model_sha256 = {q(identity['model_sha256'])}
atomic_numbers = {[29] * (4 * replica ** 3)}
composition = "total"
pbc = [true, true, true]
asap_cutoff = false
timeout_ms = 120000
transport = {q(transport)}
[integrator]
type = "velocity_verlet"
dt = 0.0001
[run]
steps = {steps}
initial_temperature = 10.0
seed = 42
neighbor_skin = 0.3
log_every = 1
log_file = "thermo.csv"
trajectory_every = {steps + 1}
restart_every = {steps}
restart_file = "final.checkpoint"
'''


def json_line(stdout, label):
    rows = [line.strip().removeprefix(label) for line in stdout.splitlines()
            if line.strip().startswith(label)]
    if len(rows) != 1:
        raise ValueError(f'expected one {label} line')
    value = strict_json(rows[0])
    if not isinstance(value, dict):
        raise ValueError(f'{label} metadata must be an object')
    return value


def verify_run(directory, stdout, identity, atoms, steps):
    completed = re.findall(r'MODENA: finished (\d+) steps', stdout)
    if completed != [str(steps)]:
        raise ValueError('completed steps do not match requested steps')
    metadata = json_line(stdout, 'run metadata = ')
    if metadata.get('provider') != identity or metadata.get('pressure_status') != 'unconstrained_provider_stress':
        raise ValueError('run metadata/provider identity mismatch')
    diagnostics = json_line(stdout, 'energy diagnostics = ')
    samples = diagnostics.get('diagnostics')
    if not isinstance(samples, dict) or samples.get('samples') != steps + 1:
        raise ValueError('full-rate diagnostic sample count mismatch')
    header, payload = (directory / 'final.checkpoint').read_text().split('\n', 1)
    checkpoint = strict_json(payload)
    if header != 'MODENA_RESTART 2' or checkpoint['step'] != steps or checkpoint['atoms'] != atoms:
        raise ValueError('checkpoint step/atom/version mismatch')
    for values in [checkpoint['positions'], checkpoint['velocities'], checkpoint['continuation']['forces']]:
        if len(values) != atoms * 3 or any(type(x) not in (int, float) or not math.isfinite(x) for x in values):
            raise ValueError('invalid checkpoint position/velocity/force data')
    for key in ['potential_energy', 'virial']:
        if not math.isfinite(checkpoint['continuation'][key]):
            raise ValueError('invalid checkpoint energy/virial')
    signature = checkpoint['continuation'].get('physical_signature')
    # Validate the pinned deck's checkpoint contract before cross-run equality:
    # the same missing physics in both transports is not evidence.
    expected = dict(schema='modena-physical-v2', initialization_seed=42,
                    masses=[1.0]*atoms, species=[0]*atoms, charges=[0.0]*atoms,
                    elements=['Cu']*atoms, barostat=None, constraints=None,
                    electrostatics=None, thermostat=None, topology=None, water_geometry=None,
                    model=dict(type='ase_emt', provider=identity, atomic_numbers=[29]*atoms,
                               pbc=[True]*3, composition='total', asap_cutoff=False))
    if not isinstance(signature, dict) or any(k not in signature or signature[k] != v
                                               for k, v in expected.items()):
        raise ValueError('invalid checkpoint physical signature')
    integ = signature.get('integrator')
    factors = signature.get('output_unit_factors')
    if (not isinstance(integ, dict) or integ.get('type') != 'velocity_verlet'
            or not isinstance(factors, list) or len(factors) != 7
            or any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0
                   for v in [integ.get('dt'), signature.get('initial_temperature'), *factors])):
        raise ValueError('invalid checkpoint integrator, temperature or unit factors')
    with (directory / 'thermo.csv').open() as stream:
        thermo = list(csv.reader(stream))
    columns = ['step', 'time', 'kinetic', 'potential', 'total', 'temperature', 'pressure', 'virial']
    if len(thermo) != steps + 2 or thermo[0] != columns:
        raise ValueError('thermo cadence/header mismatch')
    for index, row in enumerate(thermo[1:]):
        if int(row[0]) != index or len(row) != len(thermo[0]) or not all(math.isfinite(float(x)) for x in row):
            raise ValueError('invalid thermo steps/numeric values')
    # Numeric/structural equality includes full-precision final E/F and physical
    # signature, positions, velocities, and all logged thermo (not rounded stdout).
    return dict(checkpoint=checkpoint, thermo=thermo, diagnostics=diagnostics, metadata=metadata)


def run_one(out, binary, python, identity, replica, steps, repetition, mode, warmup, env, rows):
    name = f'cu{4 * replica ** 3}-' + ('warmup' if warmup else f'r{repetition}') + '-' + mode
    directory = out / name
    directory.mkdir()
    input_path = directory / 'input.toml'
    input_path.write_text(deck(identity, python, mode, replica, steps))
    command = [binary, 'run', '--diagnostics-json', '-i', str(input_path)]
    row = dict(directory=name, transport=mode, atoms=4 * replica ** 3, replica=replica,
               steps=steps, repetition=repetition, warmup=warmup, command=command,
               environment=THREADS, input_sha256=sha(input_path), parity=False,
               load_before=list(os.getloadavg()))
    start = time.perf_counter()
    try:
        process = subprocess.run(command, cwd=directory, env=env, text=True, capture_output=True)
        row.update(returncode=process.returncode)
        stdout, stderr = process.stdout, process.stderr
    except OSError as error:
        row.update(returncode=127)
        stdout, stderr = '', str(error)
    row.update(wall_seconds=time.perf_counter() - start, load_after=list(os.getloadavg()))
    (directory / 'stdout').write_text(stdout)
    (directory / 'stderr').write_text(stderr)
    rows.append(row)
    if row['returncode']:
        raise ValueError(f'{name}: command exited {row["returncode"]}')
    evidence = verify_run(directory, stdout, identity, row['atoms'], steps)
    row.update(completed_steps=steps, checkpoint_sha256=sha(directory / 'final.checkpoint'),
               thermo_sha256=sha(directory / 'thermo.csv'))
    return row, evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modena', type=Path, default=ROOT / 'target/release/modena')
    parser.add_argument('--python', required=True, help='exact optional ASE venv executable; symlink preserved')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--replicas', type=int, nargs='+', choices=[1, 3], default=[1, 3])
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--warmup-steps', type=int, default=10)
    parser.add_argument('--runs', type=int, default=5)
    args = parser.parse_args()
    if min(args.steps, args.warmup_steps, args.runs) <= 0 or len(set(args.replicas)) != len(args.replicas):
        parser.error('steps, warmup and runs must be positive; replicas must be unique')
    out = args.output.resolve()
    try:
        out.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        parser.error('output directory already exists; use a fresh path')
    artifact = dict(schema_version=1, passed=False, runs=[], summary=[],
                    scope='fresh full CLI processes including startup/setup/integration/output; OS caches may be warm',
                    warmup='separate discarded fresh process per transport and atom count',
                    command=sys.argv, created_utc=datetime.now(timezone.utc).isoformat())
    env = {**os.environ, **THREADS}
    try:
        binary = str(args.modena.resolve(strict=True))
        # abspath intentionally preserves the venv executable symlink.
        python = os.path.abspath(args.python)
        source_paths = capture(['git', '-C', str(ROOT), 'ls-files', '-z']).split('\0')
        source_files = {p: sha(ROOT / p) for p in source_paths if p and (ROOT / p).is_file()}
        artifact['provenance'] = dict(git_commit=capture(['git', '-C', str(ROOT), 'rev-parse', 'HEAD']),
            git_status=subprocess.run(['git', '-C', str(ROOT), 'status', '--porcelain'],
                                      text=True, capture_output=True, check=True).stdout,
            source_files=source_files, source_manifest_sha256=hashlib.sha256(json.dumps(source_files, sort_keys=True).encode()).hexdigest(),
            driver_sha256=sha(__file__), binary=binary, binary_sha256=sha(binary),
            rustc=capture(['rustc', '-Vv']), python=capture([python, '--version']),
            python_executable=python, driver_python=sys.version, platform=platform.platform(),
            machine=platform.machine(), logical_cpus=os.cpu_count(), environment=THREADS,
            hardware=capture(['sysctl', '-n', 'machdep.cpu.brand_string']) if platform.system() == 'Darwin' else capture(['uname', '-a']),
            load_before=list(os.getloadavg()))
        command = [binary, 'provider-info', '--python', python, '--timeout-ms', '120000']
        discovery = subprocess.run(command, env=env, text=True, capture_output=True, check=False)
        (out / 'provider-info.stdout').write_text(discovery.stdout)
        (out / 'provider-info.stderr').write_text(discovery.stderr)
        artifact['discovery'] = dict(command=command, returncode=discovery.returncode)
        if discovery.returncode:
            raise ValueError('provider discovery failed')
        identity = strict_json(discovery.stdout)
        validate_identity(identity)
        artifact['identity'] = identity
        for replica in args.replicas:
            baseline = {}
            for repetition in range(args.runs + 1):
                warmup = repetition == 0
                modes = ['one_shot', 'persistent'] if repetition % 2 else ['persistent', 'one_shot']
                for mode in modes:
                    row, evidence = run_one(out, binary, python, identity, replica,
                        args.warmup_steps if warmup else args.steps, repetition, mode, warmup, env, artifact['runs'])
                    if warmup not in baseline:
                        baseline[warmup] = evidence
                    if evidence != baseline[warmup]:
                        raise ValueError(f'{row["directory"]}: checkpoint/thermo/diagnostics/metadata parity failure')
                    row['parity'] = True
            measured = [r for r in artifact['runs'] if r['replica'] == replica and not r['warmup']]
            medians = {m: statistics.median(r['wall_seconds'] for r in measured if r['transport'] == m)
                       for m in ['one_shot', 'persistent']}
            artifact['summary'].append(dict(atoms=4 * replica ** 3, measured_runs_per_transport=args.runs,
                median_seconds=medians, one_shot_over_persistent=medians['one_shot'] / medians['persistent']))
        artifact['passed'] = True
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        artifact.update(error=str(error), summary=[])
        print(f'FAILED: {error}; retained evidence: {out}', file=sys.stderr)
    finally:
        artifact['finished_utc'] = datetime.now(timezone.utc).isoformat()
        (out / 'ase_transport.json').write_text(json.dumps(artifact, indent=2, allow_nan=False) + '\n')
    if artifact['passed']:
        print(json.dumps(artifact['summary'], indent=2))
    return 0 if artifact['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
