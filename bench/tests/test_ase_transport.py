"""Exercise benchmark failure and evidence contracts across actual processes."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

TOOL = Path(__file__).resolve().parents[1] / 'ase_transport.py'

FIXTURE = r'''
import json, os, pathlib, sys, tomllib
identity = dict(protocol="modena-ase-emt-v1", provider="ase_emt", model="EMT",
 task="energy_forces_stress", composition="total", units=dict(positions="angstrom",
 energy="eV", forces="eV/angstrom", stress="eV/angstrom^3"),
 capabilities=["energy","forces","stress","full_pbc","provider_neighbors","conservative","total_momentum"],
 ase_version="3.29.0", numpy_version="2.5.4", model_sha256="a"*64, asap_cutoff=False)
if sys.argv[1] == 'provider-info':
 print(json.dumps(identity)); sys.exit(0)
deck=pathlib.Path(sys.argv[-1]); config=tomllib.loads(deck.read_text())
run=config['run']; mode=config['model']['transport']; n=len(config['model']['atomic_numbers'])
assert all(os.environ[k]=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','RAYON_NUM_THREADS'])
assert config['model']['ase_version']=='3.29.0' and config['model']['model_sha256']=='a'*64
assert config['integrator']['dt']==0.0001 and run['seed']==42
assert run['log_every']==1 and run['restart_every']==run['steps']
fault=os.environ.get('ASE_BENCH_FAULT')
print('fixture stdout'); print('fixture stderr',file=sys.stderr)
if fault=='exit': sys.exit(17)
steps=run['steps']-(fault=='steps')
value=2.0 if fault in ['checkpoint','thermo'] and mode=='persistent' else 1.0
cp=dict(atoms=n,step=steps,time=steps*0.0001,box_edges=[3.61]*3,box_vectors=None,
 positions=[0.0]*(3*n),velocities=[0.0]*(3*n),continuation=dict(forces=[value if fault=='checkpoint' else 1.0]*(3*n),
 potential_energy=1.0,virial=1.0,physical_signature=dict(provider=identity)))
(deck.parent/run['restart_file']).write_text('MODENA_RESTART 2\n'+json.dumps(cp)+'\n')
(deck.parent/run['log_file']).write_text('step,potential\n'+''.join(f'{i},{value if fault=="thermo" else 1.0}\n' for i in range(run['steps']+1)))
if fault=='identity' and mode=='persistent': identity['numpy_version']='wrong'
print(f'MODENA: finished {steps} steps in 1s')
print('  run metadata = '+json.dumps(dict(provider=identity,pressure_status='unconstrained_provider_stress')))
print('  energy diagnostics = '+json.dumps(dict(diagnostics=dict(samples=run['steps']+1))))
'''


class AseTransportTests(unittest.TestCase):
    def invoke(self, tmp, fault=None, out=None):
        executable = Path(tmp) / 'engine'
        executable.write_text('#!' + sys.executable + '\n' + FIXTURE)
        executable.chmod(0o755)
        out = out or Path(tmp) / 'evidence'
        result = subprocess.run([sys.executable, str(TOOL), '--modena', str(executable),
            '--python', sys.executable, '--output', str(out), '--replicas', '1',
            '--steps', '3', '--warmup-steps', '1', '--runs', '2'],
            text=True, capture_output=True,
            env={**os.environ, 'ASE_BENCH_FAULT': fault or ''})
        return result, out

    def test_success_retains_alternating_full_runs_and_exact_parity(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, out = self.invoke(tmp)
            self.assertEqual(result.returncode, 0, result.stderr)
            data = json.loads((out / 'ase_transport.json').read_text())
            self.assertTrue(data['passed'])
            measured = [r for r in data['runs'] if not r['warmup']]
            self.assertEqual([r['transport'] for r in measured],
                             ['one_shot', 'persistent', 'persistent', 'one_shot'])
            self.assertEqual(len({r['directory'] for r in data['runs']}), 6)
            self.assertTrue(all(r['completed_steps'] == 3 for r in measured))
            self.assertTrue(all(r['parity'] for r in measured))
            self.assertTrue(all(r['wall_seconds'] > 0 for r in measured))
            self.assertEqual(data['summary'][0]['measured_runs_per_transport'], 2)
            self.assertGreater(data['summary'][0]['one_shot_over_persistent'], 0)
            self.assertEqual(len(data['provenance']['binary_sha256']), 64)
            self.assertTrue(data['provenance']['source_files'])
            for row in data['runs']:
                directory = out / row['directory']
                self.assertIn('fixture stdout', (directory / 'stdout').read_text())
                self.assertIn('fixture stderr', (directory / 'stderr').read_text())

    def test_command_failure_retains_exit_and_raw_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, out = self.invoke(tmp, 'exit')
            self.assertEqual(result.returncode, 1, result.stderr)
            data = json.loads((out / 'ase_transport.json').read_text())
            self.assertFalse(data['passed'])
            self.assertEqual(data['runs'][0]['returncode'], 17)
            self.assertIn('fixture stderr', (out / data['runs'][0]['directory'] / 'stderr').read_text())

    def test_physics_identity_and_incomplete_step_failures_cannot_be_benchmarks(self):
        for fault in ['checkpoint', 'thermo', 'identity', 'steps']:
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as tmp:
                result, out = self.invoke(tmp, fault)
                self.assertEqual(result.returncode, 1, result.stderr)
                data = json.loads((out / 'ase_transport.json').read_text())
                self.assertFalse(data['passed'])
                self.assertTrue(data['error'])
                self.assertTrue(list(out.glob('*/stdout')))
                self.assertFalse(data.get('summary'))

    def test_existing_directory_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'existing'
            out.mkdir()
            sentinel = out / 'sentinel'
            sentinel.write_text('preserve')
            result, _ = self.invoke(tmp, out=out)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(list(out.iterdir()), [sentinel])
            self.assertEqual(sentinel.read_text(), 'preserve')

    def test_missing_binary_retains_failure_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'evidence'
            result = subprocess.run([sys.executable, str(TOOL), '--modena', '/missing/modena',
                '--python', sys.executable, '--output', str(out)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 1, result.stderr)
            data = json.loads((out / 'ase_transport.json').read_text())
            self.assertFalse(data['passed'])
            self.assertIn('error', data)


if __name__ == '__main__':
    unittest.main()
