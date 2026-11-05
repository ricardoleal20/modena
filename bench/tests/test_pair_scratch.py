"""Exercise the real driver/child boundary without timing a simulation."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

DRIVER = Path(__file__).resolve().parents[1] / 'pair_scratch.py'
SPEC = importlib.util.spec_from_file_location('pair_scratch', DRIVER)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CHILD = '''#!{python}
import csv, json, sys, tomllib
from pathlib import Path
d=tomllib.loads(Path(sys.argv[sys.argv.index('-i')+1]).read_text())
steps=d['run']['steps']; n=4
for c in d['structure']['replicas']: n*=c
sig=dict(schema='modena-physical-v2', initialization_seed=87287, initial_temperature=1.0,
 masses=[1.0]*n,species=[0]*n,charges=[0.0]*n,elements=['Ar']*n,
 integrator=d['integrator'],thermostat=None,barostat=None,electrostatics=None,topology=None,
 constraints=None,water_geometry=None,output_unit_factors=None,
 model=dict(d['model'],mixing='lorentz_berthelot',overrides=[],types=[]))
mode={mode!r}
if mode=='exit': sys.exit(7)
if mode=='physics': sig['model']['epsilon']=2.0
if mode=='incomplete': steps-=1
state=dict(step=steps,atoms=n,time=steps*.005,box_edges=[1.6796*r for r in d['structure']['replicas']],box_vectors=None,positions=[0.0]*(3*n),velocities=[0.0]*(3*n),
 continuation=dict(forces=[0.0]*(3*n),masses=[1.0]*n,species=[0]*n,
 atomic_numbers=[0]*n,charges=[0.0]*n,elements=['Ar']*n,constraint_count=0,
 conserved_com_dof=3,potential_energy=-1.0,virial=2.0,physical_signature=sig,
 integrator=None,thermostat=None,barostat=None))
if mode=='saved_mass': state['continuation']['masses'][0]=2.0
if mode.startswith('missing:'): del state['continuation'][mode.split(':',1)[1]]
if mode=='nan': state['continuation']['forces'][0]=float('nan')
if mode=='box': state['box_edges'][0]*=2
if mode=='time': state['time']+=1
Path('final.checkpoint').write_text('MODENA_RESTART 2\\n'+json.dumps(state))
with Path('thermo.csv').open('w') as f:
 w=csv.writer(f); w.writerow(['step','time','kinetic','potential','total','temperature','pressure','virial'])
 for i in range(0,steps+1,100): w.writerow([i,i*.005,1,-1,0,1,1,2])
print(f'MODENA: finished {{steps}} steps')
'''


class DriverTests(unittest.TestCase):
    def run_driver(self, root, mode='ok', after_mode=None, env=None):
        binaries = []
        for name, behavior in [('before', mode), ('after', after_mode or mode)]:
            path = root / name
            path.write_text(CHILD.format(python=sys.executable, mode=behavior))
            path.chmod(0o755)
            binaries.append(path)
        command = [sys.executable, str(DRIVER), '--before', str(binaries[0]), '--after', str(binaries[1]),
                   '--sizes', '864', '--threads', '1', '--steps', '100', '--repetitions', '2',
                   '--output', str(root / 'evidence')]
        result = subprocess.run(command, capture_output=True, text=True, env=env)
        artifact = json.loads((root / 'evidence/pair_scratch.json').read_text())
        return result, artifact, command

    def test_success_retains_alternation_and_resources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, artifact, command = self.run_driver(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(artifact['passed'])
            runs = artifact['runs']
            self.assertEqual([r['variant'] for r in runs], ['after','before','before','after','after','before'])
            self.assertEqual(len(artifact['summary'][0]['paired_ratios']), 2)
            for row in runs:
                self.assertTrue(row['physics_verified'])
                self.assertEqual(row['completed_steps'], 50 if row['warmup'] else 100)
                self.assertGreater(row['max_rss_bytes'], 0)
                self.assertGreaterEqual(row['user_seconds'], 0)
                self.assertGreaterEqual(row['system_seconds'], 0)
                self.assertTrue((root / 'evidence' / row['directory'] / 'input.toml').is_file())
            self.assertEqual(len(artifact['provenance']['binaries']['before']['sha256']), 64)
            original = (root / 'evidence/pair_scratch.json').read_bytes()
            repeated = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual(original, (root / 'evidence/pair_scratch.json').read_bytes())

    def test_failures_are_nonzero_with_raw_evidence(self):
        for mode in ['exit', 'physics', 'incomplete', 'nan', 'box', 'time']:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result, artifact, _ = self.run_driver(root, after_mode=mode)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(artifact['passed'])
                self.assertIn('error', artifact)
                self.assertTrue(artifact['runs'])
                self.assertTrue((root / 'evidence' / artifact['runs'][0]['directory'] / 'stdout').exists())

    def test_one_thread_parity_rejects_full_precision_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Equal physics but a different trajectory must fail the 1T gate.
            result, _, command = self.run_driver(root)
            self.assertEqual(result.returncode, 0)
            after = root / 'after'
            after.write_text(after.read_text().replace('potential_energy=-1.0', 'potential_energy=-1.000000000000001'))
            command[-1] = str(root / 'different')
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            artifact = json.loads((root / 'different/pair_scratch.json').read_text())
            self.assertIn('exact', artifact['error'])

    def test_success_without_rtk_on_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directory = root / 'path-without-rtk'
            directory.mkdir()
            commands = ['git', 'rustc'] + (['sysctl'] if sys.platform == 'darwin' else [])
            for command in commands:
                executable = shutil.which(command)
                self.assertIsNotNone(executable, command)
                (directory / command).symlink_to(executable)
            env = {**os.environ, 'PATH': str(directory)}
            self.assertIsNone(shutil.which('rtk', path=env['PATH']))
            result, artifact, _ = self.run_driver(root, env=env)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(artifact['passed'])
            self.assertEqual(len(artifact['runs']), 6)
            self.assertTrue(all(r['physics_verified'] for r in artifact['runs']))

    def test_one_thread_parity_rejects_signed_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, _, command = self.run_driver(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            after = root / 'after'
            after.write_text(after.read_text().replace("if mode=='nan':", "state['continuation']['forces'][0]=-0.0\nif mode=='nan':"))
            command[-1] = str(root / 'signed-zero')
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            artifact = json.loads((root / 'signed-zero/pair_scratch.json').read_text())
            self.assertIn('exact', artifact['error'])

    def test_conflicting_saved_mass_fails_when_both_children_agree(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, artifact, _ = self.run_driver(Path(tmp), mode='saved_mass')
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(artifact['passed'])
            self.assertIn('masses', artifact['error'])

    def test_missing_saved_com_dof_fails_when_both_children_agree(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, artifact, _ = self.run_driver(Path(tmp), mode='missing:conserved_com_dof')
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(artifact['passed'])
            self.assertIn('conserved_com_dof', artifact['error'])

    def test_saved_continuation_requires_complete_fixed_deck_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, artifact, _ = self.run_driver(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            row = artifact['runs'][0]
            directory = root / 'evidence' / row['directory']
            checkpoint = directory / 'final.checkpoint'
            state = json.loads(checkpoint.read_text().split('\n', 1)[1])
            stdout = (directory / 'stdout').read_text()
            fields = ['masses', 'species', 'atomic_numbers', 'charges', 'elements',
                      'constraint_count', 'conserved_com_dof', 'integrator', 'thermostat', 'barostat']
            mutations = [(f'missing-{key}', lambda c, key=key: c.pop(key)) for key in fields]
            mutations += [
                ('conflicting-mass', lambda c: c['masses'].__setitem__(0, 2.0)),
                ('malformed-masses', lambda c: c.__setitem__('masses', 'not an atom array')),
                ('truncated-masses', lambda c: c['masses'].pop()),
                ('boolean-mass', lambda c: c['masses'].__setitem__(0, True)),
                ('conflicting-species', lambda c: c['species'].__setitem__(0, 1)),
                ('boolean-species', lambda c: c['species'].__setitem__(0, False)),
                ('inferred-atomic-number', lambda c: c['atomic_numbers'].__setitem__(0, 18)),
                ('boolean-atomic-number', lambda c: c['atomic_numbers'].__setitem__(0, False)),
                ('conflicting-charge', lambda c: c['charges'].__setitem__(0, 1.0)),
                ('boolean-charge', lambda c: c['charges'].__setitem__(0, False)),
                ('conflicting-element', lambda c: c['elements'].__setitem__(0, 'Ne')),
                ('constraint-count', lambda c: c.__setitem__('constraint_count', 1)),
                ('boolean-constraint-count', lambda c: c.__setitem__('constraint_count', False)),
                ('com-dof', lambda c: c.__setitem__('conserved_com_dof', 0)),
                ('floating-com-dof', lambda c: c.__setitem__('conserved_com_dof', 3.0)),
                ('integrator-state', lambda c: c.__setitem__('integrator', {})),
                ('thermostat-state', lambda c: c.__setitem__('thermostat', {})),
                ('barostat-state', lambda c: c.__setitem__('barostat', {})),
            ]
            for name, mutate in mutations:
                with self.subTest(name=name):
                    changed = copy.deepcopy(state)
                    mutate(changed['continuation'])
                    checkpoint.write_text('MODENA_RESTART 2\n' + json.dumps(changed))
                    with self.assertRaises((ValueError, KeyError, TypeError)):
                        MODULE.verify_run(directory, stdout, row['atoms'], row['steps'])
            for invalid in [None, [], 'invalid']:
                with self.subTest(continuation=invalid):
                    changed = copy.deepcopy(state)
                    changed['continuation'] = invalid
                    checkpoint.write_text('MODENA_RESTART 2\n' + json.dumps(changed))
                    with self.assertRaises((ValueError, KeyError, TypeError)):
                        MODULE.verify_run(directory, stdout, row['atoms'], row['steps'])


if __name__ == '__main__':
    unittest.main()
