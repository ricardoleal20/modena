"""Exercise the real driver/child boundary without timing a simulation."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

DRIVER = Path(__file__).resolve().parents[1] / 'pair_scratch.py'

CHILD = '''#!{python}
import csv, json, sys, tomllib
from pathlib import Path
d=tomllib.loads(Path(sys.argv[sys.argv.index('-i')+1]).read_text())
steps=d['run']['steps']; n=4
for c in d['structure']['replicas']: n*=c
sig=dict(schema='modena-physical-v2', initialization_seed=87287, initial_temperature=1.0,
 masses=[1.0]*n,species=[0]*n,charges=[0.0]*n,elements=['Ar']*n,
 integrator=d['integrator'],thermostat=None,barostat=None,electrostatics=None,topology=None,
 constraints=None,water_geometry=None,output_unit_factors=None,model=d['model'])
mode={mode!r}
if mode=='exit': sys.exit(7)
if mode=='physics': sig['model']['epsilon']=2.0
if mode=='incomplete': steps-=1
state=dict(step=steps,atoms=n,positions=[0.0]*(3*n),velocities=[0.0]*(3*n),
 continuation=dict(forces=[0.0]*(3*n),potential_energy=-1.0,virial=2.0,physical_signature=sig))
if mode=='nan': state['continuation']['forces'][0]=float('nan')
Path('final.checkpoint').write_text('MODENA_RESTART 2\\n'+json.dumps(state))
with Path('thermo.csv').open('w') as f:
 w=csv.writer(f); w.writerow(['step','time','kinetic','potential','total','temperature','pressure','virial'])
 for i in range(0,steps+1,100): w.writerow([i,i*.005,1,-1,0,1,1,2])
print(f'MODENA: finished {{steps}} steps')
'''


class DriverTests(unittest.TestCase):
    def run_driver(self, root, mode='ok', after_mode=None):
        binaries = []
        for name, behavior in [('before', mode), ('after', after_mode or mode)]:
            path = root / name
            path.write_text(CHILD.format(python=sys.executable, mode=behavior))
            path.chmod(0o755)
            binaries.append(path)
        command = [sys.executable, str(DRIVER), '--before', str(binaries[0]), '--after', str(binaries[1]),
                   '--sizes', '864', '--threads', '1', '--steps', '100', '--repetitions', '2',
                   '--output', str(root / 'evidence')]
        result = subprocess.run(command, capture_output=True, text=True)
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
        for mode in ['exit', 'physics', 'incomplete', 'nan']:
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


if __name__ == '__main__':
    unittest.main()
