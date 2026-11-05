import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class ValidationTests(unittest.TestCase):
    def test_scalar_json_root_retains_a_failed_validation_artifact(self):
        for root in ('null', '17', '"scalar"'):
            with self.subTest(root=root), tempfile.TemporaryDirectory() as tmp:
                executable = pathlib.Path(tmp) / 'fixture'
                executable.write_text('#!/bin/sh\nprintf \'%s\\n\' \'' + root + '\'\n')
                executable.chmod(0o755)
                out = pathlib.Path(tmp) / 'out'
                process = subprocess.run([sys.executable, str(ROOT / 'bench/validate.py'),
                    '--binary', str(executable), '--output', str(out)], capture_output=True, text=True)
                self.assertEqual(process.returncode, 1)
                self.assertTrue((out / 'validation.json').exists(), process.stderr)
                self.assertFalse(json.loads((out / 'validation.json').read_text())['passed'])

    def test_sampling_analysis_preserves_raw_series_and_reports_replica_spread(self):
        with tempfile.TemporaryDirectory() as tmp:
            executable = pathlib.Path(tmp) / 'fixture'
            payload = json.dumps({
                'observations': [{'name': 'analytic', 'observed': 1.0,
                                  'expected': 1.0, 'absolute_tolerance': 0.0}],
                'timeseries': [
                    {'name': 'kinetic', 'seed': seed, 'ensemble': 'nvt',
                     'sample_interval': 0.1, 'discard': 4,
                     'values': [99.0] * 4 + [mean - 1.0, mean + 1.0] * 100}
                    for seed, mean in [(7, 2.0), (42, 4.0)]
                ],
            })
            executable.write_text('#!/bin/sh\ncat <<\'JSON\'\n' + payload + '\nJSON\n')
            executable.chmod(0o755)
            out = pathlib.Path(tmp) / 'out'
            result = subprocess.run([sys.executable, str(ROOT / 'bench/validate.py'),
                '--binary', str(executable), '--output', str(out)], capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr.decode())
            artifact = json.loads((out / 'validation.json').read_text())
            self.assertIn('time_series_analysis', artifact)
            analysis = artifact['time_series_analysis']
            first = analysis['series'][0]['summary']
            self.assertEqual(first['discarded'], 4)
            self.assertEqual(first['samples'], 200)
            self.assertAlmostEqual(first['mean'], 2.0)
            self.assertGreater(first['effective_samples'], 0.0)
            self.assertLessEqual(first['effective_samples'], 200)
            replica = analysis['replicas'][0]
            self.assertEqual(replica['replica_mean_range'], [2.0, 4.0])
            self.assertAlmostEqual(replica['between_replica_standard_error'], 1.0)
            self.assertNotIn('timeseries', artifact, 'raw arrays belong to retained stdout')
            self.assertEqual(json.loads((out / 'validation.stdout').read_text()),
                             json.loads(payload))

    def test_invalid_sampling_series_fails_even_when_scalar_observations_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            executable = pathlib.Path(tmp) / 'fixture'
            payload = json.dumps({
                'observations': [{'name': 'analytic', 'observed': 1.0,
                                  'expected': 1.0, 'absolute_tolerance': 0.0}],
                'timeseries': [{'name': 'kinetic', 'seed': 7, 'ensemble': 'nvt',
                                'sample_interval': 0.1, 'discard': 0,
                                'values': [1.0, float('nan'), 2.0]}],
            })
            executable.write_text('#!/bin/sh\ncat <<\'JSON\'\n' + payload + '\nJSON\n')
            executable.chmod(0o755)
            out = pathlib.Path(tmp) / 'out'
            result = subprocess.run([sys.executable, str(ROOT / 'bench/validate.py'),
                '--binary', str(executable), '--output', str(out)], capture_output=True)
            self.assertEqual(result.returncode, 1)
            self.assertTrue((out / 'validation.json').exists(), result.stderr.decode())
            artifact = json.loads((out / 'validation.json').read_text())
            self.assertFalse(artifact['passed'])
            self.assertIn('finite', artifact['error'])

    def test_missing_binary_has_failure_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = pathlib.Path(tmp) / 'out'
            result = subprocess.run([sys.executable, str(ROOT / 'bench/validate.py'),
                '--binary', str(pathlib.Path(tmp) / 'missing'), '--output', str(out)], capture_output=True)
            self.assertEqual(result.returncode, 1)
            self.assertTrue((out / 'validation.json').exists())
            artifact = json.loads((out / 'validation.json').read_text())
            self.assertEqual(artifact['provenance']['returncode'], 127)
            self.assertIsNone(artifact['provenance']['binary_sha256'])

    def test_nonfinite_or_outside_tolerance_is_failure_with_raw_evidence(self):
        for value in (3.0, float('nan')):
            with tempfile.TemporaryDirectory() as tmp:
                executable = pathlib.Path(tmp) / 'fixture'
                payload = json.dumps({'observations': [{'name': 'analytic', 'observed': value,
                                      'expected': 1.0, 'absolute_tolerance': 0.1}]})
                executable.write_text('#!/bin/sh\ncat <<\'JSON\'\n' + payload + '\nJSON\n')
                executable.chmod(0o755)
                out = pathlib.Path(tmp) / 'out'
                result = subprocess.run([sys.executable, str(ROOT / 'bench/validate.py'),
                    '--binary', str(executable), '--output', str(out)], capture_output=True)
                self.assertEqual(result.returncode, 1)
                self.assertTrue((out / 'validation.stdout').exists())
                artifact = json.loads((out / 'validation.json').read_text())
                self.assertFalse(artifact['passed'])


if __name__ == '__main__':
    unittest.main()
