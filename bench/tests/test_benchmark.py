import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOL = ROOT / 'bench/benchmark.py'


class BenchmarkTests(unittest.TestCase):
    def fixture(self, exit_code):
        with tempfile.TemporaryDirectory() as tmp:
            out = pathlib.Path(tmp) / 'evidence'
            executable = pathlib.Path(tmp) / 'engine'
            executable.write_text('#!/bin/sh\necho fixture_stdout\necho fixture_stderr >&2\nexit ' + str(exit_code) + '\n')
            executable.chmod(0o755)
            result = subprocess.run([sys.executable, str(TOOL), '--modena', str(executable),
                '--ladder', '864', '--steps', '7', '--warmup-steps', '2', '--runs', '1',
                '--threads', '1', '--output', str(out)], capture_output=True, text=True)
            self.assertTrue((out / 'benchmark.json').exists(), result.stderr)
            artifact = json.loads((out / 'benchmark.json').read_text())
            self.assertEqual(result.returncode, 0 if exit_code == 0 else 1)
            self.assertEqual(artifact['runs'][0]['returncode'], exit_code)
            self.assertEqual(artifact['runs'][0]['steps'], 2)
            self.assertIn('fixture_stdout', (out / artifact['runs'][0]['stdout']).read_text())
            self.assertIn('fixture_stderr', (out / artifact['runs'][0]['stderr']).read_text())
            self.assertNotEqual(artifact['runs'][0]['rss_bytes'], 0)
            self.assertEqual(artifact['provenance']['python'], sys.version)
            if exit_code == 0:
                measured = artifact['runs'][1]
                self.assertEqual(measured['steps'], 7)
                self.assertEqual(measured['threads'], 1)
                self.assertEqual(measured['precision'], 'f64')
                self.assertEqual(measured['dof'], 2589)
                self.assertIn('steps = 7', (out / measured['input']).read_text())
                self.assertIn('864', result.stdout)
                self.assertNotIn('trajectory_file', (out / measured['input']).read_text())

    def test_failure_retains_evidence_and_fails(self):
        self.fixture(17)

    def test_success_has_actual_steps_numeric_renderer_and_provenance(self):
        self.fixture(0)

    def test_real_release_cli_uses_valid_no_file_deck_and_reports_drift(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run([sys.executable, str(TOOL), '--ladder', '864',
                '--steps', '3', '--warmup-steps', '1', '--runs', '1', '--threads', '1',
                '--output', tmp], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            artifact = json.loads((pathlib.Path(tmp) / 'benchmark.json').read_text())
            self.assertIsNotNone(artifact['runs'][1]['relative_energy_drift'])
            self.assertFalse(list(pathlib.Path(tmp).glob('*.xyz')))
            self.assertFalse(list(pathlib.Path(tmp).glob('*.csv')))


if __name__ == '__main__':
    unittest.main()
