import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class ValidationTests(unittest.TestCase):
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
