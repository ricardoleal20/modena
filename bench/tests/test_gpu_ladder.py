import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOL = ROOT / 'bench/gpu_ladder.py'


class GpuLadderTests(unittest.TestCase):
    def execute(self, executable, out, timer=None):
        command = [sys.executable, str(TOOL), '--binary', str(executable), '--ladder',
                   '864', '--engines', 'cpu', '--steps', '7', '--runs', '1',
                   '--output', str(out)]
        if timer is not None:
            command.extend(['--timer', str(timer)])
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertTrue((out / 'gpu.json').is_file(), result.stderr)
        artifact = json.loads((out / 'gpu.json').read_text())
        self.assertEqual(len(artifact['runs']), 1)
        row = artifact['runs'][0]
        self.assertEqual(row['steps'], 7)
        self.assertEqual(row['threads'], 1)
        self.assertEqual(row['dof'], 2589)
        self.assertTrue((out / row['stdout']).is_file())
        self.assertTrue((out / row['stderr']).is_file())
        return result, artifact, row

    def test_missing_executable_retains_nullable_hash_and_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            result, artifact, row = self.execute(base / 'missing', base / 'out')
            self.assertEqual(result.returncode, 1)
            self.assertEqual(row['returncode'], 127)
            self.assertIsNone(artifact['provenance']['binary_sha256'])
            self.assertTrue((base / 'out' / row['stderr']).read_text())

    def test_failed_process_without_timer_retains_raw_output_and_null_resources(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            executable = base / 'engine'
            executable.write_text('#!/bin/sh\necho fixture_stdout\necho fixture_stderr >&2\nexit 17\n')
            executable.chmod(0o755)
            result, artifact, row = self.execute(executable, base / 'out', base / 'missing_timer')
            self.assertEqual(result.returncode, 1)
            self.assertEqual(row['returncode'], 17)
            self.assertIsNotNone(artifact['provenance']['binary_sha256'])
            self.assertIsNone(row['resource_timer'])
            self.assertIsNone(row['user_cpu_seconds'])
            self.assertIsNone(row['system_cpu_seconds'])
            self.assertIsNone(row['rss_bytes'])
            self.assertIn('fixture_stdout', (base / 'out' / row['stdout']).read_text())
            self.assertIn('fixture_stderr', (base / 'out' / row['stderr']).read_text())

    def test_os_executable_launch_error_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            executable = base / 'not_executable'
            executable.write_text('#!/bin/sh\nexit 0\n')
            executable.chmod(0o644)
            result, artifact, row = self.execute(executable, base / 'out', base / 'missing_timer')
            self.assertEqual(result.returncode, 1)
            self.assertEqual(row['returncode'], 127)
            self.assertIsNotNone(artifact['provenance']['binary_sha256'])
            self.assertIn('Permission denied', (base / 'out' / row['stderr']).read_text())

    def test_os_timer_launch_error_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            executable = base / 'engine'
            executable.write_text('#!/bin/sh\nexit 0\n')
            executable.chmod(0o755)
            timer = base / 'not_executable_timer'
            timer.write_text('#!/bin/sh\nexit 0\n')
            timer.chmod(0o644)
            result, _, row = self.execute(executable, base / 'out', timer)
            self.assertEqual(result.returncode, 1)
            self.assertEqual(row['returncode'], 127)
            self.assertEqual(row['resource_timer'], str(timer))
            self.assertIn('Permission denied', (base / 'out' / row['stderr']).read_text())

    def test_success_without_timer_retains_csv_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = pathlib.Path(tmp)
            executable = base / 'engine'
            executable.write_text('#!/bin/sh\necho "run,engine,atoms,steps,wall_s,mean_ms_per_step"\n'
                                  'echo "1,cpu,864,7,0.07,10.0"\n')
            executable.chmod(0o755)
            result, _, row = self.execute(executable, base / 'out', base / 'missing_timer')
            self.assertEqual(result.returncode, 0)
            self.assertEqual(row['returncode'], 0)
            self.assertIsNone(row['resource_timer'])
            self.assertIsNone(row['rss_bytes'])
            self.assertEqual(row['integration_rows'][0]['steps'], '7')
            self.assertEqual(row['integration_rows'][0]['mean_ms_per_step'], '10.0')


if __name__ == '__main__':
    unittest.main()
