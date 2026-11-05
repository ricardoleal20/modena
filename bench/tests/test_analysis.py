import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analyze import analyze_records, analyze_series


class SamplingAnalysisTests(unittest.TestCase):
    def test_explicit_discard_and_mean_uncertainty_have_an_independent_expectation(self):
        result = analyze_series([99.0] * 4 + [1.0, 3.0] * 100,
                                discard=4, sample_interval=0.5)
        self.assertEqual(result['samples'], 200)
        self.assertEqual(result['discarded'], 4)
        self.assertEqual(result['mean'], 2.0)
        self.assertAlmostEqual(result['variance'], 200 / 199)
        self.assertEqual(result['statistical_inefficiency'], 1.0)
        self.assertEqual(result['effective_samples'], 200.0)
        self.assertAlmostEqual(result['standard_error'], 1 / math.sqrt(199))
        self.assertEqual(result['discard_policy'], 'explicit')

    def test_correlated_samples_have_fewer_effective_samples_and_a_wider_mean_band(self):
        values = ([0.0] * 20 + [2.0] * 20) * 16
        result = analyze_series(values, sample_interval=0.1)
        self.assertGreater(result['statistical_inefficiency'], 3.0)
        self.assertLess(result['effective_samples'], len(values) / 3)
        self.assertGreater(result['standard_error'], math.sqrt(result['variance'] / len(values)))
        self.assertGreater(result['autocorrelation_time'], 0.0)
        self.assertGreater(result['block_count'], 1)

    def test_an_energy_offset_changes_only_the_mean_and_location_of_the_interval(self):
        values = ([0.0] * 10 + [2.0] * 10) * 16
        base = analyze_series(values)
        shifted = analyze_series([x + 1024.0 for x in values])
        self.assertEqual(shifted['mean'] - base['mean'], 1024.0)
        for key in ('variance', 'statistical_inefficiency', 'effective_samples',
                    'standard_error', 'autocorrelation_time', 'half_mean_difference'):
            self.assertEqual(base[key], shifted[key], key)

    def test_constant_series_does_not_claim_many_independent_samples(self):
        result = analyze_series([2.0] * 64, discard='auto')
        self.assertEqual(result['effective_samples'], 1.0)
        self.assertIsNone(result['standard_error'])
        self.assertIsNone(result['mean_interval_95'])
        self.assertIn('constant_series', result['warnings'])

    def test_coarse_equilibration_heuristic_removes_a_large_known_initial_transient(self):
        result = analyze_series([10000.0] * 256 + [1.0, 3.0] * 512, discard='auto')
        self.assertGreaterEqual(result['discarded'], 256)
        self.assertLessEqual(result['discarded'], 384)
        self.assertEqual(result['mean'], 2.0)
        self.assertEqual(result['discard_policy'], 'auto_batch_means_heuristic')

    def test_invalid_inputs_are_rejected_instead_of_silently_dropping_samples(self):
        for values in ([], [1.0], [1.0, float('nan')], [1.0, float('inf')]):
            with self.assertRaises(ValueError):
                analyze_series(values)
        for kwargs in ({'discard': -1}, {'discard': 3}, {'sample_interval': 0.0}, {'max_lag': 0}):
            with self.assertRaises(ValueError):
                analyze_series([1.0, 2.0, 1.0, 2.0], **kwargs)

    def test_truncation_is_exposed_when_the_autocorrelation_window_is_too_short(self):
        result = analyze_series([float(i) for i in range(100)], max_lag=2)
        self.assertIn('autocorrelation_window_truncated', result['warnings'])

    def test_replicas_are_not_pooled_into_fictitious_independent_frames(self):
        records = [dict(name='K', seed=seed, ensemble='nvt', sample_interval=0.1,
                        discard=0, values=[mean - 1, mean + 1] * 100)
                   for seed, mean in [(7, 2.0), (42, 4.0), (91, 6.0)]]
        result = analyze_records(records)
        group = result['replicas'][0]
        self.assertEqual(group['replicas'], 3)
        self.assertEqual(group['mean_of_replica_means'], 4.0)
        self.assertEqual(group['replica_mean_range'], [2.0, 6.0])
        self.assertAlmostEqual(group['between_replica_standard_error'], 2 / math.sqrt(3))
        records[1]['seed'] = 7
        with self.assertRaisesRegex(ValueError, 'seed'):
            analyze_records(records)

    def test_csv_analysis_runs_against_an_actual_thermodynamic_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'thermo.csv'
            path.write_text('step,potential\n' + ''.join(
                f'{i},{value}\n' for i, value in enumerate([99] * 4 + [1, 3] * 100)))
            output = Path(tmp) / 'analysis.json'
            process = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / 'analyze.py'),
                str(path), '--observable', 'potential', '--discard', '4',
                '--sample-interval', '0.5', '--output', str(output)], capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            result = json.loads(output.read_text())
            self.assertEqual(result['summary']['mean'], 2.0)
            self.assertEqual(result['summary']['samples'], 200)
            self.assertEqual(result['provenance']['input'], str(path.resolve()))


if __name__ == '__main__':
    unittest.main()
