"""Synthetic durations exercise statistics only, never molecular physics."""
import copy
import importlib.util
import pathlib
import tempfile
import unittest

PATH = pathlib.Path(__file__).parents[1] / "triclinic_neighbors.py"
spec = importlib.util.spec_from_file_location("triclinic_neighbors", PATH)
campaign = importlib.util.module_from_spec(spec)
if PATH.exists():
    spec.loader.exec_module(campaign)


def table(ratio=0.8, skin=0.0):
    metadata = dict(schema_version=1, base_rev="a" * 40,
                    source_sha256="b" * 64, build_rev="a" * 40,
                    build_source_sha256="b" * 64, binary_sha256="c" * 64,
                    snapshot_sha256="d" * 64, case_id="moderate-2048",
                    family="moderate", n=2048, model="LennardJones",
                    precision="f64", cutoff_treatment="shifted_force",
                    epsilon=1.0, sigma=1.0, mass=1.0, cutoff=2.5,
                    tail_correction=False, host="test", compiler="rustc test", profile="release", features=["neighbor-diagnostics"],
                    threads=1, mode="step", skin=skin, dt=0.0005,
                    warmup_blocks=5, sample_blocks=20, block_steps=10,
                    candidate="fractional", paired=True,
                    memory_scope="requested_neighbor_workspace; retained_process_requested_live_bytes",
                    allocation_scope="process_alloc_and_realloc_calls_in_timed_kernel",
                    instrumentation="counting_global_allocator_relaxed_atomics",
                    timing_scope="complete_steps_with_per_step_diagnostic_accounting")
    records = []
    for phase, count in [("warmup", 5), ("measured", 20)]:
        for i in range(count):
            order = ["scalar", "fractional"] if i % 2 == 0 else ["fractional", "scalar"]
            for at, strategy in enumerate(order):
                records.append(dict(phase=phase, sample_index=i, order=at,
                    strategy=strategy, requested=strategy, selected="scalar" if strategy == "scalar" else "fractional",
                    fallback=None, selected_counts={"scalar" if strategy == "scalar" else "fractional": 10},
                    fallback_counts={}, block_steps=10,
                    elapsed_ns=int((1000000 + i * 1000) * (ratio if strategy != "scalar" else 1)),
                    rebuild_count=10, candidate_checks=100, n_pairs=20,
                    requested_bytes=0, retained_bytes=1000, allocation_count=0,
                    max_energy_error=0.0, max_force_error=0.0, virial_error=0.0,
                    max_position_error=0.0, max_velocity_error=0.0,
                    pair_equal=True, physics_passed=True))
    return dict(schema_version=1, status="complete", metadata=metadata,
                validation={"initial": True, "final": True, "minimum_distance": 1.0,
                            "minimum_cutoff_margin": 0.01}, records=records, error=None)


class ParserTests(unittest.TestCase):
    def test_known_paired_ratios_and_seeded_bootstrap(self):
        report = campaign.analyze_case(table())
        self.assertAlmostEqual(report["median_ratio"], 0.8, places=5)
        self.assertEqual(report["bootstrap_samples"], 10000)
        self.assertEqual(report, campaign.analyze_case(table()))
        self.assertLess(report["ci95"][1], 1.0)
        self.assertEqual(report["status"], "eligible_gain")

    def test_nonconstant_seeded_bootstrap_quantiles(self):
        # Fixed timing-ratio table; independently pinned order-statistic result.
        interval = campaign.bootstrap([0.70 + 0.02*i for i in range(20)])
        self.assertAlmostEqual(interval[0], 0.81, places=14)
        self.assertAlmostEqual(interval[1], 0.97, places=14)

    def assert_invalid(self, data):
        with self.assertRaises(campaign.DataError):
            campaign.analyze_case(data)

    def test_missing_scalar_counterpart(self):
        d = table(); d["records"].pop()
        self.assert_invalid(d)

    def test_strict_numeric_metadata(self):
        for key, bad in [("threads", True), ("n", False), ("block_steps", 9), ("dt", 0.0)]:
            d = table(); d["metadata"][key] = bad
            self.assert_invalid(d)

    def test_malformed_durations_counters_and_physics(self):
        for key, bad in [("elapsed_ns", 0), ("elapsed_ns", True), ("elapsed_ns", -1),
                         ("candidate_checks", -1), ("max_force_error", float("nan")),
                         ("max_force_error", 1e-4), ("pair_equal", False), ("physics_passed", False)]:
            d = table(); d["records"][0][key] = bad
            self.assert_invalid(d)

    def test_missing_fields_unknown_fields_order_and_mixing(self):
        for change in [lambda d: d["records"][0].pop("n_pairs"),
                       lambda d: d["records"][0].update(unrecognized=0),
                       lambda d: d["records"][0].update(order=1),
                       lambda d: d["records"][0].update(requested="auto"),
                       lambda d: d["records"][0].update(selected_counts={"fractional": 10}),
                       lambda d: d["metadata"].update(build_source_sha256="f" * 64)]:
            d = table(); change(d); self.assert_invalid(d)

    def test_child_failure_and_unavailable_sample(self):
        for status in ["invalid", "unavailable", "case_cap"]:
            d = table(); d["status"] = status; d["error"] = "child failed"
            self.assert_invalid(d)

    def test_snapshot_revision_binary_mismatch(self):
        for key in ["snapshot_sha256", "base_rev", "binary_sha256", "source_sha256"]:
            d = table()
            expected = copy.deepcopy(d["metadata"]); expected[key] = "e" * len(expected[key])
            with self.assertRaises(campaign.DataError):
                campaign.analyze_case(d, expected=expected)

    def test_skin_without_two_natural_rebuilds_is_inconclusive(self):
        d = table(skin=0.3)
        for r in d["records"]:
            r.update(rebuild_count=0, candidate_checks=0, selected_counts={}, fallback_counts={})
        self.assertEqual(campaign.analyze_case(d)["status"], "inconclusive_rebuilds")

    def test_thermal_order_effect_and_noise_are_inconclusive(self):
        d = table()
        for r in d["records"]:
            if r["strategy"] == "fractional":
                r["elapsed_ns"] = int(1000000 * (0.65 if r["sample_index"] % 2 == 0 else 1.2))
        self.assertEqual(campaign.analyze_case(d)["status"], "inconclusive_order_or_noise")

    def test_duplicate_json_keys_and_nonfinite_rejected(self):
        for raw in ['{"a":1,"a":2}', '{"a":NaN}']:
            with self.assertRaises(campaign.DataError): campaign.decode_json(raw)

    def test_deterministic_snapshots_and_explicit_elongated_matrix(self):
        d = campaign.make_snapshot(256, "elongated")
        self.assertEqual(d, campaign.make_snapshot(256, "elongated"))
        length = (256 / 0.8) ** (1 / 3)
        self.assertEqual(d["basis"], [[2 * length, 0.0, 0.0], [0.2 * length, 0.5 * length, 0.0],
                                      [0.1 * length, 0.15 * length, length]])
        self.assertLess(abs(sum(v[0] for v in d["velocities"])), 1e-14)

    def test_counterpart_pair_and_cadence_mismatch(self):
        for key, value in [("n_pairs", 21), ("rebuild_count", 9)]:
            d = table(); d["records"][0][key] = value
            self.assert_invalid(d)

    def test_rebuild_speed_cannot_be_a_fullstep_gain(self):
        d = table(); d["metadata"].update(mode="rebuild", timing_scope="forced_neighbor_builds_only_continuous_VV_advance_outside_timers")
        self.assertEqual(campaign.analyze_case(d)["status"], "rebuild_gain_only")

    def test_partial_or_false_manifest_rejected(self):
        with self.assertRaises(campaign.DataError): campaign.validate_manifest({"schema_version": True})

    def test_small_case_requires_actual_scalar_selection_for_overhead(self):
        d = table(ratio=1.0); d["metadata"]["n"] = 256
        self.assertFalse(campaign.analyze_case(d)["fallback_overhead_ok"])
        for r in d["records"]:
            r.update(selected="scalar", selected_counts={"scalar": 10})
        self.assertTrue(campaign.analyze_case(d)["fallback_overhead_ok"])


if __name__ == "__main__": unittest.main()
