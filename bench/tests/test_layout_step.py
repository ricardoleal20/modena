"""Real child processes exercise malformed output, rather than subprocess mocks."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("layout_step", Path(__file__).parents[1] / "layout_step.py")
try:
    driver = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(driver)
except FileNotFoundError:
    driver = None



def fixture():
    n = 256
    p, v, length = driver.initial(n)
    energy = dict(potential=0., kinetic=0., total=0., virial=0., temperature=0., pressure=0.)
    state = dict(n=n, logical_ids=list(range(n)), periodic=[True]*3,
        masses=[1.]*n, species=[0]*n, charges=[0.]*n, atomic_numbers=[0]*n,
        positions=p, velocities=v, forces=[0.]*(3*n),
        cell_vectors=[[length if i==j else 0. for j in range(3)] for i in range(3)],
        steps_completed=0, checksum_fnv1a64='fixture', energies=None,
        neighbor_pairs=None, pair_scales=None)
    ready = json.loads(json.dumps(state))
    ready.update(energies=energy, neighbor_pairs=[], pair_scales=[])
    final = json.loads(json.dumps(ready)); final['steps_completed']=1
    expected=dict(mode='soa', n=n, requested_threads=1, steps_requested=1,
        steps_completed=1, seed=17, dt=0.000001, rebuild_policy='no-rebuild')
    row=dict(expected, schema_version=1, actual_threads=1, box_lengths=[length]*3,
        workload=driver.WORKLOAD, potential=driver.POTENTIAL, setup_seconds=.01,
        elapsed_seconds=.01, pool_setup_seconds=.01, force_max_abs=0.,
        energies=energy, finite_efw=True, simd_claim=False, pair_execution='serial-scalar',
        comparison_class='serial-layout', setup_rebuilds=1, rebuilds=0,
        initial_neighbor_pairs=0, final_neighbor_pairs=0, neighbor_pair_visits=0,
        final_active_pairs=0, retained_packed_bytes=0, cluster_reconstructions=0,
        setup_cluster_reconstructions=0, evaluated_pairs=0,
        initial_checksum_fnv1a64='fixture', final_checksum_fnv1a64='fixture')
    snap=dict(schema_version=1, row={k:v for k,v in row.items() if k!='pool_setup_seconds'},
        initial_input=state, initial_force_ready=ready, final=final)
    return row,snap,expected

class DriverTests(unittest.TestCase):
    def test_real_process_failure_is_preserved(self):
        self.assertIsNotNone(driver, "campaign driver is not implemented")
        with tempfile.TemporaryDirectory() as tmp:
            row = driver.capture([sys.executable, "-c", "import sys; print(\"raw failure\"); sys.exit(7)"], Path(tmp) / "run")
            self.assertEqual(row["returncode"], 7)
            self.assertIn("raw failure", Path(row["stdout_path"]).read_text())
            self.assertGreater(row["wall_seconds"], 0)

    def test_malformed_real_child_output_is_rejected(self):
        self.assertIsNotNone(driver, "campaign driver is not implemented")
        # Removing an identity/physics check makes its corresponding child pass.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good, snap, expected = fixture()
            for mutation in ("mode", "n", "steps_completed", "missing_energy", "nonfinite", "mass", "missing_state", "snapshot_subset", "pool_time", "rebuilds", "seed", "actual_threads", "saved_nonfinite", "input_position"):
                r, s = json.loads(json.dumps(good)), json.loads(json.dumps(snap))
                if mutation in ("mode", "n", "steps_completed"):
                    r[mutation] = "wrong" if mutation == "mode" else 999
                elif mutation == "missing_energy":
                    del r["energies"]["potential"]
                elif mutation == "nonfinite":
                    r["force_max_abs"] = float("nan")
                elif mutation == "mass":
                    # Both compared children sharing this invalid mass must fail.
                    for state in (s["initial_input"], s["initial_force_ready"], s["final"]):
                        state["masses"][0] = 2
                elif mutation == "missing_state":
                    del s["initial_input"]
                elif mutation == "snapshot_subset":
                    s["row"]["dt"] = 0.5
                elif mutation == "pool_time":
                    r["pool_setup_seconds"] = -1
                elif mutation == "saved_nonfinite":
                    s["final"]["forces"][0] = float("inf")
                elif mutation == "input_position":
                    s["initial_input"]["positions"][0] += 0.1
                else:
                    r[mutation] += 1
                if mutation != "snapshot_subset":
                    s["row"] = {k:v for k,v in r.items() if k != "pool_setup_seconds"}
                child = root / (mutation + ".py")
                snapshot = root / (mutation + ".snapshot")
                child.write_text("import pathlib\npathlib.Path(" + repr(str(snapshot)) + ").write_text(" + repr(json.dumps(s)) + ")\nprint(" + repr(json.dumps(r)) + ")\n")
                for repeat in range(2):
                    with self.subTest(mutation=mutation, repeat=repeat):
                        result = driver.capture([sys.executable, str(child)], root / (mutation + str(repeat)))
                        with self.assertRaises((ValueError, KeyError, AssertionError)):
                            driver.validate(json.loads(Path(result["stdout_path"]).read_text()), driver.read_json(snapshot), expected)

    def test_valid_snapshot_subset_and_stale_output_rejection(self):
        self.assertIsNotNone(driver, "campaign driver is not implemented")
        good, snap, expected = fixture()
        driver.validate(good, snap, expected)
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "used"
            driver.capture([sys.executable, "-c", "print(1)"], prefix)
            with self.assertRaises(FileExistsError):
                driver.capture([sys.executable, "-c", "print(2)"], prefix)


if __name__ == "__main__":
    unittest.main()
