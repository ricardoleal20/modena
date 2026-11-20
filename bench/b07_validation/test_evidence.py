"""Receipt and cancellation controls. The only launched programs are finite nonnumerical helpers."""

import contextlib
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import campaign
import freeze


class Evidence(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name).resolve()
        self.out = self.base / "current"
        self.out.mkdir()
        (self.out / "inputs").mkdir()
        (self.out / "jobs").mkdir()
        self.exe = str(Path(sys.executable).resolve())
        self.pin = "a" * 64
        self.clock = campaign.new_clock(20)
        freeze.dump(self.out/"protocol.json", {"test":"receipt control only"})
        freeze.dump(self.out/"external-epoch.json", {"test":"synthetic external pin"})
        self.ledger = self.base/"ledger"
        self.ledger.mkdir()
        freeze.dump(self.ledger/"absolute-clock.json",self.clock)
        freeze.dump(self.out/"clock.json",dict(self.clock,
            absolute_anchor_sha256=freeze.digest(self.ledger/"absolute-clock.json"),
            internal=freeze.digest(self.out/"protocol.json"),external=freeze.digest(self.out/"external-epoch.json")))
        self.source = {
            "source_head": "current-core",
            "compiler": "test-rustc",
            "build_profile": "release",
            "features": [],
            "old_source_head": "old-core",
            "prior_projected_draws": [65536, 65536],
            "host_python": self.exe,
            "driver_source_sha256": self.pin,
            "executable_builds": {},
        }
        self.external = {
            "expected_gate_families": {
                "pair_probes": ["test_pair_energy", "test_pair_tensor"]
            },
            "external_source_revision": "native-source",
            "external_packages": ["RIGID", "KSPACE"],
            "native_build": {
                "compiler": "test-native-compiler",
                "build_profile": "cpu-double",
                "features": ["RIGID", "KSPACE"],
            },
        }
        self.jobs = {}
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.stack.enter_context(patch.object(campaign, "LEDGER", self.ledger))
        self.stack.enter_context(
            patch.object(campaign, "verify", return_value=({}, self.external))
        )
        self.stack.enter_context(patch.object(campaign, "jobs", return_value=self.jobs))

    def tearDown(self):
        evidence = os.getenv("B07_CONTROL_EVIDENCE")
        if evidence:
            shutil.copytree(self.base, Path(evidence) / self._testMethodName)

    def job(self, name, role="rust"):
        if name not in self.jobs:
            self.jobs[name] = {
                "argv": [self.exe, "-c", "pass"],
                "role": role,
                "cwd": str(self.out),
                "env": campaign.ENV,
            }
        self.source["executable_builds"][self.exe] = {
            "source_head": "old-core"
            if name.startswith("cost-old")
            else "current-core",
            "compiler": "test-rustc",
            "build_profile": "release",
            "features": [],
            "driver_sha256": self.pin,
            "driver_source_head": "current-core",
            "linked_core_source_head": "old-core"
            if name.startswith("cost-old")
            else "current-core",
            "binary_sha256": freeze.digest(self.exe),
        }
        return self.out / "jobs" / name

    def settings(self):
        freeze.dump(self.out / "inputs/settings.json", self.source)

    def run_report(self):
        campaign.report(self.out)
        run = max((self.out / "reports").iterdir())
        return freeze.load(
            run
            / (
                "gate-ledger.json"
                if (run / "gate-ledger.json").exists()
                else "physical-report.json"
            )
        ), freeze.load(run / "work.json")

    def test_denial_before_first_child_and_cross_epoch_prefixes(self):
        self.job("deterministic")
        self.settings()
        old = self.base / "old"
        folder = old / "jobs" / "ensemble"
        folder.mkdir(parents=True)
        freeze.dump(folder / "launch.json", {"role": "rust"})
        (old / "gas-progress").write_text("2048")
        ou = old / "jobs" / "ou-draws"
        ou.mkdir()
        freeze.dump(ou / "launch.json", {"role": "rust"})
        (ou / "stdout.log").write_text("conditional seed=0 draws=8192\n")
        freeze.dump(self.ledger / "epoch-01.json", {"out": str(old)})
        with (
            patch.object(campaign, "remaining", return_value=-1),
            self.assertRaises(campaign.Budget),
        ):
            campaign.run_job(self.out, "deterministic")
        self.assertTrue((self.out / "jobs/deterministic/denied.json").is_file())
        ledger, work = self.run_report()
        statuses = {g["name"]: g["status"] for g in ledger["gates"]}
        self.assertEqual(statuses["execution_deterministic"], "BudgetExceeded")
        self.assertEqual(work["all_execution_provenance"], [])
        self.assertEqual(work["current_epoch_work"]["gas_molecule_steps"], 0)
        self.assertEqual(work["cumulative_work"]["gas_molecule_steps"], 2048)
        self.assertEqual(work["cumulative_work"]["conditional_draws"], 8192)
        self.assertTrue(work["cumulative_work"]["unknown_suffixes"])

    def test_required_gates_interrupted_rung_and_wrapper_reconciliation(self):
        self.job("deterministic")
        self.job("ensemble")
        self.settings()
        # Real no-op process; an exit zero supplies no deterministic observations.
        self.assertEqual(campaign.run_job(self.out, "deterministic"), 0)
        results = [campaign.gate("test_pair_energy", [1.0], [0.0], [0.01])]
        freeze.dump(self.out / "external-gates.json", results)
        failed = campaign.gate(
            "gas-rotor_400_endpoint_xx_intercept", [1.0], [0.0], [0.01]
        )
        freeze.dump(self.out / "gas-rotor-400-gates.json", [failed])
        ledger, _ = self.run_report()
        rows = {g["name"]: g for g in ledger["gates"]}
        self.assertEqual(rows["test_pair_tensor"]["status"], "Pending")
        self.assertEqual(
            rows["required_deterministic_physical_gates"]["status"], "Pending"
        )
        self.assertEqual(rows[failed["name"]]["status"], "Failed")
        self.assertEqual(
            rows["historical_spme_outside-excluded"]["observed"],
            [0.00024207810466960128, 0.0004829819977220807],
        )
        wrapper = self.out / "jobs/deterministic/wrapper.json"
        metadata = freeze.load(wrapper)
        metadata["exit"] = 17
        wrapper.unlink()
        freeze.dump(wrapper, metadata)
        ledger, _ = self.run_report()
        self.assertEqual(
            next(
                g["status"]
                for g in ledger["gates"]
                if g["name"] == "execution_deterministic"
            ),
            "Inconclusive",
        )
        metadata["exit"] = 0
        wrapper.unlink()
        freeze.dump(wrapper, metadata)
        (self.out / "jobs/deterministic/host.time").write_text(
            "missing resource measurements"
        )
        ledger, _ = self.run_report()
        self.assertEqual(
            next(
                g["status"]
                for g in ledger["gates"]
                if g["name"] == "execution_deterministic"
            ),
            "Inconclusive",
        )

    def test_truncated_rung_and_final_preserve_other_evidence(self):
        self.settings()
        folder = self.job("ensemble")
        folder.mkdir()
        freeze.dump(folder / "launch.json", {"role": "rust"})
        (self.out / "gas-progress").write_text("96")
        failure = campaign.gate(
            "gas-rotor_400_endpoint_xx_intercept", [1.0], [0.0], [0.01]
        )
        freeze.dump(self.out / "gas-rotor-400-gates.json", [failure])
        for filename in (
            "gas-water-400-gates.json",
            "ensemble-gates.json",
            "external-gates.json",
            "ou-gates.json",
            "deterministic-gates.json.partial",
            "gas-water-800-gates.json.partial",
        ):
            with self.subTest(file=filename):
                path = self.out / filename
                path.write_bytes(b"[{")
                try:
                    ledger, work = self.run_report()
                    rows = {r["name"]: r for r in ledger["gates"]}
                    self.assertEqual(rows[failure["name"]]["status"], "Failed")
                    defect = ledger["incomplete_gate_files"][str(path)]
                    self.assertEqual(defect["sha256"], freeze.digest(path))
                    self.assertEqual(defect["bytes"], 2)
                    self.assertTrue(defect["reason"])
                    self.assertEqual(
                        Path(defect["preserved_bytes"]).read_bytes(), b"[{"
                    )
                    self.assertEqual(rows[defect["gate"]]["status"], "Inconclusive")
                    self.assertEqual(work["cumulative_work"]["gas_molecule_steps"], 96)
                    self.assertIn(
                        str(folder), work["current_epoch_work"]["unknown_suffixes"]
                    )
                    self.assertEqual(path.read_bytes(), b"[{")
                finally:
                    path.rename(path.with_suffix(".saved"))

    def test_python_gate_publication_is_atomic_and_no_clobber(self):
        path = self.out / "external-gates.json"
        with self.assertRaises(ValueError):
            freeze.dump(path, [17, float("nan")])
        self.assertFalse(path.exists())
        self.assertTrue(
            (self.out / "external-gates.json.partial").read_bytes().startswith(b"[17")
        )
        published = self.out / "valid-gates.json"
        freeze.dump(published, [1, 2])
        with self.assertRaises(FileExistsError):
            freeze.dump(published, [3])
        self.assertEqual(freeze.load(published), [1, 2])

    def test_interrupted_cost_suffix_and_completed_cost_count(self):
        old = self.job("cost-old-rotor")
        new = self.job("cost-new-rotor")
        self.jobs[old.name]["argv"] = [self.exe, "-c", "import time; time.sleep(2)"]
        self.jobs[new.name]["argv"] = [
            self.exe,
            "-c",
            "import json; print(json.dumps(dict(steps=sum(1 for _ in range(3)),fixture='rotor',elapsed_ns=3)))",
        ]
        self.settings()
        with patch.object(campaign, "remaining", return_value=0.2):
            self.assertEqual(campaign.run_job(self.out, old.name), 2)
        self.assertEqual(campaign.run_job(self.out, new.name), 0)
        ledger, work = self.run_report()
        self.assertEqual(
            next(
                g["status"]
                for g in ledger["gates"]
                if g["name"] == "execution_cost-old-rotor"
            ),
            "BudgetExceeded",
        )
        self.assertEqual(
            campaign.measured_work(
                self.out,
                new.name,
                self.jobs[new.name],
                {"exit": 0, "budget_stop": "late wall check"},
            ),
            (3, 0, False),
        )
        self.assertEqual(
            campaign.measured_work(
                self.out, new.name, self.jobs[new.name], {"exit": -9}
            ),
            (3, 0, True),
        )
        for section in ["current_epoch_work", "cumulative_work"]:
            self.assertIn(str(old), work[section]["unknown_suffixes"])
            self.assertNotIn(str(new), work[section]["unknown_suffixes"])
            self.assertEqual(work[section]["all_molecule_steps"], 3)
        # A partially written old stdout must also leave the successful comparison readable.
        freeze.dump(self.out / "costs.json", {"steps": [{"fixture": "rotor"}]})
        self.run_report()
        comparison = freeze.load(
            max((self.out / "reports").iterdir()) / "cost-comparison.json"
        )
        self.assertEqual(set(comparison[0]["common_legacy_step"]), {"new"})

    def test_old_and_native_provenance_from_registered_build(self):
        old = self.job("cost-old-rotor")
        self.jobs[old.name]["argv"] = [
            self.exe,
            "-c",
            'print(\'{"steps":0,"fixture":"rotor"}\')',
        ]
        self.job("native-helper", "reference")
        self.source["executable_builds"][self.exe]["source_head"] = "old-core"
        self.source["executable_builds"][self.exe]["linked_core_source_head"] = (
            "old-core"
        )
        self.settings()
        for name in self.jobs:
            self.assertEqual(campaign.run_job(self.out, name), 0)
        _ledger, work = self.run_report()
        records = work["all_execution_provenance"]
        self.assertEqual(
            [(p["source_head"], p["compiler"]) for p in records],
            [("old-core", "test-rustc"), ("native-source", "test-native-compiler")],
        )
        self.assertEqual(records[1]["external_packages"], ["RIGID", "KSPACE"])
        self.assertEqual(
            work["registered_builds"][old.name]["driver_source_head"], "current-core"
        )
        self.assertEqual(
            work["registered_builds"][old.name]["linked_core_source_head"], "old-core"
        )

    def test_admission_blocks_known_live_group(self):
        self.settings()
        self.job("next")
        folder = self.out / "jobs" / "unreaped"
        folder.mkdir()
        child = subprocess.Popen(["/bin/sleep", "3"], start_new_session=True)
        try:
            freeze.dump(folder / "launch.json", {"role": "rust"})
            freeze.dump(folder / "started.json", {"pid": child.pid})
            with self.assertRaises(campaign.Budget):
                campaign.run_job(self.out, "next")
            denial = freeze.load(self.out / "jobs/next/denied.json")
            self.assertIn("live process group", denial["reason"])
        finally:
            child.kill()
            child.wait()

    def test_non_rust_embedded_source_and_host_python_library(self):
        files = freeze.cargo_sources(
            [
                {
                    "reason": "compiler-artifact",
                    "manifest_path": str(freeze.OLD / "crates/modena-core/Cargo.toml"),
                }
            ]
        )
        self.assertIn(
            freeze.OLD / "crates/modena-core/src/external/ase_emt_worker.py", files
        )
        import sysconfig
        library = (Path(self.exe).parent.parent / sysconfig.get_config_var("PYTHONFRAMEWORK")
                   if sysconfig.get_config_var("PYTHONFRAMEWORK")
                   else Path(sysconfig.get_config_var("LIBDIR")) / sysconfig.get_config_var("LDLIBRARY"))
        self.assertIn(library, freeze.python_runtime_files(Path(self.exe)))

    def test_cancelled_monitor_reaps_helper_and_refuses_unresolved_prior_launch(self):
        folder = self.job("helper")
        self.settings()
        folder.mkdir()
        work = folder / "work.json"
        script = (
            "import pathlib,time; pathlib.Path("
            + repr(str(work))
            + ').write_text(\'{"system_steps":0,"molecule_steps":0,"complete":false}\'); time.sleep(5)'
        )
        freeze.dump(
            folder / "launch.json",
            dict(
                self.jobs["helper"],
                argv=[self.exe, "-c", script],
                name="helper",
                protocol_sha256=self.pin,
                **self.clock,
            ),
        )
        proc = subprocess.Popen(
            [
                self.exe,
                *campaign.PYTHON_FLAGS,
                str(Path(campaign.__file__).resolve()),
                "_monitor",
                str(self.out),
                str(folder),
            ]
        )
        pid = None
        try:
            end = time.monotonic() + 2
            while not work.exists() and time.monotonic() < end:
                time.sleep(0.005)
            self.assertTrue(work.exists())
            pid = freeze.load(folder / "started.json")["pid"]
            proc.send_signal(signal.SIGTERM)
            proc.wait(timeout=2)
            self.assertTrue(
                (folder / "child.json").exists(),
                "handled monitor cancellation must retain the child's exact wait4 receipt",
            )
            self.assertLess(freeze.load(folder / "child.json")["exit"], 0)
            with self.assertRaises(ProcessLookupError):
                os.kill(pid, 0)
            self.assertFalse(freeze.load(work)["complete"])
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
            if pid is not None:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(pid, signal.SIGKILL)
        unknown = self.out / "jobs" / "lost"
        unknown.mkdir()
        freeze.dump(unknown / "launch.json", {"role": "rust"})
        self.job("next")
        with self.assertRaises(campaign.Budget):
            campaign.run_job(self.out, "next")

    def test_synthetic_provider_worker_is_accounted(self):
        folder = self.job("ase-six-strain")
        self.settings()
        folder.mkdir()
        import ctypes
        lib=ctypes.CDLL("/usr/lib/libproc.dylib")
        path=ctypes.create_string_buffer(4096)
        self.assertGreater(lib.proc_pidpath(os.getpid(),path,len(path)),0)
        actual=str(Path(os.fsdecode(path.value)).resolve())
        # This synthetic monitor invokes monitor directly, with no epoch admission.
        (self.out/"external-epoch.json").write_text(__import__("json").dumps({"ase_runtime":{"process_executable":actual}}))
        script = "import subprocess,sys;subprocess.run([sys.executable,'-I','-S','-c','import time;time.sleep(.2)'],check=True)"
        freeze.dump(folder/"launch.json",dict(self.jobs["ase-six-strain"],
                    argv=[self.exe,*campaign.PYTHON_FLAGS,"-c",script],
                    name="ase-six-strain",protocol_sha256=self.pin,**self.clock))
        self.assertEqual(campaign.monitor(self.out,folder),0)
        child=freeze.load(folder/"child.json")
        self.assertGreaterEqual(child["elapsed_seconds"],.2)
        self.assertLess(child["elapsed_seconds"],5)
        self.assertTrue(campaign.worker_coverage(child))
        self.assertGreater(child["sampled_ase_tree_max_rss_bytes"],0)
        self.assertTrue(all(r["rss_samples"]>0 for r in child["ase_workers"].values()))
