"""No scientific campaign: launcher stubs and tiny recorded monitor controls."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
SPEC = importlib.util.spec_from_file_location("spme_campaign", Path(__file__).resolve().parents[1] / "spme_calibration_campaign.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

class LauncherTests(unittest.TestCase):
    def test_recorded_tiny_time_fixture(self):
        fixture=json.loads((Path(__file__).parent/"fixtures"/"spme-calibration-time-tiny.json").read_text())
        self.assertEqual(M.parse_rss(fixture["stderr"]),fixture["peak_rss_bytes"])
        self.assertEqual(len(fixture["source_sha256"]),64)

    def test_numeric_mac_rss(self):
        self.assertEqual(M.parse_rss("  12345678  maximum resident set size\n"), 12345678)
        self.assertIsNone(M.parse_rss("maximum resident set size: unavailable"))
        self.assertIsNone(M.parse_rss("-1 maximum resident set size"))
        self.assertFalse(M.rss_accepted(536870913))
        self.assertFalse(M.rss_accepted(None))
    def test_roster_frozen(self):
        self.assertEqual(M.FAMILIES, [("cubic",2),("ortho",2),("special",2),("cubic",3),("ortho",3)])
        self.assertEqual(M.THREADS, [1,4])
        with self.assertRaises(ValueError): M.validate_threads([4,1])
    def test_admission_replaced_not_added(self):
        first = {"reference":{"calls":7,"work":{"real_terms":70,"atom_k_terms":20},"elapsed_ns":3},"accuracy_calls":0,"candidate_real_pair_upper_bound":0,"process_requested_rust_heap_bytes":0,"peak_process_requested_rust_heap_bytes":0}
        last = json.loads(json.dumps(first)); last["reference"]["calls"] = 6
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"events.jsonl"
            path.write_text(json.dumps({"sequence":0,"manifest_sha256":"e","kind":"candidate","data":{"Admitted":{"reserved":first}}})+"\n"+json.dumps({"sequence":1,"manifest_sha256":"e","kind":"candidate","data":{"CandidateRow":{"usage":last}}})+"\n")
            state=M.recover_usage(path, M.empty_usage(), "e")
            self.assertEqual(state["candidate"]["reference"]["calls"],6)
    def test_monitor_actual_exit_independent(self):
        with tempfile.TemporaryDirectory() as tmp:
            receipt=Path(tmp)/"child.json"
            import sys
            M.monitor(receipt, 5., [sys.executable,"-c","raise SystemExit(87)"])
            row=json.loads(receipt.read_text())
            self.assertEqual(row["actual_child_exit"],87)
            self.assertFalse(row["timed_out"])
    def test_create_new(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"a.json"; M.write_new(path,{"a":1})
            with self.assertRaises(FileExistsError): M.write_new(path,{"a":2})
            self.assertEqual(json.loads(path.read_text()),{"a":1})

    def test_stub_campaign_roster_target4_and_independent_failure(self):
        calls=[]
        def stub(command, env, out, phase, remaining):
            out.mkdir(); threads=int(env["RAYON_NUM_THREADS"])
            flags=dict(zip(command[2::2],command[3::2]))
            family=out.parent.name.split("-")[0]
            calls.append((family,threads,phase,flags.copy(),remaining))
            code=0
            if phase=="prepare":
                prepared=Path(flags["--out"]); prepared.mkdir()
                M.write_new(prepared/"manifest.json",{"source_sha256":{"launcher":"fixture"}})
                M.write_new(prepared/"protocol.json",{"extended":True}); M.write_new(prepared/"snapshots.json",[])
                self.assertIn("--reference-extended",command)
            elif phase=="verify":
                prepared=Path(flags["--protocol"]).parent
                M.write_new(prepared/"verification.json",{"usage":M.empty_usage()})
            elif phase in ("sweep","confirm","profile"):
                dest=Path(flags["--out"]); dest.mkdir()
                M.write_new(dest/"summary.json",{"usage":M.empty_usage()})
                if phase=="sweep": code=4
                if family=="ortho" and threads==1 and "ortho-2-" in str(out) and phase=="sweep": code=86
                if phase=="profile": (dest/"phase-times.jsonl").write_text("{}\n")
            return {"phase":phase,"actual_child_exit":code,"monitor_exit":0,"wrapper_exit":0,"peak_rss_bytes":1000000,"active_ns":13,"timed_out":False,"command":command,"stdout":"fixture.stdout","stderr":"fixture.stderr"}
        with tempfile.TemporaryDirectory() as tmp:
            campaign=M.Campaign(Path(tmp),Path(tmp)/"fake-binary",stub)
            campaign.run()
            for threads in M.THREADS:
                calls_t=[c for c in calls if c[1]==threads]
                self.assertEqual(calls_t[0][:3],("special",threads,"prepare"))
                self.assertEqual(calls_t[1][:3],("special",threads,"verify"))
                starts=[(c[0],int(Path(c[3]["--out"]).parent.name.split("-")[1])) for c in calls_t if c[2]=="sweep"]
                self.assertEqual(starts,M.FAMILIES)
                special_sweep=next(c for c in calls_t if c[0]=="special" and c[2]=="sweep")
                self.assertEqual(int(special_sweep[3]["--elapsed-ns"]),13)
                self.assertEqual(sum(c[0]=="special" and c[2]=="verify" for c in calls_t),1)
            bad=[c[2] for c in calls if c[0]=="ortho" and c[1]==1 and "ortho-2-" in str(c[3].get("--out",""))]
            self.assertNotIn("profile",bad)
            self.assertTrue(any(c[0]=="special" and c[2]=="profile" for c in calls))
            self.assertTrue(any(c[2]=="confirm" and c[1]==4 for c in calls))
            self.assertTrue((Path(tmp)/"ortho-2-threads-1"/"sweep"/"process-failure.json").exists())

if __name__ == "__main__": unittest.main()
