"""No scientific campaign: launcher stubs and tiny recorded monitor controls."""
import importlib.util
import json
from pathlib import Path
import tempfile
import os
import subprocess
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
                previous=Path(flags["--protocol"]).parent/"verification.json" if phase=="sweep" else Path(flags["--campaign"])
                M.write_new(dest/"summary.json",{"usage":M.empty_usage(),"prior":M.artifact(previous)})
                if phase=="sweep": code=4
                if family=="ortho" and threads==1 and "ortho-2-" in str(out) and phase=="sweep": code=86
                if phase=="profile": (dest/"phase-times.jsonl").write_text("{}\n")
            return {"phase":phase,"actual_child_exit":code,"monitor_exit":0,"wrapper_exit":0,"peak_rss_bytes":1000000,"active_ns":13,"timed_out":False,"command":command,"stdout":"fixture.stdout","stderr":"fixture.stderr"}
        with tempfile.TemporaryDirectory() as tmp:
            binary=Path(tmp)/"fake-binary"; binary.write_bytes(b"registered executable")
            campaign=M.Campaign(Path(tmp),M.artifact(binary),stub)
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


class ReviewRegressionTests(unittest.TestCase):
    def drive_stub(self, tmp, hook):
        calls=[]
        binary=Path(tmp)/"fake-binary"; binary.write_bytes(b"registered executable")
        def stub(command, env, out, phase, remaining):
            out.mkdir(); flags=dict(zip(command[2::2],command[3::2])); key=out.parent.name
            calls.append((key,phase,flags,remaining))
            row={"phase":phase,"actual_child_exit":0,"monitor_exit":0,"wrapper_exit":0,"peak_rss_bytes":1000000,"active_ns":13,"timed_out":False,"command":command,"stdout":"fixture.stdout","stderr":"fixture.stderr"}
            if phase=="prepare":
                prepared=Path(flags["--out"]); prepared.mkdir(); M.write_new(prepared/"manifest.json",{})
                M.write_new(prepared/"protocol.json",{});M.write_new(prepared/"snapshots.json",[])
            elif phase=="verify":
                prepared=Path(flags["--protocol"]).parent
                M.write_new(prepared/"verification.json",{"phase":"verify","usage":M.empty_usage(),"prior":None,"inputs":{},"artifacts":{},"failure":None,"finalization_failure":None})
            elif phase in ("sweep","confirm","profile"):
                dest=Path(flags["--out"]);dest.mkdir()
                previous=Path(flags["--protocol"]).parent/"verification.json" if phase=="sweep" else Path(flags["--campaign"])
                M.write_new(dest/"summary.json",{"phase":phase,"usage":M.empty_usage(),"prior":M.artifact(previous),"inputs":{},"artifacts":{},"failure":None,"finalization_failure":None})
                if phase=="profile":(dest/"phase-times.jsonl").write_text("{}\n")
            hook(key,phase,flags,row,remaining)
            return row
        return M.Campaign(Path(tmp),M.artifact(binary),stub),calls,binary

    def test_exhausted_timeout_reports_and_continues(self):
        with tempfile.TemporaryDirectory() as tmp:
            def timeout(key,phase,flags,row,remaining):
                if key=="cubic-2-threads-1" and phase=="sweep":
                    row.update(actual_child_exit=-9,active_ns=remaining+1,timed_out=True)
            campaign,calls,_=self.drive_stub(tmp,timeout)
            campaign.run()
            report=next(c for c in calls if c[:2]==("cubic-2-threads-1","report"))
            self.assertGreater(report[3],0)
            self.assertTrue(any(c[:2]==("ortho-2-threads-1","sweep") for c in calls))
            self.assertTrue(any(c[:2]==("ortho-3-threads-4","report") for c in calls))
            self.assertEqual(report[3], M.REPORT_NS)
            self.assertGreater(campaign.active[("cubic",2,1)], M.WALL_NS)
            self.assertEqual(sum(c[1] == "report" for c in calls),10)

    def test_exhaustion_before_next_admission_reports_without_fake_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            def exhausted(key,phase,flags,row,remaining):
                if key=="cubic-2-threads-1" and phase=="verify": row["active_ns"]=remaining
            campaign,calls,_=self.drive_stub(tmp,exhausted)
            campaign.run()
            self.assertFalse(any(c[:2]==("cubic-2-threads-1","sweep") for c in calls))
            self.assertTrue(any(c[:2]==("cubic-2-threads-1","report") for c in calls))
            self.assertTrue(any(c[:2]==("ortho-2-threads-1","sweep") for c in calls))
            memory=json.loads((Path(tmp)/"cubic-2-threads-1/memory.json").read_text())
            self.assertEqual([row["phase"] for row in memory["children"]],["prepare","verify"])
            process=json.loads((Path(tmp)/"cubic-2-threads-1/sweep/process-error.json").read_text())
            self.assertIsNone(process["observation"]["actual_child_exit"])

    def test_replaced_produced_proof_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            def replace(key,phase,flags,row,remaining):
                if key=="cubic-2-threads-1" and phase=="prepare":
                    proof=Path(tmp)/"special-2-threads-1/prepared/verification.json"
                    value=json.loads(proof.read_text());value["usage"]["candidate"]["accuracy_calls"]=999
                    proof.write_text(json.dumps(value))
            campaign,calls,_=self.drive_stub(tmp,replace)
            with self.assertRaisesRegex(ValueError,"artifact.*changed"):
                campaign.run()
            self.assertFalse(any(c[:2]==("cubic-2-threads-1","verify") for c in calls))
            saved=json.loads((Path(tmp)/"special-2-threads-1/process-verify/accepted-receipt.json").read_text())
            self.assertNotEqual(saved["sha256"],M.digest(saved["path"]))

    def test_replaced_registered_binary_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign,calls,binary=self.drive_stub(tmp,lambda *args:None)
            campaign.call(("special",2,1),"report",[])
            binary.write_bytes(b"replacement executable")
            with self.assertRaisesRegex(ValueError,"artifact.*changed"):
                campaign.call(("cubic",2,1),"report",[])
            self.assertEqual(len(calls),1)

    def test_replaced_finalized_receipt_is_rejected_before_wrapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);prepared=root/"prepared";prepared.mkdir();out=root/"failed";out.mkdir()
            M.write_new(prepared/"manifest.json",{})
            summary=out/"summary.json";M.write_new(summary,{"usage":M.empty_usage()})
            saved=M.retain_output(summary,out)
            summary.write_text("{}")
            with self.assertRaisesRegex(ValueError,"artifact.*changed"):
                M.failed_receipt(prepared,out,"verify",None,{"actual_child_exit":1,"timed_out":False},1,saved)
            self.assertFalse((out/"process-failure.json").exists())

    @unittest.skipUnless(os.environ.get("MODENA_B05_REPORT_BINARY"),"requires explicit report-only helper executable")
    def test_real_report_torn_journal_with_and_without_finalized_summary(self):
        binary=os.environ["MODENA_B05_REPORT_BINARY"]
        for finalized in (False,True):
            with self.subTest(finalized=finalized), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);prepared=root/"prepared";out=root/"failed";out.mkdir()
                subprocess.run([binary,"prepare","--out",str(prepared),"--family","special","--replicas","2","--threads","1","--reference-extended"],check=True,capture_output=True)
                epoch=M.digest(prepared/"manifest.json")
                event={"sequence":0,"manifest_sha256":epoch,"kind":"note","data":{}}
                events=out/"events.jsonl";events.write_text(json.dumps(event)+"\n"+'{"sequence":1,"data":')
                detail=out/"retained.json";M.write_new(detail,{"retained":True})
                summary=out/"summary.json"
                if finalized:
                    M.write_new(summary,{"schema":"modena-spme-phase-v1","phase":"verify","prepared":str(prepared),"manifest_sha256":epoch,"prior":None,"inputs":{},"artifacts":{"events":M.artifact(events),"retained-detail":M.artifact(detail)},"usage":M.empty_usage(),"failure":{"ReferenceUnconverged":"primary reference cause"},"finalization_failure":{"Io":"original final journal cause"}})
                row={"phase":"verify","actual_child_exit":3,"monitor_exit":0,"wrapper_exit":0,"peak_rss_bytes":1000000,"active_ns":17,"timed_out":False,"command":[],"stdout":"fixture.stdout","stderr":"fixture.stderr"}
                recovered=M.failed_receipt(prepared,out,"verify",None,row,17,M.artifact(summary) if finalized else None)
                memory=root/"memory.json";M.write_new(memory,{"schema":"modena-spme-memory-v1","manifest_sha256":epoch,"children":[row]})
                result=subprocess.run([binary,"report","--campaign",str(recovered),"--expected-sha256",M.digest(recovered),"--memory",str(memory),"--memory-sha256",M.digest(memory),"--out",str(root/"report")],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                data=json.loads(recovered.read_text());self.assertIsNone(data["prior"])
                diagnostics=json.loads((root/"report/diagnostics.json").read_text())
                self.assertEqual(len(diagnostics["phases"]),1)
                self.assertEqual(events.read_text().splitlines()[-1],'{"sequence":1,"data":')
                self.assertEqual(len(Path(data["artifacts"]["events"]["path"]).read_text().splitlines()),1)
                if finalized:
                    self.assertEqual(data["failure"],{"ReferenceUnconverged":"primary reference cause"})
                    self.assertEqual(data["finalization_failure"],{"Io":"original final journal cause"})
                    self.assertIn("retained-detail",data["artifacts"])
                    self.assertIn("raw-summary",data["artifacts"])

if __name__ == "__main__": unittest.main()
