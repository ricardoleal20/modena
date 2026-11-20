"""Synthetic regressions; no native or ASE calculator is called here."""
import contextlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import b07_pressure_reference as ref
import campaign
import evidence
from b07_reference import epoch
import test_b07_pressure_reference as controls


class BehavioralControls(unittest.TestCase):
    def setUp(self):
        self.fixture = controls.ReviewedEvidence("runTest")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.out, self.x = self.fixture.out, self.fixture.x

    def test_cached_pass_cannot_survive_lost_raw_data(self):
        self.fixture.seal()
        runner=ref.Campaign(self.out,self.x,{"fixtures":[]})
        self.assertEqual(runner.family("primary",{},lambda:[campaign.gate("one",[0.],[0.],[0.])]),"Passed")
        (self.fixture.folder/"native.log").unlink()
        self.assertEqual(ref.collect(self.out,self.x)[0]["status"],"Inconclusive")

    def test_replaced_derived_bytes_and_sidecar_fail(self):
        self.fixture.seal()
        ref.derived(self.out,self.x,self.fixture.fix,"probe")
        src=self.out/"reference/probe.jsonl"
        side=Path(str(src)+".origin.json")
        binding=campaign.load(side)
        src.write_bytes(src.read_bytes().replace(b"-1.5",b"-1.6"))
        binding["sha256"]=campaign.digest(src)
        side.write_text(json.dumps(binding))
        with self.assertRaisesRegex(ValueError,"prescribed native projection"):
            ref.derived(self.out,self.x,self.fixture.fix,"probe")

    def test_projection_checks_cumulative_artifact_budget_during_append(self):
        self.fixture.seal()
        # Include another retained epoch and let one frame fit, but not the second.
        previous=tempfile.TemporaryDirectory()
        self.addCleanup(previous.cleanup)
        prior=Path(previous.name)
        (prior/"raw").write_bytes(b"x"*500)
        ledger=prior/"ledger"
        ledger.mkdir()
        campaign.dump(ledger/"epoch-prior.json",{"out":str(prior)})
        accounting=patch.object(campaign,"LEDGER",ledger)
        accounting.start()
        self.addCleanup(accounting.stop)
        initial=campaign.raw_bytes(self.out)
        size=len(next(evidence.projected_frames(self.out,self.x,self.fixture.fix,"probe")))
        cap=initial+4096+size+1
        with patch.dict(campaign.LIMITS,artifact_bytes=cap),self.assertRaisesRegex(campaign.Budget,"during"):
            ref.derived(self.out,self.x,self.fixture.fix,"probe")
        src=self.out/"reference/probe.jsonl"
        self.assertFalse(src.exists())
        self.assertEqual(Path(str(src)+".partial").stat().st_size,size)
        receipt=campaign.load(str(src)+".derivation.json")
        self.assertEqual((receipt["status"],receipt["confirmed_frames"]),("BudgetExceeded",1))
        self.assertLessEqual(campaign.raw_bytes(self.out),cap)
        with self.assertRaises(ValueError):
            evidence.validate_projection(self.out,self.x,"probe",self.fixture.fix)

    def test_native_environment_discards_ambient_loader(self):
        captured={}
        def spawn(*args,**kwargs):
            captured.update(kwargs["env"])
            return subprocess.CompletedProcess(args[0],0,b"synthetic",b"")
        ambient=dict(LAMMPS_PLUGIN_PATH="/unfrozen",DYLD_LIBRARY_PATH="/unfrozen",PYTHONPATH="/unfrozen")
        with patch.dict(os.environ,ambient),patch("subprocess.run",side_effect=spawn):
            epoch.help_output(Path("/synthetic"),self.out)
        self.assertFalse(set(ambient)&set(captured))
        self.assertEqual(captured,campaign.ENV)
        with self.assertRaises(ValueError):
            campaign.numerical_env(dict(campaign.ENV,LAMMPS_PLUGIN_PATH="/unfrozen"))

    def test_attempted_infrastructure_error_is_not_missing(self):
        (self.fixture.f/"child.json").write_text(json.dumps(dict(exit=1,budget_stop=None,interrupted=None)))
        self.fixture.seal()
        runner=ref.Campaign(self.out,self.x,{"fixtures":[]})
        self.assertEqual(runner.family("primary",{},lambda:self.fail()),"Inconclusive")
        self.assertEqual(evidence.job_state(self.out,self.x,"probe")[0],"Inconclusive")
        self.assertEqual(ref.collect(self.out,self.x)[0]["status"],"Inconclusive")

    def test_unobserved_worker_prevents_passing_receipt(self):
        folder=self.out/"jobs/ase-six-strain"
        folder.mkdir()
        child=dict(exit=0,budget_stop=None,interrupted=None,rss_samples=9,
                   ase_workers={},ase_worker_interpreter="/provider")
        campaign.dump(folder/"child.json",child)
        self.assertEqual(ref.retained_result(self.out,"ase-six-strain"),"Inconclusive")

    def test_same_count_origin_cannot_select_another_successful_producer(self):
        # Both producers have the same valid transcript/count and successful receipts.
        import shutil
        other=dict(self.fixture.s,tag="other",producer="other")
        otherdir=self.out/"reference/other"
        otherdir.mkdir()
        (otherdir/"native.log").write_text((self.fixture.folder/"native.log").read_text().replace("probe","other"))
        shutil.copyfile(self.fixture.folder/"probe.dump",otherdir/"other.dump")
        cfg=dict(self.fixture.cfg,cwd=str(otherdir),artifacts=[str(otherdir/"native.log"),str(otherdir/"other.dump")])
        self.x["sections"]["other"]=other
        self.x["jobs"]["other"]=cfg
        src=self.out/"reference/probe.jsonl"
        action=dict(role="rust",argv=["/usr/bin/true","--child","matched-evaluate",str(self.out),"rotor","zero",str(src),str(self.out/"matched.jsonl")],
                    expected_frames=2,native_origin=evidence.origin_contract(self.out,self.x,"probe"))
        self.x["jobs"]["evaluate-probe-zero"]=action
        (self.out/"external-epoch.json").write_text(json.dumps(self.x))
        self.fixture.seal()
        folder=self.out/"jobs/other"
        shutil.copytree(self.fixture.f,folder,ignore=shutil.ignore_patterns("artifacts.json"))
        (folder/"launch.json").write_text(json.dumps(dict(cfg,name="other",protocol_sha256=campaign.digest(self.out/"protocol.json"),binary_sha256=campaign.digest("/usr/bin/true"))))
        evidence.seal_job(self.out,self.x,"other")
        ref.derived(self.out,self.x,self.fixture.fix,"probe")
        self.assertTrue(evidence.validate_matched_input(self.out,self.x,"evaluate-probe-zero"))
        binding=evidence.origin_binding(self.out,self.x,"other")
        binding.update(input_path=str(src),sha256=campaign.digest(src))
        Path(str(src)+".origin.json").write_text(json.dumps(binding))
        with self.assertRaisesRegex(ValueError,"prescribed native projection"):
            evidence.validate_matched_input(self.out,self.x,"evaluate-probe-zero")

    def test_worker_monitor_live_unreadable_exited_and_measured(self):
        rows={}
        with self.assertRaisesRegex(ValueError,"RSS unreadable"):
            campaign.sample_ase_workers(rows,[101],lambda p:None,lambda p:"/provider",lambda p:True,"/provider")
        self.assertEqual(campaign.sample_ase_workers(rows,[101],lambda p:None,lambda p:None,lambda p:False,"/provider"),0)
        self.assertTrue(rows["101"]["exited"])
        self.assertFalse(campaign.worker_coverage(dict(ase_workers=rows,ase_worker_interpreter="/provider")))
        rows={}
        self.assertEqual(campaign.sample_ase_workers(rows,[102],lambda p:32768,lambda p:"/provider",lambda p:True,"/provider"),32768)
        self.assertTrue(campaign.worker_coverage(dict(ase_workers=rows,ase_worker_interpreter="/provider")))
        self.assertEqual(rows["102"]["interpreter_samples"],1)
        with self.assertRaisesRegex(ValueError,"executable unreadable"):
            campaign.sample_ase_workers({},[103],lambda p:32768,lambda p:None,lambda p:True,"/provider")

    def test_failed_and_unknown_producers_in_shared_report_and_blocked_family(self):
        # Use hand-written receipts, with metadata/build admission stubbed only for this control.
        self.x["expected_gate_families"]["downstream"]=["blocked"]
        self.x["family_inputs"]["downstream"]={"producers":[],"families":["primary"]}
        (self.out/"external-epoch.json").write_text(json.dumps(self.x))
        campaign.dump(self.out/"inputs/settings.json",{"prior_projected_draws":[]})
        build=dict(source_head="synthetic",compiler="synthetic",build_profile="synthetic",features=[],external_source_revision="synthetic",external_packages=[])
        for code,diagnostic,status in [(1,"ERROR: Lost atoms: original 2 current 1", "Failed"),
                                       (1,"ERROR: file system unavailable", "Inconclusive"),
                                       (-9,"", "BudgetExceeded"), (None,"", "Pending")]:
            with self.subTest(status=status):
                if code is None:
                    (self.fixture.f/"launch.json").unlink()
                    (self.fixture.f/"child.json").unlink()
                else:
                    child=dict(exit=code,budget_stop="wall cap" if code==-9 else None,interrupted=None,
                               started_utc="synthetic",ended_utc="synthetic",elapsed_seconds=0.1,
                               peak_rss_bytes=4096,user_cpu_seconds=0.,system_cpu_seconds=0.)
                    (self.fixture.f/"child.json").write_text(json.dumps(child))
                    (self.fixture.f/"stdout.log").write_text(diagnostic)
                    self.fixture.seal()
                runner=ref.Campaign(self.out,self.x,{"fixtures":[]})
                self.assertEqual(runner.family("primary",{},lambda:self.fail()),status)
                self.assertEqual(runner.family("downstream",{},lambda:self.fail()),"Pending")
                self.assertEqual(campaign.load(runner.results/"downstream.json")["details"]["blocked_by"]["primary"],status)
                if status == "Failed":
                    self.assertEqual(ref.collect(self.out,self.x)[0]["observed"],[])
                campaign.dump(self.out/"external-gates.json",ref.collect(self.out,self.x))
                if status == "Failed":
                    # Simulate interruption after producer receipts, before any family publication.
                    for interrupted in runner.results.iterdir(): interrupted.unlink()
                empty=dict(gas_molecule_steps=0,conditional_draws=0,all_molecule_steps=0,unknown_suffixes=[])
                with patch.object(campaign,"verify",return_value=({},self.x)),patch.object(campaign,"jobs",return_value=self.x["jobs"]),patch.object(campaign,"executable_build",return_value=build),patch.object(campaign,"measured_work",return_value=(0,0,False)),patch.object(campaign,"consumption",return_value=empty),patch.object(campaign,"external_spent",return_value=(0,0)),patch.object(campaign,"epochs",return_value=[self.out]):
                    campaign.report(self.out)
                report=max((self.out/"reports").iterdir())
                gates={g["name"]:g for g in campaign.load(report/"gate-ledger.json")["gates"]}
                self.assertEqual(gates["execution_probe"]["status"],status)
                self.assertEqual(gates["one"]["status"],status)
                self.assertEqual(gates["blocked"]["status"],"Pending")
                (self.out/"external-gates.json").unlink()
                (self.fixture.f/"artifacts.json").unlink(missing_ok=True)
                for path in runner.results.iterdir(): path.unlink()

    def test_direct_admission_rejects_replaced_prerequisite_receipt(self):
        self.fixture.seal()
        runner=ref.Campaign(self.out,self.x,{"fixtures":[]})
        self.assertEqual(runner.family("primary",{},lambda:[campaign.gate("one",[0.],[0.],[0.])]),"Passed")
        (self.fixture.f/"wrapper.json").write_text('{"replaced":true}')
        roster={"next":dict(role="rust",argv=["/synthetic"],requires_families=["primary"])}
        with patch.object(campaign,"exclusive",return_value=contextlib.nullcontext()),patch.object(campaign,"verify",return_value=({},self.x)),patch.object(campaign,"jobs",return_value=roster),patch.object(campaign,"run_wrapper",side_effect=AssertionError("must not launch")):
            with self.assertRaisesRegex(ValueError,"blocked by external prerequisite"):
                campaign.run_job(self.out,"next")
        self.assertFalse((self.out/"jobs/next").exists())

    def test_final_report_rejects_unbound_and_forged_empty_failed_rows(self):
        self.fixture.seal()
        runner=ref.Campaign(self.out,self.x,{"fixtures":[]})
        self.assertEqual(runner.family("primary",{},lambda:[campaign.gate("one",[0.],[0.],[0.])]),"Passed")
        path=runner.results/"primary.json"
        row=campaign.load(path)
        row["gates"]=[campaign.unmeasured("one","Failed")]
        row["details"]["producer_diagnostics"]=row["binding"]["dependencies"]
        path.write_text(json.dumps(row))
        campaign.dump(self.out/"external-gates.json",row["gates"])
        self.assertEqual(ref.collect(self.out,self.x)[0]["status"],"Inconclusive")
        with self.assertRaises(ValueError): campaign.checked_gates(row["gates"])
        campaign.dump(self.out/"inputs/settings.json",{"prior_projected_draws":[]})
        build=dict(source_head="synthetic",compiler="synthetic",build_profile="synthetic",features=[],external_source_revision="synthetic",external_packages=[])
        empty=dict(gas_molecule_steps=0,conditional_draws=0,all_molecule_steps=0,unknown_suffixes=[])
        child=campaign.load(self.fixture.f/"child.json")
        child.update(started_utc="synthetic",ended_utc="synthetic",elapsed_seconds=0.1,peak_rss_bytes=4096,user_cpu_seconds=0.,system_cpu_seconds=0.)
        (self.fixture.f/"child.json").write_text(json.dumps(child))
        # Recreate a consistent successful producer/family binding before the forged status.
        (self.fixture.f/"artifacts.json").unlink()
        self.fixture.seal()
        row["binding"]=evidence.family_binding(self.out,self.x,"primary")
        row["details"]["producer_diagnostics"]=row["binding"]["dependencies"]
        path.write_text(json.dumps(row))
        def shared_report():
            with patch.object(campaign,"verify",return_value=({},self.x)),patch.object(campaign,"jobs",return_value=self.x["jobs"]),patch.object(campaign,"executable_build",return_value=build),patch.object(campaign,"measured_work",return_value=(0,0,False)),patch.object(campaign,"consumption",return_value=empty),patch.object(campaign,"external_spent",return_value=(0,0)),patch.object(campaign,"epochs",return_value=[self.out]):
                campaign.report(self.out)
            return campaign.load(max((self.out/"reports").iterdir())/"gate-ledger.json")
        report=shared_report()
        gates,defects=report["gates"],report["incomplete_gate_files"]
        self.assertEqual(next(g["status"] for g in gates if g["name"]=="one"),"Inconclusive")
        self.assertIn(str(self.out/"external-gates.json"),defects)
        # A free-standing published empty Failed record cannot create a physical failure.
        path.unlink()
        gates=shared_report()["gates"]
        self.assertEqual(next(g["status"] for g in gates if g["name"]=="one"),"Pending")

    def test_failed_native_dependency_does_not_claim_unlaunched_matched_attempt(self):
        self.x["jobs"]["matched"] = dict(role="rust",argv=["/usr/bin/true"],
            cwd=str(self.out),env=campaign.ENV,artifacts=[str(self.out/"matched.jsonl")])
        self.x["family_inputs"]["primary"] = dict(producers=["probe","matched"],
            owners=["matched"],families=[])
        (self.out/"external-epoch.json").write_text(json.dumps(self.x))
        (self.fixture.f/"child.json").write_text(json.dumps(dict(exit=1,budget_stop=None,interrupted=None)))
        (self.fixture.f/"stdout.log").write_text("ERROR: Lost atoms: original 2 current 1")
        self.fixture.seal()
        self.assertEqual(ref.collect(self.out,self.x)[0]["status"],"Pending")
        family=campaign.load(self.out/"reference/comparisons/primary.json")
        self.assertEqual(family["details"]["blocked_by"]["probe"],"Failed")
        self.assertEqual(family["binding"]["dependencies"]["matched"]["status"],"Pending")
        self.assertFalse((self.out/"jobs/matched").exists())


class ControlledImports(unittest.TestCase):
    def test_actual_host_and_provider_startup_ignore_readable_timestamp_cache(self):
        interpreters = {Path(sys.executable).resolve()}
        optional = epoch.ASEROOT / "bin/python"
        if optional.is_file():
            interpreters.add(optional.resolve())
        launch=campaign.ROOT/"bench/b07_validation/launch.sh"
        host=campaign.ROOT/"bench/b07_pressure_reference.py"
        # Synthetic cache code is never scientific work. Preserve any previous cache exactly.
        poison='raise RuntimeError("SYNTHETIC_TIMESTAMP_CACHE_EXECUTED")'
        for python in interpreters:
            with self.subTest(python=str(python)), tempfile.TemporaryDirectory() as td:
                static=Path(td)
                runtime=epoch.provider_wrapper(static,python,static)
                module=static/"synthetic_provider.py"
                module.write_text('VALUE="source"\n')
                for source,entry in [(campaign.ROOT/"bench/b07_validation/freeze.py","host"),(module,"provider")]:
                    code="import importlib.util,sys;print(importlib.util.cache_from_source(sys.argv[1]))"
                    cache=Path(subprocess.check_output([str(python),*campaign.PYTHON_FLAGS,"-c",code,str(source)],env=campaign.numerical_env(campaign.ENV),text=True).strip())
                    # The tested interpreter's cache must use its normal cache location.
                    cache=source.parent/"__pycache__"/cache.name
                    existed=cache.exists()
                    previous=cache.read_bytes() if existed else None
                    cache.parent.mkdir(exist_ok=True)
                    build="import importlib.util,marshal,struct,sys,os; p,c,code=sys.argv[1:]; s=os.stat(p); open(c,'wb').write(importlib.util.MAGIC_NUMBER+struct.pack('<III',0,int(s.st_mtime),s.st_size)+marshal.dumps(compile(code,p,'exec')))"
                    try:
                        subprocess.run([str(python),*campaign.PYTHON_FLAGS,"-c",build,str(source),str(cache),poison],env=campaign.numerical_env(campaign.ENV),check=True)
                        unsafe="import sys;sys.path.insert(0,sys.argv[1]);__import__(sys.argv[2])"
                        result=subprocess.run([str(python),"-I","-S","-B","-c",unsafe,str(source.parent),source.stem],env=campaign.numerical_env(campaign.ENV),capture_output=True,text=True)
                        self.assertNotEqual(result.returncode,0)
                        self.assertIn("SYNTHETIC_TIMESTAMP_CACHE_EXECUTED",result.stderr)
                        cmd=(["/bin/sh",str(launch),str(python),str(host),"--help"] if entry=="host" else
                             [runtime["wrapper"],"-u","-c","import synthetic_provider; print(synthetic_provider.VALUE)","--persistent"])
                        result=subprocess.run(cmd,env=dict(campaign.ENV,PYTHONPATH="/unfrozen",LAMMPS_PLUGIN_PATH="/unfrozen"),capture_output=True,text=True)
                        self.assertEqual(result.returncode,0,result.stderr)
                        self.assertIn("usage:" if entry=="host" else "source",result.stdout)
                    finally:
                        if existed: cache.write_bytes(previous)
                        else: cache.unlink(missing_ok=True)


class FrozenGraph(unittest.TestCase):
    @unittest.skipUnless((epoch.V4 / "external-epoch.json").is_file(), "requires the frozen local v4 archive")
    def test_registered_graph_preserves_roster_and_exact_origins(self):
        import shutil
        archived=epoch.V4
        old=campaign.load(archived/"external-epoch.json")
        with tempfile.TemporaryDirectory() as td:
            out=Path(td).resolve()
            shutil.copytree(archived/"inputs",out/"inputs")
            shutil.copyfile(archived/"protocol.json",out/"protocol.json")
            static=out/"reference/inputs"
            static.mkdir(parents=True)
            runtime=epoch.provider_wrapper(static)
            runtime.update(process_executable="/usr/bin/true",process_executable_sha256=campaign.digest("/usr/bin/true"))
            epoch.ase_input(static)
            old_ase=campaign.load(archived/"reference/inputs/ase.json")
            self.assertEqual(campaign.load(static/"ase.json")["snapshots"],old_ase["snapshots"])
            fixtures={f["id"]:f for f in campaign.load(out/"protocol.json")["fixtures"]}
            source=campaign.inventory([Path("/usr/bin/true"),Path(runtime["wrapper"])])
            x=epoch.register(out,Path("/usr/bin/true"),fixtures,source,["synthetic"],
                dict(compiler="synthetic",build_profile="synthetic",features=[]),runtime)
            self.assertEqual(x["expected_gate_families"],old["expected_gate_families"])
            self.assertEqual(set(x["jobs"]),set(old["jobs"]))
            self.assertEqual(sum(j["role"]=="reference" for j in x["jobs"].values()),24)
            self.assertEqual(len(x["evaluations"]),68)
            for name, action in x["jobs"].items():
                self.assertTrue(action["artifacts"])
                if name in x["evaluations"]:
                    origin=evidence.origin_contract(out,x,x["evaluations"][name]["tag"])
                    self.assertEqual(action["native_origin"],origin)
                    self.assertIn(origin["producer_job"],x["family_inputs"][name]["producers"])
                    self.assertEqual(campaign.load(origin["input_path"]+".derivation.json")["status"],"Pending")
            for name,binding in x["family_inputs"].items():
                self.assertTrue(binding["producers"])
                self.assertFalse(name in binding["families"])
                self.assertTrue(set(binding["producers"])<=set(x["jobs"]))
                self.assertTrue(set(binding["families"])<=set(x["expected_gate_families"]))
            self.assertFalse((out/"clock.json").exists())
            self.assertFalse((out/"jobs").exists())
