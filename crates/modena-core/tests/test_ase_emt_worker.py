import json
import pathlib
import subprocess
import sys
import struct
import tempfile
import unittest
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT

WORKER = pathlib.Path(__file__).parents[1] / "src/external/ase_emt_worker.py"

class WorkerTests(unittest.TestCase):
    def call(self, request):
        return subprocess.run([sys.executable, str(WORKER)], input=json.dumps(request), text=True, capture_output=True)

    def test_displaced_copper_matches_independent_ase(self):
        for cell in [[[3.6,0,0],[0,3.6,0],[0,0,3.6]], [[3.6,0,0],[0.3,3.6,0],[0.2,0.1,3.6]]]:
            positions=[[0.03,0,0],[0,1.8,1.8],[1.8,0,1.8],[1.8,1.8,0]]
            request={"protocol":"modena-ase-emt-v1","action":"evaluate","atomic_numbers":[29]*4,"positions_angstrom":positions,"cell_angstrom":cell,"pbc":[True]*3,"asap_cutoff":False}
            process=self.call(request)
            self.assertEqual(process.returncode,0,process.stderr)
            response=json.loads(process.stdout)
            atoms=Atoms(numbers=[29]*4,positions=positions,cell=cell,pbc=True)
            atoms.calc=EMT(asap_cutoff=False)
            actual=response["evaluation"]
            self.assertAlmostEqual(actual["energy_ev"],atoms.get_potential_energy(),delta=1e-10)
            np.testing.assert_allclose(actual["forces_ev_per_angstrom"],atoms.get_forces(),atol=1e-10,rtol=0)
            np.testing.assert_allclose(actual["stress_ev_per_angstrom3"],atoms.get_stress(voigt=False),atol=1e-12,rtol=0)

    def test_discovery_has_model_identity(self):
        p=self.call({"protocol":"modena-ase-emt-v1","action":"info"})
        self.assertEqual(p.returncode,0,p.stderr)
        identity=json.loads(p.stdout)["identity"]
        self.assertEqual(identity["provider"],"ase_emt")
        self.assertEqual(len(identity["model_sha256"]),64)
        self.assertEqual(identity["composition"],"total")

    def test_persistent_uses_one_calculator_and_refreshes_positions_cell_numbers(self):
        snapshots=[]
        for numbers, x, cell in [([29]*4,0.03,3.6),([29]*4,0.17,3.6),([29]*4,0.17+1e-16,3.6),([29]*4,0.17,3.7),([28,29,29,29],0.17,3.7)]:
            snapshots.append({"protocol":"modena-ase-emt-v1","action":"evaluate","atomic_numbers":numbers,"positions_angstrom":[[x,0,0],[0,1.8,1.8],[1.8,0,1.8],[1.8,1.8,0]],"cell_angstrom":[[cell,0,0],[0.2,cell,0],[0.1,0.3,cell]],"pbc":[True]*3,"asap_cutoff":False})
        requests=[{"protocol":"modena-ase-emt-v1","action":"info"}]+snapshots
        data=b""
        for i, request in enumerate(requests,1):
            frame=json.dumps(dict(session_protocol="modena-ase-emt-session-v1",request_id=i,nonce=[i]*32,request=request)).encode()
            data+=struct.pack(">I",len(frame))+frame
        with tempfile.TemporaryDirectory() as td:
            counter=pathlib.Path(td)/"calculators"
            evaluations=pathlib.Path(td)/"evaluations"
            bootstrap=f'''import runpy,sys
from pathlib import Path
from ase.calculators import emt
original=emt.EMT
class CountedEMT(original):
 def __init__(self,*a,**kw):
  with Path({str(counter)!r}).open('a') as f:f.write('start\\n')
  super().__init__(*a,**kw)
 def calculate(self,*a,**kw):
  with Path({str(evaluations)!r}).open('a') as f:f.write('evaluate\\n')
  return super().calculate(*a,**kw)
emt.EMT=CountedEMT
sys.argv=[{str(WORKER)!r},"--persistent"]
runpy.run_path({str(WORKER)!r},run_name="__main__")
'''
            p=subprocess.run([sys.executable,"-c",bootstrap],input=data,capture_output=True,timeout=20)
            self.assertEqual(p.returncode,0,p.stderr)
            self.assertEqual(counter.read_text(),"start\n")
            self.assertEqual(evaluations.read_text().splitlines(),["evaluate"]*len(snapshots),"every changed snapshot, even below ASE's default cache tolerance, must recalculate")
        remaining=p.stdout
        replies=[]
        while remaining:
            size=struct.unpack(">I",remaining[:4])[0]
            replies.append(json.loads(remaining[4:4+size]))
            remaining=remaining[4+size:]
        self.assertEqual(len(replies),len(requests))
        for i, reply in enumerate(replies,1):
            self.assertEqual(reply["request_id"],i)
            self.assertEqual(reply["nonce"],[i]*32)
            self.assertEqual(reply["session_protocol"],"modena-ase-emt-session-v1")
        for request,reply in zip(snapshots,replies[1:]):
            atoms=Atoms(numbers=request["atomic_numbers"],positions=request["positions_angstrom"],cell=request["cell_angstrom"],pbc=True)
            atoms.calc=EMT(asap_cutoff=False)
            got=reply["response"]["evaluation"]
            self.assertAlmostEqual(got["energy_ev"],atoms.get_potential_energy(),delta=1e-10)
            np.testing.assert_allclose(got["forces_ev_per_angstrom"],atoms.get_forces(),atol=1e-10,rtol=0)
            np.testing.assert_allclose(got["stress_ev_per_angstrom3"],atoms.get_stress(voigt=False),atol=1e-12,rtol=0)

    def test_malformed_requests_fail_without_output(self):
        valid={"protocol":"modena-ase-emt-v1","action":"evaluate","atomic_numbers":[29],"positions_angstrom":[[0,0,0]],"cell_angstrom":[[4,0,0],[0,4,0],[0,0,4]],"pbc":[True]*3,"asap_cutoff":False}
        for key,value in [("atomic_numbers",[1]),("positions_angstrom",[[0,0]]),("positions_angstrom",[[float('nan'),0,0]]),("pbc",[True,False,True]),("cell_angstrom",[[0,0,0]]*3),("protocol","other"),("extra",1),("asap_cutoff","false")]:
            request=dict(valid);request[key]=value
            process=self.call(request)
            self.assertNotEqual(process.returncode,0,(key,process.stdout))
            self.assertEqual(process.stdout,"")

    def test_session_rejects_invalid_nonce_and_frame_before_output(self):
        valid=dict(session_protocol="modena-ase-emt-session-v1",request_id=1,
            nonce=[0]*32,request={"protocol":"modena-ase-emt-v1","action":"info"})
        frames=[]
        for nonce in [None,[],[0]*31,[0]*33,[True]*32,[-1]*32,[256]*32]:
            frame=dict(valid);frame["nonce"]=nonce;frames.append(frame)
        missing=dict(valid);del missing["nonce"];frames.append(missing)
        wrong=dict(valid);wrong["request_id"]=2;frames.append(wrong)
        for frame in frames:
            data=json.dumps(frame).encode()
            p=subprocess.run([sys.executable,str(WORKER),"--persistent"],input=struct.pack(">I",len(data))+data,capture_output=True,timeout=20)
            self.assertNotEqual(p.returncode,0,frame)
            self.assertEqual(p.stdout,b"")

if __name__ == "__main__": unittest.main()
