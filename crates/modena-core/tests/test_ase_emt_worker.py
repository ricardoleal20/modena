import json
import pathlib
import subprocess
import sys
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

    def test_malformed_requests_fail_without_output(self):
        valid={"protocol":"modena-ase-emt-v1","action":"evaluate","atomic_numbers":[29],"positions_angstrom":[[0,0,0]],"cell_angstrom":[[4,0,0],[0,4,0],[0,0,4]],"pbc":[True]*3,"asap_cutoff":False}
        for key,value in [("atomic_numbers",[1]),("positions_angstrom",[[0,0]]),("positions_angstrom",[[float('nan'),0,0]]),("pbc",[True,False,True]),("cell_angstrom",[[0,0,0]]*3),("protocol","other"),("extra",1),("asap_cutoff","false")]:
            request=dict(valid);request[key]=value
            process=self.call(request)
            self.assertNotEqual(process.returncode,0,(key,process.stdout))
            self.assertEqual(process.stdout,"")

if __name__ == "__main__": unittest.main()
