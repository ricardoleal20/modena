"""Optional ASE calculator: native scientific oracle and real CLI/cache boundaries."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
HAS_ASE = importlib.util.find_spec("ase") is not None and importlib.util.find_spec("numpy") is not None
BINARY = os.environ.get("MODENA_BINARY", str(Path(__file__).resolve().parents[3] / "target/debug/modena"))


class BaseImportTests(unittest.TestCase):
    def test_base_import_does_not_load_optional_dependencies(self):
        result = subprocess.run([sys.executable, "-c", "import modena,sys; assert 'ase' not in sys.modules; assert 'numpy' not in sys.modules"],
                                env=dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[2])), capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)


@unittest.skipUnless(HAS_ASE, "optional ASE/NumPy dependencies absent")
class CalculatorTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(importlib.util.find_spec("modena.ase"), "missing optional ModenaCalculator")
        from ase import Atoms
        from modena.ase import ModenaCalculator
        self.atoms = Atoms("Ar2", positions=[[0, 0, 0], [4, 0, 0]], cell=[20]*3, pbc=True)
        self.calc = ModenaCalculator(sigma_angstrom=3.4, epsilon_ev=0.0103,
            cutoff_angstrom=8.5, skin_angstrom=0.3, binary=BINARY)
        self.atoms.calc = self.calc

    def test_pair_energy_forces_and_finite_difference(self):
        u = lambda r: 4*0.0103*((3.4/r)**12-(3.4/r)**6)
        f = lambda r: 24*0.0103/r*(2*(3.4/r)**12-(3.4/r)**6)
        self.assertAlmostEqual(self.atoms.get_potential_energy(), u(4)-u(8.5)+(4-8.5)*f(8.5), places=14)
        force = self.atoms.get_forces()[0, 0]
        self.assertAlmostEqual(force, -f(4)+f(8.5), places=14)
        self.atoms.positions[0, 0] = 1e-5
        plus = self.atoms.get_potential_energy()
        self.atoms.positions[0, 0] = -1e-5
        minus = self.atoms.get_potential_energy()
        self.assertAlmostEqual(force, -(plus-minus)/2e-5, places=9)

    def test_changes_invalidate_and_parameter_changes_reset(self):
        from ase.calculators.calculator import CalculatorSetupError
        e = self.atoms.get_potential_energy()
        cached = self.calc.results
        self.assertEqual(self.atoms.get_potential_energy(), e)
        self.assertIs(cached, self.calc.results)
        self.atoms.positions[1, 0] = 4.5
        self.assertNotEqual(self.atoms.get_potential_energy(), e)
        cached = self.calc.results
        self.atoms.cell = [21]*3
        self.atoms.get_potential_energy()
        self.assertIsNot(cached, self.calc.results)
        cached = self.calc.results
        self.atoms.numbers[:] = 29
        self.atoms.get_potential_energy()
        self.assertIsNot(cached, self.calc.results)
        for kwargs in [dict(epsilon_ev=0.02), dict(sigma_angstrom=3.3), dict(cutoff_angstrom=8.),
                       dict(skin_angstrom=0.2), dict(timeout=10.), dict(binary="/nonexistent/modena")]:
            self.calc.set(**kwargs)
            self.assertEqual(self.calc.results, {})
            if "binary" not in kwargs:
                self.atoms.get_potential_energy()
        with self.assertRaises(Exception): self.atoms.get_potential_energy()
        self.assertEqual(self.calc.results, {})
        self.calc.set(binary=BINARY)
        self.atoms.get_potential_energy()
        self.atoms.numbers[1] = 18
        with self.assertRaises(CalculatorSetupError): self.atoms.get_potential_energy()
        self.assertEqual(self.calc.results, {})

    def test_rejects_unsupported_atoms_even_after_cached_success(self):
        from ase.calculators.calculator import CalculatorSetupError
        from ase.constraints import FixAtoms
        mutations = [lambda a: a.set_constraint(FixAtoms(indices=[0])),
            lambda a: a.set_initial_charges([1., 0.]), lambda a: a.set_initial_magnetic_moments([1., 0.]),
            lambda a: a.set_pbc([True, False, True]), lambda a: a.set_cell([[20,1,0],[0,20,0],[0,0,20]]),
            lambda a: a.set_cell([0,20,20]), lambda a: a.positions.__setitem__((0,0), float("nan")),
            lambda a: a.set_atomic_numbers([0,0]), lambda a: a.set_atomic_numbers([119,119])]
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                a = self.atoms.copy(); a.calc = self.calc
                a.get_potential_energy()
                mutate(a)
                with self.assertRaises(CalculatorSetupError): a.get_potential_energy()
                self.assertEqual(self.calc.results, {})

    def test_stress_is_explicitly_unsupported_and_clears_results(self):
        from ase.calculators.calculator import PropertyNotImplementedError
        self.atoms.get_potential_energy()
        with self.assertRaises(PropertyNotImplementedError): self.atoms.get_stress()
        self.assertEqual(self.calc.results, {})

    def test_transport_and_invalid_output_failures_clear_results(self):
        from ase.calculators.calculator import CalculationFailed
        valid = dict(protocol="modena-native-lj-v1", units=dict(length="angstrom", energy="eV", force="eV/angstrom"),
                     energy_ev=0., forces_ev_per_angstrom=[[0.,0.,0.]]*2)
        scripts = ["raise SystemExit(7)", "import time; time.sleep(2)", "print('garbage')"]
        for field, value in [("protocol", "wrong"), ("units", {}), ("energy_ev", float("nan")),
                             ("energy_ev", "0.0"), ("energy_ev", False),
                             ("forces_ev_per_angstrom", [["0",0,0]]*2),
                             ("forces_ev_per_angstrom", [[False,0,0]]*2),
                             ("forces_ev_per_angstrom", [[0.,0.,0.]]), ("forces_ev_per_angstrom", [[float("inf"),0,0]]*2)]:
            payload = dict(valid); payload[field] = value
            scripts.append("print(" + repr(json.dumps(payload)) + ")")
        with tempfile.TemporaryDirectory() as td:
            for script in scripts:
                self.calc.set(binary=BINARY, timeout=30.)
                self.atoms.get_potential_energy()
                executable = Path(td) / "cli with spaces"
                executable.write_text("#!" + sys.executable + "\nimport sys; sys.stdin.read()\n" + script + "\n")
                executable.chmod(0o755)
                self.calc.set(binary=str(executable), timeout=0.1)
                with self.assertRaises(CalculationFailed): self.atoms.get_potential_energy()
                self.assertEqual(self.calc.results, {})

    def test_invalid_parameters_fail_before_dispatch_and_clear_results(self):
        from ase.calculators.calculator import CalculatorSetupError
        for kwargs in [dict(sigma_angstrom=0), dict(epsilon_ev=-1), dict(cutoff_angstrom=10),
                       dict(skin_angstrom=-1), dict(timeout=0), dict(timeout=float("nan")), dict(unknown=1)]:
            self.atoms.get_potential_energy()
            with self.assertRaises(CalculatorSetupError):
                self.calc.set(**kwargs)
                self.atoms.get_potential_energy()
            self.assertEqual(self.calc.results, {})
            self.calc.set(sigma_angstrom=3.4,epsilon_ev=0.0103,cutoff_angstrom=8.5,skin_angstrom=0.3,timeout=30.)


if __name__ == "__main__":
    unittest.main()
