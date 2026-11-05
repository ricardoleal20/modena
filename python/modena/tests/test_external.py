"""Typed and actual CLI boundary coverage for the optional ASE EMT provider."""
import importlib
import os
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import modena
client = importlib.import_module("modena.run")
BINARY = os.environ.get("MODENA_BINARY", str(Path(__file__).resolve().parents[3] / "target/debug/modena"))

class ExternalConfigTests(unittest.TestCase):
    def test_discovery_failure_is_reported_with_the_cli_diagnostic(self):
        with self.assertRaisesRegex(RuntimeError, "failed to discover ASE EMT provider"):
            modena.provider_info("/nonexistent/modena-python", binary=BINARY, timeout_ms=100)

    def test_required_pins_and_total_composition_render_without_native_defaults(self):
        self.assertTrue(hasattr(modena, "AseEmt"), "missing public typed AseEmt model")
        with self.assertRaises(TypeError):
            modena.AseEmt(ase_version="3.26.0", model_sha256="a" * 64, atomic_numbers=[29]*4)
        model = modena.AseEmt(ase_version="3.26.0", model_sha256="a" * 64,
                              atomic_numbers=[29]*4, composition="total")
        cfg = modena.Config(modena.Fcc(3.61, [1,1,1], mass=63.546), model,
                            modena.VelocityVerlet(0.0001), units="real")
        rendered = tomllib.loads(cfg.to_toml())["model"]
        self.assertEqual(rendered, dict(type="ase_emt", python="python3", ase_version="3.26.0",
             model_sha256="a"*64, atomic_numbers=[29]*4, composition="total",
             pbc=[True]*3, asap_cutoff=False, timeout_ms=30000, transport="one_shot"))
        with tempfile.TemporaryDirectory() as td:
            cfg.model.python = "venv/bin/python"
            resolved = client._resolve_paths(cfg, Path(td))
            self.assertEqual(resolved.model.python, str((Path(td) / "venv/bin/python").resolve()))
            self.assertEqual(cfg.model.python, "venv/bin/python")
            if os.name != "nt":  # Windows symlink creation can require elevated privileges.
                link = Path(td) / "python-link"
                link.symlink_to(sys.executable)
                cfg.model.python = "./python-link"
                self.assertEqual(client._resolve_paths(cfg, Path(td)).model.python,
                                 str(Path(td).resolve() / "python-link"))
            cfg.model.python = "python3"
            self.assertEqual(client._resolve_paths(cfg, Path(td)).model.python, "python3")

    def test_explicit_persistent_transport_is_rendered(self):
        model=modena.AseEmt(ase_version="3.29.0",model_sha256="a"*64,
            atomic_numbers=[29]*4,composition="total",transport="persistent")
        cfg=modena.Config(modena.Fcc(3.61,[1,1,1],mass=63.546),model,
            modena.VelocityVerlet(0.0001),units="real")
        self.assertEqual(tomllib.loads(cfg.to_toml())["model"]["transport"],"persistent")

    def test_metadata_field_is_appended_and_optional_for_legacy_results(self):
        result = modena.RunResult([], "Healthy", None, "", {"verdict":"healthy"})
        self.assertTrue(hasattr(result, "metadata"), "missing optional RunResult metadata")
        self.assertIsNone(result.metadata)
        self.assertEqual(result.energy_diagnostics, {"verdict":"healthy"})

@unittest.skipUnless(os.environ.get("MODENA_ASE_PYTHON"), "set MODENA_ASE_PYTHON for actual ASE CLI reference")
class ActualExternalTests(unittest.TestCase):
    def test_discover_pin_check_and_run_return_identical_core_metadata(self):
        python = os.environ["MODENA_ASE_PYTHON"]
        identity = modena.provider_info(python, binary=BINARY, timeout_ms=120000)
        self.assertEqual(identity["provider"], "ase_emt")
        self.assertEqual(identity["protocol"], "modena-ase-emt-v1")
        cfg = modena.Config(modena.Fcc(3.61, [1,1,1], element="unrelated label", mass=63.546),
            modena.AseEmt(ase_version=identity["ase_version"], model_sha256=identity["model_sha256"],
                atomic_numbers=[29]*4, composition="total", python=python, timeout_ms=120000),
            modena.VelocityVerlet(0.0001), modena.Run(steps=2, initial_temperature=10.0, log_every=100), units="real")
        self.assertIn("OK: 4 atoms", modena.check(cfg, binary=BINARY))
        result = modena.run(cfg, binary=BINARY)
        self.assertEqual(result.metadata["provider"], identity)
        self.assertEqual(result.metadata["pressure_status"], "unconstrained_provider_stress")
        self.assertEqual(result.energy_diagnostics["diagnostics"]["samples"], 3)
        self.assertEqual([row["step"] for row in result.steps_rows], [0])  # sparse log, full-rate diagnostics
        cfg.model.transport="persistent"
        persistent=modena.run(cfg,binary=BINARY)
        self.assertEqual(persistent.metadata,result.metadata)
        self.assertEqual(persistent.energy_diagnostics,result.energy_diagnostics)
        self.assertEqual(persistent.steps_rows,result.steps_rows)
        cfg.model.ase_version = "0.0.0"
        with self.assertRaisesRegex(RuntimeError, "version"):
            modena.run(cfg, binary=BINARY)

if __name__ == "__main__":
    unittest.main()
