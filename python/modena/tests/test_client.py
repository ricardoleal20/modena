"""Dependency-free regressions for CLI discovery and temporary input files."""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modena import Config, Fcc, LennardJones, VelocityVerlet, check, run

client = importlib.import_module("modena.run")


class BinaryTests(unittest.TestCase):
    def test_honors_cargo_target_dir_for_a_source_build(self):
        # A custom Cargo target directory must work without changing PATH.
        with tempfile.TemporaryDirectory() as td:
            binary = Path(td) / "release" / "modena"
            binary.parent.mkdir(parents=True)
            binary.write_text("#!/bin/sh\n")
            binary.chmod(0o755)
            with patch("shutil.which", return_value=None):
                with patch.dict(os.environ, {"CARGO_TARGET_DIR": td}):
                    self.assertEqual(client._binary(), str(binary.resolve()))

    def test_missing_binary_explains_how_to_build_or_select_it(self):
        with tempfile.TemporaryDirectory() as td:
            with patch.object(client, "__file__", str(Path(td) / "python/modena/run.py")):
                with patch("shutil.which", return_value=None):
                    with patch.dict(os.environ, {"CARGO_TARGET_DIR": ""}):
                        with self.assertRaisesRegex(FileNotFoundError, "binary="):
                            client._binary()

    def test_nonexecutable_build_artifact_is_not_selected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            invalid = root / "target" / "release" / "modena"
            invalid.parent.mkdir(parents=True)
            invalid.write_text("not an executable")
            valid = root / "target" / "debug" / "modena"
            valid.parent.mkdir(parents=True)
            valid.write_text("#!/bin/sh\n")
            valid.chmod(0o755)
            with patch.object(client, "__file__", str(root / "python/modena/run.py")):
                with patch("shutil.which", return_value=None):
                    with patch.dict(os.environ, {"CARGO_TARGET_DIR": ""}):
                        self.assertEqual(client._binary(), str(valid.resolve()))


class TemporaryInputTests(unittest.TestCase):
    def test_serialization_failure_leaves_no_temporary_input(self):
        cfg = Config(
            Fcc(2.4, [1, 1, 1]), LennardJones(2.5, 1.0, 1.0),
            VelocityVerlet(0.001),
        )
        cfg.structure.element = object()  # Invalid Python value fails serialization.
        for operation in (check, run):
            with self.subTest(operation=operation.__name__):
                with tempfile.TemporaryDirectory() as td:
                    with patch.object(tempfile, "tempdir", td):
                        with self.assertRaises(TypeError):
                            operation(cfg, binary="unused")
                    self.assertEqual(list(Path(td).iterdir()), [])

    def test_cli_failure_is_reported_before_parsing_an_incomplete_log(self):
        cfg = Config(
            Fcc(2.4, [1, 1, 1]), LennardJones(2.5, 1.0, 1.0),
            VelocityVerlet(0.001),
        )

        def fail(args, **kwargs):
            parsed = tomllib.loads(Path(args[-1]).read_text())
            Path(parsed["run"]["log_file"]).write_text("step,total\n0,broken\n")
            return subprocess.CompletedProcess(args, 1, "", "simulation failed: fixture")

        with patch.object(client.subprocess, "run", side_effect=fail):
            with self.assertRaisesRegex(RuntimeError, "simulation failed: fixture"):
                run(cfg, binary="unused")


if __name__ == "__main__":
    unittest.main()
