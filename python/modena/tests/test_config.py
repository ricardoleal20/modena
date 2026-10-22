"""Regression tests for the Python-to-TOML boundary (stdlib only)."""

from __future__ import annotations

import sys
import tomllib
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from modena import Config, Fcc, LennardJones, Run, VelocityVerlet
from modena.config import Buckingham, Morse


class ConfigTests(unittest.TestCase):
    def test_strings_round_trip_quotes_backslashes_controls_and_unicode(self):
        # A raw quote/backslash used to corrupt TOML or silently change a path.
        value = 'Ar"\\\n\t\r\b\f\x00\x1f\x7f — 🧪'
        cfg = Config(
            Fcc(1.7, [2, 2, 2], element=value),
            LennardJones(2.5, sigma=1.0, epsilon=1.0),
            VelocityVerlet(0.005),
            Run(trajectory_file=value),
        )
        parsed = tomllib.loads(cfg.to_toml())
        self.assertEqual(parsed["structure"]["element"], value)
        self.assertEqual(parsed["run"]["trajectory_file"], value)

    def test_species_tables_and_overrides_belong_to_the_model(self):
        cases = (
            (LennardJones, {"sigma": 1.0, "epsilon": 1.0}, "lennard_jones"),
            (Buckingham, {"a": 10.0, "b": 1.0, "rho": 0.5}, "buckingham"),
            (Morse, {"d": 1.0, "alpha": 2.0, "r0": 1.0}, "morse"),
        )
        for cls, parameters, kind in cases:
            with self.subTest(model=kind):
                override = {"types": [0, 1], **parameters}
                cfg = Config(
                    Fcc(1.7, [2, 2, 2]),
                    cls(2.5, types=[parameters, parameters], overrides=[override]),
                    VelocityVerlet(0.005),
                )
                parsed = tomllib.loads(cfg.to_toml())
                self.assertEqual(parsed["model"]["type"], kind)
                self.assertEqual(parsed["model"]["types"], [parameters, parameters])
                self.assertEqual(parsed["model"]["overrides"], [override])
                self.assertNotIn("types", parsed)
                self.assertNotIn("overrides", parsed)


if __name__ == "__main__":
    unittest.main()
