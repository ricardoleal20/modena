"""Engineering controls for the admission/receipt layer. These tests launch no physics."""

import pathlib
import tempfile
import unittest
from unittest.mock import patch

from campaign import (
    Budget,
    admit_work,
    checked_gates,
    external_spent,
    gate,
    remaining,
    new_clock,
    verify_files,
)
from freeze import digest, dump, load


class Admission(unittest.TestCase):
    def test_changed_and_missing_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "q"
            p.write_bytes(b"fixed")
            files = {str(p): digest(p)}
            verify_files(files)
            p.write_bytes(b"changed")
            with self.assertRaises(ValueError):
                verify_files(files)
            p.unlink()
            with self.assertRaises(FileNotFoundError):
                verify_files(files)

    def test_expired_and_uncounted_jobs(self):
        self.assertLessEqual(remaining(new_clock(0)),0)
        with self.assertRaises(Budget):
            remaining({"start_monotonic":0,"deadline_monotonic":0,"deadline_utc":0})
        with tempfile.TemporaryDirectory() as d:
            out = pathlib.Path(d)
            folder = out / "jobs" / "ensemble"
            folder.mkdir(parents=True)
            with patch("campaign.LEDGER", out / "isolated-ledger"):
                with self.assertRaises(Budget):
                    admit_work(out, "ensemble")
                dump(folder / "child.json", {"exit": 0})
                (out / "gas-progress").write_text("316800000")
                with self.assertRaises(Budget):
                    admit_work(out, "ensemble")
                dump(folder / "launch.json", {"role": "reference"})
                (folder / "child.json").unlink()
                with self.assertRaises(Budget):
                    external_spent(out)

    def test_gates_and_strict_json(self):
        for invalid in [float("nan"), float("inf"), -float("inf")]:
            with self.assertRaises(ValueError):
                gate("bad", [invalid], [0], [1])
        row = gate("q", [0.1], [0.0], [0.01])
        row["status"] = "Passed"
        with self.assertRaises(ValueError):
            checked_gates([row])
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "bad"
            p.write_text('{"x":1,"x":2}')
            with self.assertRaises(ValueError):
                load(p)


if __name__ == "__main__":
    unittest.main()
