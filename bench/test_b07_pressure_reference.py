import unittest

from b07_pressure_reference import parse_pressure_line


class PressureParserTests(unittest.TestCase):
    def test_pressure_columns_preserve_shear(self):
        row = parse_pressure_line("12 0.12 1 2 3 .4 .5 .6")
        self.assertEqual(row["step"], 12)
        self.assertEqual(row["pressure"], [1.0, 2.0, 3.0, 0.4, 0.5, 0.6])


import contextlib
import json
import os
import sys
import time
import math
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch

import b07_pressure_reference as ref
from b07_reference import convergence, epoch, native, physics, records, statistics

HERE = Path(__file__).parent
SYNTHETIC_ROTOR = {
    "id": "rotor",
    "masses": [1.0, 3.0],
    "positions": [[1.0, 2.0, 3.0], [3.0, 2.0, 3.0]],
    "velocities": [[0.0, -1.5, 0.0], [0.0, 0.5, 0.0]],
    "charges": [0.0, 0.0],
    "constraints": {"waters": [], "bonds": [{"i": 0, "j": 1, "d": 2.0}]},
    "pair_scales": [],
    "cell": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
}


def output(dir, steps=(0,), pressure=(1.0, 2.0, 3.0, 0.4, 0.5, 0.6)):
    """Deliberately synthetic unwrapped native transcript; no numerical engine call."""
    folder = dir / "reference/probe"
    folder.mkdir(parents=True, exist_ok=True)
    s = dict(
        fixture="rotor",
        tag="probe",
        producer="probe",
        steps=steps[-1],
        dt=0.01,
        interval=1,
        burn_steps=0,
        molecules=1,
        mode="snapshot",
    )
    log = ["B07_BEGIN probe", " ".join(records.HEADER)]
    raw = []
    f = SYNTHETIC_ROTOR
    for step in steps:
        vals = [step, step * 0.01, 0.0, *pressure, *pressure, *([0.0] * 12)]
        log.append(" ".join(map(str, vals)))
        raw += [
            "ITEM: TIMESTEP",
            str(step),
            "ITEM: NUMBER OF ATOMS",
            "2",
            "ITEM: BOX BOUNDS pp pp pp",
            *(["0 10"] * 3),
            records.ATOMS,
        ]
        raw += [
            " ".join(map(str, [i + 1, *q, *v, 0, 0, 0]))
            for i, (q, v) in enumerate(zip(f["positions"], f["velocities"]))
        ]
    log += ["Loop time", "B07_END probe"]
    (folder / "native.log").write_text("\n".join(log) + "\n")
    (folder / "probe.dump").write_text("\n".join(raw) + "\n")
    return s, folder


class StrictRecords(unittest.TestCase):
    def test_retained_real_columns_and_units(self):
        lines = (
            (HERE / "fixtures/b07/lammps-thermo-excerpt.txt").read_text().splitlines()
        )
        self.assertEqual(lines[3].split(), ["units", "metal"])
        dt = float(lines[4].split()[1])
        self.assertEqual(
            lines[5].split(), ["Step", "Temp", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"]
        )
        # Authentic pressure in bar/time in ps. Do not feed metal-unit evidence to LJ gates.
        step, temperature, *fields = lines[7].split()
        self.assertEqual(float(temperature), 461.95579)
        row = parse_pressure_line(" ".join([step, str(int(step) * dt), *fields]))
        self.assertEqual(
            row,
            dict(
                step=100,
                time=0.2,
                pressure=[
                    11262.405,
                    9918.4702,
                    7373.1896,
                    1389.9833,
                    -165.54737,
                    -128.04989,
                ],
            ),
        )

    def test_shape_order_and_nonfinite(self):
        for bad in ["12 .1 1 2", "-1 .1 1 2 3 4 5 6", "12 .1 1 2 3 4 5 6 7"]:
            with self.assertRaises(ValueError):
                parse_pressure_line(bad)
        for field in range(1, 8):
            for token in ["nan", "inf", "1e999"]:
                fields = "12 .12 1 2 3 .4 .5 .6".split()
                fields[field] = token
                with self.assertRaises(ValueError):
                    parse_pressure_line(" ".join(fields))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            s, folder = output(out)
            a = list(records.native_rows(out, s, SYNTHETIC_ROTOR))[0]
            self.assertEqual(a["pressure"], [1.0, 2.0, 3.0, 0.4, 0.5, 0.6])
            log = folder / "native.log"
            wrong = (
                (HERE / "fixtures/b07/synthetic-swapped-pressure.txt")
                .read_text()
                .splitlines()
            )
            log.write_text(log.read_text().replace(wrong[1], wrong[2]))
            with self.assertRaises(ValueError):
                list(records.native_rows(out, s, SYNTHETIC_ROTOR))

    def test_unbalanced_dump_id_and_partial_work(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            s, folder = output(out, (0, 1))
            s["steps"] = 3
            s["molecules"] = 32
            self.assertEqual(
                records.confirmed_work(folder / "native.log", [s], False),
                dict(system_steps=1, molecule_steps=32, complete=False),
            )
            with self.assertRaises(ValueError):
                list(records.native_rows(out, s, SYNTHETIC_ROTOR))
            p = folder / "probe.dump"
            p.write_text(p.read_text().replace("1 1.0 2.0 3.0", "2 1.0 2.0 3.0"))
            with self.assertRaisesRegex(ValueError, "IDs"):
                list(records.dumps(p, SYNTHETIC_ROTOR))

    def test_derived_frames_producer_and_tampering(self):
        setup = ReviewedEvidence("runTest")
        setup.setUp()
        self.addCleanup(setup.doCleanups)
        setup.seal()
        out, x, fixtures, folder = setup.out, setup.x, setup.fix, setup.folder
        ref.derived(out, x, fixtures, "probe")
        src = out / "reference/probe.jsonl"
        origin = ref.load(Path(str(src) + ".origin.json"))
        self.assertEqual(origin["sha256"],ref.digest(src))
        self.assertEqual(origin["raw_path"],str(folder/"probe.dump"))
        self.assertEqual(origin["producer_job"],"probe")
        self.assertEqual(len(src.read_text().splitlines()), 2)
        src.write_text("{}\n")
        with self.assertRaises(ValueError):
            ref.derived(out, x, fixtures, "probe")

    def test_rotor_same_snapshot_production_reaction_not_analytic_substitute(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            s, _ = output(out)
            f = SYNTHETIC_ROTOR
            row = next(records.native_rows(out, s, f))
            force, p = physics.reaction(f, row)
            self.assertEqual(force, [1.5, 0.0, 0.0, -1.5, 0.0, 0.0])
            self.assertEqual(p, [-0.003, 0.003, 0.0, 0.0, 0.0, 0.0])
            id = "evaluate-probe-zero"
            expected = dict(
                fixture="rotor",
                backend="zero",
                potential_energy=0.0,
                forces=[0.0] * 6,
                conservative_tensor=[0.0] * 6,
                endpoint=p,
                production_reaction_forces=[force[:3], force[3:]],
                frame=row["frame"],
                analytic_force_pressure=[[force[:3], force[3:]], p],
            )
            path = out / "reference" / (id + ".jsonl")
            x = dict(
                sections={"probe": s},
                evaluations={id: dict(tag="probe", backend="zero")},
                jobs={id: dict(expected_frames=1)},
            )

            def compare(e):
                path.write_text(json.dumps(e) + "\n")
                return physics.compare(out, x, {"rotor": f}, id)

            self.assertTrue(all(r["status"] == "Passed" for r in compare(expected)))
            expected["production_reaction_forces"][0][0] *= -1
            self.assertEqual(
                [r["name"] for r in compare(expected) if r["status"] == "Failed"],
                [id + "-reaction"],
            )
            expected["forces"][2] = math.nan
            with self.assertRaises(ValueError):
                compare(expected)

    def test_water_euler_cartesian_forces_and_pressure(self):
        f = dict(SYNTHETIC_ROTOR, masses=[1.0, 1.0, 1.0])
        q = [[4.0, 4.0, 4.0], [5.0, 4.0, 4.0], [4.0, 5.0, 4.0]]
        v = [[-(r[1] - 13 / 3), r[0] - 13 / 3, 0.0] for r in q]
        raw = dict(
            frame=dict(positions=q, velocities=v),
            forces=[0.0] * 9,
            conservative_tensor=[0.0] * 6,
        )
        force, p = physics.reaction(f, raw)
        expected = [k for r in q for k in (13 / 3 - r[0], 13 / 3 - r[1], 0.0)]
        for actual, truth in zip(force, expected):
            self.assertAlmostEqual(actual, truth)
        for actual, truth in zip(p, [0.0, 0.0, 0.0, 2 / 3000, 0.0, 0.0]):
            self.assertAlmostEqual(actual, truth)

    def test_bounded_comparisons_validate_every_component(self):
        for k in range(6):
            for bad in (math.nan, math.inf):
                a = [0.0] * 6
                a[k] = bad
                with self.assertRaises(ValueError):
                    physics.bounded("negative", a, [0.0] * 6, 1.0)
                with self.assertRaises(ValueError):
                    physics.bounded("negative", [0.0] * 6, a, 1.0)
        self.assertAlmostEqual(statistics.t95(3), 3.182446305, places=7)
        self.assertAlmostEqual(statistics.t95(19), 2.093024054, places=7)
        self.assertIsNone(statistics.ips([0.0] * 40))
        self.assertAlmostEqual(
            convergence.intercept(
                [0.08, 0.04, 0.02, 0.01],
                [2 + 7 * h * h for h in (0.08, 0.04, 0.02, 0.01)],
            ),
            2.0,
        )


class SourceGates(unittest.TestCase):
    def test_real_subprocess_failure_timeout_and_missing_executable(self):
        for script, error in [("exit 7", "failed"), ("exec sleep 5", "timed out")]:
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp)
                binary = out / "synthetic-help"
                binary.write_text("#!/bin/sh\necho synthetic >&2\n" + script + "\n")
                binary.chmod(0o700)
                runner = subprocess.run

                def short_timeout(args, **kwargs):
                    kwargs["timeout"] = 0.5 if error == "timed out" else 5
                    return runner(args, **kwargs)

                with (
                    patch.object(epoch.subprocess, "run", side_effect=short_timeout),
                    self.assertRaisesRegex(ValueError, error),
                ):
                    epoch.help_output(binary, out)
                result = ref.load(out / "native-help-status.json")
                self.assertIn("synthetic", (out / "native-help.stderr").read_text())
                self.assertEqual(result["exit"], 7 if error == "failed" else None)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp).resolve()
            with patch.object(
                ref, "prepare", side_effect=AssertionError("must not prepare")
            ):
                with patch("sys.stderr"):
                    status = ref.main(
                        [
                            "prepare",
                            "--binary",
                            str(out / "absent"),
                            "--protocol",
                            str(out / "protocol.json"),
                            "--out",
                            str(out),
                        ]
                    )
            self.assertEqual(status, 4)
            self.assertEqual(
                ref.load(out / "missing-reference-binary.json")["status"], "Pending"
            )

    def test_missing_package_and_source_hash(self):
        styles = {
            "Atom": "full",
            "Fix": "rattle store/force rigid/nvt nve",
            "KSpace": "ewald",
            "Pair": "zero hybrid/overlay coul/long table",
            "Bond": "zero table",
            "Angle": "zero",
            "Compute": "pressure temp",
            "Dump": "custom",
        }
        h = "22 Jul 2025 - Update 6\nFFT precision  = double\nMPI STUBS\nCompiler: synthetic\nInstalled packages:\nRIGID KSPACE MOLECULE\nList of individual\n"
        h += (
            "\n".join("* " + g + " styles:\n" + v + "\n" for g, v in styles.items())
            + "\n* Dummy styles:"
        )
        self.assertEqual(epoch.inspect_help(h)[1], "synthetic")
        for remove in ("KSPACE", "rattle", "coul/long", "pressure"):
            with self.assertRaises(ValueError):
                epoch.inspect_help(h.replace(remove, "missing"))
        with tempfile.TemporaryDirectory() as tmp:
            sample = Path(tmp) / "binary"
            sample.write_bytes(b"synthetic")
            with self.assertRaises(ValueError):
                epoch.verify_files({str(sample): "0" * 64})

    def test_failed_denied_and_timeout_never_report_native_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            folder = out / "jobs/probe"
            folder.mkdir(parents=True)
            self.assertEqual(ref.retained_result(out, "probe"), "Pending")
            with patch("campaign.reconciled_wrapper", return_value=True):
                for exit, reason, status in [
                    (1, None, "Inconclusive"),
                    (-9, "wall cap", "BudgetExceeded"),
                    (0, None, "Passed"),
                ]:
                    (folder / "child.json").write_text(
                        json.dumps(dict(exit=exit, budget_stop=reason))
                    )
                    self.assertEqual(ref.retained_result(out, "probe"), status)
            with patch("campaign.reconciled_wrapper", return_value=False):
                self.assertEqual(ref.retained_result(out, "probe"), "Inconclusive")
            (folder / "denied.json").write_text("{}")
            self.assertEqual(ref.retained_result(out, "probe"), "BudgetExceeded")

    def test_failed_prerequisite_blocks_action_retains_named_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "reference").mkdir()
            (out / "inputs").mkdir()
            ref.dump(out / "inputs/roster.json", [])
            ref.dump(out / "clock.json", {})
            x = dict(expected_gate_families={"family": ["wanted-a", "wanted-b"]},
                     family_inputs={"family":{"producers":[],"families":[]}})
            ref.dump(out/"external-epoch.json",x)
            c = ref.Campaign(out, x, {"fixtures": []})
            with patch.object(
                ref, "run_job", side_effect=AssertionError("must not launch")
            ):
                status = c.family(
                    "family",
                    {"preflight": "Failed"},
                    lambda: self.fail("dependent comparison"),
                )
            self.assertEqual(status, "Pending")
            self.assertEqual(
                [r["name"] for r in ref.collect(out, x)], ["wanted-a", "wanted-b"]
            )
            self.assertEqual(
                ref.load(c.results / "family.json")["details"]["blocked_by"],
                {"preflight": "Failed"},
            )

    def test_frozen_native_inputs_and_ase_shape(self):
        f = SYNTHETIC_ROTOR
        base = dict(tag="synthetic", mode="nve", dt=0.01, steps=100, interval=1)
        a = native.section(f, base, Path("/synthetic"))
        self.assertLess(
            a.index("fix cons all store/force"), a.index("fix rigidconstraint")
        )
        self.assertIn("thermo_modify norm no format float %.17g", a)
        b = native.section(f, dict(base, mode="nvt"), Path("/synthetic"))
        self.assertIn("rigid/nvt molecule temp 1 1 1 tparam 3 1 3", b)
        self.assertNotIn("all nve", b)
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            epoch.ase_input(d)
            cfg = ref.load(d / "ase.json")
            self.assertEqual(len(cfg["snapshots"]), 49)
            for i in range(6):
                self.assertEqual(
                    [s["eta"] for s in cfg["snapshots"] if s["component"] == i],
                    [-1e-4, 1e-4, -3e-5, 3e-5, -1e-5, 1e-5, -3e-6, 3e-6],
                )


class SyntheticProvider(unittest.TestCase):
    def test_all_four_eps_and_request_journal(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            area = out / "reference"
            static = area / "inputs"
            static.mkdir(parents=True)
            epoch.ase_input(static)
            cfg = ref.load(static / "ase.json")
            events = [
                dict(id=i, kind=kind, stage=stage)
                for i, kind in enumerate(["info", "wire-base"] + ["reduced"] * 49, 1)
                for stage in ("attempt", "complete")
            ]

            def stream(path, rows):
                path.write_text("".join(json.dumps(row) + "\n" for row in rows))

            stream(area / "ase-requests.jsonl", events)
            tens = [1.0, 2.0, 3.0, 0.4, 0.5, 0.6]
            unit = 1.602176634e-19 * 6.02214076e23 / 4184
            stress = [[0.0] * 3 for _ in range(3)]
            for (a, b), v in zip(physics.AXES, tens):
                stress[a][b] = stress[b][a] = -v / (3.6**3 * unit)
            wire = dict(
                forces_ev_per_angstrom=[[0.0] * 3 for _ in range(4)],
                energy_ev=0.0,
                stress_ev_per_angstrom3=stress,
            )
            rows = [
                dict(
                    completed_info_requests=1,
                    identity={
                        k: cfg["provider"][k] for k in ("ase_version", "model_sha256")
                    },
                )
            ]
            rows[0]["identity"]["numpy_version"] = "2.2.6"
            for i, snap in enumerate(cfg["snapshots"]):
                rows.append(
                    dict(
                        snap,
                        index=i,
                        completed_evaluate_requests=i + 2,
                        forces=[0.0] * 12,
                        virial=tens,
                        trace=6.0,
                        volume=3.6**3,
                        energy=-snap["eta"] * tens[snap["component"] or 0],
                        wire=wire if i == 0 else None,
                    )
                )
            stream(area / "ase.jsonl", rows)
            gates, details = convergence.ase(out)
            self.assertTrue(all(r["status"] == "Passed" for r in gates))
            self.assertEqual(len(details["components"]["yz"]["derivatives"]), 4)
            rows[1]["virial"] = [1.0, 2.0, 3.0, 0.6, 0.5, 0.4]
            stream(area / "ase.jsonl", rows)
            failed = {
                r["name"] for r in convergence.ase(out)[0] if r["status"] == "Failed"
            }
            self.assertEqual(
                failed, {"ase-wire-conversion", "ase-strain-xy", "ase-strain-yz"}
            )
            stream(area / "ase-requests.jsonl", events[:-1])
            with self.assertRaisesRegex(ValueError, "journal"):
                convergence.ase(out)


class HostPrerequisites(unittest.TestCase):
    def test_direct_host_entry_rejects_failed_required_family(self):
        import contextlib
        import campaign

        with tempfile.TemporaryDirectory() as td:
            out = Path(td).resolve()
            folder = out / "reference/comparisons"
            folder.mkdir(parents=True)
            (out / "jobs").mkdir()
            g = campaign.gate("physical", [2.0], [0.0], [1.0])
            campaign.dump(folder / "table.json", dict(gates=[g]))
            x = {"expected_gate_families": {"table": ["physical"]}}
            roster = {
                "native": {
                    "argv": ["/synthetic"],
                    "role": "reference",
                    "requires_families": ["table"],
                }
            }
            with (
                patch.object(
                    campaign, "exclusive", return_value=contextlib.nullcontext()
                ),
                patch.object(campaign, "verify", return_value=(None, x)),
                patch.object(campaign, "jobs", return_value=roster),
            ):
                with self.assertRaisesRegex(
                    ValueError, "blocked by external prerequisite"
                ):
                    campaign.run_job(out, "native")
            self.assertFalse((out / "jobs/native").exists())

class ReviewedEvidence(unittest.TestCase):
    def setUp(self):
        import campaign
        self.host=campaign
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out=Path(self.tmp.name).resolve()
        out=self.out
        self.s,self.folder=output(out,(0,1))
        (out/"inputs").mkdir()
        (out/"jobs").mkdir()
        ref.dump(out/"inputs/roster.json",[])
        ref.dump(out/"clock.json",dict(start_monotonic=time.monotonic(),deadline_monotonic=time.monotonic()+600,deadline_utc=time.time()+600))
        ref.dump(out/"protocol.json",{"fixtures":[SYNTHETIC_ROTOR]})
        self.fix={"rotor":SYNTHETIC_ROTOR}
        self.cfg=dict(argv=["/usr/bin/true"],cwd=str(self.folder),env={},role="reference",
                      artifacts=[str(self.folder/"native.log"),str(self.folder/"probe.dump")])
        self.x=dict(sections={"probe":self.s},jobs={"probe":self.cfg},expected_gate_families={"primary":["one"]},
                    family_inputs={"primary":{"producers":["probe"],"families":[]}})
        ref.dump(out/"external-epoch.json",self.x)
        from evidence import derivation_receipt
        derivation_receipt(out/"reference/probe.jsonl.derivation.json",{"status":"Pending"},create=True)
        self.f=out/"jobs/probe"
        self.f.mkdir()
        ref.dump(self.f/"launch.json",dict(self.cfg,name="probe",protocol_sha256=ref.digest(out/"protocol.json"),binary_sha256=ref.digest("/usr/bin/true")))
        ref.dump(self.f/"child.json",dict(exit=0,budget_stop=None,interrupted=None))
        for name in ("started.json","wrapper.json","host.time","stdout.log","stderr.log","wrapper.stdout","wrapper.stderr"):
            (self.f/name).write_text("{}")
        lease=patch.object(ref,"exclusive",return_value=contextlib.nullcontext())
        lease.start()
        self.addCleanup(lease.stop)
        self.mock=patch.object(campaign,"reconciled_wrapper",return_value=True)
        self.mock.start()
        self.addCleanup(self.mock.stop)

    def seal(self):
        from evidence import seal_job
        seal_job(self.out,self.x,"probe")

    def test_cached_family_rechecks_changed_raw_direct_and_collect(self):
        self.seal()
        c=ref.Campaign(self.out,self.x,{"fixtures":[]})
        g=self.host.gate("one",[1.],[1.],[0.])
        self.assertEqual(c.family("primary",{},lambda:[g]),"Passed")
        (self.folder/"native.log").write_text("replaced thermo\n")
        self.assertEqual(c.family("primary",{},lambda:self.fail()),"Inconclusive")
        self.assertEqual(ref.collect(self.out,self.x)[0]["status"],"Inconclusive")
        from evidence import require_passes
        with self.assertRaises(ValueError):
            require_passes(self.out,self.x,["primary"])

    def test_deleted_matched_output_invalidates_seal(self):
        path=self.out/"matched.jsonl"
        path.write_text("{}\n")
        self.cfg["artifacts"].append(str(path))
        self.f.joinpath("launch.json").unlink()
        ref.dump(self.f/"launch.json",dict(self.cfg,name="probe",protocol_sha256=ref.digest(self.out/"protocol.json"),binary_sha256=ref.digest("/usr/bin/true")))
        self.seal()
        c=ref.Campaign(self.out,self.x,{"fixtures":[]})
        self.assertEqual(c.family("primary",{},lambda:[self.host.gate("one",[0.],[0.],[0.])]),"Passed")
        path.unlink()
        self.assertEqual(ref.collect(self.out,self.x)[0]["status"],"Inconclusive")

    def test_exact_origin_and_reprojection(self):
        self.seal()
        ref.derived(self.out,self.x,self.fix,"probe")
        src=self.out/"reference/probe.jsonl"
        dest=Path(str(src)+".origin.json")
        original=ref.load(dest)
        for change in ("producer","raw","payload"):
            obj=dict(original)
            if change=="producer": obj["producer_job"]="another"
            elif change=="raw": obj["raw_path"]=str(self.folder/"other.dump")
            else:
                src.write_bytes(src.read_bytes().replace(b"1.5",b"1.6"))
                obj["sha256"]=ref.digest(src)
            dest.write_text(json.dumps(obj))
            with self.subTest(change=change),self.assertRaises(ValueError):
                ref.derived(self.out,self.x,self.fix,"probe")

    def test_projection_budget_before_append(self):
        self.seal()
        cap=self.host.raw_bytes(self.out)+4200
        with patch.dict(self.host.LIMITS,artifact_bytes=cap),self.assertRaises(self.host.Budget):
            ref.derived(self.out,self.x,self.fix,"probe")
        self.assertFalse((self.out/"reference/probe.jsonl").exists())
        self.assertEqual(ref.load(self.out/"reference/probe.jsonl.derivation.json")["status"],"BudgetExceeded")
        self.assertLessEqual(self.host.raw_bytes(self.out),cap)

    def test_actual_producer_failure_in_both_ledgers(self):
        import evidence
        for code,diagnostic,expected in [(1,"ERROR: Lost atoms: original 96 current 90","Failed"),(1,"ERROR: could not open file","Inconclusive")]:
            (self.f/"child.json").write_text(json.dumps(dict(exit=code,budget_stop=None,interrupted=None)))
            (self.f/"stdout.log").write_text(diagnostic)
            evidence.seal_job(self.out,self.x,"probe")
            self.assertEqual(evidence.job_state(self.out,self.x,"probe")[0],expected)
            c=ref.Campaign(self.out,self.x,{"fixtures":[]})
            self.assertEqual(c.family("primary",{"probe":expected},lambda:self.fail()),expected)
            self.assertEqual(ref.collect(self.out,self.x)[0]["status"],expected)
            for path in [self.f/"artifacts.json",*c.results.iterdir()]: path.unlink()
        self.f.joinpath("launch.json").unlink()
        self.f.joinpath("child.json").unlink()
        self.assertEqual(evidence.job_state(self.out,self.x,"probe")[0],"Pending")

    def test_clean_runtime_and_worker_visibility(self):
        import campaign
        self.assertIn("pycache_prefix=/dev/null",campaign.PYTHON_FLAGS)
        with patch.dict(os.environ,LAMMPS_PLUGIN_PATH="/ambient",DYLD_LIBRARY_PATH="/ambient"):
            env=campaign.numerical_env(campaign.ENV)
            self.assertNotIn("DYLD_LIBRARY_PATH",env)
            self.assertNotIn("LAMMPS_PLUGIN_PATH",env)
        observed={}
        with self.assertRaisesRegex(ValueError,"unreadable"):
            campaign.sample_ase_workers(observed,[19],lambda pid:None,lambda pid:None,lambda pid:True,"/provider")
        self.assertEqual(campaign.sample_ase_workers(observed,[19],lambda pid:None,lambda pid:None,lambda pid:False,"/provider"),0)
        self.assertEqual(campaign.sample_ase_workers(observed,[23],lambda pid:100,lambda pid:"/provider",lambda pid:True,"/provider"),100)
        self.assertEqual(observed["23"]["interpreter_samples"],1)
