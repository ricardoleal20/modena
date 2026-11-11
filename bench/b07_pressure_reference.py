#!/usr/bin/env python3
"""Prepare or measure the frozen B07 external reference through the shared campaign ledger."""

import argparse
import json
import os
from pathlib import Path
import sys
import contextlib
import signal

sys.path.insert(0, str(Path(__file__).resolve().parent / "b07_validation"))
from campaign import (
    Budget,
    checked_gates,
    child_status,
    digest,
    dump,
    load,
    reconciled_wrapper,
    remaining,
    run_job,
    unmeasured,
    verify,
)
from b07_reference.records import (
    parse_pressure_line as parse_pressure_line,
    confirmed_work,
    native_rows,
)
from b07_reference.epoch import prepare, PAIRS, BASES
from b07_reference.physics import compare, keep, outside_table, snap_checks
from b07_reference.statistics import gas_analysis, interaction
from b07_reference.convergence import ase, nve


@contextlib.contextmanager
def reduction_deadline(clock):
    """The host-side parser/statistic belongs to the same absolute scientific clock."""
    seconds = remaining(clock)
    if seconds <= 0:
        raise Budget("combined deadline reached")

    def expired(signum, frame):
        raise Budget("combined deadline reached during external reduction")

    previous = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def retained_result(out, id):
    folder = out / "jobs" / id
    if (folder / "denied.json").exists():
        return "BudgetExceeded"
    if not (folder / "child.json").exists():
        return "Pending"
    c = load(folder / "child.json")
    if not reconciled_wrapper(folder, c):
        return "Inconclusive"
    code = child_status(c)
    return (
        "Passed"
        if code == 0
        else "BudgetExceeded"
        if c["budget_stop"]
        else "Inconclusive"
    )


def derived(out, x, fixtures, tag):
    s = x["sections"][tag]
    src = out / "reference" / (tag + ".jsonl")
    origin = Path(str(src) + ".origin.json")
    raw = out / "reference" / s["producer"] / (tag + ".dump")
    binding = dict(
        raw_path=str(raw),
        raw_sha256=digest(raw),
        producer_job=s["producer"],
        producer_receipt_sha256=digest(out / "jobs" / s["producer"] / "child.json"),
    )
    if src.exists():
        if load(origin) != dict(binding, sha256=digest(src)):
            raise ValueError("changed native-derived snapshot/receipt")
        return
    partial = Path(str(src) + ".partial")
    with partial.open("x") as f:
        for r in native_rows(out, s, fixtures[s["fixture"]]):
            frame = r["frame"]
            if keep(s, frame["step"]):
                # All q/v and step originate in this single sorted native dump.
                frame["time"] = frame["step"] * s["dt"]
                f.write(
                    json.dumps(frame, separators=(",", ":"), allow_nan=False) + "\n"
                )
        f.flush()
        os.fsync(f.fileno())
    os.link(partial, src)
    partial.unlink()
    dump(origin, dict(binding, sha256=digest(src)))


class Campaign:
    def __init__(self, out, x, p):
        self.out, self.x = out, x
        self.area = out / "reference"
        self.fixtures = {f["id"]: f for f in p["fixtures"]}
        self.roster = load(out / "inputs/roster.json")
        self.results = self.area / "comparisons"
        self.results.mkdir(exist_ok=True)
        self.clock = load(out / "clock.json")

    def launch(self, id):
        out, x = self.out, self.x
        j = x["jobs"][id]
        folder = out / "jobs" / id
        if not folder.exists():
            if remaining(self.clock) <= 0:
                return "BudgetExceeded"
            try:
                run_job(out, id, j["argv"], j["role"])
            except Budget:
                pass  # The shared monitor preserves a durable denial.
        status = retained_result(out, id)
        if (
            j["role"] == "reference"
            and (folder / "launch.json").exists()
            and not (folder / "work.json").exists()
        ):
            sections = [s for s in x["sections"].values() if s["producer"] == id]
            dump(
                folder / "work.json",
                confirmed_work(
                    Path(j["cwd"]) / "native.log", sections, status == "Passed"
                ),
            )
        return status

    def family(self, name, prereqs, fn):
        path = self.results / (name + ".json")
        expected = self.x["expected_gate_families"][name]
        if not path.exists():
            blocked = {k: v for k, v in prereqs.items() if v != "Passed"}
            info = {}
            if blocked:
                status = (
                    "BudgetExceeded"
                    if "BudgetExceeded" in blocked.values()
                    else "Pending"
                )
                rows = [unmeasured(n, status) for n in expected]
                info["blocked_by"] = blocked
            elif remaining(self.clock) <= 0:
                rows = [unmeasured(n, "BudgetExceeded") for n in expected]
                info["reason"] = "combined time exhausted before reduction"
            else:
                try:
                    with reduction_deadline(self.clock):
                        result = fn()
                    rows, info = result if isinstance(result, tuple) else (result, {})
                    checked_gates(rows)
                    if set(r["name"] for r in rows) != set(expected):
                        raise ValueError("incomplete physical gate family")
                except Budget as exc:
                    info = {"budget": str(exc)}
                    rows = [unmeasured(n, "BudgetExceeded") for n in expected]
                except (
                    KeyError,
                    ValueError,
                    TypeError,
                    ArithmeticError,
                    OSError,
                ) as exc:
                    info = {"invalid_or_missing_evidence": str(exc)}
                    rows = [unmeasured(n, "Inconclusive") for n in expected]
            dump(path, dict(gates=rows, details=info))
        rows = checked_gates(load(path)["gates"])
        return next(
            (
                s
                for s in (
                    "Failed",
                    "BudgetExceeded",
                    "Inconclusive",
                    "Unsupported",
                    "Pending",
                )
                if any(r["status"] == s for r in rows)
            ),
            "Passed",
        )

    def evaluate(self, id, dependency):
        x, out = self.x, self.out
        e = x["evaluations"][id]
        if all(v == "Passed" for v in dependency.values()):
            try:
                with reduction_deadline(self.clock):
                    derived(out, x, self.fixtures, e["tag"])
                dependency = dict(dependency, **{id: self.launch(id)})
            except Budget:
                dependency = dict(dependency, **{"snapshot_budget": "BudgetExceeded"})
            except (ValueError, OSError, KeyError) as exc:
                error_path = self.results / (id + "-input-error.json")
                if not error_path.exists():
                    dump(error_path, {"reason": str(exc)})
                dependency = dict(dependency, **{"native_input": "Inconclusive"})
        return self.family(id, dependency, lambda: compare(out, x, self.fixtures, id))

    def run(self):
        x, out, fix = self.x, self.out, self.fixtures
        # Distinct native processes, all launches including preflights enter the same ledger.
        snapshot_status = {
            id: self.launch(id) for id in ("preflight-probes", "preflight-water")
        }
        prerequisites = {}
        spme = {}
        for name in [*PAIRS, *BASES]:
            producer = "preflight-probes" if name in PAIRS else "preflight-water"
            dep = {producer: snapshot_status[producer]}
            tag = name + ("-r0" if name == "pair-oo" else "-r2")
            evaluated = {}
            for backend in ("ewald", "ewald-10", "ewald-12"):
                id = f"evaluate-{tag}-{backend}"
                evaluated[id] = self.evaluate(id, dep)
            prerequisites.update(evaluated)
            # A measured native mismatch must not suppress the independent cutoff comparisons.
            outputs = {id: retained_result(out, id) for id in evaluated}
            convergence = "converge-" + name
            prerequisites[convergence] = self.family(
                convergence,
                dict(dep, **outputs),
                lambda: snap_checks(out, x, fix, name),
            )
            if name == "pair-excluded-out":
                prerequisites["outside-table"] = self.family(
                    "outside-table", dep, lambda: outside_table(out, x, fix)
                )
        for id, e in x["evaluations"].items():
            if e["preflight"] and e["backend"] == "spme":
                spme[id] = self.evaluate(id, prerequisites)
        native_status = {}
        free_gas = {}
        for r in self.roster:
            if r["seed"]:
                continue
            tag = r["id"].replace("spme-", "ewald-")
            dep = {} if r["backend"] == "zero" else dict(prerequisites)
            if r["backend"] == "spme":
                dep.update(spme)
            self.trajectory(r, tag, dep, native_status)
        for group in ("rotor", "water-nonprincipal", "ewald", "spme"):
            cases = [
                r
                for r in self.roster
                if not r["seed"] and group in (r["fixture"], r["backend"])
            ]
            deps = {r["id"]: native_status[r["id"]] for r in cases}
            self.family(
                "nve-" + group, deps, lambda: nve(out, x, fix, self.roster, group)
            )
        self.family(
            "ase-six-strain",
            {"provider": self.launch("ase-six-strain")},
            lambda: ase(out),
        )
        for gas in ("gas-rotor", "gas-water"):
            id = "preflight-" + gas
            free_gas[id] = self.family(
                id,
                {"native": self.launch(id)},
                lambda: gas_analysis(out, x, fix[gas + "-2001"], id),
            )
        for r in self.roster:
            if r["seed"]:
                dep = dict(prerequisites, **free_gas)
                if r["backend"] == "spme":
                    dep.update(spme)
                self.trajectory(
                    r, r["id"].replace("spme-", "ewald-"), dep, native_status
                )
        for backend in ("ewald", "spme"):
            cases = [r for r in self.roster if r["seed"] and backend == r["backend"]]
            self.family(
                "nvt-" + backend,
                {r["id"]: native_status[r["id"]] for r in cases},
                lambda: interaction(out, x, fix, self.roster, backend),
            )

    def trajectory(self, r, tag, dep, statuses):
        native = "native-" + tag
        if all(v == "Passed" for v in dep.values()):
            dep = dict(dep, **{native: self.launch(native)})
        id = "evaluate-" + tag + "-" + r["backend"]
        status = self.evaluate(id, dep)
        rust = "trajectory-" + r["id"]
        statuses[r["id"]] = self.launch(rust) if status == "Passed" else status


def collect(out, x):
    rows = []
    for family, expected in x["expected_gate_families"].items():
        path = out / "reference/comparisons" / (family + ".json")
        if not path.exists():
            rows.extend(unmeasured(name) for name in expected)
            continue
        group = checked_gates(load(path)["gates"])
        if {r["name"] for r in group} != set(expected):
            raise ValueError("family differs from frozen roster")
        rows.extend(group)
    return checked_gates(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("phase", choices=("prepare", "run", "report"))
    for key in ("binary", "protocol", "out"):
        ap.add_argument("--" + key, type=Path, required=True)
    a = ap.parse_args(argv)
    if not all(p.is_absolute() for p in (a.binary, a.protocol, a.out)):
        ap.error("binary, protocol and out must be absolute paths")
    out = a.out.resolve()
    if a.protocol.resolve() != out / "protocol.json":
        ap.error("use the same campaign directory for the linked epochs")
    if a.phase == "prepare":
        if not a.binary.is_file():
            receipt = out / "missing-reference-binary.json"
            if not receipt.exists():
                dump(
                    receipt,
                    dict(status="Pending", binary=str(a.binary), launched=False),
                )
            print("Pending: reference executable missing", file=sys.stderr)
            return 4
        prepare(a.binary.resolve(), a.protocol.resolve(), out)
        return 0
    p, x = verify(out)
    if {j["argv"][0] for j in x["jobs"].values() if j["role"] == "reference"} != {
        str(a.binary.resolve())
    }:
        raise ValueError("requested native binary differs from frozen roster")
    if a.phase == "run":
        if (out / "external-gates.json").exists():
            raise ValueError("external evidence already published")
        Campaign(out, x, p).run()
    rows = collect(out, x)
    if (out / "clock.json").exists() and not (out / "external-gates.json").exists():
        dump(out / "external-gates.json", rows)
    print(
        json.dumps(
            {
                s: sum(r["status"] == s for r in rows)
                for s in sorted({r["status"] for r in rows})
            }
        )
    )
    return (
        3
        if any(r["status"] == "Failed" for r in rows)
        else 2
        if any(r["status"] == "BudgetExceeded" for r in rows)
        else 4
        if any(r["status"] != "Passed" for r in rows)
        else 0
    )


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, Budget) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
