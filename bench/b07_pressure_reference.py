#!/usr/bin/env python3
"""Prepare or measure the frozen B07 external reference through the shared campaign ledger."""

import sys
if __name__ == "__main__" and not (
    sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode
    and sys.pycache_prefix == "/dev/null"
):
    raise SystemExit("use bench/b07_validation/launch.sh INTERPRETER SCRIPT to control imports before startup")
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import argparse
import json
import os
from pathlib import Path
import sys
import contextlib
import signal
import hashlib

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
    raw_bytes,
    LIMITS,
    exclusive,
)
from evidence import (origin_binding, projected_frames, validate_projection, job_state,
    family_binding, read_family, worst, process_status, derivation_receipt)
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
        return "Pending" if not (folder/"launch.json").exists() else "Inconclusive"
    child = load(folder/"child.json")
    role = load(folder/"launch.json").get("role", "rust") if (folder/"launch.json").exists() else "rust"
    return process_status(folder, child, role)[0]


def derived(out, x, fixtures, tag):
    # Own the same lease as launches: only this writer can grow campaign artifacts.
    with exclusive(out):
        return _derived(out, x, fixtures, tag)


def _derived(out, x, fixtures, tag):
    binding = origin_binding(out, x, tag)
    src = Path(binding["input_path"])
    if src.exists():
        validate_projection(out, x, tag, fixtures)
        return
    partial = Path(str(src)+".partial")
    receipt = Path(str(src)+".derivation.json")
    if partial.exists() or load(receipt)["status"] != "Pending":
        raise Budget("native projection already interrupted")
    count = 0
    sha = hashlib.sha256()
    # Reserve enough space for publication metadata and a compact budget-stop receipt.
    reserve = 4096
    used = raw_bytes(out)
    derivation_receipt(receipt, dict(status="Running", confirmed_frames=0))
    try:
        if used+reserve > LIMITS["artifact_bytes"]:
            raise Budget("artifact cap before native projection")
        with partial.open("xb") as stream:
            for data in projected_frames(out, x, fixtures, tag):
                if used+len(data)+reserve > LIMITS["artifact_bytes"]:
                    raise Budget("artifact cap during native projection")
                stream.write(data)
                stream.flush()
                sha.update(data)
                count += 1
                used += len(data)
            os.fsync(stream.fileno())
        if origin_binding(out, x, tag) != binding:
            raise ValueError("native origin changed while deriving")
        sidecar = dict(binding, sha256=sha.hexdigest())
        metadata_bytes = len((json.dumps(sidecar,sort_keys=True,allow_nan=False,separators=(",",":"))+"\n").encode())
        if used+metadata_bytes > LIMITS["artifact_bytes"]:
            raise Budget("artifact cap before projection publication")
        # Atomic hard links do not allocate another data stream. Publish origin first.
        dump(str(src)+".origin.json", sidecar)
        os.link(partial, src)
        partial.unlink()
        derivation_receipt(receipt, dict(status="Completed", confirmed_frames=count, sha256=sha.hexdigest()))
    except BaseException as exc:
        derivation_receipt(receipt, dict(status="BudgetExceeded" if isinstance(exc,Budget) else "Inconclusive",
                                        confirmed_frames=count, reason=str(exc)[:1000]))
        raise


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
            except (ValueError, OSError, KeyError):
                if not (folder/"launch.json").exists():
                    raise
                # Preserve an attempted producer whose monitor/receipt failed.
                # job_state keeps it Inconclusive and family collection records its blocker.
        status = job_state(out, x, id)[0]
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
        if path.exists():
            return worst(r["status"] for r in read_family(self.out, self.x, name))
        info = {}
        try:
            binding = family_binding(self.out, self.x, name)
            actual = {key: value["status"] for key,value in binding["dependencies"].items()}
            blocked = {key:value for key,value in dict(prereqs,**actual).items() if value!="Passed"}
            contract = self.x["family_inputs"][name]
            own = contract.get("owners",contract["producers"])
            attempted = [actual[k] for k in own if actual[k] not in ("Passed","Pending")]
            if blocked:
                status = worst(attempted) if attempted else "Pending"
                rows = [unmeasured(n,status) for n in expected]
                # A definite producer error has no invented U/F/S measurement.
                # Collection independently verifies the seal before accepting empty Failed rows.
                info["blocked_by"] = blocked
                info["producer_diagnostics"] = binding["dependencies"]
            elif remaining(self.clock)<=0:
                rows=[unmeasured(n,"BudgetExceeded") for n in expected]
                info["reason"]="combined time exhausted before reduction"
            else:
                with reduction_deadline(self.clock):
                    result=fn()
                rows,info=result if isinstance(result,tuple) else (result,{})
                checked_gates(rows)
                if {r["name"] for r in rows} != set(expected):
                    raise ValueError("incomplete physical gate family")
            if family_binding(self.out,self.x,name)!=binding:
                raise ValueError("family inputs changed during reduction")
        except Budget as exc:
            info={"budget":str(exc)}
            rows=[unmeasured(n,"BudgetExceeded") for n in expected]
        except (KeyError,ValueError,TypeError,ArithmeticError,OSError) as exc:
            info={"invalid_or_missing_evidence":str(exc)}
            rows=[unmeasured(n,"Inconclusive") for n in expected]
            binding=None
        dump(path, dict(gates=rows,details=info,binding=binding,
                        epoch_sha256=digest(self.out/"external-epoch.json")))
        return worst(r["status"] for r in rows)

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
            outputs = {id: job_state(out, x, id)[0] for id in evaluated}
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
    return [row for name in x["expected_gate_families"]
            for row in read_family(out,x,name)]


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
