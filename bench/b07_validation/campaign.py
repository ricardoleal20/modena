"""Shared B07 campaign ledger and child monitor. All science runs are serial."""

import argparse
import contextlib
import ctypes
import datetime
import fcntl
import math
import os
import platform
import re
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path

from freeze import (
    EVIDENCE,
    LIMITS,
    MESH,
    ROOT,
    SPEC,
    digest,
    dump,
    host_identity,
    inventory,
    load,
    prepare,
)

ENV = {
    "RAYON_NUM_THREADS": "1",
    "RUST_TEST_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}
LEDGER = EVIDENCE / "b07-campaign-budget-v1"
OU = "integrator::conditional_ou::fixed_water_conditional_momentum_law"


class Budget(Exception):
    pass


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def epochs(out):
    entries = sorted(LEDGER.glob("epoch-*.json")) if LEDGER.is_dir() else []
    return list(dict.fromkeys([Path(load(p)["out"]) for p in entries] + [out]))


def raw_bytes(out):
    total = 0
    for epoch in epochs(out):
        for path in epoch.rglob("*"):
            try:
                if path.is_file():
                    total += path.stat().st_size
            except FileNotFoundError:
                continue  # Atomic progress replacement can race a directory scan.
    return total


def verify_files(files):
    if not isinstance(files, dict) or not files:
        raise ValueError("empty inventory")
    for name, sha in files.items():
        if not Path(name).is_absolute() or digest(name) != sha:
            raise ValueError("frozen file changed: " + name)


def verify(out, both=True):
    p = load(out / "protocol.json")
    pin = load(out / "internal-pin.json")["protocol_sha256"]
    if (
        digest(out / "protocol.json") != pin
        or digest(SPEC) != p["specification_sha256"]
        or p["limits"] != LIMITS
    ):
        raise ValueError("protocol changed")
    verify_files(p["source_sha256"])
    verify_files(p["snapshot_sha256"])
    if load(out / "inputs/settings.json")["host_identity"] != host_identity():
        raise ValueError("host or shared library cache changed")
    fixtures = [
        load(f)
        for f in sorted((out / "inputs").glob("*.json"))
        if f.name not in ("roster.json", "settings.json")
    ]
    if (
        p["fixtures"] != fixtures
        or len(fixtures) != 19
        or len({f["id"] for f in fixtures}) != 19
    ):
        raise ValueError("immutable fixture array mismatch")
    x = None
    if both:
        x = load(out / "external-epoch.json")
        if x["internal_protocol_sha256"] != pin:
            raise ValueError("external epoch targets another core")
        verify_files(x["source_sha256"])
        verify_files(x["snapshot_sha256"])
        if (out / "clock.json").exists():
            clock = load(out / "clock.json")
            if clock["internal"] != pin or clock["external"] != digest(
                out / "external-epoch.json"
            ):
                raise ValueError("pinned measurement epoch changed")
    return p, x


@contextlib.contextmanager
def exclusive(out):
    LEDGER.mkdir(exist_ok=True)
    with (LEDGER / "lease.lock").open("a") as lease:
        fcntl.flock(lease, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def internal_jobs(out):
    exe = str(out / "bin/stress_validation")

    def job(cmd):
        return {
            "role": "rust",
            "argv": list(map(str, cmd)),
            "cwd": str(ROOT),
            "env": ENV,
        }

    a = {
        k: job([exe, "--child", k, out])
        for k in ("deterministic", "ou", "costs", "ensemble")
    }
    a["ou-draws"] = job(
        [out / "bin/ou-libtest", OU, "--exact", "--nocapture", "--test-threads=1"]
    )
    for label in ("old", "new"):
        for fixture in ("rotor", "gas-water"):
            a[f"cost-{label}-{fixture}"] = job(
                [
                    out / "bin" / ("whole-step-" + label),
                    out / "inputs" / (fixture + ".json"),
                ]
            )
    return a


def jobs(out, x):
    result = internal_jobs(out)
    for name, j in x["jobs"].items():
        if name in result or not re.fullmatch(r"[A-Za-z0-9_-]+", name):
            raise ValueError("invalid or duplicate job ID")
        if (
            j["role"] not in ("rust", "reference")
            or not Path(j["cwd"]).is_absolute()
            or not Path(j["argv"][0]).is_absolute()
        ):
            raise ValueError("unsupported role/cwd/argv")
        if not all(j["env"].get(k) == v for k, v in ENV.items()):
            raise ValueError("numeric/thread environment changed")
        if j["role"] == "rust" and j["argv"][0] != str(out / "bin/stress_validation"):
            raise ValueError("unexpected internal numerical binary")
        if j["role"] == "reference":
            if j["argv"][0] not in x["source_sha256"] or not Path(
                j["cwd"]
            ).resolve().is_relative_to(out):
                raise ValueError("reference binary or output directory is not bound")
        else:
            argv = j["argv"]
            if len(argv) < 4 or argv[1] != "--child" or argv[3] != str(out):
                raise ValueError("invalid matched command")
            if argv[2] == "matched-evaluate" and len(argv) == 8:
                fixture, backend, src, dest = argv[4:]
                if type(j.get("expected_frames")) != int or j["expected_frames"] < 1:
                    raise ValueError("register snapshot stream frame count")
                if not (
                    out / "inputs" / (fixture + ".json")
                ).is_file() or backend not in (
                    "zero",
                    "ewald",
                    "ewald-10",
                    "ewald-12",
                    "spme",
                ):
                    raise ValueError("unregistered matched evaluation")
                if not Path(src).resolve().is_relative_to(out) or not Path(
                    dest
                ).resolve().is_relative_to(out):
                    raise ValueError("snapshot path outside bounded campaign")
            elif argv[2] == "matched-trajectory" and len(argv) == 6:
                if argv[4] not in {
                    r["id"] for r in load(out / "inputs/roster.json")
                } or not Path(argv[5]).resolve().is_relative_to(out):
                    raise ValueError("unregistered matched trajectory")
            else:
                raise ValueError("unsupported matched action")
        result[name] = j
    if sum(j["role"] == "reference" for j in result.values()) != 24:
        raise ValueError("expected 20 reference trajectories and 4 preflights")
    return result


def remaining(clock):
    # Fail closed on reboot, backward monotonic discontinuity or an expired wall deadline.
    if time.monotonic() < clock["start_monotonic"]:
        raise Budget("host clock discontinuity")
    return min(
        clock["deadline_utc"] - time.time(),
        clock["deadline_monotonic"] - time.monotonic(),
    )


def start(out):
    with exclusive(out):
        _, x = verify(out)
        jobs(out, x)
        if platform.system() != "Darwin":
            raise ValueError("Darwin wait4 RSS/host-time adapter required")
        if raw_bytes(out) >= LIMITS["artifact_bytes"]:
            raise Budget("preparation artifacts exceed cap")
        # Verify actual host monitor capability before dispatching any physical job.
        pre = subprocess.run(
            [
                "/usr/bin/time",
                "-l",
                "-o",
                str(out / "host-preflight.time"),
                "/usr/bin/true",
            ],
            check=False,
        )
        if (
            pre.returncode
            or "maximum resident set size"
            not in (out / "host-preflight.time").read_text()
        ):
            raise ValueError("host-time preflight unavailable")
        if (out / "clock.json").exists():
            raise ValueError("campaign epoch already started")
        absolute = LEDGER / "absolute-clock.json"
        if not absolute.exists():
            now = time.monotonic()
            wall = time.time()
            dump(
                absolute,
                {
                    "started_utc": utc(),
                    "start_monotonic": now,
                    "deadline_utc": wall + 7200,
                    "deadline_monotonic": now + 7200,
                },
            )
        clock = load(absolute)
        if remaining(clock) <= 0:
            raise Budget("shared campaign deadline expired")
        pins = {
            "internal": digest(out / "protocol.json"),
            "external": digest(out / "external-epoch.json"),
        }
        dump(LEDGER / f"epoch-{time.time_ns()}.json", dict(out=str(out), **pins))
        dump(out / "clock.json", dict(clock, **pins))


def monitor(out, folder):
    """Runs INSIDE /usr/bin/time. wait4 measures the direct numerical child alone."""
    launch = load(folder / "launch.json")
    cmd = launch["argv"]
    role = launch["role"]
    lib = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
    buf = ctypes.create_string_buffer(96)
    cap = LIMITS["rust_rss_bytes"] if role == "rust" else LIMITS["reference_rss_bytes"]

    def caps():
        if launch["name"] == "deterministic":
            resource.setrlimit(resource.RLIMIT_CPU, (60, 61))
        resource.setrlimit(
            resource.RLIMIT_FSIZE, (LIMITS["artifact_bytes"], LIMITS["artifact_bytes"])
        )

    env = dict(os.environ, **launch["env"])
    env["B07_PROTOCOL_SHA256"] = launch["protocol_sha256"]
    env["B07_WORK_PATH"] = str(folder / "work.json")
    env["B07_HOST_PYTHON"] = sys.executable
    started = utc()
    begin = time.monotonic()
    with (
        (folder / "stdout.log").open("xb") as so,
        (folder / "stderr.log").open("xb") as se,
    ):
        child = subprocess.Popen(
            cmd,
            cwd=launch["cwd"],
            env=env,
            stdout=so,
            stderr=se,
            start_new_session=True,
            preexec_fn=caps,  # noqa: PLW1509 -- this dedicated monitor is single-threaded
        )
        dump(
            folder / "started.json",
            {"pid": child.pid, "utc": started, "monotonic": begin, "role": role},
        )
        last_scan = 0.0
        samples = 0
        sampled = 0
        stop = None
        while True:
            pid, status, usage = os.wait4(child.pid, os.WNOHANG)
            if pid:
                break
            now = time.monotonic()
            if lib.proc_pidinfo(child.pid, 4, 0, buf, 96) == 96:
                rss = int.from_bytes(buf.raw[8:16], "little")
                sampled = max(sampled, rss)
                samples += 1
                if rss > cap:
                    stop = "RSS cap"
            if (
                now >= launch["deadline_monotonic"]
                or time.time() >= launch["deadline_utc"]
            ):
                stop = "wall cap"
            if now - last_scan >= 0.5:
                if raw_bytes(out) >= LIMITS["artifact_bytes"]:
                    stop = "artifact cap"
                last_scan = now
            if stop:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            time.sleep(0.01)
        code = os.waitstatus_to_exitcode(status)
        child.returncode = code
    elapsed = time.monotonic() - begin
    if usage.ru_maxrss > cap:
        stop = "wait4 peak RSS exceeded"
    if launch["name"] == "deterministic" and usage.ru_utime + usage.ru_stime > 60:
        stop = "deterministic CPU cap"
    if raw_bytes(out) > LIMITS["artifact_bytes"]:
        stop = "artifact cap"
    if (
        time.monotonic() > launch["deadline_monotonic"]
        or time.time() > launch["deadline_utc"]
    ):
        stop = "wall cap"
    dump(
        folder / "child.json",
        {
            "started_utc": started,
            "ended_utc": utc(),
            "elapsed_seconds": elapsed,
            "exit": code,
            "peak_rss_bytes": usage.ru_maxrss,
            "user_cpu_seconds": usage.ru_utime,
            "system_cpu_seconds": usage.ru_stime,
            "rss_samples": samples,
            "sampled_max_rss_bytes": sampled,
            "budget_stop": stop,
        },
    )
    return 2 if stop else code if code >= 0 else 1


def external_spent(out, strict=True):
    total = 0.0
    count = 0
    for folder in (p for epoch in epochs(out) for p in (epoch / "jobs").iterdir()):
        if not (folder / "launch.json").exists():
            continue
        launch = load(folder / "launch.json")
        if launch["role"] != "reference":
            continue
        if not (folder / "child.json").exists():
            if strict:
                raise Budget("unfinished reference job; elapsed work unknown")
            count += int((folder / "started.json").exists())
            continue
        total += load(folder / "child.json")["elapsed_seconds"]
        count += 1
    return total, count


def admit_work(out, name):
    if name not in ("ou-draws", "ensemble"):
        return
    gas = 0
    projected = 0
    for epoch in epochs(out):
        folder = epoch / "jobs" / name
        if not folder.exists():
            continue
        if not (folder / "child.json").exists():
            raise Budget("previous interrupted work unknown")
        child = load(folder / "child.json")
        if name == "ensemble":
            if child["exit"] not in (0, 3, 4) or not (epoch / "gas-progress").is_file():
                raise Budget("previous gas suffix unknown")
            gas += int((epoch / "gas-progress").read_text())
        else:
            counts = [
                int(v)
                for v in re.findall(
                    r"conditional seed=\d+ draws=(\d+)",
                    (folder / "stdout.log").read_text(),
                )
            ]
            if counts != [8192] * 8:
                raise Budget("previous projected-draw suffix unknown")
            projected += sum(counts)
    if (
        gas + 316800000 > LIMITS["molecule_steps"]
        or 131072 + 4 + projected + 65536 > LIMITS["projected_draws"]
    ):
        raise Budget("remaining work below complete frozen upper bound")


def run_job(out, name, argv=None, role=None):
    out = Path(out).resolve()
    with exclusive(out):
        _, external = verify(out)
        roster = jobs(out, external)
        j = roster[name]
        admit_work(out, name)
        if (argv is not None and list(map(str, argv)) != j["argv"]) or (
            role is not None and role != j["role"]
        ):
            raise ValueError("unregistered invocation")
        clock = load(out / "clock.json")
        absolute = load(LEDGER / "absolute-clock.json")
        if any(clock[k] != v for k, v in absolute.items()):
            raise ValueError("absolute clock changed")
        sec = remaining(clock)
        spent, used = external_spent(out, strict=j["role"] == "reference")
        if j["role"] == "reference":
            if used >= 24:
                raise Budget("reference jobs exhausted")
            sec = min(sec, 3600 - spent)
        if sec <= 0 or raw_bytes(out) >= LIMITS["artifact_bytes"]:
            raise Budget("combined or external budget exhausted")
        folder = out / "jobs" / name
        folder.mkdir()  # One execution per frozen action; no budget reset or implicit retry.
        input_digest = None
        if "matched-evaluate" in j["argv"]:
            src = Path(j["argv"][-2])
            input_digest = digest(src)
            with src.open("rb") as stream:
                if sum(1 for _ in stream) != j["expected_frames"]:
                    raise ValueError("unexpected actual snapshot count")
            frozen = external["snapshot_sha256"].get(str(src))
            if frozen is None:
                # The adapter binds derived q/v only after parsing an actual native dump.
                binding = load(Path(str(src) + ".origin.json"))
                native = Path(binding["raw_path"])
                producer = out / "jobs" / binding["producer_job"]
                if (
                    binding["sha256"] != input_digest
                    or not native.resolve().is_relative_to(out)
                    or digest(native) != binding["raw_sha256"]
                    or binding["producer_receipt_sha256"]
                    != digest(producer / "child.json")
                    or load(producer / "child.json")["exit"] != 0
                    or binding["producer_job"] not in external["jobs"]
                    or external["jobs"][binding["producer_job"]]["role"] != "reference"
                ):
                    raise ValueError("invalid native snapshot origin")
            elif frozen != input_digest:
                raise ValueError("matched snapshot identity changed")
        launch = dict(
            j,
            name=name,
            input_sha256=input_digest,
            protocol_sha256=digest(out / "protocol.json"),
            binary_sha256=digest(j["argv"][0]),
            deadline_monotonic=min(clock["deadline_monotonic"], time.monotonic() + sec),
            deadline_utc=min(clock["deadline_utc"], time.time() + sec),
        )
        dump(folder / "launch.json", launch)
        with (
            (folder / "wrapper.stdout").open("xb") as so,
            (folder / "wrapper.stderr").open("xb") as se,
        ):
            wrapper = [
                "/usr/bin/time",
                "-l",
                "-o",
                str(folder / "host.time"),
                sys.executable,
                __file__,
                "_monitor",
                str(out),
                str(folder),
            ]
            receipt = subprocess.run(wrapper, stdout=so, stderr=se, check=False)
        dump(folder / "wrapper.json", {"argv": wrapper, "exit": receipt.returncode})
        child = load(folder / "child.json")
        if "maximum resident set size" not in (folder / "host.time").read_text():
            raise ValueError("missing independent host-time receipt")
        status = (
            2 if child["budget_stop"] else child["exit"] if child["exit"] >= 0 else 1
        )
        if status != receipt.returncode:
            raise ValueError("child and wrapper exit disagreement")
        return status


def gate(name, val, target, limits):
    errors = [abs(v - t) for v, t in zip(val, target)]
    if (
        any(t < 0 for t in limits)
        or not val
        or len(val) != len(target)
        or len(val) != len(limits)
        or not all(
            math.isfinite(x) for seq in (val, target, limits, errors) for x in seq
        )
    ):
        raise ValueError("invalid observed gate " + name)
    return {
        "name": name,
        "observed": val,
        "expected": target,
        "absolute_error": errors,
        "tolerance": limits,
        "status": "Passed" if all(e <= t for e, t in zip(errors, limits)) else "Failed",
    }


def checked_gates(rows):
    if not rows or len({r["name"] for r in rows}) != len(rows):
        raise ValueError("empty/duplicate gate collection")
    for r in rows:
        keys = {"name", "observed", "expected", "absolute_error", "tolerance", "status"}
        if set(r) != keys or r["status"] not in (
            "Passed",
            "Failed",
            "Inconclusive",
            "Unsupported",
            "Pending",
            "BudgetExceeded",
        ):
            raise ValueError("invalid gate record")
        if not r["observed"] and r["status"] in (
            "Unsupported",
            "Pending",
            "Inconclusive",
            "BudgetExceeded",
        ):
            if r["expected"] or r["absolute_error"] or r["tolerance"]:
                raise ValueError("partial gate shape")
            continue
        rebuilt = gate(r["name"], r["observed"], r["expected"], r["tolerance"])
        if r["absolute_error"] != rebuilt["absolute_error"] or (
            r["status"] == "Passed" and rebuilt["status"] != "Passed"
        ):
            raise ValueError("gate classified inconsistent with measurements")
    return rows


def report(out):
    p, x = verify(out)
    roster = jobs(out, x)
    gates = []
    for case in load(MESH)["cases"]:
        fine = case["rows"][-1]
        if fine["mesh"] != [32] * 3:
            raise ValueError("historical mesh roster")
        gates.append(
            gate(
                "historical_spme_" + case["fixture"],
                [fine["max_component_error"]],
                [0.0],
                [2e-4],
            )
        )
    # Native solver exit status is separate from external physical comparisons.
    external_gates = out / "external-gates.json"
    if external_gates.is_file():
        gates.extend(checked_gates(load(external_gates)))
    else:
        gates.append(
            {
                "name": "required_external_physical_comparisons",
                "observed": [],
                "expected": [],
                "absolute_error": [],
                "tolerance": [],
                "status": "Pending",
            }
        )
    records = []
    unknown = []
    projected = 0
    gas = 0
    for name, j in roster.items():
        folder = out / "jobs" / name
        if not (folder / "child.json").exists():
            if (folder / "launch.json").exists():
                unknown.append(name)
            gates.append(
                {
                    "name": name,
                    "observed": [],
                    "expected": [],
                    "absolute_error": [],
                    "tolerance": [],
                    "status": "Pending",
                }
            )
            continue
        child = load(folder / "child.json")
        launch = load(folder / "launch.json")
        actual_draws = 0
        actual_steps = 0
        if name == "ou-draws":
            counts = [
                int(n)
                for n in re.findall(
                    r"conditional seed=\d+ draws=(\d+)",
                    (folder / "stdout.log").read_text(),
                )
            ]
            actual_draws = sum(counts)
            projected += actual_draws
            if counts != [8192] * 8:
                unknown.append(name)
        if (folder / "work.json").exists():
            measured = load(folder / "work.json")
            if (
                set(measured) != {"system_steps", "molecule_steps", "complete"}
                or any(
                    type(measured[k]) != int or measured[k] < 0
                    for k in ("system_steps", "molecule_steps")
                )
                or type(measured["complete"]) != bool
            ):
                raise ValueError("malformed work sidecar")
            actual_steps = measured["molecule_steps"]
            if not measured["complete"]:
                unknown.append(name)
        elif name.startswith(("cost-old-", "cost-new-")):
            if child["exit"] == 0:
                data = load(folder / "stdout.log")
                actual_steps = data["steps"] * (
                    32 if data["fixture"] == "gas-water" else 1
                )
            else:
                unknown.append(name)
        elif (
            j["role"] == "reference"
            or "matched-trajectory" in j["argv"]
            or name in ("deterministic", "costs")
        ):
            unknown.append(name)
        if name == "ensemble":
            actual_steps = (
                int((out / "gas-progress").read_text())
                if (out / "gas-progress").exists()
                else 0
            )
            gas += actual_steps
            if child["exit"] != 0 and child["exit"] != 3 and child["exit"] != 4:
                unknown.append(name)
        if child["budget_stop"]:
            status = "BudgetExceeded"
        elif child["exit"] == 0:
            status = "Passed"
        elif child["exit"] == 2:
            status = "Unsupported"
        elif child["exit"] == 3:
            status = "Failed"
        elif child["exit"] == 4:
            status = "Inconclusive"
        else:
            status = "Pending"
        if (
            not (folder / "wrapper.json").exists()
            or not (folder / "host.time").exists()
        ):
            status = "Pending"
        if (
            name in ("deterministic", "ou", "ensemble")
            and (out / (name + "-gates.json")).exists()
        ):
            gates.extend(checked_gates(load(out / (name + "-gates.json"))))
        gates.append(
            {
                "name": "execution_" + name,
                "observed": [child["exit"]],
                "expected": [0.0],
                "absolute_error": [abs(child["exit"])],
                "tolerance": [0.0],
                "status": status,
            }
        )
        records.append(
            {
                "source_head": load(out / "inputs/settings.json")["source_head"],
                "compiler": load(out / "inputs/settings.json")["compiler"],
                "build_profile": "release",
                "features": [],
                "command": launch["argv"],
                "protocol_sha256": launch["protocol_sha256"],
                "binary_sha256": launch["binary_sha256"],
                "external_source_revision": x.get("external_source_revision")
                if j["role"] == "reference"
                else None,
                "external_packages": x.get("external_packages", [])
                if j["role"] == "reference"
                else [],
                "threads": 1,
                "started_utc": child["started_utc"],
                "ended_utc": child["ended_utc"],
                "elapsed_seconds": child["elapsed_seconds"],
                "peak_rss_bytes": child["peak_rss_bytes"],
                "molecule_steps": actual_steps,
                "projected_draws": actual_draws,
                "reference_jobs": int(j["role"] == "reference"),
            }
        )
    if not records:
        raise ValueError("no new actual execution receipts")
    # Typed provenance always describes one actual executable. The index records combined consumption separately.
    dest = out / "reports"
    dest.mkdir(exist_ok=True)
    run = dest / str(time.time_ns())
    run.mkdir()
    costs = out / "costs.json"
    if costs.exists():
        summary = []
        for row in load(costs)["steps"]:
            fixture = row["fixture"]
            by_version = {}
            for version in ["old", "new"]:
                log = out / "jobs" / f"cost-{version}-{fixture}" / "stdout.log"
                if log.is_file():
                    data = load(log)
                    by_version[version] = {
                        "steps": data["steps"],
                        "elapsed_ns": data["elapsed_ns"],
                        "ns_per_step": data["elapsed_ns"] / data["steps"],
                        "sha256": digest(log),
                    }
            summary.append(
                {
                    "fixture": fixture,
                    "common_legacy_step": by_version,
                    "new_reporting": row,
                }
            )
        dump(run / "cost-comparison.json", summary)
    files = inventory(
        f for f in out.rglob("*") if f.is_file() and dest not in f.parents
    )
    worst = next(
        (
            s
            for s in [
                "Failed",
                "BudgetExceeded",
                "Unsupported",
                "Inconclusive",
                "Pending",
            ]
            if any(g["status"] == s for g in gates)
        ),
        "Passed",
    )
    dump(
        run / "physical-report.json",
        {
            "protocol": p,
            "provenance": records[-1],
            "gates": gates,
            "raw_files": files,
            "status": worst,
        },
    )
    dump(
        run / "work.json",
        {
            "validation": "Diagnostic",
            "all_execution_provenance": records,
            "unknown_or_partial_work": unknown,
            "combined_elapsed_seconds": time.monotonic()
            - load(out / "clock.json")["start_monotonic"],
            "external_consumption_known_prefix": external_spent(out, False),
            "gas_molecule_steps_observed": gas,
            "new_conditional_draws_observed": projected,
            "prior_conditional_draws": load(out / "inputs/settings.json")[
                "prior_projected_draws"
            ],
            "prior_zero_amplitude_projector_calls": 4,
            "status": worst,
        },
    )
    return {
        "Passed": 0,
        "Failed": 3,
        "BudgetExceeded": 2,
        "Unsupported": 2,
        "Inconclusive": 4,
        "Pending": 4,
    }[worst]


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_monitor":
        return monitor(Path(sys.argv[2]), Path(sys.argv[3]))
    ap = argparse.ArgumentParser()
    ap.add_argument("phase")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    out = a.out.resolve()
    if a.phase == "prepare":
        prepare(out)
        return 0
    if a.phase == "verify":
        verify(out, (out / "external-epoch.json").exists())
        return 0
    if a.phase == "start":
        start(out)
        return 0
    if a.phase == "report":
        return report(out)
    if a.phase == "ou":
        status = run_job(out, "ou-draws")
        if status not in (0, 101):
            return status
        # A libtest assertion can fail after emitting the full data. Classify its raw values.
        measured = run_job(out, "ou")
        return 3 if status == 101 and measured == 0 else measured
    return run_job(out, a.phase)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Budget as e:
        print(e, file=sys.stderr)
        sys.exit(2)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as e:
        print(e, file=sys.stderr)
        sys.exit(1)
