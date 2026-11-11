"""Shared B07 campaign ledger and child monitor. All science runs are serial."""

import argparse
import contextlib
import ctypes
import datetime
import fcntl
import hashlib
import json
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
    loads,
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
    settings = load(out / "inputs/settings.json")
    if Path(sys.executable).resolve() != Path(settings["host_python"]):
        raise ValueError("unexpected host reporter/interpreter")
    if settings["host_identity"] != host_identity():
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
    expected_external(x)
    native = x["native_build"]
    if (
        not x["external_source_revision"]
        or not x["external_packages"]
        or not native.get("compiler")
        or not native.get("build_profile")
        or not isinstance(native.get("features"), list)
    ):
        raise ValueError("incomplete native build identity")
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
    env["B07_HOST_PYTHON"] = load(out / "inputs/settings.json")["host_python"]
    started = utc()
    begin = time.monotonic()
    with (
        (folder / "stdout.log").open("xb") as so,
        (folder / "stderr.log").open("xb") as se,
    ):
        # Handlers request cleanup without interrupting the fork/started-receipt ownership window.
        cancelled = []

        def on_signal(sig, _frame):
            cancelled.append(signal.Signals(sig).name)

        saved = {
            sig: signal.signal(sig, on_signal)
            for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
        }
        child = None
        reaped = False
        samples = 0
        sampled = 0
        stop = None
        interrupted = None
        try:
            child = subprocess.Popen(
                cmd,
                cwd=launch["cwd"],
                env=env,
                stdout=so,
                stderr=se,
                start_new_session=True,
                preexec_fn=caps,  # noqa: PLW1509 -- dedicated single-threaded monitor
            )
            dump(
                folder / "started.json",
                {"pid": child.pid, "utc": started, "monotonic": begin, "role": role},
            )
            last_scan = 0.0
            while True:
                if cancelled:
                    interrupted = "monitor cancelled: " + cancelled[0]
                    break
                pid, status, usage = os.wait4(child.pid, os.WNOHANG)
                if pid:
                    reaped = True
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
                    break
                time.sleep(0.01)
        except BaseException as exc:
            if child is None:
                raise
            interrupted = f"monitor error: {type(exc).__name__}: {exc}"
        finally:
            if child is not None and not reaped:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(child.pid, signal.SIGKILL)
                _, status, usage = os.wait4(child.pid, 0)
                reaped = True
            for sig, handler in saved.items():
                signal.signal(sig, handler)
        code = os.waitstatus_to_exitcode(status)
        child.returncode = code
    elapsed = time.monotonic() - begin
    if usage.ru_maxrss > cap:
        stop = "wait4 peak RSS exceeded"
    if launch["name"] == "deterministic" and usage.ru_utime + usage.ru_stime > 60:
        stop = "deterministic CPU cap"
    if not interrupted and raw_bytes(out) > LIMITS["artifact_bytes"]:
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
            "interrupted": interrupted,
            "reaped": reaped,
        },
    )
    return 2 if stop else 1 if interrupted else code if code >= 0 else 1


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
        if not folder.exists() or (folder / "denied.json").exists():
            continue
        if not (folder / "child.json").exists():
            raise Budget("previous interrupted work unknown")
        child = load(folder / "child.json")
        if name == "ensemble":
            if (
                child.get("budget_stop")
                or child.get("interrupted")
                or child["exit"] not in (0, 3, 4)
                or not (epoch / "gas-progress").is_file()
            ):
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


def ensure_no_live_launch(out):
    for epoch in epochs(out):
        for launch in (epoch / "jobs").glob("*/launch.json"):
            folder = launch.parent
            if (folder / "child.json").is_file() and load(folder / "child.json").get(
                "reaped", True
            ):
                continue
            started = folder / "started.json"
            if not started.is_file():
                raise Budget(f"unresolved launch, unknown PID: {folder}")
            pid = load(started).get("pid")
            if type(pid) != int or pid <= 1:
                raise Budget(f"unresolved PID: {folder}")
            try:
                os.killpg(pid, 0)
            except ProcessLookupError:
                continue
            except PermissionError:
                pass
            raise Budget(f"unresolved live process group {pid}: {folder}")


def run_wrapper(argv, folder, so, se):
    """Forward cancellation to the time/monitor group and retain its independent exit."""
    signals = []

    def cancel(sig, _):
        signals.append(sig)

    saved = {
        sig: signal.signal(sig, cancel)
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)
    }
    proc = None
    try:
        proc = subprocess.Popen(argv, stdout=so, stderr=se, start_new_session=True)
        while proc.poll() is None and not signals:
            time.sleep(0.02)
    finally:
        if proc is not None and proc.poll() is None:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=2)
            except subprocess.TimeoutExpired:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
        if proc is not None:
            host = folder / "host.time"
            dump(
                folder / "wrapper.json",
                {
                    "argv": argv,
                    "exit": proc.returncode,
                    "cancellation": [signal.Signals(s).name for s in signals],
                    "host_time_sha256": digest(host) if host.exists() else None,
                },
            )
            # A killed monitor cannot wait4 its child. Terminate the live group for this launch;
            # a missing child receipt remains unknown and further admission rechecks the group.
            if not (folder / "child.json").is_file():
                started = folder / "started.json"
                if started.is_file():
                    pid = load(started)["pid"]
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(pid, signal.SIGKILL)
                dump(
                    folder / "monitor-loss.json",
                    {"utc": utc(), "wrapper_exit": proc.returncode, "wait4": "unknown"},
                )
        for sig, handler in saved.items():
            signal.signal(sig, handler)
    return proc.returncode


def run_job(out, name, argv=None, role=None):
    out = Path(out).resolve()
    with exclusive(out):
        _, external = verify(out)
        roster = jobs(out, external)
        j = roster[name]
        if (argv is not None and list(map(str, argv)) != j["argv"]) or (
            role is not None and role != j["role"]
        ):
            raise ValueError("unregistered invocation")
        folder = out / "jobs" / name
        if folder.exists():
            raise ValueError(f"action already recorded: {name}")
        try:
            clock = load(out / "clock.json")
            absolute = load(LEDGER / "absolute-clock.json")
            if any(clock[k] != v for k, v in absolute.items()):
                raise ValueError("absolute clock changed")
            sec = remaining(clock)
            if sec <= 0 or raw_bytes(out) >= LIMITS["artifact_bytes"]:
                raise Budget("combined budget exhausted")
            ensure_no_live_launch(out)
            admit_work(out, name)
            if j["role"] == "reference":
                spent, used = external_spent(out)
                if used >= 24:
                    raise Budget("reference jobs exhausted")
                sec = min(sec, 3600 - spent)
                if sec <= 0:
                    raise Budget("external wall budget exhausted")
        except Budget as exc:
            folder.mkdir()
            dump(
                folder / "denied.json",
                {
                    "name": name,
                    "utc": utc(),
                    "status": "BudgetExceeded",
                    "reason": str(exc),
                    "launched": False,
                },
            )
            raise
        folder.mkdir()
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
                load(out / "inputs/settings.json")["host_python"],
                __file__,
                "_monitor",
                str(out),
                str(folder),
            ]
            run_wrapper(wrapper, folder, so, se)
        if not (folder / "child.json").is_file():
            raise ValueError("monitor lost; direct child receipt unavailable")
        child = load(folder / "child.json")
        status = child_status(child)
        if not reconciled_wrapper(folder, child):
            raise ValueError("inconsistent child and host-time wrapper receipts")
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


def unmeasured(name, status="Pending"):
    return {
        "name": name,
        "observed": [],
        "expected": [],
        "absolute_error": [],
        "tolerance": [],
        "status": status,
    }


def expected_external(epoch):
    families = epoch.get("expected_gate_families")
    if not isinstance(families, dict) or not families:
        raise ValueError("missing frozen external gate families")
    flat = []
    for family, names in families.items():
        if (
            not isinstance(family, str)
            or not family
            or not isinstance(names, list)
            or not names
            or any(not isinstance(n, str) or not n for n in names)
        ):
            raise ValueError("invalid external gate family")
        flat.extend(names)
    if len(set(flat)) != len(flat):
        raise ValueError("duplicate external gate expectation")
    return flat


def child_status(child):
    if child["budget_stop"]:
        return 2
    return child["exit"] if child["exit"] >= 0 and not child.get("interrupted") else 1


def reconciled_wrapper(folder, child):
    if not (folder / "wrapper.json").is_file() or not (folder / "host.time").is_file():
        return False
    wrapper = load(folder / "wrapper.json")
    settings = load(folder.parent.parent / "inputs/settings.json")
    argv = [
        "/usr/bin/time",
        "-l",
        "-o",
        str(folder / "host.time"),
        settings["host_python"],
        __file__,
        "_monitor",
        str(folder.parent.parent),
        str(folder),
    ]
    if (
        wrapper.get("argv") != argv
        or wrapper.get("host_time_sha256") != digest(folder / "host.time")
        or wrapper.get("exit") != child_status(child)
        or wrapper.get("cancellation")
    ):
        return False
    txt = (folder / "host.time").read_text()
    timing = re.search(r"(?m)^\s*([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys\s*$", txt)
    peak = re.search(r"(?m)^\s*(\d+)\s+maximum resident set size\s*$", txt)
    if not timing or not peak or int(peak[1]) <= 0:
        return False
    try:
        measured = [float(v) for v in timing.groups()]
    except ValueError:
        return False
    return all(
        math.isfinite(v) and v >= 0
        for v in measured
        + [
            child["elapsed_seconds"],
            child["user_cpu_seconds"],
            child["system_cpu_seconds"],
            child["peak_rss_bytes"],
        ]
    )


def required_rows(rows, expected):
    names = {r["name"] for r in rows}
    return list(rows) + [unmeasured(n) for n in expected if n not in names]


def internal_names(phase):
    if phase == "ou":
        return [
            f"ou_{kind}_mean_{i}"
            for i in range(9)
            for kind in ("projected", "raw", "reaction")
        ] + [f"ou_cov_{i}_{j}" for i in range(9) for j in range(9)]
    if phase == "deterministic":
        names = ["vv_rotor_mixed_time_negative_control"]
        for f in (
            "rotor",
            "rotor-heavy",
            "water-principal",
            "water-principal-heavy",
            "water-nonprincipal",
            "water-nonprincipal-heavy",
        ):
            for n in (100, 200, 400, 800):
                names.extend(
                    [
                        f"{f}_{n}_finite_window_ledger",
                        f"{f}_{n}_euler_or_rotor_reaction",
                    ]
                )
            names.extend(
                f"{f}_{observable}_last_two_ratios"
                for observable in ("positions", "velocities", "interval", "endpoint")
            )
            if "nonprincipal" in f:
                names.append(f"{f}_rk4_refinement")
            if f.startswith("rotor"):
                names.extend(
                    [f"{f}_vv_maps_all_images", f"{f}_baoab_two_half_maps_all_images"]
                )
                for tilt in ("true", "false"):
                    names.extend(
                        f"{f}_forced_trap_{tilt}_{field}"
                        for field in ("pressure", "atom_0", "atom_1")
                    )
        return names
    return []


def physics_rows(out, external, defects=None, fragment_folder=None):
    gates = []
    if defects is None:
        defects = {}

    def reject(path, data, reason):
        sha = hashlib.sha256(data).hexdigest()
        marker = "incomplete_evidence_" + path.name
        detail = {"gate": marker, "reason": reason, "sha256": sha, "bytes": len(data)}
        if fragment_folder is not None:
            fragment_folder.mkdir(exist_ok=True)
            snapshot = fragment_folder / path.name
            with snapshot.open("xb") as f:
                f.write(data)
            detail["preserved_bytes"] = str(snapshot)
        defects[str(path)] = detail
        gates.append(unmeasured(marker, "Inconclusive"))

    def collect(path):
        # Read once: the digest and a report's preserved copy cover exactly the
        # bytes that failed, including pre-atomic writer output in an old epoch.
        try:
            data = path.read_bytes()
        except FileNotFoundError:
            data = None
        stage = path.with_name(path.name + ".partial")
        try:
            staged = stage.read_bytes()
        except FileNotFoundError:
            staged = None
        if staged is not None and staged != data:
            reject(stage, staged, "unpublished staged gate bytes")
        if data is None:
            return []
        try:
            decoded = loads(data)
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            reject(path, data, f"incomplete or invalid JSON: {exc}")
            return []
        except ValueError as exc:
            reject(path, data, f"invalid JSON data: {exc}")
            return []
        try:
            return checked_gates(decoded)
        except (ValueError, TypeError, KeyError) as exc:
            reject(path, data, f"invalid gate record: {exc}")
            return []

    for case in load(MESH)["cases"]:
        fine = case["rows"][-1]
        if fine["mesh"] != [32] * 3:
            raise ValueError("historical mesh roster")
        gates.append(
            gate(
                "historical_spme_" + case["fixture"],
                [fine["max_component_error"], fine["max_translation_change"]],
                [0.0, 0.0],
                [2e-4, 2e-4],
            )
        )
    fp = out / "external-gates.json"
    rows = collect(fp)
    gates.extend(required_rows(rows, expected_external(external)))
    for phase in ("deterministic", "ou"):
        path = out / (phase + "-gates.json")
        rows = collect(path)
        if rows:
            gates.extend(required_rows(rows, internal_names(phase)))
        else:
            gates.append(unmeasured(f"required_{phase}_physical_gates"))
    final = out / "ensemble-gates.json"
    final_rows = collect(final)
    if not final_rows:
        gates.append(unmeasured("required_ensemble_physical_gates"))
    latest = []
    rungs = []
    for species in ("gas-rotor", "gas-water"):
        most_recent = []
        for duration in (400, 800, 1600):
            path = out / f"{species}-{duration}-gates.json"
            observed = collect(path)
            if not observed:
                continue
            prefix = f"{species}_{duration}_"
            if any(not row["name"].startswith(prefix) for row in observed):
                raise ValueError("gas rung row from different trajectory set")
            names = {row["name"] for row in observed}
            expected = []
            pressure = [
                f"{timing}_{k}"
                for timing in ("endpoint", "interval")
                for k in ("xx", "yy", "zz", "xy", "xz", "yz", "scalar")
            ]
            energy = [
                f"{k}_{m}"
                for k in ("com", "rotation")
                for m in ("mean", "variance", "cdf10", "cdf50", "cdf90")
            ]
            for metric in pressure + energy + [f"orientation_{i}" for i in range(6)]:
                expected.extend(
                    f"{prefix}{h}_{metric}_precision"
                    for h in ("0.04", "0.02", "0.01", "0.005")
                )
                if f"{prefix}{metric}_unresolved_fit" not in names:
                    expected.append(f"{prefix}{metric}_intercept")
                    if metric in pressure:
                        expected.extend(
                            [
                                f"{prefix}{metric}_finest",
                                f"{prefix}{metric}_drop_coarse",
                            ]
                        )
            most_recent = required_rows(observed, expected)
            rungs.append(
                {
                    "species": species,
                    "measurement_duration": duration,
                    "sha256": digest(path),
                    "gates": most_recent,
                }
            )
        if not most_recent:
            gates.append(unmeasured(f"required_{species}_rung"))
        latest.extend(most_recent)
    if final_rows:
        by_name = {row["name"]: row for row in final_rows}
        for row in latest:
            other = by_name.get(row["name"])
            if other is None:
                gates.append(unmeasured("missing_final_" + row["name"]))
                by_name[row["name"]] = row
            elif other != row:
                raise ValueError("final gas record differs from completed rung")
        gates.extend(by_name.values())
    else:
        gates.extend(latest)
    return gates, rungs


def measured_work(out, name, job, child):
    folder = out / "jobs" / name
    count = 0
    draws = 0
    launched = (folder / "launch.json").exists()
    incomplete = launched and child is None
    if (folder / "work.json").exists():
        w = load(folder / "work.json")
        if (
            set(w) != {"system_steps", "molecule_steps", "complete"}
            or type(w["complete"]) != bool
            or any(
                type(w[k]) != int or w[k] < 0
                for k in ("system_steps", "molecule_steps")
            )
        ):
            raise ValueError("malformed actual work")
        count = w["molecule_steps"]
        incomplete |= not w["complete"]
    elif launched and name.startswith(("cost-old-", "cost-new-")):
        interrupted = child is None or child["exit"] != 0
        # A driver prints the measured count after its loop. Preserve an emitted
        # final record even if termination happens before a clean process exit.
        try:
            data = load(folder / "stdout.log")
        except (FileNotFoundError, json.JSONDecodeError, UnicodeDecodeError):
            if interrupted:
                return count, draws, True
            raise
        incomplete |= interrupted
        steps = data["steps"]
        if (
            type(steps) != int
            or steps < 0
            or data["fixture"] not in ("rotor", "gas-water")
            or not name.endswith("-" + data["fixture"])
        ):
            raise ValueError("invalid cost work")
        count = steps * (32 if data["fixture"] == "gas-water" else 1)
    elif launched and (
        job.get("role") == "reference"
        or "matched-trajectory" in job.get("argv", [])
        or name in ("deterministic", "costs")
    ):
        incomplete = True
    if name == "ensemble":
        progress = out / "gas-progress"
        if progress.is_file():
            count = int(progress.read_text())
            if count < 0:
                raise ValueError("negative completed gas work")
        incomplete |= launched and (
            not progress.exists()
            or not child
            or child["exit"] not in (0, 3, 4)
            or bool(child.get("budget_stop"))
        )
    if name == "ou-draws":
        path = folder / "stdout.log"
        pairs = (
            [
                (int(i), int(n))
                for i, n in re.findall(
                    r"conditional seed=(\d+) draws=(\d+)", path.read_text()
                )
            ]
            if path.exists()
            else []
        )
        if len(pairs) > 8 or pairs != [(i, 8192) for i in range(len(pairs))]:
            raise ValueError("invalid conditional seed counts")
        draws = sum(n for _, n in pairs)
        incomplete |= launched and len(pairs) != 8
    return count, draws, incomplete


def consumption(epoch):
    counts = {
        "gas_molecule_steps": 0,
        "conditional_draws": 0,
        "all_molecule_steps": 0,
        "unknown_suffixes": [],
    }
    for folder in sorted((epoch / "jobs").iterdir()):
        if not folder.is_dir():
            continue
        launch = folder / "launch.json"
        child = folder / "child.json"
        n, d, unknown = measured_work(
            epoch,
            folder.name,
            load(launch) if launch.exists() else {},
            load(child) if child.exists() else None,
        )
        counts["all_molecule_steps"] += n
        counts["conditional_draws"] += d
        if folder.name == "ensemble":
            counts["gas_molecule_steps"] += n
        if unknown:
            counts["unknown_suffixes"].append(str(folder))
    return counts


def executable_build(settings, external, job, launch):
    path = job["argv"][0]
    if launch["argv"] != job["argv"] or launch["binary_sha256"] != digest(path):
        raise ValueError("execution differs from registered executable")
    if job["role"] == "reference":
        return dict(
            external["native_build"],
            source_head=external["external_source_revision"],
            binary_sha256=launch["binary_sha256"],
            external_source_revision=external["external_source_revision"],
            external_packages=external["external_packages"],
        )
    build = settings["executable_builds"][path]
    if build["binary_sha256"] != launch["binary_sha256"]:
        raise ValueError("binary differs from registered build")
    return dict(build, external_source_revision=None, external_packages=[])


def report(out):
    p, x = verify(out)
    roster = jobs(out, x)
    dest = out / "reports"
    dest.mkdir(exist_ok=True)
    run = dest / str(time.time_ns())
    run.mkdir()
    defects = {}
    gates, rungs = physics_rows(out, x, defects, run / "incomplete-gates")
    settings = load(out / "inputs/settings.json")
    records = []
    builds = {}
    receipt_checks = {}
    denials = {}
    for name, job in roster.items():
        folder = out / "jobs" / name
        if (folder / "denied.json").exists():
            denial = load(folder / "denied.json")
            if (
                denial["name"] != name
                or denial["status"] != "BudgetExceeded"
                or denial["launched"] is not False
            ):
                raise ValueError("malformed denial")
            denials[name] = denial
            gates.append(unmeasured("execution_" + name, "BudgetExceeded"))
            continue
        if not (folder / "child.json").is_file():
            gates.append(unmeasured("execution_" + name))
            continue
        child = load(folder / "child.json")
        launch = load(folder / "launch.json")
        confirmed = reconciled_wrapper(folder, child)
        receipt_checks[name] = {
            "matched_host_wrapper": confirmed,
            "child_exit": child["exit"],
            "budget_stop": child["budget_stop"],
            "interruption": child.get("interrupted"),
        }
        status = {
            0: "Passed",
            2: "BudgetExceeded" if child["budget_stop"] else "Unsupported",
            3: "Failed",
            4: "Inconclusive",
        }.get(child_status(child), "Pending")
        if not confirmed:
            status = "Pending"
        execution = gate("execution_" + name, [child["exit"]], [0.0], [0.0])
        execution["status"] = status
        gates.append(execution)
        n, d, _ = measured_work(out, name, job, child)
        build = executable_build(settings, x, job, launch)
        builds[name] = build
        records.append(
            {
                k: build[k]
                for k in (
                    "source_head",
                    "compiler",
                    "build_profile",
                    "features",
                    "external_source_revision",
                    "external_packages",
                )
            }
        )
        records[-1].update(
            command=launch["argv"],
            protocol_sha256=launch["protocol_sha256"],
            binary_sha256=launch["binary_sha256"],
            threads=1,
            started_utc=child["started_utc"],
            ended_utc=child["ended_utc"],
            elapsed_seconds=child["elapsed_seconds"],
            peak_rss_bytes=child["peak_rss_bytes"],
            molecule_steps=n,
            projected_draws=d,
            reference_jobs=int(job["role"] == "reference"),
        )
    # Typed provenance always describes one actual executable. The index records combined consumption separately.
    costs = out / "costs.json"
    if costs.exists():
        summary = []
        for row in load(costs)["steps"]:
            fixture = row["fixture"]
            by_version = {}
            for version in ["old", "new"]:
                log = out / "jobs" / f"cost-{version}-{fixture}" / "stdout.log"
                child_file = log.parent / "child.json"
                if (
                    child_file.is_file()
                    and child_status(load(child_file)) == 0
                    and log.is_file()
                ):
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
        run / "gate-ledger.json",
        {
            "status": worst,
            "gates": gates,
            "completed_gas_rungs": rungs,
            "incomplete_gate_files": defects,
            "raw_files": files,
        },
    )
    if records:
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
    # With no measured child there is no valid §12 ExecutionProvenance to serialize.
    by_epoch = {str(e): consumption(e) for e in epochs(out)}
    totals = {
        k: sum(v[k] for v in by_epoch.values())
        for k in ("gas_molecule_steps", "conditional_draws", "all_molecule_steps")
    }
    totals["unknown_suffixes"] = [
        item for v in by_epoch.values() for item in v["unknown_suffixes"]
    ]
    dump(
        run / "work.json",
        {
            "status": worst,
            "validation": "Diagnostic",
            "registered_builds": builds,
            "all_execution_provenance": records,
            "receipt_checks": receipt_checks,
            "budget_denials": denials,
            "current_epoch_work": by_epoch[str(out)],
            "work_by_epoch": by_epoch,
            "cumulative_work": totals,
            "combined_elapsed_seconds": time.monotonic()
            - load(out / "clock.json")["start_monotonic"],
            "external_consumption_known_prefix": external_spent(out, False),
            "prior_conditional_draws": settings["prior_projected_draws"],
            "prior_zero_amplitude_projector_calls": 4,
            "charged_conditional_calls_known_prefix": sum(
                settings["prior_projected_draws"]
            )
            + 4
            + totals["conditional_draws"],
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
