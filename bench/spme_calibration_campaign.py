#!/usr/bin/env python3
"""Sequential immutable B05 campaigns; standard library, no shell execution.

Builds release before the measurement lease. A wait4 monitor retains the actual
Rust process status independently of the /usr/bin/time wrapper's exit status.
"""
import argparse
import copy
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

FAMILIES = [("cubic", 2), ("ortho", 2), ("special", 2), ("cubic", 3), ("ortho", 3)]
THREADS = [1, 4]
ACTIVE_PHASES = {"verify", "sweep", "confirm", "profile"}
WALL_NS = 1_800_000_000_000
REPORT_NS = 120_000_000_000  # Bounded noncompute rendering; never charged to family work.
RSS_CAP = 536_870_912
REPORT_FILES = ("summary.json", "diagnostics.json", "ordinary-costs.csv",
                "errors.csv", "targets.csv", "cost-error.svg", "held-out-errors.csv")
ROOT = Path(__file__).resolve().parents[1]


def validate_threads(threads):
    if threads != THREADS:
        raise ValueError("frozen campaign requires --threads 1 4, in that order")


def digest(path):
    with Path(path).open("rb") as stream:
        value = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
        return value.hexdigest()


def source_inventory():
    paths = [ROOT / name for name in ("Cargo.toml", "Cargo.lock", "crates/modena-core/Cargo.toml", "crates/modena-core/examples/spme_convergence.rs", "bench/spme_calibration_campaign.py", "docs/scientific-hardening/B05-SPME-PROTOCOL.md")]
    paths += sorted((ROOT / "crates/modena-core/src").rglob("*.rs"))
    paths += sorted((ROOT / "crates/modena-core/examples/spme_calibration").rglob("*.rs"))
    return {str(path.relative_to(ROOT)): digest(path) for path in paths}


def artifact(path):
    return {"path": str(Path(path).resolve()), "sha256": digest(path)}


def checked_artifact(reference):
    """Compare imported bytes against the saved producer identity, never refresh it."""
    path = Path(reference["path"])
    if not path.is_file() or digest(path) != reference["sha256"]:
        raise ValueError("artifact changed since production: " + str(path))
    return path


def retain_output(path, directory, name="receipt"):
    """Called once at owned-child completion (or our own new-file creation)."""
    reference = artifact(path)
    write_new(Path(directory) / ("accepted-" + name + ".json"), reference)
    return reference


def write_new(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def parse_rss(stderr):
    matches = re.findall(r"^\s*([0-9]+)\s+maximum resident set size\s*$", stderr, re.M)
    return int(matches[-1]) if len(matches) == 1 and int(matches[0]) > 0 else None


def rss_accepted(value):
    return isinstance(value, int) and 0 < value <= RSS_CAP


def monitor(receipt, timeout, command):
    """Wait directly for the Rust child. Never substitute Cargo/time's status."""
    started = time.monotonic_ns()
    child = subprocess.Popen(command, start_new_session=True)
    write_new(str(receipt) + ".pid.json", {"pid": child.pid})
    expired = False
    while True:
        pid, status, usage = os.wait4(child.pid, os.WNOHANG)
        if pid:
            break
        if time.monotonic_ns() - started >= int(timeout * 1e9):
            expired = True
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            _, status, usage = os.wait4(child.pid, 0)
            break
        time.sleep(0.01)
    child.returncode = os.waitstatus_to_exitcode(status)
    write_new(receipt, {"actual_child_exit": child.returncode,
                        "wait_status": status, "monitor_exit": 0, "timed_out": expired,
                        "active_ns": time.monotonic_ns() - started,
                        "child_peak_rss_bytes": int(usage.ru_maxrss) if sys.platform == "darwin" else int(usage.ru_maxrss) * 1024,
                        "command": command})
    return 0


def execute(command, env, out, phase, remaining_ns):
    """One child at a time, distinct durable stdout/stderr and process receipt."""
    out = Path(out)
    out.mkdir()
    child_receipt = out / "child.json"
    argv = ["/usr/bin/time", "-l", sys.executable, str(Path(__file__).resolve()),
            "--_monitor", str(child_receipt), str(remaining_ns / 1e9), "--", *command]
    started = time.monotonic_ns()
    expired = False
    with (out / "stdout.log").open("xb") as stdout, (out / "stderr.log").open("xb") as stderr:
        wrapper = subprocess.Popen(argv, env=env, cwd=ROOT, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            wrapper.wait(timeout=remaining_ns / 1e9 + 10)
        except subprocess.TimeoutExpired:
            expired = True
            pidfile = Path(str(child_receipt) + ".pid.json")
            if pidfile.exists():
                try:
                    os.killpg(json.loads(pidfile.read_text())["pid"], signal.SIGKILL)
                except ProcessLookupError:
                    pass
            os.killpg(wrapper.pid, signal.SIGKILL)
            wrapper.wait(timeout=5)
        stdout.flush(); os.fsync(stdout.fileno())
        stderr.flush(); os.fsync(stderr.fileno())
    elapsed = time.monotonic_ns() - started
    child = json.loads(child_receipt.read_text()) if child_receipt.exists() else {}
    # A missing receipt is an unknown child exit, never a fabricated successful 0.
    row = {"phase": phase, "actual_child_exit": child.get("actual_child_exit"),
           "monitor_exit": child.get("monitor_exit"), "wrapper_exit": wrapper.returncode,
           "peak_rss_bytes": parse_rss((out / "stderr.log").read_text()),
           "active_ns": elapsed, "timed_out": expired or child.get("timed_out", False),
           "command": command, "stdout": str(out / "stdout.log"), "stderr": str(out / "stderr.log")}
    write_new(out / "memory.json", row)
    return row


def empty_usage():
    return {"candidate": {"accuracy_calls": 0, "candidate_real_pair_upper_bound": 0,
                          "reference": {"work": {"real_terms": 0, "atom_k_terms": 0}, "calls": 0, "elapsed_ns": 0},
                          "process_requested_rust_heap_bytes": 0, "peak_process_requested_rust_heap_bytes": 0},
            "timing_calls": 0, "profile_calls": 0}


def recover_usage(path, initial, epoch):
    """Admission is a reserved cumulative state; completion replaces it.

    A torn last line is retained raw and ignored. No in-flight evaluation is made
    up, and no reservation is added to the completed actual counters.
    """
    state = copy.deepcopy(initial)

    def visit(value):
        if not isinstance(value, dict):
            return
        for key in ("reserved", "usage"):
            usage = value.get(key)
            if isinstance(usage, dict):
                if "accuracy_calls" in usage:
                    state["candidate"] = copy.deepcopy(usage)
                elif "work" in usage and "calls" in usage:
                    state["candidate"]["reference"] = copy.deepcopy(usage)
        if "timing_calls" in value:
            state["timing_calls"] = value["timing_calls"]
        if "profile_calls" in value:
            state["profile_calls"] = value["profile_calls"]
        for key, child in value.items():
            if key not in ("reserved", "usage", "evaluation", "record", "rows", "raw"):
                visit(child)
        # Physical completed records carry the actual replacement usage.
        if isinstance(value.get("record"), dict):
            visit({"usage": value["record"].get("usage")})

    if Path(path).exists():
        with Path(path).open("rb") as stream:
            for sequence, line in enumerate(stream):
                try:
                    event = json.loads(line)
                except (json.JSONDecodeError, UnicodeDecodeError):
                    break
                if event.get("sequence") != sequence or event.get("manifest_sha256") != epoch:
                    raise ValueError("journal source/sequence mismatch")
                visit(event["data"])
    return state


def failed_receipt(prepared, out, phase, previous, row, active_ns, finalized=None):
    """One recovered view of the failed phase; originals remain raw artifacts.

    `previous` is the saved earlier phase, never the failed phase's own summary.
    `finalized` is the producer-pinned optional summary from this failed child.
    """
    prepared, out = Path(prepared), Path(out)
    out.mkdir(exist_ok=True)
    epoch = digest(prepared / "manifest.json")
    prior = json.loads(checked_artifact(previous).read_text()) if previous else None
    original = None
    summary_error = None
    if finalized:
        path = checked_artifact(finalized)
        try:
            original = json.loads(path.read_text())
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            summary_error = {"Io": "incomplete finalized summary: " + str(error)}
        if original:
            if original.get("prior") != previous:
                raise ValueError("finalized phase differs from saved predecessor")
            for reference in original.get("artifacts", {}).values():
                checked_artifact(reference)
            for reference in original.get("inputs", {}).values():
                checked_artifact(reference)
    usage = recover_usage(out / "events.jsonl", prior["usage"] if prior else empty_usage(), epoch)
    # A finalized outcome is the producer's last actual continuing usage; an
    # interrupted child without one retains the last durable reservation/event.
    if original:
        usage = copy.deepcopy(original["usage"])
    usage["candidate"]["reference"]["elapsed_ns"] = max(active_ns, usage["candidate"]["reference"]["elapsed_ns"])
    code = row["actual_child_exit"]
    category = {86: "requested Rust heap cap", 87: "primary allocator callback", 88: "System allocator returned null"}.get(code, "process/monitor/resource failure")
    process_failure = {"BudgetExceeded": f"{category}; actual child exit {code}; timed_out={row['timed_out']}"}
    artifacts = copy.deepcopy(original.get("artifacts", {})) if original else {}
    if finalized:
        artifacts["raw-summary"] = copy.deepcopy(finalized)
    events = out / "events.jsonl"
    raw_events = original.get("artifacts", {}).get("events") if original else None
    if raw_events and checked_artifact(raw_events).resolve() != events.resolve():
        raise ValueError("finalized journal differs from owned phase path")
    if events.exists():
        prefix = out / "completed-events.jsonl"
        with events.open("rb") as source, prefix.open("xb") as dest:
            for line in source:
                try:
                    json.loads(line)
                except (json.JSONDecodeError, UnicodeDecodeError):
                    break
                dest.write(line)
            dest.flush()
            os.fsync(dest.fileno())
        artifacts["events"] = artifact(prefix)
        # The torn original is hash-checked evidence, never a traversed stream.
        artifacts["raw-events"] = copy.deepcopy(raw_events) if raw_events else artifact(events)
    process_path = out / "process-error.json"
    write_new(process_path, {"observation": row, "failure": process_failure,
                             "summary_parse_failure": summary_error})
    artifacts["process-error"] = artifact(process_path)
    target = out / "process-failure.json"
    write_new(target, {"schema": "modena-spme-phase-v1", "phase": phase,
                       "prepared": str(prepared.resolve()), "manifest_sha256": epoch,
                       "prior": copy.deepcopy(previous),
                       "inputs": copy.deepcopy(original.get("inputs", {})) if original else {},
                       "artifacts": artifacts, "usage": usage,
                       "failure": (original.get("failure") if original else None) or process_failure,
                       "finalization_failure": (original.get("finalization_failure") if original else None) or summary_error})
    return target


class Campaign:
    def __init__(self, out, registered_binary, executor=execute, registered_sources=None):
        self.out = Path(out)
        self.binary_artifact = copy.deepcopy(registered_binary)
        self.binary = str(checked_artifact(self.binary_artifact))
        self.executor = executor
        self.source_hashes = copy.deepcopy(registered_sources) if registered_sources is not None else source_inventory()
        self.prepared_artifacts = {}
        self.rows = {}
        self.active = {}

    def call(self, key, phase, argv):
        if source_inventory() != self.source_hashes:
            raise ValueError("registered source changed during campaign")
        checked_artifact(self.binary_artifact)
        if key in self.prepared_artifacts:
            checked_artifact(self.prepared_artifacts[key])
        family, replicas, threads = key
        folder = self.out / f"{family}-{replicas}-threads-{threads}"
        folder.mkdir(exist_ok=True)
        env = dict(os.environ, RAYON_NUM_THREADS=str(threads))
        remaining = REPORT_NS if phase == "report" else WALL_NS - self.active.get(key, 0)
        if remaining <= 0:
            raise ValueError("family active deadline exhausted before child admission")
        row = self.executor([self.binary, phase, *map(str, argv)], env, folder / f"process-{phase}", phase, remaining)
        self.rows.setdefault(key, []).append(row)
        if phase in ACTIVE_PHASES:
            self.active[key] = self.active.get(key, 0) + row["active_ns"]
        return row

    @staticmethod
    def usable(row):
        return row["actual_child_exit"] in (0, 4) and row["monitor_exit"] == 0 and row["wrapper_exit"] == 0 and not row["timed_out"] and rss_accepted(row["peak_rss_bytes"])

    @classmethod
    def operationally_complete(cls, row):
        # Target exit 4 is meaningful only for scientific phases.
        return row is not None and row["actual_child_exit"] == 0 and cls.usable(row)

    def prepare_slot(self, key):
        slot = self.completion[key]
        try:
            prepared, row = self.prepare(key)
            missing = [name for name in ("manifest.json", "protocol.json", "snapshots.json")
                       if not (prepared / name).is_file()]
            slot["preparation"] = {"process": row, "missing_files": missing,
                                   "complete": self.operationally_complete(row) and not missing}
            return prepared, row
        except OSError as error:
            slot["preparation"] = {"complete": False, "error": str(error)}
            return None, None

    def render_slot(self, key, prepared, receipt, profiles=None):
        slot = self.completion[key]
        try:
            row = self.render(key, prepared, receipt, profiles)
            report = prepared.parent / "report"
            missing = [name for name in REPORT_FILES if not (report / name).is_file()]
            retained = {name: artifact(report / name) for name in REPORT_FILES if name not in missing}
            slot["report"] = {"attempted": True, "process": row, "artifacts": retained,
                              "missing_files": missing,
                              "complete": self.operationally_complete(row) and not missing}
        except OSError as error:
            slot["report"] = {"attempted": True, "complete": False, "error": str(error)}

    def recover(self, key, prepared, out, phase, previous, row, finalized=None):
        checked_artifact(self.prepared_artifacts[key])
        path = failed_receipt(prepared, out, phase, previous, row, self.active.get(key, 0), finalized)
        return retain_output(path, out, "recovery")

    def prepare(self, key):
        family, replicas, threads = key
        prepared = self.out / f"{family}-{replicas}-threads-{threads}" / "prepared"
        row = self.call(key, "prepare", ["--out", prepared, "--family", family, "--replicas", replicas, "--threads", threads, "--reference-extended"])
        if (prepared / "manifest.json").exists():
            self.prepared_artifacts[key] = retain_output(prepared / "manifest.json", prepared.parent / "process-prepare", "manifest")
        return prepared, row

    def verify(self, key, prepared, special=None, madelung=None):
        argv = ["--protocol", prepared / "protocol.json", "--snapshots", prepared / "snapshots.json", "--elapsed-ns", self.active.get(key, 0)]
        for name, receipt in (("special", special), ("madelung", madelung)):
            if receipt:
                argv += [f"--{name}", checked_artifact(receipt), f"--{name}-sha256", receipt["sha256"]]
        row = self.call(key, "verify", argv)
        path = prepared / "verification.json"
        receipt = retain_output(path, prepared.parent / "process-verify") if path.exists() else None
        if receipt is None or not self.usable(row):
            receipt = self.recover(key, prepared, prepared / "verification", "verify", None, row, receipt)
        return receipt, row

    def render(self, key, prepared, receipt, profiles=None):
        folder = prepared.parent
        checked_artifact(self.prepared_artifacts[key])
        campaign_path = checked_artifact(receipt)
        memory = folder / "memory.json"
        write_new(memory, {"schema": "modena-spme-memory-v1", "manifest_sha256": self.prepared_artifacts[key]["sha256"], "children": self.rows[key]})
        memory_ref = retain_output(memory, folder, "memory")
        inputs = {"campaign": receipt, "memory": memory_ref}
        argv = ["--campaign", campaign_path, "--expected-sha256", receipt["sha256"], "--memory", memory, "--memory-sha256", memory_ref["sha256"], "--out", folder / "report"]
        if profiles:
            data, summary = profiles
            argv += ["--profiles", checked_artifact(data), "--profiles-sha256", data["sha256"], "--profile-receipt-sha256", summary["sha256"]]
            checked_artifact(summary)
            inputs.update(profiles=data, profile_receipt=summary)
        write_new(folder / "report-inputs.json", inputs)
        return self.call(key, "report", argv)

    def run(self):
        self.completion = {(family, replicas, threads): {
            "family": family, "replicas": replicas, "threads": threads,
            "preparation": {"complete": False, "reason": "not attempted"},
            "report": {"attempted": False, "complete": False, "reason": "no report"}}
            for threads in THREADS for family, replicas in FAMILIES}
        fatal = None
        try:
            self.run_roster()
        except Exception as error:
            # Integrity violations still abort; never refresh a changed proof.
            fatal = f"{type(error).__name__}: {error}"
            raise
        finally:
            slots = []
            for key, slot in self.completion.items():
                slot["processes"] = self.rows.get(key, [])
                slot["complete"] = slot["preparation"]["complete"] and slot["report"]["complete"]
                slots.append(slot)
            success = fatal is None and all(slot["complete"] for slot in slots)
            write_new(self.out / "outcome.json", {
                "schema": "modena-spme-launcher-outcome-v1", "success": success,
                "exit_code": 0 if success else 1, "fatal_error": fatal,
                "scientific_status_separate": True, "slots": slots})
        return 0 if success else 1

    def run_roster(self):
        for threads in THREADS:
            special_key = ("special", 2, threads)
            special_prepared, preparation = self.prepare_slot(special_key)
            special = None
            special_row = preparation
            if self.completion[special_key]["preparation"]["complete"]:
                try:
                    special, special_row = self.verify(special_key, special_prepared)
                except OSError as error:
                    self.completion[special_key]["error"] = str(error)
            cubic64 = None
            for family, replicas in FAMILIES:
                key = (family, replicas, threads)
                try:
                    if key == special_key:
                        prepared, receipt, row = special_prepared, special, special_row
                    else:
                        prepared, row = self.prepare_slot(key)
                        if not self.completion[key]["preparation"]["complete"]:
                            continue
                        if special is None or not self.usable(special_row) or (family == "cubic" and replicas == 3 and cubic64 is None):
                            row = dict(row, actual_child_exit=None)
                            receipt = self.recover(key, prepared, prepared / "verification", "verify", None, row)
                        else:
                            receipt, row = self.verify(key, prepared, special, cubic64 if family == "cubic" and replicas == 3 else None)
                    if receipt is None:
                        continue
                    if family == "cubic" and replicas == 2 and self.usable(row):
                        cubic64 = receipt
                    if not self.usable(row):
                        self.render_slot(key, prepared, receipt)
                        continue
                    profiles = None
                    for phase in ("sweep", "confirm", "profile"):
                        out = prepared.parent / phase
                        checked_artifact(receipt)
                        if self.active.get(key, 0) >= WALL_NS:
                            # No process was admitted: do not invent an exit code or
                            # add a nonexistent child to memory.json.
                            not_started = {"actual_child_exit": None, "timed_out": False,
                                           "reason": "family active deadline exhausted before admission"}
                            receipt = self.recover(key, prepared, out, phase, receipt, not_started)
                            break
                        if phase == "sweep":
                            argv = ["--protocol", prepared / "protocol.json", "--snapshots", prepared / "snapshots.json", "--verification-sha256", receipt["sha256"]]
                        else:
                            argv = ["--campaign", receipt["path"], "--expected-sha256", receipt["sha256"]]
                        row = self.call(key, phase, [*argv, "--out", out, "--elapsed-ns", self.active.get(key, 0)])
                        summary_path = out / "summary.json"
                        summary = retain_output(summary_path, prepared.parent / f"process-{phase}") if summary_path.exists() else None
                        profile_data = None
                        if phase == "profile" and (out / "phase-times.jsonl").exists():
                            profile_data = retain_output(out / "phase-times.jsonl", prepared.parent / "process-profile", "profile-data")
                        if not self.usable(row) or summary is None or (phase == "profile" and profile_data is None):
                            receipt = self.recover(key, prepared, out, phase, receipt, row, summary)
                            break
                        if phase == "profile":
                            profiles = (profile_data, summary)
                        else:
                            receipt = summary
                    self.render_slot(key, prepared, receipt, profiles)
                except OSError as error:
                    self.completion[key]["error"] = str(error)



def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "--_monitor":
        return monitor(Path(argv[1]), float(argv[2]), argv[4:])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, nargs="+", default=THREADS)
    args = parser.parse_args(argv)
    validate_threads(args.threads)
    if sys.platform != "darwin":
        parser.error("frozen host requires numeric macOS /usr/bin/time -l RSS")
    args.out = args.out.resolve()
    args.out.mkdir()
    target = ROOT / "target/b05-build"
    env = dict(os.environ, CARGO_TARGET_DIR=str(target))
    if (target / "spme-calibration.lease").exists():
        raise FileExistsError("another calibration lease exists")
    frozen_sources = source_inventory()
    build = ["cargo", "build", "--release", "-p", "modena-core", "--example", "spme_convergence"]
    with (args.out / "build.stdout").open("xb") as stdout, (args.out / "build.stderr").open("xb") as stderr:
        subprocess.run(build, cwd=ROOT, env=env, stdout=stdout, stderr=stderr, check=True)
    if frozen_sources != source_inventory():
        raise ValueError("source changed during release build")
    binary = target / "release/examples/spme_convergence"
    registered_binary = artifact(binary)
    write_new(args.out / "registration.json", {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "build": build,
              "binary": registered_binary, "launcher": artifact(__file__), "protocol": artifact(ROOT / "docs/scientific-hardening/B05-SPME-PROTOCOL.md"),
              "source_sha256": frozen_sources, "families": FAMILIES, "threads": THREADS, "extended": True, "requested_heap_bytes": 268435456, "rss_cap_bytes": RSS_CAP, "family_active_ns": WALL_NS})
    # Exclusive campaign lease lives outside measured family time.
    lease = target / "spme-calibration.lease"
    fd = os.open(lease, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try:
        os.write(fd, str(os.getpid()).encode()); os.fsync(fd)
        result = Campaign(args.out, registered_binary, registered_sources=frozen_sources).run()
    finally:
        os.close(fd)
        lease.unlink()
    return result


if __name__ == "__main__":
    raise SystemExit(main())
