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
RSS_CAP = 536_870_912
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
        with Path(path).open() as stream:
            for sequence, line in enumerate(stream):
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    break
                if event.get("sequence") != sequence or event.get("manifest_sha256") != epoch:
                    raise ValueError("journal source/sequence mismatch")
                visit(event["data"])
    return state


def failed_receipt(prepared, out, phase, previous, row, active_ns):
    """Parent-authorized process failure envelope; it can never mint a proof."""
    prepared = Path(prepared)
    out = Path(out)
    out.mkdir(exist_ok=True)
    epoch = digest(prepared / "manifest.json")
    prior = json.loads(Path(previous).read_text()) if previous else None
    usage = recover_usage(out / "events.jsonl", prior["usage"] if prior else empty_usage(), epoch)
    usage["candidate"]["reference"]["elapsed_ns"] = max(active_ns, usage["candidate"]["reference"]["elapsed_ns"])
    code = row["actual_child_exit"]
    category = {86: "requested Rust heap cap", 87: "primary allocator callback", 88: "System allocator returned null"}.get(code, "process/monitor/resource failure")
    # Raw partial files stay retained. Only a complete JSONL prefix can be read
    # as report events; a torn tail remains separately inspectable.
    artifacts = {}
    events = out / "events.jsonl"
    if events.exists():
        prefix = out / "completed-events.jsonl"
        with events.open("rb") as source, prefix.open("xb") as dest:
            for line in source:
                try: json.loads(line)
                except json.JSONDecodeError: break
                dest.write(line)
            dest.flush(); os.fsync(dest.fileno())
        artifacts["events"] = artifact(prefix)
        artifacts["raw-events"] = artifact(events)
    target = out / "process-failure.json"
    write_new(target, {"schema": "modena-spme-phase-v1", "phase": phase,
                       "prepared": str(prepared.resolve()), "manifest_sha256": epoch,
                       "prior": artifact(previous) if previous else None, "inputs": {},
                       "artifacts": artifacts, "usage": usage,
                       "failure": {"BudgetExceeded": f"{category}; actual child exit {code}; timed_out={row['timed_out']}"},
                       "finalization_failure": None})
    return target


class Campaign:
    def __init__(self, out, binary, executor=execute):
        self.out = Path(out)
        self.binary = str(Path(binary).resolve())
        self.executor = executor
        self.source_hashes = source_inventory()
        self.rows = {}
        self.active = {}

    def call(self, key, phase, argv):
        if source_inventory() != self.source_hashes:
            raise ValueError("registered source changed during campaign")
        family, replicas, threads = key
        folder = self.out / f"{family}-{replicas}-threads-{threads}"
        folder.mkdir(exist_ok=True)
        env = dict(os.environ, RAYON_NUM_THREADS=str(threads))
        remaining = WALL_NS - self.active.get(key, 0)
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

    def prepare(self, key):
        family, replicas, threads = key
        prepared = self.out / f"{family}-{replicas}-threads-{threads}" / "prepared"
        row = self.call(key, "prepare", ["--out", prepared, "--family", family, "--replicas", replicas, "--threads", threads, "--reference-extended"])
        return prepared, row

    def verify(self, key, prepared, special=None, madelung=None):
        argv = ["--protocol", prepared / "protocol.json", "--snapshots", prepared / "snapshots.json", "--elapsed-ns", self.active.get(key, 0)]
        for name, receipt in (("special", special), ("madelung", madelung)):
            if receipt: argv += [f"--{name}", receipt, f"--{name}-sha256", digest(receipt)]
        row = self.call(key, "verify", argv)
        receipt = prepared / "verification.json"
        if not receipt.exists() or not self.usable(row):
            receipt = failed_receipt(prepared, prepared / "verification", "verify", receipt if receipt.exists() else None, row, self.active.get(key, 0))
        return receipt, row

    def render(self, key, prepared, receipt, profiles=None):
        folder = prepared.parent
        memory = folder / "memory.json"
        write_new(memory, {"schema": "modena-spme-memory-v1", "manifest_sha256": digest(prepared / "manifest.json"), "children": self.rows[key]})
        argv = ["--campaign", receipt, "--expected-sha256", digest(receipt), "--memory", memory, "--memory-sha256", digest(memory), "--out", folder / "report"]
        if profiles:
            argv += ["--profiles", profiles / "phase-times.jsonl", "--profiles-sha256", digest(profiles / "phase-times.jsonl"), "--profile-receipt-sha256", digest(profiles / "summary.json")]
        return self.call(key, "report", argv)

    def run(self):
        for threads in THREADS:
            special_key = ("special", 2, threads)
            special_prepared, preparation = self.prepare(special_key)
            special = None
            if self.usable(preparation):
                special, special_row = self.verify(special_key, special_prepared)
            else:
                special_row = preparation
            cubic64 = None
            for family, replicas in FAMILIES:
                key = (family, replicas, threads)
                if key == special_key:
                    prepared, receipt, row = special_prepared, special, special_row
                else:
                    prepared, row = self.prepare(key)
                    if not self.usable(row):
                        continue
                    if special is None or not self.usable(special_row) or (family == "cubic" and replicas == 3 and cubic64 is None):
                        row = dict(row, actual_child_exit=None)
                        receipt = failed_receipt(prepared, prepared / "verification", "verify", None, row, 0)
                    else:
                        receipt, row = self.verify(key, prepared, special, cubic64 if family == "cubic" and replicas == 3 else None)
                if receipt is None:
                    continue
                if family == "cubic" and replicas == 2 and self.usable(row):
                    cubic64 = receipt
                if not self.usable(row):
                    self.render(key, prepared, receipt)
                    continue
                profiles = None
                for phase in ("sweep", "confirm", "profile"):
                    out = prepared.parent / phase
                    if phase == "sweep":
                        argv = ["--protocol", prepared / "protocol.json", "--snapshots", prepared / "snapshots.json", "--verification-sha256", digest(receipt)]
                    else:
                        argv = ["--campaign", receipt, "--expected-sha256", digest(receipt)]
                    row = self.call(key, phase, [*argv, "--out", out, "--elapsed-ns", self.active.get(key, 0)])
                    summary = out / "summary.json"
                    if not self.usable(row) or not summary.exists():
                        failure = failed_receipt(prepared, out, phase, summary if summary.exists() else receipt, row, self.active.get(key, 0))
                        receipt = failure
                        break
                    if phase == "profile":
                        profiles = out
                    else:
                        receipt = summary
                self.render(key, prepared, receipt, profiles)


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
    write_new(args.out / "registration.json", {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "build": build,
              "binary": artifact(binary), "launcher": artifact(__file__), "protocol": artifact(ROOT / "docs/scientific-hardening/B05-SPME-PROTOCOL.md"),
              "source_sha256": frozen_sources, "families": FAMILIES, "threads": THREADS, "extended": True, "requested_heap_bytes": 268435456, "rss_cap_bytes": RSS_CAP, "family_active_ns": WALL_NS})
    # Exclusive campaign lease lives outside measured family time.
    lease = target / "spme-calibration.lease"
    fd = os.open(lease, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try:
        os.write(fd, str(os.getpid()).encode()); os.fsync(fd)
        Campaign(args.out, binary).run()
    finally:
        os.close(fd)
        lease.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
