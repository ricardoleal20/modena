"""Bind measured files and prerequisite receipts to physical classifications."""
import hashlib
import json
from pathlib import Path
import re

import campaign as host
from freeze import digest, dump, load

RECEIPTS = ("launch.json", "child.json", "started.json", "wrapper.json", "host.time",
            "stdout.log", "stderr.log", "wrapper.stdout", "wrapper.stderr", "denied.json",
            "monitor-loss.json", "work.json")
ORDER = ("Failed", "BudgetExceeded", "Inconclusive", "Unsupported", "Pending", "Passed")


def worst(statuses):
    statuses = set(statuses)
    return next((s for s in ORDER if s in statuses), "Pending")


def files_state(paths):
    return {str(p): digest(p) if Path(p).is_file() else None for p in paths}


def verify_state(state):
    if files_state(map(Path, state)) != state:
        raise ValueError("measured evidence changed or disappeared")


def job_paths(out, epoch, name):
    job = epoch["jobs"][name]
    paths = [out/"jobs"/name/p for p in RECEIPTS]
    paths += [Path(p) for p in job["artifacts"]]
    return sorted(set(paths))


def launch_matches(out, epoch, name):
    job = epoch["jobs"][name]
    launch = load(out/"jobs"/name/"launch.json")
    if (any(launch.get(k) != v for k, v in job.items())
            or launch.get("name") != name
            or launch.get("protocol_sha256") != digest(out/"protocol.json")
            or launch.get("binary_sha256") != digest(job["argv"][0])):
        raise ValueError("producer launch differs from registered action")


def seal_job(out, epoch, name):
    """Called once by the host after it reaps the process and writes wrapper receipts."""
    folder = out/"jobs"/name
    path = folder/"artifacts.json"
    if path.exists():
        return verify_job(out, epoch, name)
    launch_matches(out, epoch, name)
    seal = {"action": name, "epoch_sha256": digest(out/"external-epoch.json"),
            "files": files_state(job_paths(out, epoch, name))}
    dump(path, seal)
    return seal


def verify_job(out, epoch, name):
    seal = load(out/"jobs"/name/"artifacts.json")
    if (seal["action"] != name or seal["epoch_sha256"] != digest(out/"external-epoch.json")
            or set(seal["files"]) != {str(p) for p in job_paths(out, epoch, name)}):
        raise ValueError("unbound producer artifact inventory")
    verify_state(seal["files"])
    launch_matches(out, epoch, name)
    return seal


def process_status(folder, child, role):
    """Only narrow native numerical diagnostics and the Rust status contract imply failure."""
    if child.get("budget_stop"):
        return "BudgetExceeded", {"reason": child["budget_stop"]}
    if child.get("interrupted") or not host.reconciled_wrapper(folder, child):
        return "Inconclusive", {"reason": "interrupted or unreconciled producer"}
    if folder.name == "ase-six-strain" and not host.worker_coverage(child):
        return "Inconclusive", {"reason": "ASE worker RSS coverage missing"}
    code = child["exit"]
    if code == 0:
        return "Passed", {}
    if role == "rust" and code in (2, 3, 4):
        return {2: "Unsupported", 3: "Failed", 4: "Inconclusive"}[code], {"exit": code}
    if role == "reference":
        pattern = re.compile(r"^ERROR(?: on proc \d+)?: (Lost atoms:|Shake determinant < 0|"
                             r"Non-numeric atom coords|Non-numeric pressure|Bond atoms .* missing|"
                             r"Shake atoms .* missing)")
        for name in ("stdout.log", "stderr.log"):
            path = folder/name
            if path.is_file():
                with path.open(errors="replace") as stream:
                    for number, line in enumerate(stream, 1):
                        if pattern.search(line):
                            return "Failed", {"exit": code, "diagnostic_path": str(path),
                                              "line": number, "diagnostic": line.rstrip()[:1024]}
    return "Inconclusive", {"exit": code, "reason": "unclassified execution error"}


def job_state(out, epoch, name):
    folder = out/"jobs"/name
    if (folder/"denied.json").exists():
        return "BudgetExceeded", {"denial": str(folder/"denied.json")}
    if not (folder/"launch.json").exists():
        return "Pending", {"reason": "action never launched"}
    try:
        verify_job(out, epoch, name)
        child = load(folder/"child.json")
        if name == "ase-six-strain" and child.get("ase_worker_interpreter") != epoch["ase_runtime"]["process_executable"]:
            raise ValueError("worker coverage targets an unfrozen executable")
        status, details = process_status(folder, child, epoch["jobs"][name]["role"])
        return status, dict(details, receipt=str(folder/"child.json"),
                            receipt_sha256=digest(folder/"child.json"),
                            artifacts_sha256=digest(folder/"artifacts.json"))
    except (OSError, KeyError, ValueError) as exc:
        return "Inconclusive", {"reason": str(exc), "receipt": str(folder/"child.json")}


def family_binding(out, epoch, name, seen=None):
    """Recursively verifies only the producer/dependency graph frozen before outcomes."""
    seen = set() if seen is None else set(seen)
    if name in seen:
        raise ValueError("cyclic family dependency")
    seen.add(name)
    contract = epoch["family_inputs"][name]
    paths = set()
    states = {}
    for producer in contract["producers"]:
        state, details = job_state(out, epoch, producer)
        states[producer] = {"status": state, "details": details}
        paths.update(job_paths(out, epoch, producer))
        paths.add(out/"jobs"/producer/"artifacts.json")
    for family in contract["families"]:
        rows = read_family(out, epoch, family, seen)
        states[family] = {"status": worst(r["status"] for r in rows)}
        paths.add(out/"reference/comparisons"/(family+".json"))
    return {"files": files_state(sorted(paths)), "dependencies": states}


def verify_family(out, epoch, name, seen=None):
    row = load(out/"reference/comparisons"/(name+".json"))
    if (row.get("epoch_sha256") != digest(out/"external-epoch.json")
            or row.get("binding") != family_binding(out, epoch, name, seen)):
        raise ValueError("physical family evidence or receipt changed")
    rows = row["gates"]
    if (len(rows) != len(epoch["expected_gate_families"][name])
            or {r["name"] for r in rows} != set(epoch["expected_gate_families"][name])):
        raise ValueError("physical family roster changed")
    for gate in rows:
        if gate.get("status") == "Failed" and not gate.get("observed"):
            # This exception is evidence-bound here; checked_gates stays strict everywhere else.
            own = epoch["family_inputs"][name].get("owners", epoch["family_inputs"][name]["producers"])
            deps = row["binding"]["dependencies"]
            if (gate != host.unmeasured(gate["name"], "Failed")
                    or not any(deps[p]["status"] == "Failed" for p in own)
                    or row["details"].get("producer_diagnostics") != deps):
                raise ValueError("empty Failed row lacks sealed numerical producer failure")
        else:
            host.checked_gates([gate])
    if any(r["status"] == "Passed" for r in rows):
        if any(v["status"] != "Passed" for v in row["binding"]["dependencies"].values()):
            raise ValueError("physical pass has unsuccessful dependency")
    return rows


def read_family(out, epoch, name, seen=None):
    path = out/"reference/comparisons"/(name+".json")
    expected = epoch["expected_gate_families"][name]
    try:
        if not path.exists():
            if not (out/"clock.json").exists():
                return [host.unmeasured(n) for n in expected]
            # An interrupted adapter may have recorded the process but no comparison.
            # Preserve its verified failure/blocker without performing a physical reduction.
            binding = family_binding(out, epoch, name, seen)
            deps = binding["dependencies"]
            blocked = {k:v["status"] for k,v in deps.items() if v["status"] != "Passed"}
            if not any(v not in ("Passed", "Pending") for v in blocked.values()):
                return [host.unmeasured(n) for n in expected]
            contract = epoch["family_inputs"][name]
            own = contract.get("owners", contract["producers"])
            attempted = [deps[k]["status"] for k in own if deps[k]["status"] not in ("Passed", "Pending")]
            status = worst(attempted) if attempted else "Pending"
            row = dict(gates=[host.unmeasured(n,status) for n in expected], binding=binding,
                       epoch_sha256=digest(out/"external-epoch.json"),
                       details=dict(blocked_by=blocked,producer_diagnostics=deps,
                                    collection="retained attempted producer; reduction not completed"))
            if family_binding(out, epoch, name, seen) != binding:
                raise ValueError("producer changed during failure collection")
            path.parent.mkdir(exist_ok=True)
            if not path.exists():
                dump(path,row)
        return verify_family(out, epoch, name, seen)
    except (OSError, KeyError, ValueError, TypeError):
        return [host.unmeasured(n, "Inconclusive") for n in expected]


def require_passes(out, epoch, names):
    for name in names:
        try:
            rows = verify_family(out, epoch, name)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise ValueError("blocked by external prerequisite " + name + ": " + str(exc)) from exc
        if any(r["status"] != "Passed" for r in rows):
            raise ValueError("blocked by external prerequisite " + name)


def origin_contract(out, epoch, tag):
    section = epoch["sections"][tag]
    producer = section["producer"]
    job = epoch["jobs"][producer]
    if job["role"] != "reference":
        raise ValueError("origin must be a native producer")
    return dict(tag=tag, producer_job=producer,
                raw_path=str(out/"reference"/producer/(tag+".dump")),
                log_path=str(out/"reference"/producer/"native.log"),
                input_path=str(out/"reference"/(tag+".jsonl")),
                fixture=section["fixture"])


def projected_frames(out, epoch, fixtures, tag):
    from b07_reference.records import native_rows
    section = epoch["sections"][tag]
    for raw in native_rows(out, section, fixtures[section["fixture"]]):
        frame = raw["frame"]
        if frame["step"] == 0 or frame["step"] > section["burn_steps"]:
            frame["time"] = frame["step"] * section["dt"]
            yield (json.dumps(frame, separators=(",", ":"), allow_nan=False)+"\n").encode()


def origin_binding(out, epoch, tag):
    contract = origin_contract(out, epoch, tag)
    state, details = job_state(out, epoch, contract["producer_job"])
    if state != "Passed":
        raise ValueError("native origin producer unavailable: " + state + " " + str(details))
    return dict(contract, raw_sha256=digest(contract["raw_path"]),
                log_sha256=digest(contract["log_path"]),
                producer_receipt_sha256=digest(out/"jobs"/contract["producer_job"]/"child.json"),
                producer_artifacts_sha256=digest(out/"jobs"/contract["producer_job"]/"artifacts.json"))


def validate_projection(out, epoch, tag, fixtures):
    binding = origin_binding(out, epoch, tag)
    src = Path(binding["input_path"])
    if load(str(src)+".derivation.json").get("status") != "Completed" or Path(str(src)+".partial").exists():
        raise ValueError("native projection is incomplete")
    expected = hashlib.sha256()
    count = 0
    for data in projected_frames(out, epoch, fixtures, tag):
        expected.update(data)
        count += 1
    sha = expected.hexdigest()
    if digest(src) != sha or load(str(src)+".origin.json") != dict(binding, sha256=sha):
        raise ValueError("derived snapshot differs from prescribed native projection")
    # Detect a producer changed during a long projection.
    if origin_binding(out, epoch, tag) != binding:
        raise ValueError("native producer changed during projection")
    return sha, count


def validate_matched_input(out, epoch, name):
    action = epoch["jobs"][name]
    contract = action["native_origin"]
    if contract != origin_contract(out, epoch, contract["tag"]):
        raise ValueError("frozen native origin does not match section")
    if action["argv"][-2] != contract["input_path"] or action["argv"][4] != contract["fixture"]:
        raise ValueError("matched invocation differs from frozen native origin")
    fixtures = {f["id"]: f for f in load(out/"protocol.json")["fixtures"]}
    sha, count = validate_projection(out, epoch, contract["tag"], fixtures)
    if count != action["expected_frames"]:
        raise ValueError("native projection frame count")
    return sha


DERIVATION_RECEIPT_BYTES = 2048


def derivation_receipt(path, receipt, create=False):
    """A reserved slot can record a budget stop without growing the campaign."""
    encoded = (json.dumps(receipt, separators=(",", ":"))+"\n").encode()
    if len(encoded) > DERIVATION_RECEIPT_BYTES:
        raise ValueError("derivation receipt exceeds reserved slot")
    with Path(path).open("xb" if create else "r+b") as stream:
        if not create and Path(path).stat().st_size != DERIVATION_RECEIPT_BYTES:
            raise ValueError("missing reserved derivation receipt")
        stream.write(encoded.ljust(DERIVATION_RECEIPT_BYTES, b" "))
        stream.flush()
        import os
        os.fsync(stream.fileno())
