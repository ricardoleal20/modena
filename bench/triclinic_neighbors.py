"""Dedicated B04 cost campaign. stdlib only; durations are never synthesized here.

Generation freezes input bytes. Run requires an explicitly acquired CPU lease.
Analysis rejects malformed/partial records and retains every warmup and failure.
"""
import argparse
import hashlib
import json
import math
import os
import pathlib
import platform
import random
import statistics
import subprocess
import sys
import time

SEED = 20261106
ROOT = pathlib.Path(__file__).resolve().parents[1]
PHYSICS = dict(model="LennardJones", precision="f64", sigma=1.0, epsilon=1.0,
               mass=1.0, cutoff=2.5, cutoff_treatment="shifted_force", tail_correction=False)
PROTOCOL = dict(dt=0.0005, warmup_blocks=5, sample_blocks=20, block_steps=10)


class DataError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise DataError(message)


def decode_json(raw):
    def pairs(items):
        d = {}
        for key, value in items:
            require(key not in d, "duplicate JSON key: " + key)
            d[key] = value
        return d
    try:
        return json.loads(raw, object_pairs_hook=pairs,
                          parse_constant=lambda x: (_ for _ in ()).throw(DataError("nonfinite JSON: " + x)))
    except (ValueError, TypeError) as exc:
        raise DataError(str(exc)) from exc


def integer(x, lower=0):
    return type(x) is int and x >= lower


def number(x, lower=0):
    return type(x) in (int, float) and math.isfinite(x) and x >= lower


def exact_keys(d, keys, where):
    require(type(d) is dict and set(d) == set(keys.split()), "invalid fields in " + where)


def bootstrap(ratios):
    rng = random.Random(SEED)
    draws = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(10000))
    return [draws[249], draws[9749]]


def analyze_case(data, expected=None):
    exact_keys(data, "schema_version status metadata validation records error", "case")
    require(type(data["schema_version"]) is int and data["schema_version"] == 1, "schema")
    require(data["status"] == "complete" and data["error"] is None, "incomplete or failed child")
    m = data["metadata"]
    exact_keys(m, "schema_version base_rev source_sha256 build_rev build_source_sha256 binary_sha256 snapshot_sha256 case_id family n model precision cutoff_treatment epsilon sigma mass cutoff tail_correction host compiler profile features threads mode skin dt warmup_blocks sample_blocks block_steps candidate paired memory_scope allocation_scope instrumentation timing_scope", "metadata")
    for key, width in [("base_rev", 40), ("build_rev", 40), ("source_sha256", 64),
                       ("build_source_sha256", 64), ("binary_sha256", 64), ("snapshot_sha256", 64)]:
        require(type(m[key]) is str and len(m[key]) == width and all(c in "0123456789abcdef" for c in m[key]), key)
    require(m["base_rev"] == m["build_rev"] and m["source_sha256"] == m["build_source_sha256"], "build/runtime source mismatch")
    if expected:
        for key, value in expected.items():
            require(key in m and type(m[key]) is type(value) and m[key] == value, "identity mismatch: " + key)
    require(type(m["schema_version"]) is int and m["schema_version"] == 1, "metadata schema")
    for key in ["case_id", "family", "host", "compiler", "memory_scope", "allocation_scope", "instrumentation", "timing_scope"]:
        require(type(m[key]) is str and bool(m[key]), key)
    require(integer(m["n"], 2) and integer(m["threads"], 1) and m["threads"] in (1, 4), "N/threads")
    require(m["mode"] in ("rebuild", "step") and number(m["skin"]) and m["skin"] in (0, 0.3), "mode/skin")
    require(m["candidate"] in ("fractional", "auto") and m["paired"] is True, "paired strategies required")
    require(m["profile"] == "release", "cost evidence requires release profile")
    require(m["features"] == ["neighbor-diagnostics"], "feature contract")
    require(m["memory_scope"] == "requested_neighbor_workspace; retained_process_requested_live_bytes" and
            m["allocation_scope"] == "process_alloc_and_realloc_calls_in_timed_kernel" and
            m["instrumentation"] == "counting_global_allocator_relaxed_atomics", "counter scope")
    require(m["timing_scope"] == ("complete_steps_with_per_step_diagnostic_accounting" if m["mode"] == "step" else
            "forced_neighbor_builds_only_continuous_VV_advance_outside_timers"), "timing scope")
    for key, value in PHYSICS.items():
        require(type(m[key]) is type(value) and m[key] == value, "physical setting: " + key)
    for key, value in PROTOCOL.items():
        require(type(m[key]) is type(value) and m[key] == value, "protocol: " + key)
    v = data["validation"]
    exact_keys(v, "initial final minimum_distance minimum_cutoff_margin", "validation")
    require(v["initial"] is True and v["final"] is True and number(v["minimum_distance"]) and v["minimum_distance"] > 0 and number(v["minimum_cutoff_margin"]) and v["minimum_cutoff_margin"] > 1e-8, "physical preflight/final failure")
    require(type(data["records"]) is list and len(data["records"]) == 50, "incomplete paired roster")
    ratios = []; counts = {"scalar": 0, m["candidate"]: 0}; measured = []
    cursor = 0
    errors = "max_energy_error max_force_error virial_error max_position_error max_velocity_error".split()
    for phase, count in [("warmup", 5), ("measured", 20)]:
        for i in range(count):
            order = ["scalar", m["candidate"]] if i % 2 == 0 else [m["candidate"], "scalar"]
            pair = {}
            for at, strategy in enumerate(order):
                r = data["records"][cursor]; cursor += 1
                exact_keys(r, "phase sample_index order strategy requested selected fallback selected_counts fallback_counts block_steps elapsed_ns rebuild_count candidate_checks n_pairs requested_bytes retained_bytes allocation_count max_energy_error max_force_error virial_error max_position_error max_velocity_error pair_equal physics_passed", "record")
                require(r["phase"] == phase and integer(r["sample_index"]) and r["sample_index"] == i and integer(r["order"]) and r["order"] == at, "phase/index/ABBA order")
                require(r["strategy"] == r["requested"] == strategy, "strategy mixing")
                require(r["selected"] in ("scalar", "fractional") and (strategy != "scalar" or r["selected"] == "scalar"), "selected strategy")
                reasons = {"uncertified_geometry", "fractional_boundary", "image_domain", "grid_cap", "stencil_cap", "work_cap", "memory_cap", "all_bins", "candidate_cost"}
                require(r["fallback"] is None or r["fallback"] in reasons, "fallback")
                require(r["fallback"] is None or r["selected"] == "scalar", "fallback selected fractional")
                for key in ["rebuild_count", "candidate_checks", "n_pairs", "requested_bytes", "retained_bytes", "allocation_count", "block_steps", "elapsed_ns"]:
                    require(integer(r[key], 1 if key == "elapsed_ns" else 0), key)
                require(r["block_steps"] == 10 and r["n_pairs"] > 0, "steps/pairs")
                for key, allowed in [("selected_counts", {"scalar", "fractional"}), ("fallback_counts", reasons)]:
                    require(type(r[key]) is dict and set(r[key]) <= allowed and all(integer(x, 1) for x in r[key].values()), key)
                require(sum(r["selected_counts"].values()) == r["rebuild_count"] and sum(r["fallback_counts"].values()) <= r["rebuild_count"], "build accounting")
                if r["rebuild_count"] > 0:
                    require(r["candidate_checks"] > 0 and r["selected"] in r["selected_counts"] and
                            sum(r["fallback_counts"].values()) <= r["selected_counts"].get("scalar", 0), "new build diagnostics")
                if strategy == "scalar":
                    require(set(r["selected_counts"]) <= {"scalar"} and not r["fallback_counts"] and r["fallback"] is None, "scalar report")
                if r["rebuild_count"] == 0:
                    require(r["candidate_checks"] == 0 and r["requested_bytes"] == 0, "cached report counted twice")
                require(r["pair_equal"] is True and r["physics_passed"] is True, "pair/physics parity")
                # Binary enforces atol+rtol for each actual component; these are
                # maxima of |a-b|/(1+max(|a|,|b|)), so this independent guard is meaningful.
                require(all(number(r[key]) and r[key] <= 1e-10 for key in errors), "E/F/W/state parity")
                pair[strategy] = r
                if phase == "measured": counts[strategy] += r["rebuild_count"]
            require(pair["scalar"]["n_pairs"] == pair[m["candidate"]]["n_pairs"] and
                    pair["scalar"]["rebuild_count"] == pair[m["candidate"]]["rebuild_count"], "paired pair/cadence mismatch")
            if phase == "measured":
                ratios.append(pair[m["candidate"]]["elapsed_ns"] / pair["scalar"]["elapsed_ns"])
                measured.append(pair)
    median = statistics.median(ratios); interval = bootstrap(ratios)
    even = statistics.median(ratios[::2]); odd = statistics.median(ratios[1::2])
    drift = abs(statistics.median(ratios[:10]) / statistics.median(ratios[10:]) - 1)
    status = "no_gain"
    if m["mode"] == "step" and m["skin"] == 0.3 and min(counts.values()) < 2:
        status = "inconclusive_rebuilds"
    elif abs(even / odd - 1) > 0.10 or drift > 0.10 or interval[1] - interval[0] > 0.10 * median:
        status = "inconclusive_order_or_noise"
    elif median <= 0.95 and interval[1] < 1.0:
        status = "eligible_gain" if m["mode"] == "step" else "rebuild_gain_only"
    fallback_regime = m["n"] == 256 or m["family"] in ("grid_cap", "boundary")
    fallback_overhead_ok = not status.startswith("inconclusive") and all(r["selected"] == "scalar" and not r["selected_counts"].get("fractional", 0)
                               for p in measured for r in [p[m["candidate"]]]) and median <= 1.05
    return dict(status=status, mode=m["mode"], median_ratio=median, ci95=interval,
                bootstrap_samples=10000, bootstrap_seed=SEED, measured_rebuilds=counts,
                order_medians=[even, odd], first_last_ratio_drift=drift,
                fallback_overhead_required=fallback_regime, fallback_overhead_ok=fallback_overhead_ok,
                absolute_ns={s: [p[s]["elapsed_ns"] for p in measured] for s in counts})


def make_snapshot(n, family):
    length = (n / 0.8) ** (1 / 3)
    matrices = {
        "moderate": [[length, 0.0, 0.0], [0.2*length, length, 0.0], [0.1*length, 0.15*length, length]],
        "substantial": [[length, 0.0, 0.0], [0.7*length, length, 0.0], [0.2*length, 0.4*length, length]],
        "unreduced": [[length, 0.0, 0.0], [4.9*length, length, 0.0], [0.0, 0.0, length]],
        "elongated": [[2*length, 0.0, 0.0], [0.2*length, 0.5*length, 0.0], [0.1*length, 0.15*length, length]],
    }
    basis = matrices[family]
    side = math.ceil(n ** (1/3))
    while side ** 3 < n: side += 1
    positions = []; velocities = []
    for i in range(n):
        indices = [i % side, (i // side) % side, i // (side * side)]
        frac = [(indices[k] + 0.31 + 0.06 * (((i*1103515245 + (k+1)*12345 + SEED) % 65536) / 65535 - 0.5)) / side for k in range(3)]
        positions.append([(basis[0][k]*frac[0] + basis[1][k]*frac[1]) + basis[2][k]*frac[2] for k in range(3)])
        velocities.append([0.05 * (2 * ((i*1664525 + (k+1)*1013904223 + SEED) % 65536) / 65535 - 1) for k in range(3)])
    means = [sum(v[k] for v in velocities) / n for k in range(3)]
    velocities = [[v[k] - means[k] for k in range(3)] for v in velocities]
    return dict(schema_version=1, case_id=f"{family}-{n}", family=family, n=n, basis=basis,
                positions=positions, velocities=velocities, physical=PHYSICS.copy(),
                protocol=PROTOCOL.copy(), cell_after=None, generation=dict(seed=SEED, modulo=65536,
                integer_arithmetic="u64_no_overflow_then_modulo", density=0.8, side=side))


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def source_identity():
    # Actual bytes, including uncommitted edits. Docs/evidence do not change code identity.
    names = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    names += ["bench/triclinic_neighbors.py", "crates/modena-bench/src/bin/triclinic_neighbors.rs"]
    names = sorted({n for n in names if n and (n.endswith(".rs") or n.endswith("Cargo.toml") or n == "Cargo.lock" or n == "bench/triclinic_neighbors.py")})
    h = hashlib.sha256()
    for name in names:
        h.update(name.encode()); h.update(b"\0"); h.update((ROOT/name).read_bytes()); h.update(b"\0")
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), h.hexdigest()


def write_new(path, value):
    with pathlib.Path(path).open("x") as f:
        json.dump(value, f, sort_keys=True, indent=2, allow_nan=False); f.write("\n")


def generate(output, binary):
    require(binary is not None and pathlib.Path(binary).is_file(), "--binary must identify the frozen release executable")
    output.mkdir(parents=True, exist_ok=True)
    require(not (output/"manifest.json").exists(), "refusing to overwrite historical manifest")
    rev, source = source_identity()
    cases = []
    for n in [256, 2048, 8192]:
        for family in ["moderate", "substantial", "unreduced", "elongated"]:
            cases.append(make_snapshot(n, family))
    # Separate fallback and full-basis cache-invalidation controls. No selector tuning.
    for family in ["grid_cap", "boundary", "basis_change"]:
        d = make_snapshot(8 if family != "basis_change" else 256, "moderate")
        d.update(case_id=family, family=family)
        if family == "grid_cap":
            d["basis"] = [[1000., 0., 0.], [0.2, 1000., 0.], [0.1, 0.15, 1000.]]
            d["generation"]["density"] = 8 / 1000**3
            # Keep the serialized positions: sparse control with interacting pairs.
        elif family == "boundary":
            d["positions"][0] = [0., 0., 0.]
        else:
            d["cell_after"] = [row[:] for row in d["basis"]]
            d["cell_after"][1][0] = 0.21 * d["basis"][0][0]
        cases.append(d)
    snapshots = []
    for d in cases:
        path = output/(d["case_id"] + ".snapshot.json")
        write_new(path, d)
        snapshots.append(dict(case_id=d["case_id"], family=d["family"], n=d["n"], file=path.name, sha256=sha(path)))
    manifest = dict(schema_version=1, status="generated_not_evaluated", base_rev=rev,
                    source_sha256=source, binary_sha256=sha(binary) if binary else None,
                    binary=str(pathlib.Path(binary).resolve()) if binary else None,
                    host=platform.platform(), protocol=PROTOCOL, physics=PHYSICS, snapshots=snapshots,
                    skins=[0., 0.3], threads=[1, 4], modes=["rebuild", "step"],
                    candidate="fractional", case_cap_seconds=600, campaign_cap_seconds=7200,
                    bootstrap_samples=10000, bootstrap_seed=SEED,
                    scope="G4 preparation only; no model/tool-engine parity claim")
    write_new(output/"manifest.json", manifest)
    return manifest


def roster(manifest):
    return [(s, mode, skin, threads) for s in manifest["snapshots"]
            for mode in manifest["modes"] for skin in manifest["skins"] for threads in manifest["threads"]
            if s["family"] != "basis_change" or mode == "rebuild"]


def validate_manifest(m):
    exact_keys(m, "schema_version status base_rev source_sha256 binary_sha256 binary host protocol physics snapshots skins threads modes candidate case_cap_seconds campaign_cap_seconds bootstrap_samples bootstrap_seed scope", "manifest")
    require(type(m["schema_version"]) is int and m["schema_version"] == 1 and m["status"] == "generated_not_evaluated", "manifest schema/status")
    for key, width in [("base_rev", 40), ("source_sha256", 64), ("binary_sha256", 64)]:
        require(type(m[key]) is str and len(m[key]) == width and all(c in "0123456789abcdef" for c in m[key]), "manifest " + key)
    require(type(m["binary"]) is str and m["binary"] and type(m["host"]) is str and m["host"], "manifest binary/host")
    require(m["protocol"] == PROTOCOL and all(type(m["protocol"][k]) is type(v) for k, v in PROTOCOL.items()), "manifest protocol")
    require(m["physics"] == PHYSICS and all(type(m["physics"][k]) is type(v) for k, v in PHYSICS.items()), "manifest physics")
    require(m["skins"] == [0., 0.3] and all(type(x) is float for x in m["skins"]) and m["threads"] == [1, 4] and all(type(x) is int for x in m["threads"]), "manifest regimes")
    require(m["modes"] == ["rebuild", "step"] and m["candidate"] in ("fractional", "auto"), "manifest modes/candidate")
    for key, value in [("case_cap_seconds", 600), ("campaign_cap_seconds", 7200), ("bootstrap_samples", 10000), ("bootstrap_seed", SEED)]:
        require(type(m[key]) is int and m[key] == value, "manifest cap/bootstrap")
    expected = [(f"{f}-{n}", f, n) for n in [256, 2048, 8192] for f in ["moderate", "substantial", "unreduced", "elongated"]]
    expected += [("grid_cap", "grid_cap", 8), ("boundary", "boundary", 8), ("basis_change", "basis_change", 256)]
    require(type(m["snapshots"]) is list and len(m["snapshots"]) == len(expected), "snapshot roster")
    for s, (case_id, family, n) in zip(m["snapshots"], expected):
        exact_keys(s, "case_id family n file sha256", "snapshot manifest")
        require(s["case_id"] == case_id and s["family"] == family and type(s["n"]) is int and s["n"] == n and s["file"] == case_id+".snapshot.json", "snapshot roster/order")
        require(type(s["sha256"]) is str and len(s["sha256"]) == 64 and all(c in "0123456789abcdef" for c in s["sha256"]), "snapshot hash")
    return m


def case_name(s, mode, skin, threads):
    return f'{s["case_id"]}-{mode}-skin{skin:g}-t{threads}'


def run(output, binary):
    require(os.environ.get("MODENA_B04_CPU_LEASE") == "exclusive", "root must reserve exclusive team CPU/write pause and set MODENA_B04_CPU_LEASE=exclusive")
    manifest = validate_manifest(decode_json((output/"manifest.json").read_text()))
    rev, source = source_identity()
    require(rev == manifest["base_rev"] and source == manifest["source_sha256"] and sha(binary) == manifest["binary_sha256"], "frozen source/binary mismatch")
    start = time.monotonic()
    for s, mode, skin, threads in roster(manifest):
        name = case_name(s, mode, skin, threads); result = output/(name+".json")
        require(not result.exists() and not (output/(name+".child.json")).exists(), "refusing result reuse")
        require(sha(output/s["file"]) == s["sha256"], "snapshot changed")
        remaining = 7200 - (time.monotonic() - start)
        if remaining <= 0:
            write_new(output/(name+".child.json"), dict(status="campaign_cap", error="7200-second campaign cap")); break
        command = [str(binary), "--snapshot", str(output/s["file"]), "--strategy", "paired",
                   "--candidate", manifest["candidate"], "--mode", mode, "--skin", str(skin),
                   "--threads", str(threads), "--output", str(result)]
        for key, value in PROTOCOL.items(): command += ["--"+key.replace("_", "-"), str(value)]
        try:
            p = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=min(600, remaining))
            receipt = dict(command=command, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr,
                           campaign_elapsed_seconds=time.monotonic()-start)
        except subprocess.TimeoutExpired as exc:
            receipt = dict(command=command, status="case_or_campaign_cap", error=str(exc),
                           stdout=(exc.stdout or b"").decode() if isinstance(exc.stdout, bytes) else exc.stdout,
                           stderr=(exc.stderr or b"").decode() if isinstance(exc.stderr, bytes) else exc.stderr)
        except OSError as exc:
            receipt = dict(command=command, status="child_error", error=str(exc))
        write_new(output/(name+".child.json"), receipt)
        # Fail closed, keep partial output and all child diagnostics.
        if receipt.get("returncode") != 0: break


def analyze(output):
    manifest = validate_manifest(decode_json((output/"manifest.json").read_text()))
    results = []; invalid = False
    for s, mode, skin, threads in roster(manifest):
        name = case_name(s, mode, skin, threads)
        try:
            require(sha(output/s["file"]) == s["sha256"], "snapshot bytes changed")
            child = decode_json((output/(name+".child.json")).read_text())
            require(type(child.get("returncode")) is int and child["returncode"] == 0, "child failure")
            expected = dict(base_rev=manifest["base_rev"], source_sha256=manifest["source_sha256"],
                            binary_sha256=manifest["binary_sha256"], snapshot_sha256=s["sha256"],
                            case_id=s["case_id"], family=s["family"], n=s["n"], mode=mode, skin=skin, threads=threads,
                            candidate=manifest["candidate"])
            summary = analyze_case(decode_json((output/(name+".json")).read_text()), expected)
        except (DataError, OSError) as exc:
            summary = dict(status="invalid", error=str(exc)); invalid = True
        results.append(dict(case=name, **summary))
    # A harness report never silently promotes production Auto.
    intended = [r for r in results if r["case"].startswith(("moderate-2048-step-", "moderate-8192-step-"))]
    fallback = [r for r in results if r.get("mode") == "step" and r.get("fallback_overhead_required")]
    report = dict(schema_version=1, status="invalid_or_incomplete" if invalid else "complete_requires_gate_review",
                  g4="pending", auto="scalar", intended_step_gains=all(r["status"] == "eligible_gain" for r in intended),
                  small_fallback_overhead_ok=not invalid and all(r["fallback_overhead_ok"] for r in fallback),
                  instrumentation_limit="instrumented binary only; allocation atomic overhead unquantified; production throughput not established",
                  results=results)
    write_new(output/"analysis.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=pathlib.Path)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--generate-only", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    try:
        if args.generate_only: generate(args.output_dir, args.binary)
        elif args.run:
            require(args.binary is not None, "--binary is required for --run")
            run(args.output_dir.resolve(), args.binary.resolve())
            report = analyze(args.output_dir)
            if report["status"] == "invalid_or_incomplete": return 1
        else:
            report = analyze(args.output_dir)
            if report["status"] == "invalid_or_incomplete": return 1
    except (DataError, OSError) as exc:
        print(str(exc), file=sys.stderr); return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
