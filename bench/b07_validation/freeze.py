"""Build and bind the internal B07 epoch. Preparation performs no integration or force evaluation."""

import ctypes
import hashlib
import json
import os
import pathlib
import platform
import shutil
import subprocess
import sys

Path = pathlib.Path
ROOT = Path(__file__).resolve().parents[2]
SPEC = ROOT / "docs/superpowers/specs/2026-11-06-b07-tensor-pressure-design.md"
OLD = ROOT / "target/b07-tools/b05-baseline"
EVIDENCE = ROOT / "target/b07-evidence"
BASE = EVIDENCE / "baseline-source-preparation-20261110"
BASE_MANIFEST = "7067b18d389451510b84d5acee30afb722ed58b1236b51245e107c8c289de45e"
LIMITS = {
    "wall_seconds": 7200,
    "max_workers": 4,
    "rust_rss_bytes": 512 << 20,
    "reference_rss_bytes": 1 << 30,
    "artifact_bytes": 5 << 30,
    "molecule_steps": 400000000,
    "projected_draws": 1000000,
    "reference_jobs": 24,
}
OU_LOGS = [
    EVIDENCE / "task8-baoab-20261110/ou-02-modena_core.stdout",
    EVIDENCE / "task8-fix1-20261110/green-modena_core.stdout",
]
MESH = EVIDENCE / "task4-spme-tensor-20261110/scientific-result.json"
DFT = EVIDENCE / "task4-fix1-nonfinite-20261110/oracle-work-counts.json"
BUILD = ROOT / "target/b07-tools/whole-step-build"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def dump(path, obj):
    path = Path(path)
    if path.exists():
        raise FileExistsError(str(path))
    partial = path.with_name(path.name + ".partial")
    with partial.open("x") as f:
        json.dump(obj, f, sort_keys=True, allow_nan=False, separators=(",", ":"))
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.link(partial, path)
    partial.unlink()


def loads(data):
    def invalid(x):
        raise ValueError("nonfinite JSON " + x)

    def unique(items):
        d = {}
        for k, v in items:
            if k in d:
                raise ValueError("duplicate JSON field " + k)
            d[k] = v
        return d

    return json.loads(
        data.decode("utf-8"), parse_constant=invalid, object_pairs_hook=unique
    )


def load(p):
    return loads(Path(p).read_bytes())


def capture(cmd):
    return subprocess.check_output(list(map(str, cmd)), cwd=ROOT, text=True).strip()


def build(out, label, cmd):
    cmd = list(map(str, cmd))
    env = dict(os.environ)
    # No inherited rustflags or alternate target/target-feature selection.
    for name in (
        "RUSTFLAGS",
        "CARGO_ENCODED_RUSTFLAGS",
        "RUSTC_WRAPPER",
        "RUSTC_WORKSPACE_WRAPPER",
        "CARGO_TARGET_DIR",
        "RUSTC",
    ):
        if env.get(name):
            raise ValueError(f"unsupported compiler environment {name}")
    if Path(cmd[0]).name == "cargo":
        # Build scripts inherit the same absolute compiler, even from the archive's older toolchain file.
        env["RUSTC"] = str(Path(cmd[0]).with_name("rustc"))
    with (
        (out / f"{label}.stdout").open("xb") as stdout,
        (out / f"{label}.stderr").open("xb") as stderr,
    ):
        finished = subprocess.run(
            cmd, cwd=ROOT, env=env, stdout=stdout, stderr=stderr, check=False
        )
    dump(out / f"{label}.json", {"command": cmd, "exit": finished.returncode})
    finished.check_returncode()


def identities(outs):
    result = []
    for out in outs:
        for line in out.read_text().splitlines():
            row = json.loads(line)
            if row["reason"] in ("compiler-artifact", "build-script-executed"):
                result.append(row)
    return result


def cargo_sources(rows):
    files = set()
    for r in rows:
        if r["reason"] == "build-script-executed":
            files.update(x for x in Path(r["out_dir"]).rglob("*") if x.is_file())
            continue
        manifest = Path(r["manifest_path"])
        if manifest.is_relative_to(ROOT) and not manifest.is_relative_to(OLD):
            files.update((ROOT / "crates/modena-core/src").rglob("*.rs"))
            files.add(manifest)
        elif manifest.is_relative_to(OLD):
            files.update(
                p for p in (OLD / "crates/modena-core/src").rglob("*") if p.is_file()
            )
            files.update([manifest, OLD / "Cargo.toml", OLD / "Cargo.lock"])
        else:
            # Registry sources include build scripts, features and data included as bytes.
            files.update(p for p in manifest.parent.rglob("*") if p.is_file())
    return files


def binary(rows, name, test=False):
    matches = {
        Path(r["executable"])
        for r in rows
        if r.get("executable")
        and r["target"]["name"] == name
        and r["profile"]["test"] == test
    }
    if len(matches) != 1:
        raise ValueError(f"cannot identify {name}: {matches}")
    return matches.pop()


def rustlib(rows, name):
    matches = {
        p
        for r in rows
        if r.get("target", {}).get("name") == name
        for p in r["filenames"]
        if p.endswith(".rlib")
    }
    if len(matches) != 1:
        raise ValueError(f"cannot identify rlib {name}: {matches}")
    return matches.pop()


def host_identity():
    byte_uuid = ctypes.create_string_buffer(16)
    if not ctypes.CDLL(None)._dyld_get_shared_cache_uuid(byte_uuid):
        raise ValueError("cannot identify the macOS dyld shared cache")
    return {
        "dyld_shared_cache_uuid": byte_uuid.raw.hex(),
        "macos_build": capture(["/usr/bin/sw_vers", "-buildVersion"]),
        "machine": platform.machine(),
    }


def old_source():
    path = BASE / "source-manifest.json"
    if digest(path) != BASE_MANIFEST:
        raise ValueError("unrecognized pinned B05 archive manifest")
    info = load(path)
    if digest(BASE / "baseline.tar") != info["archive_sha256"]:
        raise ValueError("baseline archive changed")
    for name, item in info["files"].items():
        if item["kind"] == "file" and digest(OLD / name) != item["sha256"]:
            raise ValueError(f"old source changed {name}")
    return info


def inventory(paths):
    return {
        str(p.absolute()): digest(p)
        for p in sorted(set(map(Path, paths)))
        if p.is_file()
    }


def python_runtime_files(executable):
    """Loaded host interpreter/framework and non-shared-cache extension dependencies."""
    if executable.resolve() != Path(sys.executable).resolve():
        raise ValueError("inventory must run under the interpreter being frozen")
    dyld = ctypes.CDLL(None)
    name = dyld._dyld_get_image_name
    name.argtypes = [ctypes.c_uint32]
    name.restype = ctypes.c_char_p
    files = {executable.resolve()}
    for i in range(dyld._dyld_image_count()):
        path = Path(os.fsdecode(name(i)))
        if path.is_file():
            files.add(path.resolve())
        elif not str(path).startswith(("/usr/lib/", "/System/Library/")):
            raise ValueError(f"unresolved loaded host library {path}")
    return files


def prepare(out):
    out = Path(out).resolve()
    out.mkdir()
    inputs = out / "inputs"
    inputs.mkdir()
    bins = out / "bin"
    bins.mkdir()
    logs = out / "build"
    logs.mkdir()
    (out / "jobs").mkdir()
    baseline = old_source()
    rustc = Path(capture([str(Path.home() / ".cargo/bin/rustup"), "which", "rustc"]))
    cargo = Path(capture([str(Path.home() / ".cargo/bin/rustup"), "which", "cargo"]))
    files = {
        SPEC,
        MESH,
        DFT,
        *OU_LOGS,
        BASE / "source-manifest.json",
        BASE / "baseline.tar",
        ROOT / "Cargo.toml",
        ROOT / "Cargo.lock",
        ROOT / "rust-toolchain.toml",
        rustc,
        cargo,
        Path(sys.executable).resolve(),
        Path("/usr/bin/time"),
        Path("/usr/bin/shasum"),
    }
    # Every consumed Rust crate (including private lib-test support) and all workspace manifests.
    files.update(
        ROOT / name
        for name in capture(["git", "ls-files", "crates/modena-core"]).splitlines()
    )
    files.update((ROOT / "crates/modena-core/tests/support").rglob("*.rs"))
    files.update((ROOT / "crates/modena-core/examples/stress_validation").rglob("*.rs"))
    files.add(ROOT / "crates/modena-core/examples/stress_validation.rs")
    files.update(ROOT.glob("crates/*/Cargo.toml"))
    files.update((ROOT / "bench/b07_validation").glob("*.py"))
    files.update((ROOT / "bench/b07_validation").glob("*.sh"))
    files.add(Path("/bin/sh"))
    files.add(Path("/usr/bin/env"))
    files.update((ROOT / "bench/b07_validation").glob("*.rs"))
    prereg = out / "preregistered-protocol.md"
    shutil.copyfile(
        ROOT / "docs/scientific-hardening/B07-PRESSURE-VALIDATION.md", prereg
    )
    files.add(prereg)
    files.update(python_runtime_files(Path(sys.executable)))
    files.update(
        x for x in rustc.parent.parent.joinpath("lib").glob("*.dylib") if x.is_file()
    )
    syslib = Path(capture([rustc, "--print", "target-libdir"]))
    files.update(syslib.glob("*.rlib"))
    for d in [ROOT, *ROOT.parents, Path.home()]:
        for n in ("config", "config.toml"):
            if (d / ".cargo" / n).is_file():
                files.add(d / ".cargo" / n)
    common = [
        "--release",
        "--locked",
        "--offline",
        "-j",
        "4",
        "-p",
        "modena-core",
        "--message-format=json",
    ]
    build(logs, "example", [cargo, "build", *common, "--example", "stress_validation"])
    build(logs, "libtest", [cargo, "test", *common, "--lib", "--no-run"])
    build(
        logs,
        "oldlib",
        [
            cargo,
            "build",
            *common,
            "--lib",
            "--manifest-path",
            OLD / "Cargo.toml",
            "--target-dir",
            BUILD,
        ],
    )
    rows = identities(
        [logs / "example.stdout", logs / "libtest.stdout", logs / "oldlib.stdout"]
    )
    files.update(cargo_sources(rows))
    files.update(
        Path(v)
        for row in rows
        if row["reason"] == "compiler-artifact"
        for v in row["filenames"]
    )
    for tool in ("clang", "ld", "otool"):
        files.add(Path(capture(["xcrun", "--find", tool])))
    sdk = Path(capture(["xcrun", "--show-sdk-path"]))
    files.update(sdk.glob("usr/lib/libSystem*.tbd"))
    files.add(sdk / "SDKSettings.json")
    for source, label in [
        (binary(rows, "stress_validation"), "stress_validation"),
        (binary(rows, "modena_core", True), "ou-libtest"),
    ]:
        shutil.copy2(source, bins / label)
    for label, log in [("old", "oldlib"), ("new", "example")]:
        built = identities([logs / f"{log}.stdout"])
        core = Path(rustlib(built, "modena_core"))
        # Some cargo library targets emit a public top-level rlib; hashed dependencies are in deps.
        deps = core.parent if core.parent.name == "deps" else core.parent / "deps"
        src = ROOT / "bench/b07_validation/whole_step.rs"
        build(
            logs,
            "whole-step-" + label,
            [
                rustc,
                "--edition=2021",
                "-C",
                "opt-level=3",
                "-C",
                "lto=thin",
                "-C",
                "codegen-units=1",
                src,
                "-L",
                "dependency=" + str(deps),
                "--extern",
                "modena_core=" + str(core),
                "--extern",
                "serde_json=" + rustlib(built, "serde_json"),
                "-o",
                bins / ("whole-step-" + label),
            ],
        )
    build(
        logs,
        "freeze-inputs",
        [bins / "stress_validation", "--child", "freeze-inputs", inputs],
    )
    solver = {"tolerance": 1e-13, "max_iterations": 100}
    prior = []
    for fp in OU_LOGS:
        records = [
            json.loads(l.split("conditional_statistics ", 1)[1])
            for l in fp.read_text().splitlines()
            if "conditional_statistics " in l
        ]
        if (
            len(records) != 1
            or records[0]["total"] != sum(records[0]["draws_by_seed"])
            or records[0]["draws_by_seed"] != [8192] * 8
        ):
            raise ValueError("invalid historical draw accounting")
        prior.append(records[0]["total"])
    upper = {
        str(d): 2 * 8 * 32 * sum((50 + d) * i for i in (25, 50, 100, 200))
        for d in (400, 800, 1600)
    }
    if (
        upper["1600"] > LIMITS["molecule_steps"]
        or sum(prior) + 4 + 65536 > LIMITS["projected_draws"]
    ):
        raise ValueError("work cap")
    facts = {
        "source_head": capture(["git", "rev-parse", "HEAD"]),
        "source_tree": capture(["git", "rev-parse", "HEAD^{tree}"]),
        "old_source_head": baseline["source_head"],
        "old_source_tree": baseline["source_tree"],
        "compiler": capture([rustc, "-vV"]),
        "cargo": capture([cargo, "-Vv"]),
        "build_profile": "release",
        "features": [],
        "workers": 1,
        "platform": platform.platform(),
        "python": sys.version,
        "wall_seconds": 7200,
        "external_wall_seconds": 3600,
        "gas_work_upper": upper,
        "prior_projected_draws": prior,
        "prior_zero_amplitude_projector_calls": 4,
        "new_fixed_conditional_draw_upper": 65536,
        "conditional_cap_scope": "fixed-q conditional experiments; trajectory O stages separately counted as MD steps",
        "gas_sample_index": "k=1..10*duration after burn; dt .04 uses (5*k+1)//2; otherwise k*(5,10,20); actual times in rows; nominal .1 index clock",
        "nvt_streams": [2001, 2002, 2003, 2004],
        "nve_stream": 0,
        "orientation_seed": 901,
        "maxwell_seed": 902,
        "direct_snapshot_cutoff_ladder": [8, 10, 12],
        "spme_actual_mesh": [32, 32, 32],
        "table_cutoff": 2.5,
        "kappa": 0.8,
        "real_cutoff": 3.5,
        "benchmarks": {
            "integrator": "legacy/common default VV and new reporting VV",
            "dt": 0.005,
            "steps": {"rotor": 128, "gas-water": 32},
            "hamiltonian": "free",
            "old_and_new_workers": 1,
        },
        "external_roster": {"preflight": 4, "nve": 12, "nvt": 8, "total": 24},
    }
    facts["host_identity"] = host_identity()
    facts["binary_load_commands"] = {
        b.name: capture(["otool", "-L", b]) for b in bins.iterdir()
    }
    facts["compiled_rustc"] = str(rustc)
    facts["matched_rust_system_steps_upper"] = sum(
        j["steps"] for j in load(inputs / "roster.json")
    )
    facts["matched_rust_molecule_steps_upper"] = sum(
        j["steps"] * (32 if j["backend"] != "zero" else 1)
        for j in load(inputs / "roster.json")
    )
    from timing import legacy_snapshot, sample
    facts["clock_observation"] = sample()
    facts["clock_identity"] = facts["clock_observation"]["identity"]
    facts["legacy_clock"] = legacy_snapshot(EVIDENCE/"b07-campaign-budget-v1")
    facts["host_python"] = str(Path(sys.executable).resolve())
    from runtime import PYTHON_FLAGS, SAFE_ENV, THREAD_ENV
    facts["host_python_flags"] = PYTHON_FLAGS
    facts["numerical_environment"] = dict(SAFE_ENV, **THREAD_ENV)
    facts["host_startup"] = str(ROOT/"bench/b07_validation/launch.sh")
    facts["checksum_utility"] = "/usr/bin/shasum"
    driver = digest(ROOT / "bench/b07_validation/whole_step.rs")
    facts["driver_source_sha256"] = driver
    builds = {}
    for executable in sorted(bins.iterdir()):
        source_head = (
            baseline["source_head"]
            if executable.name == "whole-step-old"
            else facts["source_head"]
        )
        metadata = {
            "source_head": source_head,
            "compiler": facts["compiler"],
            "build_profile": "release",
            "features": [],
            "binary_sha256": digest(executable),
        }
        if executable.name.startswith("whole-step-"):
            metadata.update(
                driver_source_head=facts["source_head"],
                driver_sha256=driver,
                linked_core_source_head=source_head,
            )
        builds[str(executable)] = metadata
    facts["executable_builds"] = builds
    dump(inputs / "settings.json", facts)
    files.update(bins.iterdir())
    files.update(logs.iterdir())
    # The actual host stdlib modules loaded by the campaign/monitor are part of the launcher.
    files.update(
        Path(v.__file__)
        for v in list(sys.modules.values())
        if getattr(v, "__file__", None) and Path(v.__file__).is_file()
    )
    source = inventory(files)
    snapshots = inventory(inputs.iterdir())
    fixtures = [
        load(f)
        for f in sorted(inputs.glob("*.json"))
        if f.name not in ("roster.json", "settings.json")
    ]
    protocol = {
        "version": 1,
        "solver": solver,
        "strain_eps": [1e-4, 3e-5, 1e-5, 3e-6],
        "gas_seeds": list(range(1001, 1009)),
        "gas_dt": [0.04, 0.02, 0.01, 0.005],
        "burn_time": 50.0,
        "measure_times": [400.0, 800.0, 1600.0],
        "sample_interval": 0.1,
        "fixtures": fixtures,
        "limits": LIMITS,
        "specification_sha256": digest(SPEC),
        "source_sha256": source,
        "snapshot_sha256": snapshots,
    }
    dump(out / "protocol.json", protocol)
    pin = digest(out / "protocol.json")
    if (
        capture(["/usr/bin/shasum", "-a", "256", out / "protocol.json"]).split()[0]
        != pin
    ):
        raise ValueError("host SHA256 mismatch")
    dump(out / "internal-pin.json", {"protocol_sha256": pin})
    print(
        json.dumps(
            {
                "protocol_sha256": pin,
                "source_files": len(source),
                "fixtures": len(fixtures),
                "upper": upper,
                "binaries": inventory(bins.iterdir()),
            },
            indent=2,
        )
    )
