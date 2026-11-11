"""Metadata-only freeze of the installed reference build, inputs, jobs and physical gates."""

import os
import re
import subprocess
from pathlib import Path

from campaign import (
    ENV,
    ROOT,
    digest,
    dump,
    inventory,
    jobs as validate_jobs,
    load,
    verify,
    verify_files,
)
from .native import COMP, LADDER, NATIVE_BUILD, data_file, section, write, write_table

V3 = Path("/private/tmp/modena-b07-validation-v3")
V3_SHA = "1e2a8db5c20cd932f7aed61008addead9e05255a14eee09c94c8a2e82e879a67"
REV = "9c5ab448c78a14fd534619622162ba418d6a1fb1"
EVIDENCE = ROOT / "target/b07-evidence"
ASEROOT = ROOT / "target/b07-tools/ase-3.26"
DOCS = (
    "fix_shake",
    "fix_rigid",
    "fix_store_force",
    "compute_pressure",
    "thermo_modify",
    "dump_modify",
    "pair_coul",
    "kspace_modify",
    "pair_table",
    "bond_table",
    "bond_zero",
    "angle_zero",
    "special_bonds",
)
PAIRS = ("pair-oo", "pair-charge", "pair-excluded-in", "pair-excluded-out")
BASES = ["interacting32"] + [f"interacting32-{s}" for s in range(2001, 2005)]
ENERGY = [
    f"{energy}-{metric}"
    for energy in ("com", "rotation")
    for metric in ("mean", "variance", "cdf10", "cdf50", "cdf90")
]
GAS_METRICS = list(COMP) + ["scalar"] + ENERGY + ["orientation-" + x for x in COMP]
ASE_ENV = {
    "MODENA_ASE_PYTHON": str(ASEROOT / "bin/python"),
    "PYTHONPATH": "",
    "PYTHONNOUSERSITE": "1",
    "PYTHONDONTWRITEBYTECODE": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "PYTHONSAFEPATH": "1",
}


def inspect_help(s):
    if (
        "22 Jul 2025 - Update 6" not in s
        or "FFT precision  = double" not in s
        or "MPI STUBS" not in s
    ):
        raise ValueError("unexpected LAMMPS version, precision or MPI variant")
    packages = (
        s.split("Installed packages:", 1)[1].split("List of individual", 1)[0].split()
    )
    if not {"RIGID", "KSPACE", "MOLECULE"}.issubset(packages):
        raise ValueError("required LAMMPS packages absent")
    for group, required in {
        "Atom": ["full"],
        "Fix": ["store/force", "rigid/nvt", "nve", "rattle"],
        "KSpace": ["ewald"],
        "Pair": ["zero", "hybrid/overlay", "coul/long", "table"],
        "Bond": ["zero", "table"],
        "Angle": ["zero"],
        "Compute": ["pressure", "temp"],
        "Dump": ["custom"],
    }.items():
        styles = re.search(r"\* " + group + r" styles:?\s+(.*?)(?=\n\*)", s, re.S)
        if not styles or not set(required).issubset(styles[1].split()):
            raise ValueError("unsupported native " + group + " styles")
    return packages, re.search(r"^Compiler: (.+)$", s, re.M)[1]


def linked(binary):
    todo, names, listing = [binary], set(), {}
    while todo:
        file = todo.pop()
        if file in names:
            continue
        names.add(file)
        data = subprocess.check_output(
            ["/usr/bin/otool", "-L", str(file)], text=True, timeout=20
        )
        listing[str(file)] = data
        for name in re.findall(r"^\s*(\S+) \(compatibility", data, re.M):
            if name.startswith(("/usr/lib/", "/System/Library/")):
                continue  # verified macOS dyld cache UUID in internal protocol
            if name.startswith("@rpath/"):
                commands = subprocess.check_output(
                    ["/usr/bin/otool", "-l", str(file)], text=True, timeout=20
                )
                rpaths = re.findall(
                    r"cmd LC_RPATH\s+cmdsize \d+\s+path (.+?) \(offset", commands
                )
                matches = [
                    Path(d.replace("@loader_path", str(file.resolve().parent)))
                    / name[7:]
                    for d in rpaths
                ]
                matches = [p.resolve() for p in matches if p.is_file()]
                if len(set(matches)) != 1:
                    raise ValueError("ambiguous native LC_RPATH " + name)
                name = str(matches[0])
            if name.startswith("@loader_path/"):
                name = str(file.resolve().parent / name[13:])
            dep = Path(name)
            if not dep.is_absolute() or not dep.is_file():
                raise ValueError("unresolved linked native dependency " + name)
            if dep not in names:
                todo.append(dep)
    return names, listing


def help_output(binary, static):
    try:
        proc = subprocess.run(
            [str(binary), "-h"],
            env=dict(os.environ, **ENV),
            timeout=20,
            capture_output=True,
        )
    except subprocess.TimeoutExpired as exc:
        dump(static / "native-help-status.json", dict(exit=None, timeout_seconds=20))
        (static / "native-help.stdout").write_bytes(exc.stdout or b"")
        (static / "native-help.stderr").write_bytes(exc.stderr or b"")
        raise ValueError("metadata help timed out") from exc
    (static / "native-help.stdout").write_bytes(proc.stdout)
    (static / "native-help.stderr").write_bytes(proc.stderr)
    dump(
        static / "native-help-status.json",
        dict(exit=proc.returncode, timeout_seconds=None),
    )
    if proc.returncode:
        raise ValueError("metadata help failed")
    return proc.stdout.decode()


def provenance(binary, static):
    adopt = EVIDENCE / "adoption-20261110"
    release = load(adopt / "lammps-source-candidate.json")
    discovery = load(adopt / "lammps-discovery.json")
    prefix = binary.parent.parent
    if (
        binary != Path(discovery["resolved_binary"])
        or release["official_revision"] != REV
        or release["official_tag"] != "stable_22Jul2025_update6"
    ):
        raise ValueError("unexpected native release/binary")
    source = EVIDENCE / "task12-source-20261110/pair_coul_long.cpp"
    verify_files(
        {
            str(binary): discovery["binary_sha256"],
            release["linked_lammps_library"]: release["linked_lammps_library_sha256"],
            str(prefix / ".brew/lammps.rb"): release["formula_sha256"],
            str(
                source
            ): "bd732a11df7a3691934e2d372723d011e1132b093d1b784cd9b5a1b120ca1661",
        }
    )
    helptext = help_output(binary, static)
    packages, compiler = inspect_help(helptext)
    libs, listing = linked(binary)
    dump(static / "native-links.json", listing)
    files = libs | {
        source,
        *adopt.glob("lammps-*.json"),
        prefix / ".brew/lammps.rb",
        prefix / "INSTALL_RECEIPT.json",
    }
    docs = [prefix / f"share/lammps/doc/src/{name}.rst" for name in DOCS]
    if not all(f.is_file() for f in docs):
        raise ValueError("missing native primary documentation")
    ase_e = EVIDENCE / "ase-preparation-20261110"
    provider = load(ase_e / "environment.json")
    verify_files(provider["package_source_sha256"])
    verify_files(
        {
            str(ase_e / "requirements.lock.txt"): provider["lock_sha256"],
            str(ase_e / "pip-report.json"): provider["pip_report_sha256"],
        }
    )
    py = Path(provider["python"])
    if py != ASEROOT / "bin/python":
        raise ValueError("ASE interpreter drift")
    # Import-only metadata, no calculator construction, force call or integration.
    script = """import ase, numpy, scipy, ase.calculators.emt, ctypes, json, os, sys, sysconfig
d=ctypes.CDLL(None); g=d._dyld_get_image_name; g.argtypes=[ctypes.c_uint32]; g.restype=ctypes.c_char_p
print(json.dumps(dict(versions=[ase.__version__,numpy.__version__,scipy.__version__],
 stdlib=sysconfig.get_path('stdlib'),libraries=[os.fsdecode(g(i)) for i in range(d._dyld_image_count())],
 python=sys.version)))"""
    import json

    metadata = json.loads(
        subprocess.check_output(
            [str(py), "-c", script],
            text=True,
            cwd=static,
            env=dict(os.environ, **dict(ENV, **ASE_ENV)),
            timeout=40,
        )
    )
    if metadata["versions"] != ["3.26.0", "2.2.6", "1.15.3"]:
        raise ValueError("ASE runtime version drift")
    dump(static / "ase-runtime.json", metadata)
    libs = [
        Path(p)
        for p in metadata["libraries"]
        if not p.startswith(("/usr/lib/", "/System/Library/"))
    ]
    if not all(p.is_file() for p in libs):
        raise ValueError("missing ASE runtime library")
    site = ASEROOT / "lib/python3.12/site-packages"
    std = Path(metadata["stdlib"])
    files.update(
        [py, py.resolve(), ASEROOT / "pyvenv.cfg", *libs, *ase_e.iterdir(), *docs]
    )
    # Include importable bytecode as well as source, since a timestamp-valid cache can be loaded.
    files.update(p for p in site.rglob("*") if p.is_file())
    files.update(
        p
        for p in std.rglob("*")
        if p.is_file() and "site-packages" not in p.relative_to(std).parts
    )
    files.update((ROOT / "bench/b07_reference").glob("*.py"))
    files.update(
        [
            ROOT / "bench/b07_pressure_reference.py",
            ROOT / "crates/modena-core/src/external/ase_emt_worker.py",
        ]
    )
    return (
        inventory(files),
        packages,
        dict(compiler=compiler, build_profile=NATIVE_BUILD, features=packages),
    )


def ase_input(static):
    pos = [[0.03, 0, 0], [0, 1.8, 1.8], [1.8, 0, 1.8], [1.8, 1.8, 0]]
    cell = [[3.6, 0, 0], [0.3, 3.6, 0], [0.2, 0.1, 3.6]]
    base = dict(component=None, eta=0.0, positions=pos, cell=cell)
    snapshots = [base]
    for k, (a, b) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))):
        for eps in (1e-4, 3e-5, 1e-5, 3e-6):
            for eta in (-eps, eps):

                def deform(r):
                    v = list(r)
                    if a == b:
                        v[a] += eta * r[a]
                    else:
                        v[a] += 0.5 * eta * r[b]
                        v[b] += 0.5 * eta * r[a]
                    return v

                snapshots.append(
                    dict(
                        component=k,
                        eta=eta,
                        positions=list(map(deform, pos)),
                        cell=list(map(deform, cell)),
                    )
                )
    dump(
        static / "ase.json",
        dict(
            snapshots=snapshots,
            mass_amu=63.546,
            units=[1.0, 1.0, 1.0],
            provider=dict(
                python=ASE_ENV["MODENA_ASE_PYTHON"],
                ase_version="3.26.0",
                model_sha256="b8384ffed58deb1182134a0e5088260b7118971488c35cdb6b29e51f6e604f9d",
                atomic_numbers=[29] * 4,
                composition="total",
                pbc=[True] * 3,
                asap_cutoff=False,
                timeout_ms=30000,
                transport="persistent",
            ),
        ),
    )


def prepare(binary, protocol, out):
    if out != protocol.parent or protocol.name != "protocol.json":
        raise ValueError(
            "--out must be the registered campaign root containing --protocol"
        )
    p, _ = verify(out, False)
    if (out / "external-epoch.json").exists() or (out / "clock.json").exists():
        raise ValueError("source epoch already linked/started")
    if digest(V3 / "protocol.json") != V3_SHA:
        raise ValueError("unexpected previous core epoch")
    frozen = list((V3 / "inputs").glob("*.json"))
    for path in frozen:
        if path.name != "settings.json" and digest(
            out / "inputs" / path.name
        ) != digest(path):
            raise ValueError("fixture/roster drift: " + path.name)
    area = out / "reference"
    area.mkdir()
    static = area / "inputs"
    static.mkdir()
    source, packages, build = provenance(binary, static)
    fixtures = {f["id"]: f for f in p["fixtures"]}
    for name, f in fixtures.items():
        data_file(static / (name + ".data"), f)
    write_table(static / "oo.table")
    with (static / "outside.table").open("x") as f:
        f.write(
            "# Isolated outside-exclusion bare subtraction, q_i*q_j=-0.32\n\nBARE\nN 2049\n\n"
        )
        for i in range(2049):
            r = 3.6 + i / 2048
            f.write(f"{i + 1} {r:.17g} {0.32 / r:.17g} {0.32 / r**2:.17g}\n")
    ase_input(static)
    import hashlib, json

    def array_hash(value):
        return hashlib.sha256(
            json.dumps(value, separators=(",", ":"), allow_nan=False).encode()
        ).hexdigest()

    states = {
        k: {
            field: array_hash(fixtures[k][field])
            for field in ("positions", "velocities")
        }
        for k in BASES
    }
    if (
        len({r["positions"] for r in states.values()}) != 1
        or len({r["velocities"] for r in states.values()}) != 5
    ):
        raise ValueError("replica stream states are not distinct at fixed q")
    dump(static / "replica-identities.json", states)
    return register(out, binary, fixtures, source, packages, build)


def register(out, binary, fixtures, source, packages, build):
    area, static = out / "reference", out / "reference/inputs"
    exe = out / "bin/stress_validation"
    actions, sections, evaluations, families = {}, {}, {}, {}

    def job(name, argv, cwd, role="rust", **kw):
        actions[name] = dict(
            role=role, argv=list(map(str, argv)), cwd=str(cwd), env=dict(ENV), **kw
        )

    def native(id, rows):
        cwd = area / id
        cwd.mkdir()
        for s in rows:
            s["producer"] = id
            sections[s["tag"]] = s
        script = static / (id + ".in")
        write(
            script, "\n".join(section(fixtures[s["fixture"]], s, static) for s in rows)
        )
        job(
            id,
            [binary, "-echo", "none", "-log", "native.log", "-in", script],
            cwd,
            "reference",
        )

    def snapshot(name, tag, accuracy=None, **kw):
        return dict(
            fixture=name,
            tag=tag,
            accuracy=accuracy,
            mode="snapshot",
            interval=1,
            dt=0.005,
            steps=0,
            burn_steps=0,
            molecules=0,
            **kw,
        )

    def bind(tag, backend, preflight=False):
        s = sections[tag]
        id = f"evaluate-{tag}-{backend}"
        count = 1 + (s["steps"] - s["burn_steps"]) // s["interval"]
        src = area / (tag + ".jsonl")
        job(
            id,
            [
                exe,
                "--child",
                "matched-evaluate",
                out,
                s["fixture"],
                backend,
                src,
                area / (id + ".jsonl"),
            ],
            ROOT,
            expected_frames=count,
        )
        evaluations[id] = dict(tag=tag, backend=backend, preflight=preflight)
        names = ["U", "F", "S"]
        if backend == "zero":
            names += ["reaction", "endpoint"]
        families[id] = [id + "-" + n for n in names]

    for group, names in [("probes", PAIRS), ("water", BASES)]:
        rows = []
        for name in names:
            for i, acc in enumerate(LADDER if name != "pair-oo" else [LADDER[-1]]):
                rows.append(snapshot(name, f"{name}-r{i}", acc))
        if group == "probes":
            rows.append(
                snapshot(
                    "pair-excluded-out",
                    "outside-full-control",
                    LADDER[-1],
                    bare_control=True,
                )
            )
        native("preflight-" + group, rows)
        for name in names:
            tag = name + ("-r0" if name == "pair-oo" else "-r2")
            for backend in ("ewald", "ewald-10", "ewald-12", "spme"):
                bind(tag, backend, True)
            families["converge-" + name] = [
                name + "-" + metric
                for metric in (
                    ["lj-U", "lj-F", "lj-S"]
                    if name == "pair-oo"
                    else ["native-1-2-F", "native-1-2-S"]
                )
                + [
                    f"direct-{lo}-{hi}-{v}"
                    for lo, hi in [(8, 10), (10, 12)]
                    for v in ("F", "S")
                ]
            ]
        if group == "probes":
            families["outside-table"] = ["outside-bare-" + x for x in ("U", "F", "S")]
    for species in ("gas-rotor", "gas-water"):
        id = "preflight-" + species
        native(
            id,
            [
                dict(
                    fixture=species + "-2001",
                    tag=id,
                    mode="nvt",
                    dt=0.005,
                    steps=210000,
                    interval=20,
                    burn_steps=10000,
                    accuracy=None,
                    molecules=32,
                )
            ],
        )
        families[id] = [id + "-" + m for m in GAS_METRICS]
    roster = load(out / "inputs/roster.json")
    for r in roster:
        tag = r["id"].replace("spme-", "ewald-")
        id = "native-" + tag
        if id not in actions:
            native(
                id,
                [
                    dict(
                        r,
                        tag=tag,
                        mode="nvt" if r["seed"] else "nve",
                        accuracy=LADDER[-1],
                        molecules=(1 if r["backend"] == "zero" else 32),
                    )
                ],
            )
        bind(tag, r["backend"])
        name = "trajectory-" + r["id"]
        job(
            name,
            [
                exe,
                "--child",
                "matched-trajectory",
                out,
                r["id"],
                area / (name + ".jsonl"),
            ],
            ROOT,
        )
    for group in ("rotor", "water-nonprincipal", "ewald", "spme"):
        families["nve-" + group] = [
            f"nve-{group}-{key}"
            for key in (
                "endpoint",
                "interval",
                "native-position-ratios",
                "native-velocity-ratios",
                "rust-position-ratios",
                "rust-velocity-ratios",
            )
        ]
    for backend in ("ewald", "spme"):
        families["nvt-" + backend] = [
            f"nvt-{backend}-{t}-{k}" for t in ("endpoint", "interval") for k in COMP
        ]
    job("ase-six-strain", [exe, "--child", "ase-six-strain", out], area)
    actions["ase-six-strain"]["env"].update(ASE_ENV)
    families["ase-six-strain"] = ["ase-wire-conversion", "ase-six-nonzero"] + [
        "ase-strain-" + k for k in COMP
    ]
    x = dict(
        internal_protocol_sha256=digest(out / "protocol.json"),
        previous_internal_protocol_sha256=V3_SHA,
        source_sha256=source,
        snapshot_sha256=inventory(static.iterdir()),
        external_source_revision=REV,
        external_packages=packages,
        native_build=build,
        jobs=actions,
        sections=sections,
        evaluations=evaluations,
        expected_gate_families=families,
        reference_pressure=dict(
            units="lj",
            sign="compressive",
            order=COMP,
            kinetic="LAMMPS endpoint velocities, no COM subtraction",
            constraint="LAMMPS fix native step virial",
            ewald_accuracy=LADDER,
            ewald_gewald=0.8,
            rc=3.5,
            thermostat="rigid/nvt molecule temp 1 1 1 tparam 3 1 3",
        ),
    )
    # Protect the public run_job entry point as well as the adapter scheduler.
    primary = [
        id for id, e in evaluations.items() if e["preflight"] and e["backend"] != "spme"
    ]
    primary += [f"converge-{name}" for name in [*PAIRS, *BASES]] + ["outside-table"]
    mesh = [
        id for id, e in evaluations.items() if e["preflight"] and e["backend"] == "spme"
    ]
    gas = ["preflight-gas-rotor", "preflight-gas-water"]
    for id in mesh:
        actions[id]["requires_families"] = primary
    for r in roster:
        tag = r["id"].replace("spme-", "ewald-")
        native_id = "native-" + tag
        base = primary + (gas if r["seed"] else []) if r["backend"] != "zero" else []
        actions[native_id]["requires_families"] = base
        ev = "evaluate-" + tag + "-" + r["backend"]
        actions[ev]["requires_families"] = base + (
            mesh if r["backend"] == "spme" else []
        )
        actions["trajectory-" + r["id"]]["requires_families"] = [ev]
    validate_jobs(out, x)
    dump(out / "external-epoch.json", x)
    return x
