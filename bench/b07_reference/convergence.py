"""Deterministic h² continuum windows and independent provider strain differences."""

import math
from itertools import zip_longest
from statistics import fmean

from campaign import gate, load, unmeasured
from .native import COMP
from .physics import AXES, bounded, flat, jsonl
from .records import finite, native_rows
from .statistics import validate_trajectory


def intercept(h, y):
    z = [v * v for v in h]
    mu = fmean(z)
    mean = fmean(finite(y, len(z)))
    slope = sum((a - mu) * (b - mean) for a, b in zip(z, y)) / sum(
        (a - mu) ** 2 for a in z
    )
    return finite([mean - slope * mu])[0]


def nve(out, x, fixtures, roster, group):
    jobs = [
        r
        for r in roster
        if not r["seed"] and (r["backend"] == group or r["fixture"] == group)
    ]
    jobs.sort(key=lambda r: -r["dt"])
    if len(jobs) != 4:
        raise ValueError("NVE ladder")
    all_steps = []
    for r in jobs:
        f = fixtures[r["fixture"]]
        tag = r["id"].replace("spme-", "ewald-")
        s = x["sections"][tag]
        means = {k: [0.0] * 6 for k in ("native", "endpoint", "interval")}
        final = {}
        for i, (a, b) in enumerate(
            zip_longest(
                native_rows(out, s, f),
                jsonl(out / "reference" / ("trajectory-" + r["id"] + ".jsonl")),
            )
        ):
            if (
                not a
                or not b
                or a["frame"]["step"] != i
                or b["frame"]["step"] != i
                or b["job"] != r["id"]
            ):
                raise ValueError("NVE grid")
            validate_trajectory(b, f, r["dt"])
            w = 0.5 if i in (0, r["steps"]) else 1.0
            for k in range(6):
                means["native"][k] += a["pressure"][k] * w / r["steps"]
                means["endpoint"][k] += b["endpoint"][k] * w / r["steps"]
                if i:
                    means["interval"][k] += b["interval"][k] / r["steps"]
            if i == r["steps"]:
                final = dict(native=a["frame"], rust=b["frame"])
        if not final:
            raise ValueError("missing NVE terminal state")
        all_steps.append(dict(h=r["dt"], mean=means, final=final))
    h = [s["h"] for s in all_steps]
    zero = {
        kind: [intercept(h, [s["mean"][kind][k] for s in all_steps]) for k in range(6)]
        for kind in ("native", "endpoint", "interval")
    }
    name = "nve-" + group + "-"
    tol = 2e-5 * (1 + math.sqrt(sum(v * v for v in zero["native"])))
    rows = [
        gate(name + k, zero[k], zero["native"], [tol] * 6)
        for k in ("endpoint", "interval")
    ]
    for engine in ("native", "rust"):
        for key, label in (("positions", "position"), ("velocities", "velocity")):
            diff = []
            for lo, hi in zip(all_steps, all_steps[1:]):
                errors = []
                for a, b in zip(lo["final"][engine][key], hi["final"][engine][key]):
                    for k, (aa, bb) in enumerate(zip(a, b)):
                        v = aa - bb
                        if key == "positions":
                            l = f["cell"][k][k]
                            v -= l * round(v / l)
                        errors.append(v)
                diff.append(max(map(abs, finite(errors))))
            ratios = name + engine + "-" + label + "-ratios"
            if min(diff) <= 256 * math.ulp(1.0):
                rows.append(unmeasured(ratios, "Inconclusive"))
            else:
                rows.append(
                    gate(
                        ratios,
                        [diff[0] / diff[1], diff[1] / diff[2]],
                        [4.0, 4.0],
                        [1.0, 1.0],
                    )
                )
    return rows, dict(rungs=all_steps, intercepts=zero)


def ase(out):
    cfg = load(out / "reference/inputs/ase.json")
    attempts = list(jsonl(out / "reference/ase-requests.jsonl"))
    expected = [
        dict(id=i, kind=kind, stage=stage)
        for i, kind in enumerate(["info", "wire-base"] + ["reduced"] * 49, 1)
        for stage in ("attempt", "complete")
    ]
    if attempts != expected:
        raise ValueError("incomplete ASE request journal")
    readings = list(jsonl(out / "reference/ase.jsonl"))
    if len(readings) != 50 or readings[0]["completed_info_requests"] != 1:
        raise ValueError("incomplete ASE requests")
    identity = readings[0]["identity"]
    if (
        identity["ase_version"] != cfg["provider"]["ase_version"]
        or identity["numpy_version"] != "2.2.6"
        or identity["model_sha256"] != cfg["provider"]["model_sha256"]
    ):
        raise ValueError("provider identity")
    energies = []
    for i, (snapshot, row) in enumerate(zip(cfg["snapshots"], readings[1:])):
        if (
            row["index"] != i
            or row["completed_evaluate_requests"] != i + 2
            or any(row[k] != snapshot[k] for k in ("component", "eta"))
            or any(
                max(
                    abs(a - b)
                    for a, b in zip(
                        flat(row[k], len(snapshot[k])),
                        flat(snapshot[k], len(snapshot[k])),
                    )
                )
                > 1e-14
                for k in ("positions", "cell")
            )
        ):
            raise ValueError("provider snapshot/request stream differs")
        finite(row["forces"], 12)
        finite(row["virial"], 6)
        finite([row["energy"], row["trace"], row["volume"]], 3)
        energies.append(row["energy"])
    base = readings[1]
    s = base["virial"]
    wire = base["wire"]
    stress = wire["stress_ev_per_angstrom3"]
    flat(stress, 3)
    unit = 1.602176634e-19 * 6.02214076e23 / 4184
    expected = [-base["volume"] * stress[a][b] * unit for a, b in AXES]
    rows = [
        bounded(
            "ase-wire-conversion",
            s + [base["energy"], base["trace"]] + base["forces"],
            expected
            + [wire["energy_ev"] * unit, sum(expected[:3])]
            + [v * unit for v in flat(wire["forces_ev_per_angstrom"], 4)],
            5e-12,
        ),
        gate(
            "ase-six-nonzero", [float(abs(v) > 1e-8) for v in s], [1.0] * 6, [0.0] * 6
        ),
    ]
    evidence = {"identity": identity, "components": {}}
    for k, label in enumerate(COMP):
        values = []
        for j, eps in enumerate((1e-4, 3e-5, 1e-5, 3e-6)):
            i = 1 + 8 * k + 2 * j
            values.append(-(energies[i + 1] - energies[i]) / (2 * eps))
        finite(values, 4)
        tol = 2e-7 * (1 + abs(s[k]))
        # Any adjacent pair must meet both individual bounds and the mutual-difference bound.
        scores = [
            max(abs(a - s[k]) / tol, abs(b - s[k]) / tol, abs(a - b) / (2 * tol))
            for a, b in zip(values, values[1:])
        ]
        rows.append(gate("ase-strain-" + label, [min(finite(scores, 3))], [0.0], [1.0]))
        evidence["components"][label] = dict(
            derivatives=values,
            expected=s[k],
            tolerance=tol,
            adjacent_normalized_errors=scores,
        )
    return rows, evidence
