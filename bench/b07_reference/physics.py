"""Independent rigid-body and table oracles plus strictly paired raw comparison."""

import math
from itertools import zip_longest
from campaign import gate, loads, load
from .records import finite, native_rows, reciprocal_settings

AXES = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def flat(rows, n):
    if len(rows) != n:
        raise ValueError("wrong number of atom triples")
    return [v for row in rows for v in finite(row, 3)]


def vectors(a, b, op):
    return [op(x, y) for x, y in zip(finite(a), finite(b, len(a)))]


def sub(a, b):
    return vectors(a, b, lambda x, y: x - y)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def solve(a, b):
    # Reciprocal basis of the symmetric 3x3 inertia tensor (Cramer's rule).
    cof = [cross(a[(i + 1) % 3], a[(i + 2) % 3]) for i in range(3)]
    det = dot(a[0], cof[0])
    if not math.isfinite(det) or det <= 0:
        raise ValueError("singular rigid-body inertia")
    return finite([dot(c, b) / det for c in cof], 3)


def reaction(f, raw):
    m = f["masses"]
    n = len(m)
    q = raw["frame"]["positions"]
    v = raw["frame"]["velocities"]
    flat(q, n)
    flat(v, n)
    force = finite(raw["forces"], 3 * n)
    force = [force[3 * i : 3 * i + 3] for i in range(n)]
    # Native dump exports unwrapped atom coordinates. Unwrap each molecule about its first atom.
    edges = [f["cell"][k][k] for k in range(3)]
    r = [[d - l * round(d / l) for d, l in zip(sub(row, q[0]), edges)] for row in q]
    center = [sum(mi * a[k] for mi, a in zip(m, r)) / sum(m) for k in range(3)]
    r = [sub(a, center) for a in r]
    if n == 2:
        d = sub(r[1], r[0])
        dv = sub(v[1], v[0])
        da = [force[1][k] / m[1] - force[0][k] / m[0] for k in range(3)]
        lam = (dot(dv, dv) + dot(d, da)) / (dot(d, d) * (1 / m[0] + 1 / m[1]))
        fc = [[d[k] * lam * sign for k in range(3)] for sign in (1, -1)]
    elif n == 3:
        inertia = [
            [
                sum(
                    mi * ((dot(a, a) if i == j else 0) - a[i] * a[j])
                    for mi, a in zip(m, r)
                )
                for j in range(3)
            ]
            for i in range(3)
        ]
        angular = [
            sum(mi * cross(a, b)[k] for mi, a, b in zip(m, r, v)) for k in range(3)
        ]
        omega = solve(inertia, angular)
        torque = [sum(cross(a, b)[k] for a, b in zip(r, force)) for k in range(3)]
        alpha = solve(inertia, sub(torque, cross(omega, angular)))
        acm = [sum(row[k] for row in force) / sum(m) for k in range(3)]
        fc = [
            [
                mi * (acm[k] + cross(alpha, a)[k] + cross(omega, cross(omega, a))[k])
                - external[k]
                for k in range(3)
            ]
            for mi, a, external in zip(m, r, force)
        ]
    else:
        raise ValueError("single rotor/water required")
    volume = math.prod(edges)
    p = [
        (
            s
            + sum(
                mi * vel[i] * vel[j] + 0.5 * (a[i] * b[j] + a[j] * b[i])
                for mi, vel, a, b in zip(m, v, r, fc)
            )
        )
        / volume
        for (i, j), s in zip(AXES, finite(raw["conservative_tensor"], 6))
    ]
    return flat(fc, n), finite(p, 6)


def jsonl(path):
    with path.open("rb") as file:
        for line in file:
            yield loads(line)


def residuals(actual, expected, n, scaled=True):
    a, b = finite(actual, n), finite(expected, n)
    return finite([abs(x - y) / (1 + abs(y) if scaled else 1) for x, y in zip(a, b)], n)


def compare(out, epoch, fixtures, id):
    e = epoch["evaluations"][id]
    s = epoch["sections"][e["tag"]]
    f = fixtures[s["fixture"]]
    n = len(f["masses"])
    keys = {"U": 1, "F": 3 * n, "S": 6}
    if e["backend"] == "zero":
        keys.update(reaction=3 * n, endpoint=6)
    maxima = {k: [0.0] * size for k, size in keys.items()}
    src = (r for r in native_rows(out, s, f) if keep(s, r["frame"]["step"]))
    count = 0
    if s["mode"] != "snapshot" and any(f["charges"]):
        log = out / "reference" / s["producer"] / "native.log"
        registered = (
            out / "reference/comparisons" / ("converge-" + s["fixture"] + ".json")
        )
        settings = load(registered)["details"]["native_reciprocal"][-1]
        if reciprocal_settings(log, s["tag"]) != settings:
            raise ValueError("native dynamics differs from preflight Ewald settings")
    for native, r in zip_longest(src, jsonl(out / "reference" / (id + ".jsonl"))):
        if (
            not native
            or not r
            or r["fixture"] != f["id"]
            or r["backend"] != e["backend"]
        ):
            raise ValueError("incorrect matched dataset")
        if e["backend"] == "spme" and r.get("spme_mesh") != [32, 32, 32]:
            raise ValueError("unexpected matched SPME mesh")
        if e["backend"] != "zero":
            cutoff = {"ewald": 8.0, "ewald-10": 10.0, "ewald-12": 12.0, "spme": 8.0}[
                e["backend"]
            ]
            if r.get("reciprocal_cutoff") != cutoff:
                raise ValueError("unexpected matched reciprocal cutoff")
        ensure_frame(r["frame"], native["frame"], n)
        pairs = {
            "U": ([r["potential_energy"]], [native["potential_energy"]]),
            "F": (r["forces"], native["forces"]),
            "S": (r["conservative_tensor"], native["conservative_tensor"]),
        }
        finite(r["endpoint"], 6)
        if e["backend"] == "zero":
            fc, p = reaction(f, native)
            pairs.update(
                reaction=(flat(r["production_reaction_forces"], n), fc),
                endpoint=(r["endpoint"], p),
            )
            oracle_f, oracle_p = r["analytic_force_pressure"]
            flat(oracle_f, n)
            finite(oracle_p, 6)
        for k, (a, b) in pairs.items():
            errs = residuals(a, b, keys[k])
            maxima[k] = list(map(max, maxima[k], errs))
        count += 1
    if count != epoch["jobs"][id]["expected_frames"]:
        raise ValueError("wrong number of paired snapshots")
    return [
        gate(
            id + "-" + k,
            m,
            [0.0] * len(m),
            [2e-9 if k in ("endpoint", "reaction") else 2e-6] * len(m),
        )
        for k, m in maxima.items()
    ]


def keep(s, step):
    return step == 0 or step > s["burn_steps"]


def ensure_frame(a, b, n):
    if a["step"] != b["step"] or abs(finite([a["time"]])[0] - b["time"]) > 1e-10 * (
        1 + abs(b["time"])
    ):
        raise ValueError("matched frame time")
    for k in ("positions", "velocities"):
        if max(residuals(flat(a[k], n), flat(b[k], n), 3 * n)) > 1e-13:
            raise ValueError("matched frame geometry " + k)


def bounded(name, a, b, tol, scaled=True):
    finite(a, len(b))
    finite(b, len(a))
    return gate(name, a, b, [tol * (1 + abs(x)) if scaled else tol for x in b])


def snap_checks(out, x, fix, name):
    def row(tag):
        rows = list(native_rows(out, x["sections"][tag], fix[name]))
        if len(rows) != 1:
            raise ValueError("snapshot count")
        return rows[0]

    fine = name + ("-r0" if name == "pair-oo" else "-r2")
    f = row(fine)
    result = []
    if name == "pair-oo":
        d = sub(f["frame"]["positions"][0], f["frame"]["positions"][1])
        r = math.sqrt(dot(d, d))

        def energy(a):
            return 0.4 * (a**-12 - a**-6)

        def deriv(a):
            return 0.4 * (-12 * a**-13 + 6 * a**-7)

        u = energy(r) - energy(2.5) - (r - 2.5) * deriv(2.5)
        fi = [-v * (deriv(r) - deriv(2.5)) / r for v in d]
        result += [
            bounded(name + "-lj-U", [f["potential_energy"]], [u], 2e-7, False),
            bounded(name + "-lj-F", f["forces"], fi + [-v for v in fi], 2e-7, False),
            bounded(
                name + "-lj-S",
                f["conservative_tensor"],
                [d[i] * fi[j] for i, j in AXES],
                2e-7,
                False,
            ),
        ]
    else:
        # Fully parse all three native rungs, including the coarse rung.
        row(name + "-r0")
        m = row(name + "-r1")
        result += [
            bounded(name + "-native-1-2-F", m["forces"], f["forces"], 2e-6, False),
            bounded(
                name + "-native-1-2-S",
                m["conservative_tensor"],
                f["conservative_tensor"],
                2e-6,
            ),
        ]
    for l, h in [(8, 10), (10, 12)]:

        def direct(k):
            label = "ewald" if k == 8 else "ewald-" + str(k)
            rows = list(jsonl(out / "reference" / f"evaluate-{fine}-{label}.jsonl"))
            if len(rows) != 1 or rows[0]["reciprocal_cutoff"] != k:
                raise ValueError("direct kcut mismatch")
            return rows[0]

        lo, hi = direct(l), direct(h)
        for k, field, scaled in [
            ("F", "forces", False),
            ("S", "conservative_tensor", True),
        ]:
            result.append(
                bounded(
                    f"{name}-direct-{l}-{h}-{k}", lo[field], hi[field], 2e-6, scaled
                )
            )
    details = {}
    if name != "pair-oo":
        details["native_reciprocal"] = [
            reciprocal_settings(
                out
                / "reference"
                / x["sections"][f"{name}-r{i}"]["producer"]
                / "native.log",
                f"{name}-r{i}",
            )
            for i in range(3)
        ]
    return result, details


def outside_table(out, x, fixtures):
    f = fixtures["pair-excluded-out"]

    def row(tag):
        return list(native_rows(out, x["sections"][tag], f))[0]

    corr, base = row("pair-excluded-out-r2"), row("outside-full-control")
    d = sub(f["positions"][0], f["positions"][1])
    r = math.sqrt(dot(d, d))
    u = -math.prod(f["charges"]) / r
    fi = [u * v / r**2 for v in d]
    s = [d[i] * fi[j] for i, j in AXES]
    return [
        bounded("outside-bare-" + k, sub(a, b), c, 2e-7, False)
        for k, a, b, c in [
            ("U", [corr["potential_energy"]], [base["potential_energy"]], [u]),
            ("F", corr["forces"], base["forces"], fi + [-v for v in fi]),
            ("S", corr["conservative_tensor"], base["conservative_tensor"], s),
        ]
    ]
