"""Frozen native single-run block inference and four-independent-start comparisons."""

import math
import random
from statistics import fmean, variance, stdev
from campaign import gate, unmeasured
from .epoch import COMP, GAS_METRICS
from .records import finite, native_rows
from .physics import AXES, sub, cross, dot, flat, jsonl


def ips(seq):
    seq = finite(seq)
    avg = fmean(seq)
    x = [v - avg for v in seq]
    n = len(x)
    v = dot(x, x) / n
    if n < 2 or v <= 0:
        return None
    rho = 0.0
    for lag in range(0, n - 1, 2):
        pair = (
            sum(sum(x[i] * x[i + k] for i in range(n - k)) / n for k in (lag, lag + 1))
            / v
        )
        if not math.isfinite(pair):
            raise ValueError("nonfinite autocorrelation")
        if pair <= 0:
            break
        rho += pair
    return max(0.5, rho - 0.5)


def block_layout(series):
    counts = {len(v) for v in series.values()}
    if len(counts) != 1 or min(counts) < 4:
        raise ValueError("unbalanced series")
    half = counts.pop() // 2
    tau = {k: ips(v[:half]) for k, v in series.items()}
    usable = [t for t in tau.values() if t is not None]
    size = max(1, math.ceil(10 * max(usable, default=0.5)))
    count = (len(next(iter(series.values()))) - half) // size
    return (
        half,
        size,
        count,
        {k: count * size / (2 * t) if t is not None else 0.0 for k, t in tau.items()},
    )


def t95(df):
    if df < 1:
        raise ValueError("degrees of freedom")
    # t=sqrt(df)*tan(theta) maps its density to a power of cos(theta).
    norm = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(
        math.pi
    )
    lo, hi = 0.0, 16.0
    for _ in range(50):
        t = (hi + lo) / 2
        theta = math.atan(t / math.sqrt(df))
        m = df - 1
        integral = theta if m % 2 == 0 else math.sin(theta)
        c, s = math.cos(theta), math.sin(theta)
        for p in range(2 if m % 2 == 0 else 3, m + 1, 2):
            integral = (s * c ** (p - 1) + (p - 1) * integral) / p
        if 0.5 + norm * integral < 0.975:
            lo = t
        else:
            hi = t
    return (hi + lo) / 2


def chi_quantile(df, prob):
    lo, hi = 0.0, 4.0 * df
    for _ in range(64):
        x = (lo + hi) / 4  # x = chi-square argument / 2
        term = 1.0
        total = term
        for i in range(1, df // 2):
            term *= x / i
            total += term
        if 1 - math.exp(-x) * total < prob:
            lo = 2 * x
        else:
            hi = 2 * x
    return (lo + hi) / 2


def gas_series(out, x, fixture, id):
    s = x["sections"][id]
    metrics = {k: [] for k in GAS_METRICS}
    size = 3 if "water" in id else 2
    qs = {
        name: [chi_quantile(df, p) for p in (0.1, 0.5, 0.9)]
        for name, df in [("com", 96), ("rotation", 32 * size)]
    }
    for row in native_rows(out, s, fixture):
        if row["frame"]["step"] <= s["burn_steps"]:
            continue
        q, v = row["frame"]["positions"], row["frame"]["velocities"]
        p = finite(row["pressure"], 6)
        values = dict(zip(COMP, p), scalar=sum(p[:3]) / 3)
        ec = er = 0.0
        orientation = [0.0] * 6
        for start in range(0, len(v), size):
            m = fixture["masses"][start : start + size]
            vel = v[start : start + size]
            cm = [
                sum(a[i] * mass for a, mass in zip(vel, m)) / sum(m) for i in range(3)
            ]
            ec += sum(m) * dot(cm, cm)
            er += sum(mass * dot(sub(a, cm), sub(a, cm)) for a, mass in zip(vel, m))
            e = sub(q[start + 1], q[start])
            if size == 3:
                e = cross(e, sub(q[start + 2], q[start]))
            norm2 = dot(e, e)
            if norm2 <= 0:
                raise ValueError("zero orientation")
            for k, (a, b) in enumerate(AXES):
                orientation[k] += e[a] * e[b] / (norm2 * 32)
        for k, vv in zip(COMP, orientation):
            values["orientation-" + k] = vv
        for label, energy in [("com", ec), ("rotation", er)]:
            values[label + "-mean"] = energy
            values[label + "-variance"] = energy**2
            for qval, pct in zip(qs[label], (10, 50, 90)):
                values[f"{label}-cdf{pct}"] = float(energy <= qval)
        finite(values.values(), len(metrics))
        for k in metrics:
            metrics[k].append(values[k])
    return metrics


def gas_analysis(out, x, f, id):
    data = gas_series(out, x, f, id)
    half, size, n, neff = block_layout(data)
    rng = random.Random(7001)
    critical = t95(n - 1) if n > 1 else None
    evidence = dict(
        block_samples=size, blocks=n, first_half=half, neff=neff, metrics={}
    )
    output = []
    for k, series in data.items():
        name = id + "-" + k
        target = 0.0
        precision = None
        if k in COMP or k == "scalar":
            target = 1 / 54 if k in ("xx", "yy", "zz", "scalar") else 0.0
            tol = 0.01 / 54
            precision = 0.005 / 54
        elif k.startswith("orientation"):
            target = 1 / 3 if k.split("-")[1] in COMP[:3] else 0.0
            tol = 0.02
        else:
            en, kind = k.split("-")
            df = 64 if en == "rotation" and "rotor" in id else 96
            if kind == "variance":
                target = 2.0 * df
                tol = 0.03 * target
            elif kind == "mean":
                target = float(df)
                tol = 0.01 * target
            else:
                target = int(kind[3:]) / 100
                tol = 0.02
        if n < 2:
            output.append(unmeasured(name, "Inconclusive"))
            continue
        retained = series[half : half + n * size]
        blocks = [fmean(retained[i : i + size]) for i in range(0, len(retained), size)]
        isvar = k.endswith("-variance")
        if isvar:
            raw = data[k.replace("variance", "mean")][half : half + n * size]
            first = [fmean(raw[i : i + size]) for i in range(0, len(raw), size)]
            avg = variance(raw)
            first_mean = fmean(first)
            influence = [(v - 2 * first_mean * b) for v, b in zip(blocks, first)]
        else:
            avg = fmean(blocks)
            influence = blocks
        se = stdev(influence) / math.sqrt(n)
        # Resample whole retained blocks. Variance uses both empirical raw moments.
        boot = []
        for _ in range(2000):
            inds = rng.choices(range(n), k=n)
            stat = fmean(blocks[i] for i in inds)
            if isvar:
                stat = (
                    (stat - fmean(first[i] for i in inds) ** 2)
                    * len(retained)
                    / (len(retained) - 1)
                )
            boot.append(stat)
        finite(boot, 2000)
        boot.sort()
        width = critical * se
        finite([avg, se, width])
        evidence["metrics"][k] = dict(
            mean=avg,
            se=se,
            ci95=[avg - width, avg + width],
            bootstrap95=[boot[49], boot[1949]],
            target=target,
        )
        row = gate(name, [avg], [target], [min(tol, 3 * se)])
        if n < 20 or neff[k] < 1000 or (precision is not None and width > precision):
            row["status"] = "Inconclusive"
        output.append(row)
    return output, evidence


def nvt_series(out, x, f, r, backend):
    id = r["id"].replace("spme-", "ewald-")
    s = x["sections"][id]
    native = native_rows(out, s, f)
    if backend == "native":
        return pressure_columns(
            (a["pressure"] for a in native if a["frame"]["step"] > s["burn_steps"]),
            1000,
        )
    path = out / "reference" / ("trajectory-" + r["id"] + ".jsonl")
    steps = [0] + list(
        range(r["burn_steps"] + r["interval"], r["steps"] + 1, r["interval"])
    )
    from itertools import zip_longest

    endpoint = []
    interval = []
    for step, row in zip_longest(steps, jsonl(path)):
        if (
            step is None
            or not row
            or row["frame"]["step"] != step
            or row["job"] != r["id"]
        ):
            raise ValueError("unbalanced NVT trajectory")
        validate_trajectory(row, f, r["dt"])
        if step:
            if row["window_steps"] != r["interval"]:
                raise ValueError("invalid NVT interval exposure")
            endpoint.append(finite(row["endpoint"], 6))
            interval.append(finite(row["window_interval"], 6))
    return dict(
        endpoint=pressure_columns(endpoint, 1000),
        interval=pressure_columns(interval, 1000),
    )


def validate_trajectory(row, f, dt):
    n = len(f["masses"])
    flat(row["frame"]["positions"], n)
    flat(row["frame"]["velocities"], n)
    finite(row["forces"], 3 * n)
    finite(row["conservative_tensor"], 6)
    finite(row["endpoint"], 6)
    finite([row["potential_energy"], row["frame"]["time"]], 2)
    step = row["frame"]["step"]
    if (
        type(step) != int
        or step < 0
        or abs(row["frame"]["time"] - dt * step) > 1e-10 * (1 + abs(dt * step))
    ):
        raise ValueError("invalid trajectory timestamp")
    if step:
        finite(row["interval"], 6)
        finite(row["window_interval"], 6)
        sample = row.get("sample")
        if n == 96:
            finite(
                [sample[k] for k in ("com_energy", "rotation_energy", "time", "dt")], 4
            )
            for key in ("endpoint", "interval", "orientation_second"):
                finite(sample[key], 6)


def pressure_columns(rows, n):
    columns = {k: [] for k in COMP}
    for row in rows:
        for k, v in zip(COMP, finite(row, 6)):
            columns[k].append(v)
    if len(columns["xx"]) != n:
        raise ValueError("wrong pressure sample count")
    return columns


def interaction(out, x, fixtures, roster, backend):
    slots = {}
    details = {}
    for r in (j for j in roster if j["backend"] == backend and j["seed"]):
        f = fixtures[r["fixture"]]
        series = nvt_series(out, x, f, r, backend)
        series["native"] = nvt_series(out, x, f, r, "native")
        sample = {t + "-" + k: v for t, d in series.items() for k, v in d.items()}
        half, size, n, neff = block_layout(sample)
        means = (
            {k: fmean(v[half : half + n * size]) for k, v in sample.items()}
            if n
            else {}
        )
        slots[(r["seed"], r["dt"])] = (means, neff, n)
        details[r["id"]] = dict(
            block_samples=size, retained_blocks=n, neff=neff, means=means
        )
    rows = []
    for timing in ("endpoint", "interval"):
        for k in COMP:
            name = f"nvt-{backend}-{timing}-{k}"
            if any(slot[2] < 2 for slot in slots.values()):
                rows.append(unmeasured(name, "Inconclusive"))
                continue

            def intercept(seed, method):
                c = slots[(seed, 0.01)][0][method + "-" + k]
                f = slots[(seed, 0.005)][0][method + "-" + k]
                return (4 * f - c) / 3

            rust = [intercept(s, timing) for s in range(2001, 2005)]
            native = [intercept(s, "native") for s in range(2001, 2005)]
            diff = sub(rust, native)
            se = (
                stdev(diff) / 2
            )  # shared initial conditions and paired h, four independent starts
            a, b = fmean(rust), fmean(native)
            width = t95(3) * se
            native_width = t95(3) * stdev(native) / 2
            rust_width = t95(3) * stdev(rust) / 2
            tol = 0.01 * (max(1 / 54, abs(b)) if k in COMP[:3] else 1 / 54)
            row = gate(name, [a], [b], [min(tol, 3 * se)])
            pooled = min(
                sum(slots[(s, h)][1][m + "-" + k] for s in range(2001, 2005))
                for h in (0.01, 0.005)
                for m in (timing, "native")
            )
            if (
                min(v[2] for v in slots.values()) < 20
                or pooled < 500
                or max(width, rust_width, native_width) > 0.01 / 54
            ):
                row["status"] = "Inconclusive"
            details[name] = dict(
                rust=rust,
                native=native,
                paired_se=se,
                halfwidth=width,
                native_halfwidth=native_width,
                rust_halfwidth=rust_width,
                minimum_pooled_neff=pooled,
                estimated_precision_work_multiplier=max(
                    1.0,
                    500 / pooled if pooled else 1.0,
                    (max(width, native_width, rust_width) / (0.01 / 54)) ** 2,
                )
                if pooled
                else None,
            )
            rows.append(row)
    return rows, details
