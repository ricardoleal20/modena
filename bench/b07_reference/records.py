"""Strict readers for actual native and matched outputs; fail closed on incomplete data."""

import math
import re
from itertools import zip_longest


def finite(values, size=None):
    v = list(values)
    if (size is not None and len(v) != size) or any(
        type(x) not in (float, int) or not math.isfinite(x) for x in v
    ):
        raise ValueError("invalid numeric component/shape")
    return v


def numbers(fields, size):
    return finite(map(float, fields), size)


def parse_pressure_line(line):
    fields = line.split()
    if len(fields) != 8 or not fields[0].isdigit():
        raise ValueError("expected step time xx yy zz xy xz yz")
    t, *p = numbers(fields[1:], 7)
    return {"step": int(fields[0]), "time": t, "pressure": p}


HEADER = ["Step", "Time", "PotEng"] + [
    f"c_{kind}[{k}]" for kind in "PKCR" for k in range(1, 7)
]
ATOMS = "ITEM: ATOMS id xu yu zu vx vy vz f_cons[1] f_cons[2] f_cons[3]"


def thermo(path, tag, partial=False):
    active = reading = ended = False
    with path.open() as f:
        for line in f:
            s = line.strip()
            if s == "B07_BEGIN " + tag:
                if active or ended:
                    raise ValueError("duplicate native section")
                active = True
            elif active and s == "B07_END " + tag:
                ended, active = True, False
            elif active and s.split() == HEADER:
                if reading:
                    raise ValueError("duplicate thermo header")
                reading = True
            elif active and reading and s.startswith("Loop time"):
                reading = False
            elif active and reading and s and s.split()[0].isdigit():
                fields = s.split()
                if partial and not line.endswith("\n"):
                    break
                vals = numbers(fields[1:], 26)
                canonical = parse_pressure_line(" ".join(fields[:2] + fields[3:9]))
                yield dict(
                    **canonical,
                    potential_energy=vals[1],
                    kinetic=vals[8:14],
                    conservative=vals[14:20],
                    constraint=vals[20:26],
                )
    if not ended and not partial:
        raise ValueError("missing native section end " + tag)


def dumps(path, fixture):
    n = len(fixture["masses"])
    with path.open() as f:
        while True:
            s = f.readline()
            if not s:
                break
            if s.strip() != "ITEM: TIMESTEP":
                raise ValueError("dump timestep header")
            step = int(f.readline())
            if (
                f.readline().strip() != "ITEM: NUMBER OF ATOMS"
                or int(f.readline()) != n
            ):
                raise ValueError("dump atom count")
            if f.readline().strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError("dump orthorhombic PBC")
            for k in range(3):
                if numbers(f.readline().split(), 2) != [0.0, fixture["cell"][k][k]]:
                    raise ValueError("dump box differs")
            if f.readline().strip() != ATOMS:
                raise ValueError("dump atom columns")
            q, v, force = [], [], []
            for i in range(n):
                row = f.readline().split()
                if not row or row[0] != str(i + 1):
                    raise ValueError("dump IDs/sort")
                a = numbers(row[1:], 9)
                q.append(a[:3])
                v.append(a[3:6])
                force.extend(a[6:])
            yield dict(step=step, positions=q, velocities=v, forces=force)


def native_rows(out, section, fixture):
    folder = out / "reference" / section["producer"]
    grid = range(0, section["steps"] + 1, section["interval"])
    volume = math.prod(fixture["cell"][i][i] for i in range(3))
    for step, t, a in zip_longest(
        grid,
        thermo(folder / "native.log", section["tag"]),
        dumps(folder / (section["tag"] + ".dump"), fixture),
    ):
        if step is None or not t or not a or step != t["step"] or step != a["step"]:
            raise ValueError("native frame count/order differs from frozen grid")
        if abs(t["time"] - step * section["dt"]) > 1e-10 * (1 + abs(t["time"])):
            raise ValueError("native time differs from grid")
        if step == 0:
            for key in ("positions", "velocities"):
                for row, ref in zip(a[key], fixture[key]):
                    if any(abs(x - y) > 1e-10 for x, y in zip(row, ref)):
                        raise ValueError("native initial state differs: " + key)
        if any(
            abs(p - k - c - r) > 1e-10 * (1 + abs(p) + abs(k) + abs(c) + abs(r))
            for p, k, c, r in zip(
                t["pressure"], t["kinetic"], t["conservative"], t["constraint"]
            )
        ):
            raise ValueError("native pressure split")
        yield dict(
            frame=dict(
                step=step,
                time=t["time"],
                positions=a["positions"],
                velocities=a["velocities"],
            ),
            forces=a["forces"],
            potential_energy=t["potential_energy"],
            pressure=t["pressure"],
            conservative_tensor=finite([x * volume for x in t["conservative"]], 6),
        )


def confirmed_work(log, sections, clean):
    """A native thermo step is evidence even when its section later aborts."""
    steps = molecules = 0
    complete = bool(clean)
    for section in sections:
        last = -1
        if log.exists():
            try:
                for row in thermo(log, section["tag"], partial=True):
                    step = row["step"]
                    if (
                        step != (0 if last == -1 else last + section["interval"])
                        or step > section["steps"]
                    ):
                        raise ValueError("invalid native step prefix")
                    last = step
            except ValueError:
                complete = False
        steps += max(0, last)
        molecules += max(0, last) * section["molecules"]
        complete &= last == section["steps"]
    return dict(system_steps=steps, molecule_steps=molecules, complete=complete)


def reciprocal_settings(path, tag):
    values = {}
    active = False
    with path.open() as file:
        for line in file:
            line = line.strip()
            if line == "B07_BEGIN " + tag:
                active = True
            elif line == "B07_END " + tag:
                break
            elif active and "=" in line:
                key, value = (s.strip() for s in line.split("=", 1))
                if key == "G vector (1/distance)":
                    values["gewald"] = numbers(value.split(), 1)[0]
                if key in ("KSpace vectors: actual max1d max3d", "kxmax kymax kzmax"):
                    if not re.fullmatch(r"[0-9]+\s+[0-9]+\s+[0-9]+", value):
                        raise ValueError("invalid Ewald vector counts")
                    values[key] = [int(v) for v in value.split()]
    if (
        len(values) != 3
        or values["gewald"] != 0.8
        or min(values["kxmax kymax kzmax"]) <= 0
    ):
        raise ValueError("unconfirmed native Ewald settings")
    return values
