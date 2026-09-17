#!/usr/bin/env python3
"""Bond-length residual analysis for the SHAKE validation decks.

Reads a trajectory (LAMMPS custom dump or MODENA XYZ trajectory), computes
the three constraint distances per trimer (2 O-H bonds + the H-H 1-3
distance) with minimum-image wrapping in the 10x10x10 box, and reports the
RMS and max deviation |r - d| from the constraint targets.

Usage:
  python3 analyze_shake.py lammps results/lammps_trimers.lammpstrj
  python3 analyze_shake.py modena results/modena_traj.xyz
"""

import math
import sys

BOX = 10.0
HH = math.sqrt(8.0 / 3.0)  # 1.6329931618554521
# (i, j, d) with 0-based indices: O-H1, O-H2, H1-H2 per trimer.
CONSTRAINTS = []
for m in range(16):
    o, h1, h2 = 3 * m, 3 * m + 1, 3 * m + 2
    CONSTRAINTS.append((o, h1, 1.0))
    CONSTRAINTS.append((o, h2, 1.0))
    CONSTRAINTS.append((h1, h2, HH))


def min_image(d):
    return d - BOX * round(d / BOX)


def frames_lammps(path):
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(lines[i + 1])
            n = int(lines[i + 3])
            # ITEM: ATOMS id x y z
            pos = {}
            for j in range(n):
                parts = lines[i + 9 + j].split()
                pos[int(parts[0]) - 1] = (float(parts[1]), float(parts[2]), float(parts[3]))
            yield step, pos
            i += 9 + n
        else:
            i += 1


def frames_modena(path):
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i])
        step = int(lines[i + 1].split("step=")[1].split()[0])
        pos = {}
        for j in range(n):
            parts = lines[i + 2 + j].split()
            pos[j] = (float(parts[1]), float(parts[2]), float(parts[3]))
        yield step, pos
        i += 2 + n


def residuals(pos):
    worst = 0.0
    for (a, b, d) in CONSTRAINTS:
        pa, pb = pos[a], pos[b]
        dx = min_image(pa[0] - pb[0])
        dy = min_image(pa[1] - pb[1])
        dz = min_image(pa[2] - pb[2])
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        worst = max(worst, abs(r - d))
    return worst


def main():
    kind, path = sys.argv[1], sys.argv[2]
    it = frames_lammps(path) if kind == "lammps" else frames_modena(path)
    n_samples = 0
    sum_sq = 0.0
    worst = 0.0
    nframes = 0
    for step, pos in it:
        # Skip step 0 (pre-dynamics, constraints exactly satisfied by
        # construction); measure the post-dynamics frames.
        if step == 0:
            continue
        nframes += 1
        for (a, b, d) in CONSTRAINTS:
            pa, pb = pos[a], pos[b]
            dx = min_image(pa[0] - pb[0])
            dy = min_image(pa[1] - pb[1])
            dz = min_image(pa[2] - pb[2])
            r = math.sqrt(dx * dx + dy * dy + dz * dz)
            dev = abs(r - d)
            sum_sq += dev * dev
            n_samples += 1
            worst = max(worst, dev)
    rms = math.sqrt(sum_sq / n_samples)
    print(f"frames: {nframes}  bond samples: {n_samples}")
    print(f"bond |r-d| RMS : {rms:.6e}")
    print(f"bond |r-d| max : {worst:.6e}")


if __name__ == "__main__":
    main()
