"""Generate LAMMPS input from the immutable MODENA arrays (no native evaluation)."""

import math
from pathlib import Path

COMP = ("xx", "yy", "zz", "xy", "xz", "yz")
LADDER = ("1e-7", "1e-9", "1e-11")
N_TABLE = 32769
NATIVE_BUILD = "Homebrew serial bottle, CMake all_on/nolib, BUILD_MPI=no BUILD_OMP=no BUILD_SHARED_LIBS=yes FFT=FFTW3"


def write(path, data):
    with Path(path).open("x") as f:
        f.write(data)


def lj(r):
    def raw(d):
        z = d**-6
        return 0.4 * z * (z - 1), 2.4 * z * (2 * z - 1) / d

    u, force = raw(r)
    uc, fc = raw(2.5)
    return u - uc + (r - 2.5) * fc, force - fc


def write_table(path):
    with path.open("x") as f:
        f.write(
            f"# Reduced shifted-force O-O LJ. epsilon=.1, sigma=1\n\nOO\nN {N_TABLE} RSQ 0.5 2.5\n\n"
        )
        for i in range(N_TABLE):
            r = math.sqrt(0.25 + 6 * i / (N_TABLE - 1))
            u, force = lj(r)
            f.write(f"{i + 1} {r:.17g} {u:.17g} {force:.17g}\n")


def data_file(path, fixture):
    masses = fixture["masses"]
    types = list(dict.fromkeys(masses))
    bonds, angles = [], []
    mol = list(range(1, len(masses) + 1))
    for b in fixture["constraints"]["bonds"]:
        i, j = b["i"], b["j"]
        bonds.append((i + 1, j + 1))
        mol[j] = mol[i]
    for water in fixture["constraints"]["waters"]:
        i = water["first"]
        bonds.extend([(i + 1, i + 2), (i + 1, i + 3)])
        angles.append((i + 2, i + 1, i + 3))
        mol[i + 1] = mol[i + 2] = mol[i]
    if fixture["id"].startswith("pair-excluded"):
        bonds.append((1, 2))
        mol[1] = mol[0]
    lines = [
        "Frozen " + fixture["id"],
        "",
        f"{len(masses)} atoms",
        f"{len(bonds)} bonds",
        f"{len(angles)} angles",
        f"{len(types)} atom types",
        f"{int(bool(bonds))} bond types",
        f"{int(bool(angles))} angle types",
        "",
    ]
    for k, a in enumerate("xyz"):
        lines.append(f"0 {fixture['cell'][k][k]:.17g} {a}lo {a}hi")
    lines.extend(["", "Masses", ""])
    lines.extend(f"{i + 1} {m:.17g}" for i, m in enumerate(types))
    lines.extend(["", "Atoms # full", ""])
    lengths = [fixture["cell"][k][k] for k in range(3)]
    for i, pos in enumerate(fixture["positions"]):
        images = [math.floor(v / l) for v, l in zip(pos, lengths)]
        wrapped = [v - img * l for v, img, l in zip(pos, images, lengths)]
        lines.append(
            f"{i + 1} {mol[i]} {1 + types.index(masses[i])} {fixture['charges'][i]:.17g} "
            + " ".join(format(v, ".17g") for v in wrapped)
            + " "
            + " ".join(map(str, images))
        )
    lines.extend(["", "Velocities", ""])
    lines.extend(
        f"{i + 1} " + " ".join(format(v, ".17g") for v in vel)
        for i, vel in enumerate(fixture["velocities"])
    )
    for name, rows in (("Bonds", bonds), ("Angles", angles)):
        if rows:
            lines.extend(["", name, ""])
            lines.extend(
                f"{i + 1} 1 " + " ".join(map(str, row)) for i, row in enumerate(rows)
            )
    write(path, "\n".join(lines) + "\n")


def section(f, config, static):
    """A clear separates native preflight snapshots within the same admitted process."""
    name, tag = f["id"], config["tag"]
    waters = f["constraints"]["waters"]
    bonds = bool(waters or f["constraints"]["bonds"] or f["pair_scales"])
    charged = any(f["charges"])
    commands = [
        "clear",
        f'print "B07_BEGIN {tag}"',
        "units lj",
        "atom_style full",
        "boundary p p p",
        "bond_style zero",
        "angle_style zero",
        "special_bonds lj/coul 0 0 0",
        f'read_data "{static}/{name}.data"',
    ]
    if name == "pair-excluded-out":
        commands[commands.index("bond_style zero")] = "bond_style table linear 2049"
        commands[commands.index("special_bonds lj/coul 0 0 0")] = (
            "special_bonds lj 0 0 0 coul 1 1 1"
        )
        if config.get("bare_control"):
            commands[commands.index("bond_style table linear 2049")] = "bond_style zero"
            commands.append("bond_coeff 1 4.1")
        else:
            commands.append(f'bond_coeff 1 "{static}/outside.table" BARE')
    elif bonds:
        length = (
            waters[0]["d_oh"]
            if waters
            else (f["constraints"]["bonds"] or [{"d": 1.0}])[0]["d"]
        )
        commands.append(f"bond_coeff 1 {length:.17g}")
    if waters:
        w = waters[0]
        angle = math.degrees(2 * math.asin(w["d_hh"] / (2 * w["d_oh"])))
        commands.append(f"angle_coeff 1 {angle:.17g}")
    if charged:
        commands += [
            f"pair_style hybrid/overlay coul/long 3.5 table linear {N_TABLE}",
            "pair_coeff * * coul/long",
            f'pair_coeff 1 1 table "{static}/oo.table" OO 2.5',
            "pair_modify pair coul/long table 0",
            f"kspace_style ewald {config['accuracy']}",
            "kspace_modify gewald 0.8",
        ]
    elif name == "pair-oo":
        commands += [
            f"pair_style table linear {N_TABLE}",
            f'pair_coeff * * "{static}/oo.table" OO 2.5',
        ]
    else:
        commands += ["pair_style zero 3.5", "pair_coeff * *"]
    commands += [
        "neighbor 0.3 bin",
        "neigh_modify every 1 delay 0 check yes",
        "comm_modify cutoff 5.0",
        "fix cons all store/force",
    ]
    mode = config["mode"]
    if mode == "nve":
        commands += [
            "fix integrate all nve",
            "fix rigidconstraint all rattle 1e-13 100 0 b 1"
            + (" a 1" if waters else ""),
        ]
    elif mode == "nvt":
        commands.append("fix integrate all rigid/nvt molecule temp 1 1 1 tparam 3 1 3")
    commands += [
        "compute T all temp",
        "compute_modify T extra/dof 0",
        "compute K all pressure T ke",
        "compute C all pressure NULL pair bond angle kspace",
        "compute R all pressure NULL fix",
        "compute P all pressure T",
    ]
    fields = " ".join(f"c_{kind}[{k}]" for kind in "PKCR" for k in range(1, 7))
    stride = config["interval"]
    commands += [
        f"thermo {stride}",
        f"thermo_style custom step time pe {fields}",
        "thermo_modify norm no format float %.17g flush yes",
        f"dump atoms all custom {stride} {tag}.dump id xu yu zu vx vy vz f_cons[1] f_cons[2] f_cons[3]",
        "dump_modify atoms sort id format float %.17g flush yes",
        f"timestep {config['dt']:.17g}",
        f"run {config['steps']}",
        f'print "B07_END {tag}"',
    ]
    return "\n".join(commands) + "\n"
