"""Typed mirrors of MODENA's TOML schema (`InputConfig`).

`Config.to_toml()` renders input for the Rust `InputConfig` validator, which
is the source of truth for validation and unit conversion.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from typing import Any, Optional, Union


def _scalar(v: Any) -> str:
    import json
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v).lower()
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False).replace(chr(127), "\\u007f")
    raise TypeError(f"unhandled scalar {type(v)}: {v!r}")


def _val(v: Any) -> str:
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_val(x) for x in v) + "]"
    return _scalar(v)


def _mapping(obj: Any) -> dict:
    if is_dataclass(obj):
        return {f.name: getattr(obj, f.name) for f in fields(obj)}
    return obj


def _render_section(obj: Any, path: str) -> list[str]:
    import json
    leaves, nested = [], []
    for name, v in _mapping(obj).items():
        if v is None or v == [] or v == {}:
            continue
        child = path + "." + name
        if isinstance(v, list) and all(is_dataclass(x) or isinstance(x, dict) for x in v):
            for item in v:
                nested += [f"[[{child}]]"] + _render_section(item, child)
        elif is_dataclass(v) or isinstance(v, dict):
            nested += [f"[{child}]"] + _render_section(v, child)
        else:
            leaves.append(f"{json.dumps(name)} = {_val(v)}" if not (name.isascii() and name.isidentifier()) else f"{name} = {_val(v)}")
    return leaves + nested


############################### structure ###############################

@dataclass
class Fcc:
    lattice_constant: float
    replicas: list[int]
    element: str = "Ar"
    type: str = "fcc"
    mass: Optional[float] = None
    charge: Optional[float] = None


@dataclass
class Xyz:
    file: str
    box_size: list[float]
    element: str = "X"
    type: str = "xyz"
    mass: Optional[float] = None
    charge: Optional[float] = None


############################### model ###############################

@dataclass
class LennardJones:
    cutoff: float
    sigma: Optional[float] = None
    epsilon: Optional[float] = None
    types: list[dict] = field(default_factory=list)
    overrides: list[dict] = field(default_factory=list)
    type: str = "lennard_jones"
    mixing: str = "lorentz_berthelot"


@dataclass
class Buckingham:
    cutoff: float
    a: Optional[float] = None
    b: Optional[float] = None
    rho: Optional[float] = None
    types: list[dict] = field(default_factory=list)
    overrides: list[dict] = field(default_factory=list)
    type: str = "buckingham"


@dataclass
class Morse:
    cutoff: float
    d: Optional[float] = None
    alpha: Optional[float] = None
    r0: Optional[float] = None
    types: list[dict] = field(default_factory=list)
    overrides: list[dict] = field(default_factory=list)
    type: str = "morse"


########################################################################

@dataclass
class VelocityVerlet:
    dt: float
    type: str = "velocity_verlet"


@dataclass
class Langevin:
    dt: float
    gamma: float
    target_temperature: float
    seed: int = 42
    type: str = "langevin_baoab"


@dataclass
class Run:
    steps: int = 100
    initial_temperature: float = 0.8
    seed: int = 42
    neighbor_skin: float = 0.3
    trajectory_file: Optional[str] = None
    trajectory_every: int = 100
    log_file: Optional[str] = None
    log_every: int = 10
    restore_file: Optional[str] = None
    restart_mode: str = "continue"
    restart_file: Optional[str] = None
    restart_every: int = 500


@dataclass
class RockSalt:
    lattice_constant: float
    replicas: list[int]
    cation: str
    anion: str
    cation_charge: float
    anion_charge: float
    cation_mass: Optional[float] = None
    anion_mass: Optional[float] = None
    type: str = "rock_salt"


@dataclass
class Water:
    n_waters: int
    density: float = 0.997
    d_oh: float = 1.0
    hoh_angle: float = 109.47
    q_o: float = -0.8476
    q_h: float = 0.4238
    m_o: float = 15.9994
    m_h: float = 1.008
    seed: int = 42
    type: str = "water"


@dataclass
class RandomAlloy:
    lattice_constant: float
    replicas: list[int]
    species: list[dict]
    seed: int = 42
    type: str = "random_alloy"


@dataclass
class Tabulated:
    file: str
    interpolation: str = "linear"
    type: str = "tabulated"


@dataclass
class SoftSphere:
    cutoff: float
    sigma: Optional[float] = None
    epsilon: Optional[float] = None
    n: Optional[int] = None
    types: list[dict] = field(default_factory=list)
    overrides: list[dict] = field(default_factory=list)
    type: str = "soft_sphere"


@dataclass
class NoseHooverChain:
    dt: float
    target_temperature: float
    tau: float = 0.5
    chain_len: int = 4
    type: str = "nose_hoover_chain"


@dataclass
class Csvr:
    target_temperature: float
    tau: float = 0.5
    seed: int = 42
    type: str = "velocity_rescaling"


@dataclass
class Berendsen:
    target_temperature: float
    tau: float = 0.5
    type: str = "berendsen"


@dataclass
class NoThermostat:
    type: str = "none"


@dataclass
class Ewald:
    kappa: float
    real_cutoff: float
    kspace_cutoff: float
    method: str = "ewald"


@dataclass
class Bond:
    atoms: list[int]
    r0: float
    k: float


@dataclass
class Angle:
    atoms: list[int]
    theta0: float
    k: float


@dataclass
class DihedralTerm:
    n: int
    k: float
    delta: float


@dataclass
class Dihedral:
    atoms: list[int]
    terms: list[DihedralTerm]


@dataclass
class PairException:
    atoms: list[int]
    lj: float
    coulomb: float


@dataclass
class Topology:
    bonds: list[Bond] = field(default_factory=list)
    angles: list[Angle] = field(default_factory=list)
    dihedrals: list[Dihedral] = field(default_factory=list)
    exclude_bonded: bool = False
    lj_scales: Optional[list[float]] = None
    coulomb_scales: Optional[list[float]] = None
    exceptions: list[PairException] = field(default_factory=list)


@dataclass
class Constraint:
    atoms: list[int]
    length: float


@dataclass
class Constraints:
    bonds: list[Constraint] = field(default_factory=list)


@dataclass
class MonteCarlo:
    target_pressure: float
    target_temperature: float
    max_log_volume: float = 0.05
    interval: int = 10
    seed: int = 42
    type: str = "monte_carlo"


@dataclass
class BerendsenBarostat:
    target_pressure: float
    tau: float = 0.5
    compressibility: float = 1.0
    anisotropic: bool = False
    type: str = "berendsen"


@dataclass
class Config:
    """Input units are reduced or real (Å, ps, amu, e, kcal/mol, K, atm)."""

    structure: Union[Fcc, Xyz, Water, RockSalt, RandomAlloy]
    model: Union[LennardJones, Buckingham, Morse, SoftSphere, Tabulated]
    integrator: Union[VelocityVerlet, Langevin, NoseHooverChain]
    run: Run = field(default_factory=Run)

    units: str = "reduced"
    electrostatics: Optional[Ewald] = None
    topology: Optional[Topology] = None
    constraints: Optional[Constraints] = None
    thermostat: Optional[Union[Csvr, Berendsen, NoThermostat]] = None
    barostat: Optional[Union[MonteCarlo, BerendsenBarostat]] = None

    def to_toml(self) -> str:
        out = [f"units = {_scalar(self.units)}"]
        for name in ("structure", "model", "integrator", "run", "electrostatics", "topology", "constraints", "thermostat", "barostat"):
            obj = getattr(self, name)
            if obj is not None:
                out.append(f"[{name}]")
                out.extend(_render_section(obj, name))
        return "\n".join(out) + "\n"
