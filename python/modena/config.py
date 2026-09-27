"""Typed mirrors of MODENA's TOML schema (`InputConfig`).

`Config.to_toml()` renders field names/values exactly like the Rust serde
schema, so the Rust `InputConfig` validator is the single source of truth —
Python never re-implements validation or unit conversion (the parity tests
render the shipped validation decks and diff them byte-for-byte).
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from typing import Any, Optional, Union


def _scalar(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, str):
        return f'"{v}"'
    raise TypeError(f"unhandled scalar {type(v)}: {v!r}")


def _val(v: Any) -> str:
    if isinstance(v, list):
        return "[" + ", ".join(_scalar(x) for x in v) + "]"
    return _scalar(v)


def _render_section(obj: Any, indent: int = 0) -> str:
    """Render a dataclass as TOML key = value lines (leaf lists inline)."""
    pad = "  " * indent
    lines = []
    for f in fields(obj):  # type: ignore[arg-type]
        v = getattr(obj, f.name)
        if v is None or v == [] or v == {}:
            continue
        if isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
            # list-of-tables: [[name]] nested tables
            for item in v:
                lines.append(f"{pad}[[{f.name}]]")
                lines.append(_render_dict_table(item, indent + 1))
        elif is_dataclass(v) and not isinstance(v, type):
            lines.append(f"{pad}{f.name} = {{ ... }}")  # nested table (unused today)
        elif not is_dataclass(v):
            lines.append(f"{pad}{f.name} = {_val(v)}")
    return "\n".join(lines)


def _render_dict_table(d: dict[str, Any], indent: int) -> str:
    pad = "  " * indent
    return "\n".join(
        f"{pad}{k} = {_val(v)}" for k, v in d.items() if v is not None
    )


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


@dataclass
class Config:
    """Renderable MODENA config (reduced units). `to_toml()` feeds the CLI."""

    structure: Union[Fcc, Xyz]
    model: Union[LennardJones, Buckingham, Morse]
    integrator: Union[VelocityVerlet, Langevin]
    run: Run = field(default_factory=Run)

    def to_toml(self) -> str:
        out = []
        for name in ("structure", "model", "integrator", "run"):
            out.append(f"[{name}]")
            out.append(_render_section(getattr(self, name)))
        return "\n".join(out) + "\n"