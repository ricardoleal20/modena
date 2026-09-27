"""modena — configure and drive MODENA molecular-dynamics runs from Python.

Transport A (this package, v0.1): renders a typed config as TOML and drives
the `modena` CLI via subprocess — zero dependencies. Transport B (v0.2,
`modena.mcp_client`) drives the `modena-mcp` server for streaming/agent
workflows. The Rust `InputConfig` validator is the single source of truth.
"""

from . import config, run
from .config import Config, Fcc, Langevin, LennardJones, Run, VelocityVerlet, Xyz
from .run import RunResult, check, run

__all__ = [
    "Config", "Fcc", "Xyz", "Langevin", "LennardJones", "Run", "VelocityVerlet",
    "RunResult", "check", "run",
]
__version__ = "0.1.0"