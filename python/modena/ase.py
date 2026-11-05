"""Optional ASE energy/forces calculator for native shifted-force LJ.

Install ``modena[ase]`` and import explicitly from this module. Lengths are
Angstrom, energies eV; sigma/epsilon describe the one homogeneous species.
"""
from __future__ import annotations

import json
import math
import os
import subprocess

import numpy as np
from ase.calculators.calculator import (
    Calculator, CalculatorSetupError, CalculationFailed,
    PropertyNotImplementedError, all_changes,
)

from .run import _binary

_PROTOCOL = "modena-native-lj-v1"
_UNITS = dict(length="angstrom", energy="eV", force="eV/angstrom")
_MODEL_KEYS = ("sigma_angstrom", "epsilon_ev", "cutoff_angstrom", "skin_angstrom")


class ModenaCalculator(Calculator):
    """Single-species native shifted-force LJ via one CLI process per snapshot.

    ``binary`` is one executable path or PATH name, never a shell command.
    ``timeout`` is a positive finite subprocess deadline in seconds. Only
    energy and forces are available; full periodic diagonal cells are required.
    Supported atomic numbers are 1..118; numbers select identity, not parameters.
    """

    implemented_properties = ["energy", "forces"]
    discard_results_on_any_change = True

    def __init__(self, *, sigma_angstrom, epsilon_ev, cutoff_angstrom,
                 skin_angstrom=0.3, binary=None, timeout=30.0):
        super().__init__(sigma_angstrom=sigma_angstrom, epsilon_ev=epsilon_ev,
                         cutoff_angstrom=cutoff_angstrom, skin_angstrom=skin_angstrom,
                         binary=binary, timeout=timeout)

    def set(self, **kwargs):
        unknown = set(kwargs) - set(_MODEL_KEYS) - {"binary", "timeout"}
        if unknown:
            self.reset()
            raise CalculatorSetupError(f"unknown ModenaCalculator parameters: {sorted(unknown)}")
        try:
            return super().set(**kwargs)
        except Exception:
            self.reset()
            raise

    def _request(self, atoms):
        if atoms is None:
            raise CalculatorSetupError("atoms are required")
        p = self.parameters
        try:
            values = {k: float(p[k]) for k in (*_MODEL_KEYS, "timeout")}
        except (TypeError, ValueError, OverflowError) as exc:
            raise CalculatorSetupError("LJ parameters and timeout must be finite numbers") from exc
        if (not all(math.isfinite(v) for v in values.values())
                or values["sigma_angstrom"] <= 0 or values["epsilon_ev"] < 0
                or values["cutoff_angstrom"] <= 0 or values["skin_angstrom"] < 0
                or values["timeout"] <= 0):
            raise CalculatorSetupError("invalid finite LJ parameters, skin or timeout")
        if p["binary"] is not None:
            try:
                executable = os.fspath(p["binary"])
            except TypeError as exc:
                raise CalculatorSetupError("binary must be one executable path or PATH name") from exc
            if not isinstance(executable, str) or not executable or "\0" in executable:
                raise CalculatorSetupError("binary must be one executable path or PATH name")
        n = len(atoms)
        numbers = atoms.get_atomic_numbers()
        if (not 1 <= n <= 100000 or not np.all((numbers >= 1) & (numbers <= 118))
                or not np.all(numbers == numbers[0])):
            raise CalculatorSetupError("requires 1..100000 homogeneous atomic numbers in 1..118")
        if atoms.constraints:
            raise CalculatorSetupError("constraints are unsupported")
        if (np.any(atoms.get_initial_charges() != 0)
                or np.any(atoms.get_initial_magnetic_moments() != 0)):
            raise CalculatorSetupError("initial charge and magnetic moments must be zero")
        cell = np.asarray(atoms.cell)
        if (not np.all(atoms.pbc) or not np.all(np.isfinite(cell))
                or not np.all(np.diag(cell) > 0)
                or np.any(cell != np.diag(np.diag(cell)))):
            raise CalculatorSetupError("requires full PBC and a positive finite diagonal orthorhombic cell")
        if not np.all(np.isfinite(atoms.positions)):
            raise CalculatorSetupError("positions must be finite")
        radius = values["cutoff_angstrom"] + values["skin_angstrom"]
        if not math.isfinite(radius) or np.any(radius > np.diag(cell) / 2):
            raise CalculatorSetupError("cutoff plus skin must be at most half each cell length")
        return dict(protocol=_PROTOCOL, units=_UNITS, atomic_numbers=numbers.tolist(),
                    positions_angstrom=atoms.positions.tolist(), cell_angstrom=cell.tolist(),
                    pbc=atoms.pbc.tolist(), **{k: values[k] for k in _MODEL_KEYS})

    def get_property(self, name, atoms=None, allow_calculation=True):
        # ASE's ordinary cache compares positions/numbers/cell/PBC/charge/spin,
        # but not constraints. Validate on every access, including cache hits.
        try:
            if name not in self.implemented_properties:
                raise PropertyNotImplementedError(f"{name} is not available from ModenaCalculator")
            self._request(atoms if atoms is not None else self.atoms)
            return super().get_property(name, atoms, allow_calculation)
        except Exception:
            self.reset()
            raise

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        self.results = {}
        try:
            if any(name not in self.implemented_properties for name in properties):
                raise PropertyNotImplementedError("only energy and forces are available")
            atoms = atoms if atoms is not None else self.atoms
            request = self._request(atoms)
            encoded = json.dumps(request, allow_nan=False)
            if len(encoded.encode("utf-8")) > 8 * 1024 * 1024:
                raise CalculatorSetupError("snapshot input exceeds 8 MiB")
            super().calculate(atoms, properties, system_changes)
            try:
                result = subprocess.run([self.parameters["binary"] or _binary(), "snapshot"],
                    input=encoded, capture_output=True, text=True, timeout=float(self.parameters["timeout"]),
                    shell=False)
            except subprocess.TimeoutExpired as exc:
                raise CalculationFailed("MODENA snapshot timed out") from exc
            except OSError as exc:
                raise CalculationFailed(f"cannot execute MODENA snapshot: {exc}") from exc
            if result.returncode != 0:
                raise CalculationFailed(f"MODENA snapshot failed (exit {result.returncode}): {result.stderr[-4096:]}")
            try:
                response = json.loads(result.stdout)
                if (set(response) != {"protocol", "units", "energy_ev", "forces_ev_per_angstrom"}
                        or response["protocol"] != _PROTOCOL or response["units"] != _UNITS
                        or type(response["energy_ev"]) not in (int, float)):
                    raise ValueError("invalid snapshot schema, protocol or units")
                rows = response["forces_ev_per_angstrom"]
                if (not isinstance(rows, list) or len(rows) != len(atoms)
                        or any(not isinstance(row, list) or len(row) != 3
                               or any(type(v) not in (int, float) for v in row) for row in rows)):
                    raise ValueError("invalid numeric forces shape")
                energy = float(response["energy_ev"])
                forces = np.asarray(rows, dtype=float)
                if (not math.isfinite(energy) or forces.shape != (len(atoms), 3)
                        or not np.all(np.isfinite(forces))):
                    raise ValueError("invalid snapshot energy or forces")
            except (ValueError, KeyError, TypeError, OverflowError) as exc:
                raise CalculationFailed(f"invalid MODENA snapshot response: {exc}") from exc
            self.results = dict(energy=energy, forces=forces)
        except Exception:
            self.reset()
            raise
