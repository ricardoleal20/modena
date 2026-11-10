"""Strict one-shot and opt-in framed sessions for ASE's fixed EMT calculator."""
import hashlib
import inspect
import json
import math
import pathlib
import sys
import struct

PROTOCOL = "modena-ase-emt-v1"
SUPPORTED = {13, 28, 29, 46, 47, 78, 79}
SESSION_PROTOCOL = "modena-ase-emt-session-v1"
MAX_FRAME = 16 * 1024 * 1024


def identity():
    import ase
    import numpy
    from ase.calculators import emt
    return {
        "protocol": PROTOCOL,
        "provider": "ase_emt",
        "model": "EMT",
        "task": "energy_forces_stress",
        "composition": "total",
        "units": {"positions": "angstrom", "energy": "eV", "forces": "eV/angstrom", "stress": "eV/angstrom^3"},
        "capabilities": ["energy", "forces", "stress", "full_pbc", "provider_neighbors", "conservative", "total_momentum"],
        "ase_version": ase.__version__,
        "numpy_version": numpy.__version__,
        "model_sha256": hashlib.sha256(pathlib.Path(inspect.getfile(emt)).read_bytes()).hexdigest(),
        "asap_cutoff": False,
    }


def matrix(value, rows, name):
    if not isinstance(value, list) or len(value) != rows:
        raise ValueError(name + " must have the required row count")
    for row in value:
        if not isinstance(row, list) or len(row) != 3 or any(type(x) not in (int, float) or not math.isfinite(x) for x in row):
            raise ValueError(name + " must contain finite triples")


def handle(request, calculator=None, provider_identity=None):
    if not isinstance(request, dict) or request.get("protocol") != PROTOCOL:
        raise ValueError("unsupported protocol")
    action = request.get("action")
    keys = {"protocol", "action"}
    if action == "evaluate":
        keys |= {"atomic_numbers", "positions_angstrom", "cell_angstrom", "pbc", "asap_cutoff"}
    elif action != "info":
        raise ValueError("unsupported action")
    if set(request) != keys:
        raise ValueError("unexpected or missing request fields")
    if action == "info":
        return {"protocol": PROTOCOL, "action": action, "identity": provider_identity or identity()}
    numbers = request["atomic_numbers"]
    if not isinstance(numbers, list) or not numbers or any(type(z) is not int or z not in SUPPORTED for z in numbers):
        raise ValueError("unsupported atomic_numbers")
    if request["pbc"] != [True, True, True] or any(type(x) is not bool for x in request["pbc"]):
        raise ValueError("full pbc required")
    if type(request["asap_cutoff"]) is not bool or request["asap_cutoff"]:
        raise ValueError("asap_cutoff must be false")
    matrix(request["positions_angstrom"], len(numbers), "positions_angstrom")
    matrix(request["cell_angstrom"], 3, "cell_angstrom")
    import numpy as np
    from ase import Atoms
    from ase.calculators.emt import EMT
    det = np.linalg.det(request["cell_angstrom"])
    if not math.isfinite(det) or det <= 0:
        raise ValueError("cell must have a positive finite determinant")
    atoms = Atoms(numbers=numbers, positions=request["positions_angstrom"], cell=request["cell_angstrom"], pbc=True)
    # ASE normally tolerates tiny coordinate changes in its property cache.
    # Full snapshots require exact invalidation, including sub-tolerance edits.
    if calculator is not None and calculator.check_state(atoms, tol=0.0):
        calculator.reset()
    atoms.calc = calculator if calculator is not None else EMT(asap_cutoff=False)
    evaluation = {"energy_ev": float(atoms.get_potential_energy()), "forces_ev_per_angstrom": atoms.get_forces().tolist(), "stress_ev_per_angstrom3": atoms.get_stress(voigt=False).tolist()}
    return {"protocol": PROTOCOL, "action": action, "identity": provider_identity or identity(), "evaluation": evaluation}


def session():
    from ase.calculators.emt import EMT
    calculator = EMT(asap_cutoff=False)
    provider_identity = identity()
    expected_id = 1
    while True:
        header = sys.stdin.buffer.read(4)
        if not header:
            return 0
        if len(header) != 4:
            raise ValueError("truncated session header")
        size = struct.unpack(">I", header)[0]
        if size == 0 or size > MAX_FRAME:
            raise ValueError("session input frame limit exceeded")
        data = sys.stdin.buffer.read(size)
        if len(data) != size:
            raise ValueError("truncated session frame")
        frame = json.loads(data)
        if (not isinstance(frame, dict) or set(frame) != {"session_protocol", "request_id", "nonce", "request"}
                or frame["session_protocol"] != SESSION_PROTOCOL
                or type(frame["request_id"]) is not int or frame["request_id"] != expected_id
                or not isinstance(frame["nonce"], list) or len(frame["nonce"]) != 32
                or any(type(x) is not int or not 0 <= x < 256 for x in frame["nonce"])):
            raise ValueError("invalid session protocol, request ID or nonce")
        response = handle(frame["request"], calculator, provider_identity)
        output = json.dumps({"session_protocol": SESSION_PROTOCOL, "request_id": expected_id, "nonce": frame["nonce"], "response": response}, allow_nan=False, separators=(",", ":")).encode()
        if len(output) > MAX_FRAME:
            raise ValueError("session output frame limit exceeded")
        sys.stdout.buffer.write(struct.pack(">I", len(output)) + output)
        sys.stdout.buffer.flush()
        expected_id += 1


def main():
    try:
        if sys.argv[1:] == ["--persistent"]:
            return session()
        request = json.load(sys.stdin)
        response = handle(request)
        output = json.dumps(response, allow_nan=False, separators=(",", ":"))
    except Exception as error:
        print((type(error).__name__ + ": " + str(error))[:4096], file=sys.stderr)
        return 1
    print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
