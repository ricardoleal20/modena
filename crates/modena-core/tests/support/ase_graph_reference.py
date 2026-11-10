"""Pinned official ASE all-images comparison against actual Rust graph edges."""
import hashlib
import inspect
import json
import platform
from pathlib import Path
import sys

import ase
import numpy as np
from ase import Atoms, neighborlist

raw_request = sys.stdin.buffer.read()
request = json.loads(raw_request)
f = request["fixture"]
if ase.__version__ != "3.29.0":
    raise RuntimeError(f"ASE version mismatch: {ase.__version__}")
if np.__version__ != "2.5.4":
    raise RuntimeError(f"NumPy version mismatch: {np.__version__}")
a = Atoms(numbers=f["atomic_numbers"], positions=np.asarray(f["positions"]).reshape(-1, 3),
          cell=f["cell"], pbc=f["pbc"])
i, j, shifts, vectors = neighborlist.neighbor_list(
    "ijSD", a, f["cutoff"], self_interaction=False)
records = sorted(((int(x), int(y), tuple(map(int, shift))), vector)
                 for x, y, shift, vector in zip(i, j, shifts, vectors)
                 if f["periodic_self"] or x != y)
reference = [edge for edge, vector in records]
actual = sorted((e["source"], e["destination"], tuple(e["edge_cell_shift"]))
                for e in request["edges"])
if actual != reference:
    raise ValueError(f"integer-image tuple mismatch: fixture={f['id']!r} "
                     f"actual={actual!r} reference={reference!r}")
source = inspect.getsourcefile(neighborlist)
with open(source, "rb") as stream:
    digest = hashlib.sha256(stream.read()).hexdigest()


def vector(edge):
    x, y, shift = edge
    value = np.asarray(f["positions"][y]) - np.asarray(f["positions"][x])
    for axis in range(3):
        value = value + np.asarray(f["cell"][axis]) * shift[axis]
    return value


ref_vectors = np.asarray([v for edge, v in records]).reshape(-1, 3)
act_vectors = np.asarray([vector(e) for e in actual]).reshape(-1, 3)
np.testing.assert_allclose(act_vectors, ref_vectors, rtol=2e-12, atol=2e-12)
vector_error = float(np.max(np.abs(act_vectors - ref_vectors), initial=0.0))
cases_path = Path(__file__).parent.parent / "fixtures" / "periodic_graph" / "cases.json"
json.dump({"fixture_id": f["id"], "edge_count": len(reference),
           "versions": {"ase": ase.__version__, "numpy": np.__version__,
                        "python": platform.python_version()},
           "source_sha256": digest, "max_vector_error": vector_error,
           "source_path": source, "interpreter": sys.executable,
           "platform": platform.platform(), "machine": platform.machine(),
           "request_sha256": hashlib.sha256(raw_request).hexdigest(),
           "fixture_sha256": hashlib.sha256(json.dumps(f, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
           "cases_sha256": hashlib.sha256(cases_path.read_bytes()).hexdigest(),
           "helper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, sys.stdout)
