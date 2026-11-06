"""G2: initialized official MACE software parity, never chemical accuracy.

All acceptance checks deliberately survive Python -O. Only graph fields may
change between official and native batches; each forward receives fresh tensors.
"""
import contextlib
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import re
import sys
import uuid

ROOT = Path(__file__).resolve().parents[4]
PINS = {"torch": "2.6.0", "mace-torch": "0.3.13", "e3nn": "0.4.4",
        "numpy": "1.26.4", "ase": "3.26.0", "matscipy": "1.1.1"}
RECIPE = {"seed": 20261106, "r_max": 3.1, "num_bessel": 4,
          "num_polynomial_cutoff": 5, "max_ell": 1, "num_interactions": 2,
          "num_elements": 2, "hidden_irreps": "8x0e + 8x1o", "MLP_irreps": "8x0e",
          "atomic_energies": [0.0, 0.0], "avg_num_neighbors": 6.0,
          "atomic_numbers": [6, 8], "correlation": 2, "gate": "silu",
          "radial_MLP": [8, 8], "dtype": "float64", "device": "cpu", "threads": 1,
          "interaction_cls": "RealAgnosticResidualInteractionBlock",
          "interaction_cls_first": "RealAgnosticInteractionBlock"}
GRAPH_FIELDS = {"edge_index", "unit_shifts", "shifts"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def file_hash(path):
    return digest(Path(path).read_bytes())


def initialize():
    require(platform.python_version() == "3.12.8", "Python pin mismatch")
    for name, expected in PINS.items():
        require(metadata.version(name) == expected, f"package pin mismatch: {name}")
    lock = (ROOT / "bench/requirements-b04-mace-py312.lock").read_text()
    for name, version in re.findall(r"^([A-Za-z0-9_.-]+)==([^\s\\]+)", lock, re.MULTILINE):
        require(metadata.version(name) == version, f"transitive package pin mismatch: {name}")
    # Controller-approved narrow compatibility scope for bundled e3nn constants.
    # Keep torch's weights-only loader; never enable arbitrary pickle loading.
    require(os.environ.get("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "").lower()
            not in {"1", "y", "yes", "true"}, "unsafe loader environment override")
    constants = metadata.distribution("e3nn").locate_file("e3nn/o3/constants.pt")
    require(file_hash(constants) == "279fb1f7844b2b85df092550bf7854a70fe3b427b312672f6a1e1fef0a786876",
            "e3nn packaged constants hash mismatch")
    import torch
    prior_globals = torch.serialization.get_safe_globals()
    with torch.serialization.safe_globals([slice]):
        from e3nn import o3
    require(torch.serialization.get_safe_globals() == prior_globals, "safe globals scope leaked")
    from mace import modules
    import numpy as np
    from mace.data import AtomicData, Configuration
    from mace.tools import AtomicNumberTable
    from mace.tools.torch_geometric.batch import Batch
    globals().update(torch=torch, np=np, AtomicData=AtomicData,
                     Configuration=Configuration, AtomicNumberTable=AtomicNumberTable, Batch=Batch)
    torch.set_default_dtype(torch.float64)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.manual_seed(20261106)
    np.random.seed(20261106)
    model = modules.MACE(
        r_max=3.1, num_bessel=4, num_polynomial_cutoff=5, max_ell=1,
        interaction_cls=modules.RealAgnosticResidualInteractionBlock,
        interaction_cls_first=modules.RealAgnosticInteractionBlock,
        num_interactions=2, num_elements=2,
        hidden_irreps=o3.Irreps('8x0e + 8x1o'),
        MLP_irreps=o3.Irreps('8x0e'), atomic_energies=np.zeros(2),
        avg_num_neighbors=6.0, atomic_numbers=[6, 8], correlation=2,
        gate=torch.nn.functional.silu, radial_MLP=[8, 8],
    ).double().cpu().eval()
    require(torch.get_num_threads() == 1 and torch.get_default_dtype() == torch.float64,
            "arithmetic configuration mismatch")
    require(all(p.device.type == "cpu" and p.dtype == torch.float64 for p in model.parameters()),
            "model arithmetic mismatch")
    return model


def state_identity(model):
    sha = hashlib.sha256()
    description = []
    for name, tensor in sorted(model.state_dict().items()):
        tensor = tensor.detach().cpu().contiguous()
        record = {"name": name, "dtype": str(tensor.dtype), "shape": list(tensor.shape)}
        header = encoded(record)
        raw = tensor.numpy().tobytes(order="C")
        sha.update(len(header).to_bytes(8, "big")); sha.update(header)
        sha.update(len(raw).to_bytes(8, "big")); sha.update(raw)
        description.append(record)
    return sha.hexdigest(), description


def tuples(edges, atoms):
    result = []
    for e in edges:
        require(set(e) == {"source", "destination", "edge_cell_shift"}, "edge schema mismatch")
        i, j, shift = e["source"], e["destination"], e["edge_cell_shift"]
        require(type(i) is int and type(j) is int and 0 <= i < atoms and 0 <= j < atoms,
                "edge index mismatch")
        require(len(shift) == 3 and all(type(v) is int and abs(v) <= 2**53 for v in shift),
                "exact integer shift required")
        result.append((i, j, *shift))
    return result


def batch_tuples(batch):
    shifts = batch.unit_shifts.detach().cpu()
    require(torch.isfinite(shifts).all().item() and torch.equal(shifts, shifts.round()),
            "official shifts are not integral")
    return [(int(i), int(j), *(int(v) for v in s)) for (i, j), s in
            zip(batch.edge_index.T.tolist(), shifts.tolist())]


def nongraph_equal(reference, alternate):
    a, b = reference.to_dict(), alternate.to_dict()
    require(set(a) == set(b), "non-graph field keys mismatch")
    for key in set(a) - GRAPH_FIELDS:
        if torch.is_tensor(a[key]):
            require(torch.is_tensor(b[key]) and a[key].dtype == b[key].dtype
                    and a[key].device == b[key].device and torch.equal(a[key], b[key]),
                    f"non-graph field mismatch: {key}")
        else:
            require(a[key] == b[key], f"non-graph field mismatch: {key}")


def replace_graph(reference, edges, cell):
    alternate = reference.clone()
    alternate.edge_index = torch.tensor([[e[0], e[1]] for e in edges], dtype=torch.long).reshape(-1, 2).T.contiguous()
    alternate.unit_shifts = torch.tensor([e[2:] for e in edges], dtype=torch.float64).reshape(-1, 3)
    alternate.shifts = alternate.unit_shifts @ torch.tensor(cell, dtype=torch.float64)
    nongraph_equal(reference, alternate)
    return alternate


def prepare(case):
    f = case["fixture"]
    require(f["cutoff"] == 3.1 and f["pbc"] == [True]*3 and f["periodic_self"] is True,
            "fixture model contract mismatch")
    require(case["margin"] >= 1e-8, "fixture cutoff margin missing")
    config = Configuration(atomic_numbers=np.asarray(f["atomic_numbers"]),
        positions=np.asarray(f["positions"]), properties={}, property_weights={},
        cell=np.asarray(f["cell"]), pbc=(True, True, True))
    reference = Batch.from_data_list([AtomicData.from_config(
        config, z_table=AtomicNumberTable([6, 8]), cutoff=3.1)])
    native = tuples(case["edges"], len(f["positions"]))
    official = batch_tuples(reference)
    require(native == sorted(official), "integer-image tuple mismatch")
    require(len(native) == len(set(native)) and len(native) > 0, "duplicate/empty graph")
    alternate = replace_graph(reference, native, f["cell"])
    canonical = replace_graph(reference, sorted(official), f["cell"])
    return reference, alternate, canonical, native


def evaluate(model, batch):
    fresh = batch.clone()
    require(fresh.positions.data_ptr() != batch.positions.data_ptr(), "positions must be fresh")
    require(fresh.node_attrs.data_ptr() != batch.node_attrs.data_ptr(), "node tensors must be fresh")
    output = model(fresh.to_dict(), compute_force=True, training=False)
    energy, forces = output["energy"], output["forces"]
    require(energy is not None and forces is not None, "missing model outputs")
    require(torch.isfinite(energy).all().item() and torch.isfinite(forces).all().item(), "nonfinite model output")
    require(energy.abs().max().item() > 1e-10 and forces.abs().max().item() > 1e-10,
            "zero/trivial model output")
    require(energy.dtype == forces.dtype == torch.float64, "output dtype mismatch")
    return {"energy": energy.detach(), "forces": forces.detach()}


def errors(a, b):
    return {key: (a[key] - b[key]).abs().max().item() for key in ["energy", "forces"]}


def close(a, b, tolerance=1e-10):
    for key in ["energy", "forces"]:
        torch.testing.assert_close(a[key], b[key], atol=tolerance, rtol=tolerance)
    return errors(a, b)


def require_mutation_sensitivity(values):
    require(set(values) == {"drop_nonzero", "drop_self", "change_shift"}, "mutation schema mismatch")
    for name, value in values.items():
        require(np.isfinite(value) and value > 1e-10, f"mutation suite inadequate: {name}")


def main(raw):
    request = json.loads(raw)
    require(set(request) == {"schema", "cases", "finite_differences", "finite_difference_h"}
            and request["schema"] == 1, "request schema mismatch")
    require([x["fixture"]["id"] for x in request["cases"]] == ["mace3", "gauge", "deform", "displaced"], "fixture suite mismatch")
    require(request["finite_difference_h"] == 1e-5 and len(request["finite_differences"]) == 18,
            "finite difference schema mismatch")
    model = initialize()
    state_sha, state_description = state_identity(model)
    cases_path = ROOT / "crates/modena-core/tests/fixtures/periodic_graph/cases.json"
    frozen = {f["id"]: f for f in json.loads(cases_path.read_bytes())}
    for case in request["cases"]:
        require(case["fixture"] == frozen[case["fixture"]["id"]], "frozen fixture mismatch")
    manifest = {"status": "preflight", "recipe": RECIPE, "model_state_sha256": state_sha,
                "parameter_count": sum(p.numel() for p in model.parameters()),
                "state_description": state_description, "request_sha256": digest(raw),
                "cases_sha256": file_hash(cases_path), "helper_sha256": file_hash(__file__),
                "lock_sha256": file_hash(ROOT / "bench/requirements-b04-mace-py312.lock"),
                "python": sys.version, "interpreter": sys.executable,
                "platform": platform.platform(), "machine": platform.machine(),
                "python_optimization": sys.flags.optimize,
                "mutation_recipe": "first nonzero edge; drop it / drop all self images / first nonzero shift component increases in absolute value by one",
                "versions": {d.metadata["Name"]: d.version for d in metadata.distributions()}}
    site = Path(metadata.distribution("mace-torch").locate_file(""))
    manifest["sources"] = {name: file_hash(site/name) for name in [
        "mace/modules/models.py", "mace/modules/utils.py", "mace/data/atomic_data.py",
        "mace/data/neighborhood.py", "e3nn/o3/_wigner.py", "e3nn/o3/constants.pt"]}
    # Emit input/model/source identity before tuple checks and any numerical evaluation.
    print("MACE_PREFLIGHT " + json.dumps(manifest, sort_keys=True), file=sys.stderr, flush=True)
    evidence = os.environ.get("MODENA_MACE_EVIDENCE")
    if evidence:
        provenance = Path(evidence) / "dependency-provenance.json"
        proof = json.loads(provenance.read_bytes())
        require(proof["lock_sha256"] == manifest["lock_sha256"], "dependency provenance lock mismatch")
        manifest["dependency_provenance_sha256"] = file_hash(provenance)
        manifest["package_artifacts"] = proof["packages"]
        path = Path(evidence) / ("preflight-" + uuid.uuid4().hex + ".json")
        path.write_bytes(encoded(manifest) + b"\n")
        path.with_suffix(".request.json").write_bytes(raw)
    prepared = [prepare(c) for c in request["cases"]]
    fd_prepared = [prepare(c) for c in request["finite_differences"]]
    # Explicit negative guard self-checks also execute under -O.
    changed = prepared[0][1].clone(); changed.positions[0, 0] += 1.0
    try:
        nongraph_equal(prepared[0][0], changed)
    except ValueError as exc:
        require("non-graph field mismatch" in str(exc), "wrong non-graph guard error")
    else:
        raise ValueError("non-graph guard ignored a position change")
    try:
        require_mutation_sensitivity(dict.fromkeys(["drop_nonzero", "drop_self", "change_shift"], 0.0))
    except ValueError:
        pass
    else:
        raise ValueError("mutation guard accepted an ignored-graph path")
    results, outputs = [], {}
    adequacy = dict.fromkeys(["drop_nonzero", "drop_self", "change_shift"], 0.0)
    for case, (reference, alternate, canonical, edges) in zip(request["cases"], prepared):
        original = evaluate(model, reference)
        ordered = evaluate(model, canonical)
        actual = evaluate(model, alternate)
        parity = close(actual, original)
        order_error = close(ordered, original)
        repeated = close(evaluate(model, alternate), actual, 0.0)
        index = next(i for i, e in enumerate(edges) if any(e[2:]))
        drop = edges[:index] + edges[index+1:]
        no_self = [e for e in edges if e[0] != e[1]]
        shifted = list(edges)
        e = list(shifted[index])
        # Fixed before the rerun: change the first nonzero component away from
        # zero, preserving a nonzero image even for a periodic self edge.
        component = next(k for k in range(2, 5) if e[k] != 0)
        e[component] += 1 if e[component] > 0 else -1
        shifted[index] = tuple(e)
        controls = {}
        for name, mutated in [("drop_nonzero", drop), ("drop_self", no_self), ("change_shift", shifted)]:
            require(mutated != edges, f"mutation inapplicable: {name}")
            control = evaluate(model, replace_graph(reference, mutated, case["fixture"]["cell"]))
            controls[name] = errors(actual, control)
            adequacy[name] = max(adequacy[name], *controls[name].values())
        fid = case["fixture"]["id"]
        outputs[fid] = actual
        results.append({"fixture_id": fid, "edge_count": len(edges), "margin": case["margin"],
                        "energy": actual["energy"].tolist(), "forces": actual["forces"].tolist(),
                        "parity_error": parity, "order_error": order_error,
                        "repeat_error": repeated, "mutations": controls})
    require_mutation_sensitivity(adequacy)
    close(outputs["gauge"], outputs["mace3"])
    fd_values = {}
    base = frozen["displaced"]
    for case, (reference, alternate, canonical, edges) in zip(request["finite_differences"], fd_prepared):
        atom, axis, sign = case["atom"], case["axis"], case["sign"]
        require(type(atom) is int and 0 <= atom < 3 and type(axis) is int and 0 <= axis < 3
                and sign in [-1.0, 1.0], "finite difference coordinate mismatch")
        expected = json.loads(json.dumps(base))
        expected["id"] = case["fixture"]["id"]
        expected["positions"][atom][axis] += sign * 1e-5
        require(case["fixture"] == expected, "finite difference fixture mismatch")
        require(edges == prepared[3][3], "finite difference neighbor crossing")
        key = (atom, axis, sign)
        require(key not in fd_values, "duplicate finite difference coordinate")
        out = evaluate(model, alternate)
        close(out, evaluate(model, reference))
        fd_values[key] = out["energy"].item()
    finite = torch.empty((3, 3), dtype=torch.float64)
    for atom in range(3):
        for axis in range(3):
            finite[atom, axis] = -(fd_values[atom, axis, 1.0] - fd_values[atom, axis, -1.0]) / (2e-5)
    torch.testing.assert_close(finite, outputs["displaced"]["forces"], atol=1e-6, rtol=1e-6)
    require(state_identity(model)[0] == state_sha, "model mutated during inference")
    manifest.update(status="passed", cases=results, mutation_suite_sensitivity=adequacy,
                    finite_difference_components=9, finite_difference_h=1e-5,
                    finite_difference_max_error=(finite-outputs["displaced"]["forces"]).abs().max().item(),
                    finite_difference_forces=finite.tolist())
    return manifest


if __name__ == "__main__":
    raw = sys.stdin.buffer.read()
    # Official imports can write optional-acceleration notices to stdout.
    with contextlib.redirect_stdout(sys.stderr):
        result = main(raw)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
