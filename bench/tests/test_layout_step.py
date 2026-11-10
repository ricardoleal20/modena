"""Real child processes exercise malformed output, rather than subprocess mocks."""
import importlib.util
from contextlib import contextmanager
import json
import math
import os
import struct
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("layout_step", Path(__file__).parents[1] / "layout_step.py")
try:
    driver = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(driver)
except FileNotFoundError:
    driver = None



def fixture():
    n = 256
    p, v, length = driver.initial(n)
    energy = dict(potential=0., kinetic=0., total=0., virial=0., temperature=0., pressure=0.)
    state = dict(n=n, logical_ids=list(range(n)), periodic=[True]*3,
        masses=[1.]*n, species=[0]*n, charges=[0.]*n, atomic_numbers=[0]*n,
        positions=p, velocities=v, forces=[0.]*(3*n),
        cell_vectors=[[length if i==j else 0. for j in range(3)] for i in range(3)],
        steps_completed=0, checksum_fnv1a64='fixture', energies=None,
        neighbor_pairs=None, pair_scales=None)
    ready = json.loads(json.dumps(state))
    ready.update(energies=energy, neighbor_pairs=[], pair_scales=[])
    final = json.loads(json.dumps(ready)); final['steps_completed']=1
    expected=dict(mode='soa', n=n, requested_threads=1, steps_requested=1,
        steps_completed=1, seed=17, dt=0.000001, rebuild_policy='no-rebuild')
    row=dict(expected, schema_version=1, actual_threads=1, box_lengths=[length]*3,
        workload=driver.WORKLOAD, potential=driver.POTENTIAL, setup_seconds=.01,
        elapsed_seconds=.01, pool_setup_seconds=.01, force_max_abs=0.,
        energies=energy, finite_efw=True, simd_claim=False, pair_execution='serial-scalar',
        comparison_class='serial-layout', setup_rebuilds=1, rebuilds=0,
        initial_neighbor_pairs=0, final_neighbor_pairs=0, neighbor_pair_visits=0,
        final_active_pairs=0, retained_packed_bytes=0, cluster_reconstructions=0,
        setup_cluster_reconstructions=0, evaluated_pairs=0,
        initial_checksum_fnv1a64='fixture', final_checksum_fnv1a64='fixture')
    snap=dict(schema_version=1, row={k:v for k,v in row.items() if k!='pool_setup_seconds'},
        initial_input=state, initial_force_ready=ready, final=final)
    return row,snap,expected
def independent_fnv(state):
    values = [state['cell_vectors'][k][k] for k in range(3)]
    for key in ('positions', 'velocities', 'forces', 'masses'):
        values.extend(state[key])
    h = 14695981039346656037
    for value in values:
        for b in struct.pack('<d', value):
            h = ((h ^ b) * 1099511628211) % (2**64)
    return f'{h:016x}'


def physical_fixture():
    row, snap, expected = fixture()
    state = snap['initial_force_ready']; n=state['n']; length=state['cell_vectors'][0][0]
    p=state['positions']; forces=[0.]*(3*n); pe=virial=0.; pairs=[]; active=0
    uc=4*(2.5**-12-2.5**-6); fc=24*(2*2.5**-13-2.5**-7)
    for i in range(n):
        for j in range(i+1,n):
            d=[p[3*i+k]-p[3*j+k] for k in range(3)]
            d=[x-length*math.floor(x/length+.5) for x in d]
            r2=sum(x*x for x in d)
            if r2<2.8**2: pairs.append([i,j])
            if not 1e-14<=r2<6.25: continue
            active+=1; r=math.sqrt(r2)
            pe+=4*(r**-12-r**-6)-uc+(r-2.5)*fc
            scalar=(24*(2*r**-13-r**-7)-fc)/r
            for k in range(3):
                f=d[k]*scalar; forces[3*i+k]+=f; forces[3*j+k]-=f; virial+=d[k]*f
    ke=.5*sum(v*v for v in state['velocities'])
    energies=dict(potential=pe, kinetic=ke, total=pe+ke, virial=virial,
        temperature=2*ke/(3*n), pressure=(2*ke+virial)/(3*length**3))
    for name in ('initial_force_ready','final'):
        snap[name].update(forces=forces.copy(), energies=energies.copy(), neighbor_pairs=pairs.copy())
    row.update(energies=energies.copy(), force_max_abs=max(map(abs,forces)),
        initial_neighbor_pairs=len(pairs), final_neighbor_pairs=len(pairs),
        final_active_pairs=active, neighbor_pair_visits=len(pairs), evaluated_pairs=active)
    for name in ('initial_input','initial_force_ready','final'):
        snap[name]['checksum_fnv1a64']=independent_fnv(snap[name])
    row['initial_checksum_fnv1a64']=snap['initial_input']['checksum_fnv1a64']
    row['final_checksum_fnv1a64']=snap['final']['checksum_fnv1a64']
    snap['row']={k:v for k,v in row.items() if k!='pool_setup_seconds'}
    return row,snap,expected


@contextmanager
def retained_test_directory(name):
    """Keep new regression subprocess bytes when explicitly requested."""
    evidence = os.environ.get("B03_TEST_EVIDENCE_DIR")
    if evidence:
        path = Path(evidence).expanduser().resolve() / name
        path.mkdir(parents=True, exist_ok=False)
        yield path
    else:
        with tempfile.TemporaryDirectory(prefix="b03-"+name+"-") as tmp:
            yield Path(tmp)


def contract_child(out, label, row, snapshot):
    """Synthetic malformed-contract child; never used as the Rust gate."""
    child = out / (label+".py")
    child.write_text("#!"+sys.executable+"\nimport pathlib,sys\n"
        "pathlib.Path(sys.argv[sys.argv.index('--snapshot')+1]).write_text("
        +repr(json.dumps(snapshot))+")\nprint("+repr(json.dumps(row))+")\n")
    child.chmod(0o755)
    return child


def resolved_rust_example(test):
    override = os.environ.get("B03_LAYOUT_STEP_BINARY")
    if override:
        manifest = os.environ.get("B03_LAYOUT_STEP_PROVENANCE")
        test.assertTrue(manifest,
            "B03_LAYOUT_STEP_BINARY requires B03_LAYOUT_STEP_PROVENANCE; see docs/B03_PACKED_LAYOUT.md")
        binary = Path(override).expanduser().resolve()
        provenance_path = Path(manifest).expanduser().resolve()
    else:
        binary = (driver.ROOT/"target/b03-evidence/task2-campaign-20261106/layout_step").resolve()
        provenance_path = binary.parent/"provenance.json"
    prerequisite = ("Required real Rust layout_step contract gate has no frozen executable/provenance. "
        "Build and freeze the example, set B03_LAYOUT_STEP_BINARY and "
        "B03_LAYOUT_STEP_PROVENANCE as documented in docs/B03_PACKED_LAYOUT.md; this gate never skips.")
    test.assertTrue(binary.is_file(),prerequisite)
    test.assertTrue(provenance_path.is_file(),prerequisite)
    provenance = driver.read_json(provenance_path)
    test.assertEqual(driver.sha(binary),provenance["binary_sha256"],"Rust executable provenance SHA mismatch")
    test.assertTrue(provenance.get("rustc","").startswith("rustc "),"compiler provenance is required")
    source_files = ["Cargo.toml","Cargo.lock","crates/modena-core/Cargo.toml",
        "crates/modena-bench/Cargo.toml","crates/modena-bench/examples/layout_step.rs"]
    for crate in ("modena-core","modena-bench"):
        source_files += [str(path.relative_to(driver.ROOT)) for path in
            sorted((driver.ROOT/"crates"/crate/"src").rglob("*.rs"))]
        build = driver.ROOT/"crates"/crate/"build.rs"
        if build.exists(): source_files.append(str(build.relative_to(driver.ROOT)))
    source_hashes = provenance.get("source_sha256",{})
    for relative in source_files:
        test.assertEqual(source_hashes.get(relative),driver.sha(driver.ROOT/relative),
            "Rust source provenance missing/stale for "+relative+"; freeze a matching real build")
    return binary,provenance_path,provenance


class DriverTests(unittest.TestCase):
    def test_real_process_failure_is_preserved(self):
        self.assertIsNotNone(driver, "campaign driver is not implemented")
        with tempfile.TemporaryDirectory() as tmp:
            row = driver.capture([sys.executable, "-c", "import sys; print(\"raw failure\"); sys.exit(7)"], Path(tmp) / "run")
            self.assertEqual(row["returncode"], 7)
            self.assertIn("raw failure", Path(row["stdout_path"]).read_text())
            self.assertGreater(row["wall_seconds"], 0)

    def test_malformed_real_child_output_is_rejected(self):
        self.assertIsNotNone(driver, "campaign driver is not implemented")
        # Removing an identity/physics check makes its corresponding child pass.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            good, snap, expected = fixture()
            for mutation in ("mode", "n", "steps_completed", "missing_energy", "nonfinite", "mass", "missing_state", "snapshot_subset", "pool_time", "rebuilds", "seed", "actual_threads", "saved_nonfinite", "input_position", "mass_bool", "id_bool", "periodic_integer", "saved_step_bool"):
                r, s = json.loads(json.dumps(good)), json.loads(json.dumps(snap))
                if mutation in ("mode", "n", "steps_completed"):
                    r[mutation] = "wrong" if mutation == "mode" else 999
                elif mutation == "missing_energy":
                    del r["energies"]["potential"]
                elif mutation == "nonfinite":
                    r["force_max_abs"] = float("nan")
                elif mutation == "mass":
                    # Both compared children sharing this invalid mass must fail.
                    for state in (s["initial_input"], s["initial_force_ready"], s["final"]):
                        state["masses"][0] = 2
                elif mutation == "missing_state":
                    del s["initial_input"]
                elif mutation == "snapshot_subset":
                    s["row"]["dt"] = 0.5
                elif mutation == "pool_time":
                    r["pool_setup_seconds"] = -1
                elif mutation == "saved_nonfinite":
                    s["final"]["forces"][0] = float("inf")
                elif mutation == "input_position":
                    s["initial_input"]["positions"][0] += 0.1
                elif mutation == "mass_bool":
                    for name in ("initial_input","initial_force_ready","final"): s[name]["masses"][0] = True
                elif mutation == "id_bool":
                    s["initial_input"]["logical_ids"][0] = False
                elif mutation == "periodic_integer":
                    s["initial_input"]["periodic"] = [1,1,1]
                elif mutation == "saved_step_bool":
                    s["initial_input"]["steps_completed"] = False
                else:
                    r[mutation] += 1
                if mutation != "snapshot_subset":
                    s["row"] = {k:v for k,v in r.items() if k != "pool_setup_seconds"}
                child = root / (mutation + ".py")
                snapshot = root / (mutation + ".snapshot")
                child.write_text("import pathlib\npathlib.Path(" + repr(str(snapshot)) + ").write_text(" + repr(json.dumps(s)) + ")\nprint(" + repr(json.dumps(r)) + ")\n")
                for repeat in range(2):
                    with self.subTest(mutation=mutation, repeat=repeat):
                        result = driver.capture([sys.executable, str(child)], root / (mutation + str(repeat)))
                        with self.assertRaises((ValueError, KeyError, AssertionError)):
                            driver.validate(json.loads(Path(result["stdout_path"]).read_text()), driver.read_json(snapshot), expected)

    def test_valid_snapshot_subset_and_stale_output_rejection(self):
        self.assertIsNotNone(driver, "campaign driver is not implemented")
        good, snap, expected = fixture()
        driver.validate(good, snap, expected)
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "used"
            driver.capture([sys.executable, "-c", "print(1)"], prefix)
            with self.assertRaises(FileExistsError):
                driver.capture([sys.executable, "-c", "print(2)"], prefix)


class DirectAuditTests(unittest.TestCase):
    def test_shared_biased_forces_reject_against_oracle_not_tolerant_canonical(self):
        good, original, expected = physical_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp); records=[]
            for number, bias in enumerate((.75,1.5,1.5)):
                row=json.loads(json.dumps(good)); snap=json.loads(json.dumps(original))
                forces=snap['final']['forces']; index=max(range(len(forces)),key=lambda k:abs(forces[k]))
                tolerance=1e-10+1e-12*abs(forces[index])
                forces[index]+=math.copysign(bias*tolerance,forces[index])
                snap['final']['checksum_fnv1a64']=independent_fnv(snap['final'])
                row['final_checksum_fnv1a64']=snap['final']['checksum_fnv1a64']
                row['force_max_abs']=max(map(abs,forces))
                snap['row']={k:v for k,v in row.items() if k!='pool_setup_seconds'}
                child=out/f'fixture{number}.py'
                child.write_text('#!'+sys.executable+'\nimport sys,pathlib\npathlib.Path(sys.argv[sys.argv.index("--snapshot")+1]).write_text('+repr(json.dumps(snap))+')\nprint('+repr(json.dumps(row))+')\n')
                child.chmod(0o755)
                record=driver.run_one(child,out,f'biased{number}','soa',256,1,1,'no-rebuild',.000001)
                self.assertTrue(record['valid'],record.get('error'))
                record.update(phase='warmup',pair=number,role='B',candidate='soa')
                records.append(record)
            (out/'index.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in records))
            (out/'provenance.json').write_text(json.dumps({'binary_sha256':'fixture','loaded_driver_sha256':driver.sha(driver.__file__)}))
            result=subprocess.run([sys.executable,driver.__file__,'audit','--out',str(out)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            audit=json.loads((out/'audit.json').read_text())['rows']
            self.assertTrue(audit[0]['valid'],audit[0])
            for row in audit[1:]:
                self.assertFalse(row['valid'],'two-hop tolerance accepted physically wrong force')
                self.assertIn('force',row['error'].lower())

    def test_fnv_is_recomputed_from_actual_saved_fields(self):
        self.assertTrue(hasattr(driver,'state_fnv'),'independent saved-state FNV is not implemented')
        _,snap,_=physical_fixture()
        for state in (snap['initial_input'],snap['initial_force_ready'],snap['final']):
            self.assertEqual(driver.state_fnv(state),state['checksum_fnv1a64'])
        state=snap['final']; state['forces'][0]+=.001
        self.assertNotEqual(driver.state_fnv(state),state['checksum_fnv1a64'])


class SecondaryTests(unittest.TestCase):
    def test_real_no_snapshot_outputs_validate_all_available_fields(self):
        self.assertTrue(hasattr(driver,'validate_secondary'),'secondary validator is not implemented')
        good,snap,_=physical_fixture()
        group=dict(initial_fnv=independent_fnv(snap['initial_input']),fingerprints={
            independent_fnv(snap['final']):dict(row=good,snapshot_path='audited-reference')})
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            for mutation in ('valid','mode','steps_completed','missing_energy','nonfinite','force_max_abs',
                             'fingerprint','counter','threads_type','pool_setup_seconds'):
                row=json.loads(json.dumps(good))
                if mutation=='valid': row.update(elapsed_seconds=.5,setup_seconds=.02,pool_setup_seconds=.003)
                elif mutation=='mode': row['mode']='cluster4'
                elif mutation=='steps_completed': row['steps_completed']=999
                elif mutation=='missing_energy': del row['energies']['potential']
                elif mutation=='nonfinite': row['energies']['virial']=float('inf')
                elif mutation=='force_max_abs': row['force_max_abs']+=.001
                elif mutation=='fingerprint': row['final_checksum_fnv1a64']='0000000000000000'
                elif mutation=='counter': row['evaluated_pairs']=None
                elif mutation=='threads_type': row['requested_threads']=True
                else: row['pool_setup_seconds']=-1
                child=out/(mutation+'.py'); child.write_text('print('+repr(json.dumps(row))+')\n')
                record=driver.capture([sys.executable,str(child)],out/mutation)
                self.assertEqual(record['returncode'],0)
                actual=driver.read_json(record['stdout_path'])
                if mutation=='valid': driver.validate_secondary(actual,group)
                else:
                    with self.subTest(mutation=mutation),self.assertRaises((ValueError,KeyError)):
                        driver.validate_secondary(actual,group)

    def test_secondary_real_cli_rejects_a_different_binary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); primary=root/'primary'; primary.mkdir()
            (primary/'provenance.json').write_text(json.dumps({'binary_sha256':'0'*64}))
            (primary/'secondary-references.json').write_text(json.dumps({'primary_binary_sha256':'0'*64,'all_primary_rows_valid':True}))
            result=subprocess.run([sys.executable,driver.__file__,'secondary','--primary',str(primary),
                '--out',str(root/'secondary'),'--binary',sys.executable],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('binary SHA',result.stdout+result.stderr)


class ReviewRegressionTests(unittest.TestCase):
    def test_strict_metadata_rejects_retained_malformed_children(self):
        good, original, _ = fixture()
        mutations = ("actual_threads_bool","actual_threads_float","requested_threads_bool",
            "stdout_schema_bool","stdout_schema_float","snapshot_schema_bool","snapshot_schema_float",
            "shared_requested_threads_bool","shared_actual_threads_bool","shared_schema_bool",
            "shared_nested_energy_bool")
        with retained_test_directory("strict-metadata") as out:
            records=[]
            for mutation in mutations:
                row=json.loads(json.dumps(good)); snap=json.loads(json.dumps(original))
                if mutation=="actual_threads_bool": row["actual_threads"]=True
                elif mutation=="actual_threads_float": row["actual_threads"]=1.0
                elif mutation=="requested_threads_bool": row["requested_threads"]=True
                elif mutation=="stdout_schema_bool": row["schema_version"]=True
                elif mutation=="stdout_schema_float": row["schema_version"]=1.0
                elif mutation=="snapshot_schema_bool": snap["schema_version"]=True
                elif mutation=="snapshot_schema_float": snap["schema_version"]=1.0
                snap["row"]={k:v for k,v in row.items() if k!="pool_setup_seconds"}
                if mutation=="shared_requested_threads_bool": snap["row"]["requested_threads"]=True
                elif mutation=="shared_actual_threads_bool": snap["row"]["actual_threads"]=True
                elif mutation=="shared_schema_bool": snap["row"]["schema_version"]=True
                elif mutation=="shared_nested_energy_bool":
                    snap["row"]["energies"]=dict(row["energies"],temperature=False)
                child=contract_child(out,mutation,row,snap)
                record=driver.run_one(child,out,mutation,"soa",256,1,1,"no-rebuild",.000001)
                records.append(record)
                with self.subTest(mutation=mutation):
                    self.assertFalse(record["valid"],"malformed metadata child was accepted")
                    self.assertEqual(driver.read_json(record["stdout_path"]),row)
                    self.assertTrue(Path(record["snapshot_path"]).is_file())
            self.assertEqual(len((out/"rows.jsonl").read_text().splitlines()),len(mutations))

    def test_zero_elapsed_child_is_retained_and_cannot_crash_summary(self):
        good, original, _ = fixture()
        with retained_test_directory("zero-elapsed") as out:
            records=[]
            for pair in range(5):
                for role in ("A","B"):
                    row=json.loads(json.dumps(good)); snap=json.loads(json.dumps(original))
                    mode="interleaved-serial" if role=="A" else "soa"
                    row.update(mode=mode,evaluated_pairs=None if role=="A" else 0)
                    if pair==0 and role=="A": row["elapsed_seconds"]=0.0
                    snap["row"]={k:v for k,v in row.items() if k!="pool_setup_seconds"}
                    label=f"pair{pair}-{role}"
                    child=contract_child(out,label,row,snap)
                    record=driver.run_one(child,out,label,mode,256,1,1,"no-rebuild",.000001)
                    record.update(phase="measured",pair=pair,role=role,candidate="soa")
                    records.append(record)
            bad=records[0]
            self.assertEqual(driver.read_json(bad["stdout_path"])["elapsed_seconds"],0.0)
            self.assertEqual(driver.sha(bad["stdout_path"]),bad["stdout_sha256"])
            self.assertTrue(Path(bad["snapshot_path"]).is_file())
            self.assertEqual(len((out/"rows.jsonl").read_text().splitlines()),10)
            try:
                result=driver.summary(records)
            except ZeroDivisionError as error:
                self.fail("malformed child crashed summary instead of becoming an invalid row: "+str(error))
            self.assertFalse(bad["valid"])
            self.assertEqual(len(result),1)
            self.assertFalse(result[0]["valid"])
            self.assertFalse(result[0]["timing_pass"])
            self.assertEqual(len(result[0]["improvements"]),4)

    def test_actual_rust_example_stdout_snapshot_subset_contract(self):
        binary,provenance_path,provenance=resolved_rust_example(self)
        with retained_test_directory("actual-rust-contract") as out:
            driver.write_json(out/"binary-resolution.json",dict(binary=str(binary),
                provenance_path=str(provenance_path),binary_sha256=driver.sha(binary),
                compiler=provenance["rustc"],source_sha256=provenance["source_sha256"]))
            records=[]
            for mode in ("interleaved-serial","cluster4"):
                snapshot=out/(mode+".snapshot.json")
                command=[str(binary),"--mode",mode,"--n","256","--threads","1","--steps","1",
                    "--seed","17","--dt","0.000001","--rebuild-policy","no-rebuild",
                    "--snapshot",str(snapshot)]
                result=driver.capture(command,out/mode)
                records.append(result)
                self.assertEqual(result["returncode"],0,Path(result["stderr_path"]).read_text())
                row=driver.read_json(result["stdout_path"]); snap=driver.read_json(snapshot)
                self.assertEqual(set(row)-set(snap["row"]),{"pool_setup_seconds"})
                self.assertEqual(set(snap["row"])-set(row),set())
                shared={k:v for k,v in row.items() if k!="pool_setup_seconds"}
                # Canonical JSON exposes bool/int/float differences that == hides.
                self.assertEqual(json.dumps(shared,sort_keys=True,allow_nan=False),
                    json.dumps(snap["row"],sort_keys=True,allow_nan=False))
                self.assertIn(type(row["pool_setup_seconds"]),(int,float))
                self.assertTrue(math.isfinite(row["pool_setup_seconds"]))
                self.assertGreaterEqual(row["pool_setup_seconds"],0)
                self.assertGreater(row["elapsed_seconds"],0)
                driver.validate(row,snap,dict(mode=mode,n=256,requested_threads=1,
                    steps_requested=1,steps_completed=1,seed=17,dt=.000001,
                    rebuild_policy="no-rebuild"))
            driver.write_json(out/"actual-rust-invocations.json",records)


if __name__ == "__main__":
    unittest.main()
