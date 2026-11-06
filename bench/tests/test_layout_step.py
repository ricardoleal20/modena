"""Real child processes exercise malformed output, rather than subprocess mocks."""
import importlib.util
import json
import math
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


if __name__ == "__main__":
    unittest.main()
