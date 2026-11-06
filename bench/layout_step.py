#!/usr/bin/env python3
"""Frozen B03 campaign: direct argv children; standard library only.

Snapshot row is stdout minus pool_setup_seconds (Task1 schema-1 contract).
Each raw snapshot is retained byte-for-byte in gzip with both SHA256 hashes.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import resource
import shutil
import statistics
import subprocess
import time

MODES = ("soa", "aosoa4", "cluster4")
WORKLOAD = "truncated-fcc-a1.6796-splitmix64-jitter0.002-velocity0.1-com-removed-v1"
POTENTIAL = "shifted-force-lj-sigma1-epsilon1-mass1-cutoff2.5-skin0.3"
ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    path = Path(path)
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def finite(value, name):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError("nonfinite or nonnumeric " + name)
    return value


def near(a, b, tol=1e-10):
    finite(a, "actual"); finite(b, "expected")
    if abs(a-b) > tol + 1e-12*abs(b):
        raise ValueError(f"numerical mismatch {a} vs {b}, tolerance {tol}")


def initial(n, seed=17):
    """Independently reconstruct the documented integer RNG/FCC input."""
    mask = (1 << 64)-1
    def sample():
        nonlocal seed
        seed = (seed + 0x9e3779b97f4a7c15) & mask
        z = seed
        z = ((z ^ (z >> 30))*0xbf58476d1ce4e5b9) & mask
        z = ((z ^ (z >> 27))*0x94d049bb133111eb) & mask
        z ^= z >> 31
        return (z >> 11)/(1 << 53)-0.5
    cells = math.ceil((n/4)**(1/3)-1e-12)
    length = 1.6796*cells
    basis = ((0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0))
    p, v = [], []
    for i in range(n):
        cell = i//4
        xyz = (cell//(cells*cells), (cell//cells)%cells, cell%cells)
        for k in range(3):
            p.append((1.6796*(xyz[k]+basis[i%4][k])+.002*sample()) % length)
            v.append(.1*sample())
    com = [sum(v[k::3])/n for k in range(3)]
    v = [x-com[k%3] for k,x in enumerate(v)]
    return p, v, length


def capture(command, prefix):
    """wait4 returns the actual child CPU and RSS, never cumulative maxima."""
    prefix = Path(prefix)
    stdout, stderr = Path(str(prefix)+".stdout"), Path(str(prefix)+".stderr")
    if stdout.exists() or stderr.exists():
        raise FileExistsError(prefix)
    load_before = os.getloadavg()
    start = time.monotonic()
    with stdout.open("xb") as out, stderr.open("xb") as err:
        process = subprocess.Popen(command, stdout=out, stderr=err)
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    return {"command":command, "returncode":process.returncode,
            "stdout_path":str(stdout), "stderr_path":str(stderr),
            "stdout_sha256":sha(stdout), "stderr_sha256":sha(stderr),
            "wall_seconds":time.monotonic()-start,
            "cpu_user_seconds":usage.ru_utime, "cpu_system_seconds":usage.ru_stime,
            "max_rss_bytes":usage.ru_maxrss if os.uname().sysname == "Darwin" else usage.ru_maxrss*1024,
            "load_before":load_before, "load_after":os.getloadavg()}


def validate(row, snap, expected):
    for key, value in expected.items():
        if type(row.get(key)) is not type(value) or row.get(key) != value:
            raise ValueError(f"wrong {key}: {row.get(key)!r} != {value!r}")
    if row["schema_version"] != 1 or snap["schema_version"] != 1:
        raise ValueError("schema")
    if set(row)-set(snap["row"]) != {"pool_setup_seconds"} or set(snap["row"])-set(row):
        raise ValueError("snapshot.row must omit exactly pool_setup_seconds")
    if any(row[k] != v for k,v in snap["row"].items()):
        raise ValueError("snapshot/stdout mismatch")
    for key in ("setup_seconds", "elapsed_seconds", "pool_setup_seconds", "force_max_abs"):
        if finite(row[key], key) < 0: raise ValueError(key)
    for value in row["energies"].values(): finite(value, "energy")
    for key in ("potential", "kinetic", "total", "virial", "temperature", "pressure"):
        finite(row["energies"][key], key)
    if row["finite_efw"] is not True or row["simd_claim"] is not False:
        raise ValueError("physics/SIMD declaration")
    n, steps = row["n"], row["steps_requested"]
    p, v, length = initial(n, row["seed"])
    for a,b in zip(row["box_lengths"], [length]*3): near(a,b)
    if len(row["box_lengths"]) != 3: raise ValueError("box shape")
    if row["workload"] != WORKLOAD or row["potential"] != POTENTIAL: raise ValueError("domain")
    if row["actual_threads"] < 1 or (row["requested_threads"] and row["actual_threads"] != row["requested_threads"]): raise ValueError("threads")
    auto = row["mode"] == "interleaved"
    if row["pair_execution"] != ("production-auto" if auto else "serial-scalar"): raise ValueError("execution")
    if row["comparison_class"] != ("production-threading-control" if auto and row["actual_threads"]>1 else "serial-layout"): raise ValueError("comparison")
    for key in ("setup_rebuilds", "rebuilds", "initial_neighbor_pairs", "final_neighbor_pairs", "neighbor_pair_visits", "final_active_pairs", "retained_packed_bytes", "cluster_reconstructions", "setup_cluster_reconstructions"):
        if type(row[key]) is not int or row[key] < 0: raise ValueError(key)
    if row["setup_rebuilds"] != 1: raise ValueError("setup rebuild")
    forced = row["rebuild_policy"] == "forced"
    if row["rebuilds"] != (steps if forced else 0): raise ValueError("timed rebuild")
    cluster = row["mode"] == "cluster4"
    if row["cluster_reconstructions"] != (steps if forced and cluster else 0): raise ValueError("cluster epoch")
    if row["setup_cluster_reconstructions"] != int(cluster): raise ValueError("setup cluster")
    if not forced and row["neighbor_pair_visits"] != steps*row["initial_neighbor_pairs"]: raise ValueError("visits")
    if row["neighbor_pair_visits"] < row["final_active_pairs"]*steps: raise ValueError("visits too few")
    if row["mode"].startswith("interleaved"):
        if row["evaluated_pairs"] is not None: raise ValueError("unavailable evaluated count")
    elif type(row["evaluated_pairs"]) is not int or not 0 <= row["evaluated_pairs"] <= row["neighbor_pair_visits"]: raise ValueError("evaluated count")
    for name in ("initial_input", "initial_force_ready", "final"):
        state = snap[name]
        if state["n"] != n or state["logical_ids"] != list(range(n)) or state["periodic"] != [True]*3: raise ValueError("logical identity")
        for key, expected_values in (("masses",[1.]*n),("species",[0]*n),("charges",[0.]*n),("atomic_numbers",[0]*n)):
            if state[key] != expected_values: raise ValueError(key)
        if state["steps_completed"] != (steps if name == "final" else 0): raise ValueError("state steps")
        if state["cell_vectors"] != [[length if i==j else 0. for j in range(3)] for i in range(3)]: raise ValueError("cell")
        for key in ("positions", "velocities", "forces"):
            if len(state[key]) != 3*n: raise ValueError("state shape")
            for x in state[key]: finite(x,key)
        if name != "final":
            for a,b in zip(state["positions"], p): near(a,b)
            for a,b in zip(state["velocities"], v): near(a,b)
        if name == "initial_input":
            if any(state["forces"]) or state["energies"] is not None or state["neighbor_pairs"] is not None or state["pair_scales"] is not None: raise ValueError("input state")
        else:
            for key in ("potential", "kinetic", "total", "virial", "temperature", "pressure"): finite(state["energies"][key], key)
            pairs = state["neighbor_pairs"]
            if len(pairs) != len(set(map(tuple,pairs))) or any(type(i) is not int or type(j) is not int or not 0 <= i < j < n for i,j in pairs): raise ValueError("half pair identity")
            if state["pair_scales"] != []: raise ValueError("unexpected scales")
            count_key = "initial_neighbor_pairs" if name == "initial_force_ready" else "final_neighbor_pairs"
            if len(pairs) != row[count_key]: raise ValueError("pair count")
    if snap["final"]["energies"] != row["energies"]: raise ValueError("final energies")
    near(row["force_max_abs"], max(map(abs,snap["final"]["forces"])))
    for name,key in (("initial_input","initial_checksum_fnv1a64"),("final","final_checksum_fnv1a64")):
        if snap[name]["checksum_fnv1a64"] != row[key]: raise ValueError("checksum identity")


def compare(a, b):
    """Compare every physics field; only explicit checksums may differ."""
    for name in ("initial_input", "initial_force_ready", "final"):
        x,y = a[name],b[name]
        if set(x) != set(y): raise ValueError("unknown state fields differ")
        for key in x:
            if key == "checksum_fnv1a64": continue
            if key in ("forces", "positions", "velocities"):
                for u,v in zip(x[key],y[key]): near(u,v)
            elif key == "energies" and x[key] is not None:
                if set(x[key]) != set(y[key]): raise ValueError("energy fields")
                for k in x[key]: near(x[key][k],y[key][k],1e-8)
            elif x[key] != y[key]: raise ValueError("state mismatch "+key)


def audit(state, exact_neighbors):
    """Independent linked-cell half-pair oracle, no production helpers.

    Enumerates every neighboring cell, checks exact skin pairs at build states,
    and active-pair coverage for reused lists (which include valid old skin pairs).
    """
    n=state["n"]; p=state["positions"]; length=state["cell_vectors"][0][0]
    nc=max(1,int(length/2.8)); width=length/nc
    buckets={}
    for i in range(n):
        cell=tuple(int((p[3*i+k]%length)/width)%nc for k in range(3))
        buckets.setdefault(cell,[]).append(i)
    expected=set(); active=set(); f=[0.]*(3*n); e=w=0.
    rc=2.5; uc=4*(rc**-12-rc**-6); fc=24*(2*rc**-13-rc**-7)
    for cell, ids in buckets.items():
        nearby={tuple((cell[k]+d[k])%nc for k in range(3)) for d in itertools.product((-1,0,1),repeat=3)}
        for i in ids:
            for other in nearby:
                for j in buckets.get(other,()):
                    if j<=i: continue
                    d=[p[3*i+k]-p[3*j+k] for k in range(3)]
                    d=[x-length*math.floor(x/length+.5) for x in d]
                    r2=sum(x*x for x in d)
                    if r2<2.8**2: expected.add((i,j))
                    if not 1e-14<=r2<rc*rc: continue
                    active.add((i,j)); r=math.sqrt(r2)
                    e+=4*(r**-12-r**-6)-uc+(r-rc)*fc
                    factor=(24*(2*r**-13-r**-7)-fc)/r
                    for k in range(3):
                        fv=d[k]*factor; f[3*i+k]+=fv; f[3*j+k]-=fv; w+=d[k]*fv
    actual=set(map(tuple,state["neighbor_pairs"]))
    if (exact_neighbors and actual != expected) or not active <= actual: raise ValueError("independent neighbor set")
    near(state["energies"]["potential"],e,1e-8); near(state["energies"]["virial"],w,1e-8)
    for x,y in zip(state["forces"],f): near(x,y)
    kinetic=.5*sum(v*v for v in state["velocities"])
    near(state["energies"]["kinetic"],kinetic,1e-8)
    near(state["energies"]["total"],kinetic+e,1e-8)
    near(state["energies"]["temperature"],2*kinetic/(3*n),1e-8)
    near(state["energies"]["pressure"],(2*kinetic+w)/(3*length**3),1e-8)
    return {"active_pairs":len(active), "skin_pairs":len(expected), "max_force_error":max(abs(x-y) for x,y in zip(state["forces"],f)), "energy_error":abs(state["energies"]["potential"]-e), "virial_error":abs(state["energies"]["virial"]-w)}


def run_one(binary, out, label, mode, n, threads, steps, policy, dt):
    snapshot=out/(label+".snapshot.json")
    if snapshot.exists() or Path(str(snapshot)+".gz").exists(): raise FileExistsError(snapshot)
    command=[str(binary),"--mode",mode,"--n",str(n),"--threads",str(threads),"--steps",str(steps),"--seed","17","--dt",str(dt),"--rebuild-policy",policy,"--snapshot",str(snapshot)]
    record=capture(command,out/label)
    record.update(label=label, mode=mode, n=n, threads=threads, steps=steps, policy=policy, dt=dt)
    try:
        if record["returncode"]: raise ValueError("child failure")
        row=read_json(record["stdout_path"])
        snap=read_json(snapshot)
        validate(row,snap,{"mode":mode,"n":n,"requested_threads":threads,"steps_requested":steps,"steps_completed":steps,"seed":17,"dt":dt,"rebuild_policy":policy})
        record["row"]=row; record["valid"]=True
    except Exception as error:
        record["valid"]=False; record["error"]=repr(error)
    if snapshot.exists():
        raw=snapshot.read_bytes(); compressed=gzip.compress(raw,compresslevel=1,mtime=0)
        packed=Path(str(snapshot)+".gz"); packed.write_bytes(compressed)
        if gzip.decompress(packed.read_bytes()) != raw: raise ValueError("gzip roundtrip")
        record.update(snapshot_path=str(packed), snapshot_raw_sha256=hashlib.sha256(raw).hexdigest(), snapshot_gzip_sha256=sha(packed), snapshot_raw_bytes=len(raw))
        snapshot.unlink()  # Exact bytes remain in gzip, checked before removing duplicate.
    with (out/"rows.jsonl").open("a") as file: file.write(json.dumps(record,allow_nan=False)+"\n")
    print(label,record["valid"],record.get("row",{}).get("elapsed_seconds"),flush=True)
    return record


def summary(rows):
    groups={}
    for r in rows:
        if r.get("phase") not in ("measured","reproduction"): continue
        key=(r["n"],r["threads"],r["policy"],r["candidate"],r["phase"])
        groups.setdefault(key,{}).setdefault(r["pair"],{})[r["role"]]=r
    result=[]
    for key,pairs in groups.items():
        ratios=[]; baseline=[]; candidate=[]; valid=len(pairs)==5
        for pair in pairs.values():
            valid &= set(pair)=={"A","B"} and all(r["valid"] for r in pair.values())
            if set(pair)=={"A","B"} and all(r["valid"] for r in pair.values()): 
                baseline.append(pair["A"]["row"]["elapsed_seconds"]); candidate.append(pair["B"]["row"]["elapsed_seconds"])
                ratios.append(1-candidate[-1]/baseline[-1])
        med=(1-statistics.median(candidate)/statistics.median(baseline)) if ratios else None
        result.append(dict(zip(("n","threads","policy","candidate","phase"),key), valid=bool(valid), improvements=ratios, median_improvement=med, median_paired_improvement=statistics.median(ratios) if ratios else None, median_baseline_seconds=statistics.median(baseline) if baseline else None, median_candidate_seconds=statistics.median(candidate) if candidate else None, favorable=sum(x>0 for x in ratios), timing_pass=bool(valid and med>=.05 and sum(x>0 for x in ratios)>=4)))
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("action", choices=("pilot","campaign","audit"))
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--binary",type=Path)
    parser.add_argument("--config",type=Path)
    args=parser.parse_args(); out=args.out.resolve()
    if args.action=="audit":
        rows=[json.loads(x) for x in (out/"index.jsonl").read_text().splitlines()]
        canonical={}; oracle={}; audits=[]
        for r in rows:
            try:
                if not r["valid"]: raise ValueError(r.get("error"))
                snap=read_json(r["snapshot_path"])
                for kind in ("stdout", "stderr"):
                    if sha(r[kind+"_path"]) != r[kind+"_sha256"]: raise ValueError("raw output SHA")
                stdout=read_json(r["stdout_path"])
                if stdout != r["row"]: raise ValueError("recorded stdout")
                validate(stdout,snap,{"mode":r["mode"],"n":r["n"],"requested_threads":r["threads"],"steps_requested":r["steps"],"steps_completed":r["steps"],"seed":17,"dt":r["dt"],"rebuild_policy":r["policy"]})
                raw=gzip.decompress(Path(r["snapshot_path"]).read_bytes())
                if hashlib.sha256(raw).hexdigest()!=r["snapshot_raw_sha256"] or sha(r["snapshot_path"])!=r["snapshot_gzip_sha256"]: raise ValueError("snapshot SHA")
                key=(r["n"],r["steps"],r["policy"],r["dt"])
                if key not in canonical:
                    canonical[key]=r["snapshot_path"]
                    oracle[key]={"initial":audit(snap["initial_force_ready"],True),"final":audit(snap["final"],r["policy"]=="forced")}
                else: compare(snap,read_json(canonical[key]))
                near(r["row"]["final_active_pairs"],oracle[key]["final"]["active_pairs"],0)
                audits.append({"label":r["label"],"valid":True,"oracle_reference":canonical[key]})
            except Exception as error:
                r["valid"]=False; audits.append({"label":r["label"],"valid":False,"error":repr(error)})
            if len(audits)%25==0: print("audited",len(audits),flush=True)
        write_json(out/"audit.json",{"rows":audits,"oracles":[{"case":k,"result":v} for k,v in oracle.items()]})
        write_json(out/"summary-audited.json",summary(rows))
        return
    out.mkdir(parents=True,exist_ok=False)
    binary=out/"layout_step"; shutil.copy2(args.binary,binary)
    driver=out/"layout_step.py"; shutil.copy2(__file__,driver)
    sources=subprocess.check_output(["git","ls-files"],cwd=ROOT,text=True).splitlines()
    provenance={"binary_sha256":sha(binary),"driver_sha256":sha(driver),"loaded_driver_sha256":sha(__file__),"source_sha256":{p:sha(ROOT/p) for p in sources if (ROOT/p).is_file()},"head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),"rustc":subprocess.check_output(["rustc","-vV"],text=True),"platform":os.uname()._asdict() if hasattr(os.uname(),"_asdict") else list(os.uname()),"cpu":subprocess.check_output(["sysctl","-n","machdep.cpu.brand_string"],text=True).strip(),"free_disk_bytes":shutil.disk_usage(out).free}
    write_json(out/"provenance.json",provenance)
    if args.action=="pilot":
        for n in (864,2048,5000,10000):
            for policy in ("no-rebuild","forced"):
                for mode in ("interleaved-serial",*MODES):
                    run_one(binary,out,f"pilot-{n}-{policy}-{mode}",mode,n,1,20,policy,.000001)
        return
    config=read_json(args.config); write_json(out/"config.json",config)
    collected=[]
    def invocation(n,t,policy,candidate,phase,pair,role,mode):
        steps=config["steps"][str(n)][policy]
        label=f"{phase}-{n}-t{t}-{policy}-{candidate}-p{pair}-{role}"
        r=run_one(binary,out,label,mode,n,t,steps,policy,config["dt"])
        r.update(phase=phase,pair=pair,role=role,candidate=candidate)
        # Metadata is separate, preserving the raw record unchanged.
        with (out/"index.jsonl").open("a") as file: file.write(json.dumps(r)+"\n")
        collected.append(r)
    for n,t,policy in itertools.product((864,2048,5000,10000),(1,4,8,0),("no-rebuild","forced")):
        for candidate in MODES:
            for role,mode in (("A","interleaved-serial"),("B",candidate)):
                invocation(n,t,policy,candidate,"warmup",0,role,mode)
            for pair in range(5):
                roles=(("A","interleaved-serial"),("B",candidate))
                for role,mode in (roles if pair%2==0 else reversed(roles)):
                    invocation(n,t,policy,candidate,"measured",pair,role,mode)
        for repeat in range(6):
            invocation(n,t,policy,"interleaved","production-warmup" if repeat==0 else "production-control",repeat,"control","interleaved")
    write_json(out/"summary-preaudit.json",summary(collected))
    for case in list(summary(collected)):
        if not case["timing_pass"]: continue
        n,t,policy,candidate=(case[k] for k in ("n","threads","policy","candidate"))
        for role,mode in (("A","interleaved-serial"),("B",candidate)):
            invocation(n,t,policy,candidate,"reproduction-warmup",0,role,mode)
        for pair in range(5):
            roles=(("A","interleaved-serial"),("B",candidate))
            for role,mode in (reversed(roles) if pair%2==0 else roles):
                invocation(n,t,policy,candidate,"reproduction",pair,role,mode)
    write_json(out/"summary-preaudit.json",summary(collected))
    # Audit consumes the metadata-bearing index, raw rows always remain intact.
    write_json(out/"campaign-integrity.json",{"binary_unchanged":sha(binary)==provenance["binary_sha256"],"driver_unchanged":sha(__file__)==provenance["loaded_driver_sha256"],"sources_unchanged":all(sha(ROOT/p)==h for p,h in provenance["source_sha256"].items()),"snapshot_count":len(collected)})


if __name__=="__main__": main()
