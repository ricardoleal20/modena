#!/usr/bin/env python3
"""Retain experimental host-offload ladder logs, precision and timing scopes."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

from benchmark import ROOT, capture, resources, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', default=str(ROOT / 'target/release/examples/gpu_ladder'))
    parser.add_argument('--timer', default='/usr/bin/time', help='resource timer path; absent timer runs without CPU/RSS measurements')
    parser.add_argument('--ladder', nargs='+', type=int, default=[864, 2048, 5000, 10000, 20000])
    parser.add_argument('--engines', nargs='+', choices=['cpu', 'gpu', 'gpu-csr', 'gather'], default=['cpu', 'gpu', 'gpu-csr', 'gather'])
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if min(args.steps,args.runs,args.threads) <= 0:
        parser.error('steps, runs, threads must be positive')
    out = args.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    if (out / 'gpu.json').exists(): parser.error('use a fresh output directory')
    binary = str(Path(args.binary).resolve())
    try:
        binary_hash, hash_error = sha(binary), None
    except OSError as error:
        binary_hash, hash_error = None, str(error)
    artifact = dict(provenance=dict(binary=binary, binary_sha256=binary_hash, binary_hash_error=hash_error, git_commit=capture(['git','-C',str(ROOT),'rev-parse','HEAD']),
        rustc=capture(['rustc','--version']),python=sys.version,platform=platform.platform(),logical_cpus=os.cpu_count(),
        hardware=capture(['sysctl','-n','machdep.cpu.brand_string']) if platform.system()=='Darwin' else capture(['lscpu'])),
        scope='experimental host-offload; row integration wall excludes init and initial force; outer process resources include them',
        warmup='one fresh trajectory in the same process, discarded; not continued sampling',
        potential='shifted-force LJ sigma=epsilon=mass=1 rc=2.5 (gather rc=.01 is a different Hamiltonian)',
        seed=87287,dt=.005,skin=.3,runs=[])
    (out/'gpu.json').write_text(json.dumps(artifact,indent=2)+'\n')
    for atoms in args.ladder:
        for engine in args.engines:
            command=[binary,engine,str(atoms),str(args.steps),str(args.runs)]
            timer = args.timer if Path(args.timer).is_file() else None
            timed=[timer,'-l' if platform.system()=='Darwin' else '-v',*command] if timer else command
            start=time.perf_counter()
            try:
                result=subprocess.run(timed,env={**os.environ,'RAYON_NUM_THREADS':str(args.threads)},text=True,capture_output=True)
            except OSError as error:
                result=subprocess.CompletedProcess(timed,127,'',str(error))
            label=f'{engine}-{atoms}'
            (out/(label+'.stdout')).write_text(result.stdout);(out/(label+'.stderr')).write_text(result.stderr)
            user,system,rss=resources(result.stderr) if timer else (None,None,None)
            lines=[line for line in result.stdout.splitlines() if line.startswith('run,') or line[:1].isdigit()]
            rows=list(csv.DictReader(io.StringIO('\n'.join(lines))))
            artifact['runs'].append(dict(engine=engine,atoms=atoms,steps=args.steps,repetitions=args.runs,threads=args.threads,
                dof=3*atoms-3,conserved_com_dof=3,precision='f64' if engine=='cpu' else 'f32 kernel/Kahan; f64 host integration',
                command=command,launch_command=timed,resource_timer=timer,returncode=result.returncode,outer_wall_seconds=time.perf_counter()-start,
                user_cpu_seconds=user,system_cpu_seconds=system,rss_bytes=rss,
                stdout=label+'.stdout',stderr=label+'.stderr',integration_rows=rows))
            (out/'gpu.json').write_text(json.dumps(artifact,indent=2)+'\n')
            if result.returncode or len(rows)!=args.runs:
                print(f'FAILED {label}: exit={result.returncode}, rows={len(rows)}; logs retained',file=sys.stderr)
                return 1
    print(out/'gpu.json')
    return 0


if __name__=='__main__':
    sys.exit(main())
