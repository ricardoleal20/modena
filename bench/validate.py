#!/usr/bin/env python3
"""Run prebuilt scientific_validation example and retain numerical evidence."""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys

from benchmark import ROOT, capture, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', default=str(ROOT / 'target/release/examples/scientific_validation'))
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    command = [str(Path(args.binary).resolve())] + (['--extended'] if args.extended else [])
    try:
        result = subprocess.run(command, text=True, capture_output=True, env={**os.environ, 'RAYON_NUM_THREADS': '1'})
    except OSError as error:
        result = subprocess.CompletedProcess(command, 127, '', str(error))
    (out / 'validation.stdout').write_text(result.stdout)
    (out / 'validation.stderr').write_text(result.stderr)
    try:
        data = json.loads(result.stdout)
        observations = data['observations']
        if not observations:
            raise ValueError('no numerical observations')
        for row in observations:
            values = [row['observed'], row['expected'], row['absolute_tolerance']]
            finite = all(isinstance(x, (int, float)) and math.isfinite(x) for x in values)
            row['passed'] = finite and row['absolute_tolerance'] >= 0 and abs(row['observed'] - row['expected']) <= row['absolute_tolerance']
            if not finite:
                for key in ('observed', 'expected', 'absolute_tolerance'):
                    if not math.isfinite(row[key]):
                        row[key] = None
        data['passed'] = result.returncode == 0 and all(r['passed'] for r in observations)
    except (ValueError, KeyError, TypeError) as error:
        data = {'passed': False, 'error': str(error), 'observations': []}
    data['provenance'] = dict(created_utc=datetime.now(timezone.utc).isoformat(), command=command, returncode=result.returncode, binary_sha256=sha(command[0]) if Path(command[0]).is_file() else None,
        git_commit=capture(['git', '-C', str(ROOT), 'rev-parse', 'HEAD']),
        git_status=capture(['git', '-C', str(ROOT), 'status', '--porcelain']),
        rustc=capture(['rustc', '--version']), python=sys.version, python_executable=sys.executable,
        platform=platform.platform(), threads=1)
    (out / 'validation.json').write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    print(f"{len(data['observations'])} observations; passed={data['passed']}; {out / 'validation.json'}")
    return 0 if data['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
