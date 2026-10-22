#!/usr/bin/env bash
# Failure-sensitive whole-process runner. PYTHON must be >=3.11; no dependencies.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec "${PYTHON:-python3}" "$ROOT/bench/benchmark.py" \
  --modena "${MODENA:-$ROOT/target/release/modena}" --lammps "${LMP:-lmp_serial}" \
  --steps "${TIMEOUT_STEPS:-1000}" --warmup-steps "${WARMUP_STEPS:-200}" \
  --runs "${RUNS:-3}" --ladder ${LADDER:-864 2048 20000} \
  --threads ${THREADS:-1} --output "${BENCH_OUT:-$ROOT/bench/results/$(date +%Y%m%d_%H%M%S)}"
