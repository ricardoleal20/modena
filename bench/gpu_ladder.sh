#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec "${PYTHON:-python3}" "$ROOT/bench/gpu_ladder.py" --steps "${STEPS:-100}" \
  --runs "${RUNS:-3}" --ladder ${LADDER:-864 2048 5000 10000 20000} \
  --threads "${THREADS:-1}" --output "${BENCH_OUT:-$ROOT/bench/results/gpu-$(date +%Y%m%d_%H%M%S)}"
