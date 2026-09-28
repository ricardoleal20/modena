#!/usr/bin/env bash
# GPU-d decision ladder (committee gate G3, docs/COMMITTEE_ACTA.md).
#
# Three engines at five LJ-NVE sizes (864, 2048, 5000, 10000, 20000 atoms):
#   - lmp_serial      (serial, matched decks in bench/lammps/)
#   - MODENA CPU      (release CLI on matched decks in bench/modena/, default
#                      threads — the shipped configuration)
#   - MODENA GPU      (examples/gpu_ladder.rs: `gpu` = WGSL all-pairs kernel,
#                      `gpu-csr` = cell-list gather kernel; both fence-synced:
#                      every step blocks on the staging-buffer map)
# plus the `gather` probe (all-pairs kernel with a 0.01 cutoff — the j-loop
# gather + minimum-image + distance cost with the LJ math skipped), which
# bounds what a cell-list gather can save.
#
# Median (+min) over RUNS measured runs after warm-up discard, per the
# Hoefler&Belli protocol (as in bench.sh). The machine may be load-shared;
# the ladder reports medians and minima for that reason.
#
# Prereqs: cargo build --release --features gpu --example gpu_ladder
#          cargo build --release
#          lmp_serial on PATH.
#
# Usage:  bash bench/gpu_ladder.sh            # RUNS=5 default
#         RUNS=3 LADDER="864 2048" bash bench/gpu_ladder.sh

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODENA="$ROOT/target/release/modena"
LADDER_EX="$ROOT/target/release/examples/gpu_ladder"
OUT="${BENCH_OUT:-$ROOT/bench/results}"
RUNS="${RUNS:-5}"
LADDER="${LADDER:-864 2048 5000 10000 20000}"
mkdir -p "$OUT"

GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
RUSTC="$(rustc --version 2>/dev/null || echo unknown)"
LMP_VER="$(command -v lmp_serial >/dev/null && lmp_serial -h 2>&1 | grep -m1 "version" || echo unknown)"

CSV="$OUT/gpu_ladder_$(date +%Y%m%d_%H%M%S).csv"
echo "engine,atoms,steps,run,wall_s,ms_per_step" > "$CSV"
echo "rig: sha=$GIT_SHA rustc=$RUSTC lmp=$LMP_VER runs=$RUNS load='$(sysctl -n vm.loadavg 2>/dev/null || echo n/a)'"

# Per-size step counts for the in-process example engines (the O(N^2)
# all-pairs kernel is the budget constraint at the top of the ladder; the
# CLI/LAMMPS decks fix their own step counts).
steps_for() {
    case "$1" in
        864) echo 1000 ;;
        2048) echo 1000 ;;
        5000) echo 500 ;;
        10000) echo 250 ;;
        20000) echo 100 ;;
        *) echo 200 ;;
    esac
}

# Deck step counts (for ms/step of the CLI / LAMMPS rows).
deck_steps() {
    case "$1" in
        864|2048) echo 5000 ;;
        *) echo 2000 ;;
    esac
}

for atoms in $LADDER; do
    # --- lmp_serial (fresh process per run; first run = warm-up, discarded) ---
    tmpd="$(mktemp -d)"
    for run in $(seq 0 "$RUNS"); do
        wall="$(cd "$tmpd" && /usr/bin/time -p lmp_serial -in "$ROOT/bench/lammps/in.nve_$atoms" 2>&1 >/dev/null | grep -m1 '^real' | awk '{print $2}')"
        [ "$run" -gt 0 ] && echo "lmp,$atoms,$(deck_steps "$atoms"),$run,$wall,$(awk -v w="$wall" -v s="$(deck_steps "$atoms")" 'BEGIN{printf "%.4f", w*1000/s}')" >> "$CSV"
    done
    rm -rf "$tmpd"

    # --- MODENA CPU (release CLI, default threads; first run = warm-up) ---
    for run in $(seq 0 "$RUNS"); do
        wall="$(/usr/bin/time -p "$MODENA" run -i "$ROOT/bench/modena/nve_$atoms.toml" 2>&1 >/dev/null | grep -m1 '^real' | awk '{print $2}')"
        [ "$run" -gt 0 ] && echo "modena_cpu,$atoms,$(deck_steps "$atoms"),$run,$wall,$(awk -v w="$wall" -v s="$(deck_steps "$atoms")" 'BEGIN{printf "%.4f", w*1000/s}')" >> "$CSV"
    done

    # --- MODENA GPU engines (in-process warm-up + runs; fence-synced) ---
    steps="$(steps_for "$atoms")"
    for eng in gpu gpu-csr gather; do
        # Per-run data lines start with the run number; $6 is the run-mean
        # ms/step, $10 the run drift (kept out of this timing CSV).
        "$LADDER_EX" "$eng" "$atoms" "$steps" "$RUNS" 2>/dev/null | grep -E '^[0-9]+,' | awk -v e="$eng" -F, '{printf "%s,%s,%s,%s,%s,%s\n", e, $3, $4, $1, $5, $6}' >> "$CSV"
    done
    echo "size $atoms done"
done

echo "--- GPU-d ladder (median / min ms per step over $RUNS runs) ---"
python3 - "$CSV" <<'EOF'
import csv, sys, statistics, collections
rows = list(csv.DictReader(open(sys.argv[1])))
g = collections.defaultdict(list)
for r in rows:
    g[(r['engine'], int(r['atoms']))].append(float(r['ms_per_step']))
print(f"{'engine':12s} {'atoms':>6s} {'n':>2s} {'median_ms':>10s} {'min_ms':>10s} {'max_ms':>10s}")
for (e, a), ms in sorted(g.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    print(f"{e:12s} {a:6d} {len(ms):2d} {statistics.median(ms):10.3f} {min(ms):10.3f} {max(ms):10.3f}")
print("\nCSV:", sys.argv[1])
EOF
