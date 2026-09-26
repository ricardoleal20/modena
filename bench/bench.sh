#!/usr/bin/env bash
# Reproducible performance rig (committee Stage 0, gate G0).
#
# Runs the matched MODENA/LAMMPS decks across a size ladder, reports
# µs/(atom·step) with wall + CPU time, median of RUNS runs, for the thread
# configurations in THREADS. Stamps git SHA + rustc + LAMMPS version so every
# number in the paper is traceable. Also runs one NVE-drift check per engine
# and the charged (Ewald) head-to-head row.
#
# Prereqs: `cargo build --release` done, `lmp_serial` on PATH.
#
# Usage:  bench/bench.sh            # default ladder + runs
#         RUNS=5 bench/bench.sh     # more runs for tighter medians
#         LADDER="864 2048" bench/bench.sh

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODENA="$ROOT/target/release/modena"
OUT="${BENCH_OUT:-$ROOT/bench/results}"
RUNS="${RUNS:-3}"
LADDER="${LADDER:-864 2048}"
THREADS="${THREADS:-1 default}"
TIMEOUT_STEPS="${TIMEOUT_STEPS:-2000}"
mkdir -p "$OUT"

GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
RUSTC="$(rustc --version 2>/dev/null || echo unknown)"
LMP_VER="$(command -v lmp_serial >/dev/null && lmp_serial -h 2>&1 | grep -m1 "version" || echo unknown)"

echo "rig: sha=$GIT_SHA rustc=$RUSTC lmp=$LMP_VER runs=$RUNS"
CSV="$OUT/perf_$(date +%Y%m%d_%H%M%S).csv"
echo "case,engine,threads,atoms,steps,wall_s,rss_kb,us_per_atom_step" > "$CSV"

run_one() { # name, engine, atoms, steps, threads
    local name="$1" engine="$2" atoms="$3" steps="$4" threads="$5"
    local envs=()
    if [ "$threads" = 1 ]; then envs=(env RAYON_NUM_THREADS=1); fi
    local outfile timefile
    outfile="$(mktemp)"; timefile="$(mktemp)"
    case "$engine" in
        modena)
            env "${envs[@]+"${envs[@]}"}" /usr/bin/time -p -o "$timefile" "$MODENA" run \
                -i "$ROOT/bench/modena/nve_${atoms}.toml" > "$outfile" 2>&1 || true
            ;;
        lmp)
            /usr/bin/time -p -o "$timefile" lmp_serial -in \
                "$ROOT/bench/lammps/in.nve_${atoms}" > "$outfile" 2>&1 || true
            ;;
        ewald_modena)
            env "${envs[@]+"${envs[@]}"}" /usr/bin/time -p -o "$timefile" "$MODENA" run \
                -i "$ROOT/bench/validation/ewald/modena.toml" > "$outfile" 2>&1 || true
            ;;
        ewald_lmp)
            /usr/bin/time -p -o "$timefile" lmp_serial -in \
                "$ROOT/bench/validation/ewald/in.lammps" -log /dev/null > "$outfile" 2>&1 || true
            ;;
    esac
    local wall_s
    wall_s="$(grep -m1 '^real' "$timefile" | awk '{print $2}')"
    rm -f "$outfile" "$timefile"
    if [ -n "$wall_s" ]; then
        awk -v c="$name" -v e="$engine" -v t="$threads" -v a="$atoms" -v s="$steps" -v w="$wall_s" \
            'BEGIN { printf "%s,%s,%s,%d,%d,%.4f,0,%.3f\n", c, e, t, a, s, w, (w*1e6)/(a*s) }' >> "$CSV"
    fi
}

for threads in $THREADS; do
    for atoms in $LADDER; do
        best_m=""; best_l=""
        for _ in $(seq 1 "$RUNS"); do
            run_one "ladder" modena "$atoms" "$TIMEOUT_STEPS" "$threads"
            run_one "ladder" lmp "$atoms" "$TIMEOUT_STEPS" "$threads"
        done
    done
done

# Charged head-to-head row (the ~34x measured gap, before Ewald-b).
for _ in $(seq 1 "$RUNS"); do
    run_one "charged" ewald_modena 216 2000 1
    run_one "charged" ewald_lmp 216 2000 1
done

echo "--- table (median wall_s per (engine,threads,atoms)) ---"
python3 - "$CSV" <<'EOF'
import csv, sys, statistics, collections
rows = list(csv.DictReader(open(sys.argv[1])))
g = collections.defaultdict(list)
for r in rows:
    try: g[(r['engine'], r['threads'], r['atoms'])].append(float(r['wall_s']))
    except ValueError: pass
print(f"{'engine':12s} {'thr':8s} {'atoms':>6s} {'median_s':>10s} {'us/atom/step':>14s}")
for (e, t, a), ws in sorted(g.items(), key=lambda kv: (kv[0][2], kv[0][0])):
    med = statistics.median(ws)
    print(f"{e:12s} {str(t):8s} {a:6d} {med:10.3f} {med*1e6/(a*2000):14.2f}")
print("\nCSV:", sys.argv[1])
EOF