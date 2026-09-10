#!/usr/bin/env bash
# MODENA-vs-LAMMPS comparison harness.
#
# Runs each matched (MODENA toml, LAMMPS in.*) pair, capturing wall time and
# peak RSS via /usr/bin/time, then prints a comparison table. Outputs a CSV at
# bench/results/comparison.csv.
#
# Usage:  bench/run_comparison.sh
#
# Prereqs: `lmp_serial` on PATH (Homebrew), `cargo build --release` done.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LMP="${LMP:-lmp_serial}"
MODENA="$ROOT/target/release/modena"
RESULTS="$ROOT/bench/results"
mkdir -p "$RESULTS"

if ! command -v "$LMP" >/dev/null 2>&1; then
  echo "ERROR: '$LMP' not found. Install LAMMPS or set LMP=/path/to/lmp." >&2
  exit 1
fi
if [[ ! -x "$MODENA" ]]; then
  echo "ERROR: $MODENA not built. Run 'cargo build --release' first." >&2
  exit 1
fi

# macOS /usr/bin/time -l prints peak RSS in bytes; GNU time -v prints in KB.
# Key off the OS directly — parsing `time` output under `set -o pipefail` is
# fragile (the grep pipeline can race the resource-list flush).
RSS_OS="$(uname)"
case "$RSS_OS" in
  Darwin) RSS_OS="mac" ;;
  *)      RSS_OS="linux" ;;
esac

run_one() {
  # Echoes "wall_s,rss_kb" to stdout (captured by caller). Progress to stderr.
  local label="$1" engine="$2" cmd="$3" workdir="$4"
  local tmp
  tmp="$RESULTS/${label}_${engine}.time"
  if [[ "$RSS_OS" == "mac" ]]; then
    # macOS /usr/bin/time -l: prints "real/user/sys" then a resource list with
    # "maximum resident set size" in BYTES.
    ( cd "$workdir" && /usr/bin/time -l $cmd ) >"$RESULTS/${label}_${engine}.out" 2>"$tmp" || {
      echo "  [$engine] FAILED — see $tmp" >&2; return 1
    }
    local real rss
    real=$(awk '/real/{print $1; exit}' "$tmp")
    rss=$(awk '/maximum resident set size/ {print $1; exit}' "$tmp")
    rss=$(( rss / 1024 ))   # bytes → KB
    echo "  [$engine] ${real}s  (peak RSS ${rss} KB)" >&2
    echo "${real},${rss}"
  else
    # GNU /usr/bin/time -v: "Elapsed (wall) clock" + "Maximum resident set" in KB.
    ( cd "$workdir" && /usr/bin/time -v $cmd ) >"$RESULTS/${label}_${engine}.out" 2>"$tmp" || {
      echo "  [$engine] FAILED — see $tmp" >&2; return 1
    }
    local real rss
    real=$(awk '/Elapsed/{print $NF}' "$tmp")
    rss=$(awk '/Maximum resident set/{print $NF}' "$tmp")
    echo "  [$engine] ${real}s  (peak RSS ${rss} KB)" >&2
    echo "${real},${rss}"
  fi
}

echo "=== MODENA vs LAMMPS comparison ==="
echo "MODENA: $($MODENA --version 2>/dev/null || echo 'modena')"
echo "LAMMPS: $($LMP -h 2>/dev/null | head -1)"
echo

csv="$RESULTS/comparison.csv"
echo "case,engine,atoms,steps,wall_s,rss_kb" > "$csv"

# ---- Case 1: NVE 864 ----
echo "Case 1: NVE LJ liquid, 864 atoms, 5000 steps"
r_mod=$(run_one "nve_864" "modena" "$MODENA run -i $ROOT/bench/modena/nve_864.toml" "$RESULTS")
r_lmp=$(run_one "nve_864" "lammps" "$LMP -in $ROOT/bench/lammps/in.nve_864" "$RESULTS")
echo "nve_864,modena,864,5000,$r_mod" >> "$csv"
echo "nve_864,lammps,864,5000,$r_lmp" >> "$csv"
echo

# ---- Case 2: NVE 2048 ----
echo "Case 2: NVE LJ liquid, 2048 atoms, 5000 steps"
r_mod=$(run_one "nve_2048" "modena" "$MODENA run -i $ROOT/bench/modena/nve_2048.toml" "$RESULTS")
r_lmp=$(run_one "nve_2048" "lammps" "$LMP -in $ROOT/bench/lammps/in.nve_2048" "$RESULTS")
echo "nve_2048,modena,2048,5000,$r_mod" >> "$csv"
echo "nve_2048,lammps,2048,5000,$r_lmp" >> "$csv"
echo

# ---- Case 3: scaling sweep ----
for n in 4 6 8 10; do
  atoms=$((4*n*n*n))
  echo "Case 3 (scaling): NVE, $atoms atoms (n=$n), 2000 steps"
  r_mod=$(run_one "scale_$n" "modena" "$MODENA run -i $ROOT/bench/modena/nve_scaling_$n.toml" "$RESULTS")
  r_lmp=$(run_one "scale_$n" "lammps" "$LMP -var n $n -in $ROOT/bench/lammps/in.nve_scaling" "$RESULTS")
  echo "scale_$n,modena,$atoms,2000,$r_mod" >> "$csv"
  echo "scale_$n,lammps,$atoms,2000,$r_lmp" >> "$csv"
  echo
done

echo "=== Summary (raw CSV: $csv) ==="
column -t -s, "$csv"
echo
echo "steps/s = steps / wall_s  |  speedup = modena_wall / lammps_wall"
