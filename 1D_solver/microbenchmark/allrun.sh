#!/usr/bin/env bash
# One-shot driver: run the full microbenchmark sweep, then post-process the
# same result directory with the static SASS/PTXAS summary and the fixed-cost
# markdown summary. Solver-side profiling is opt-in.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT=""
GPU_CC="${CASE_GPU_CC:-native}"
NX="${NX:-4194304}"
NLAUNCH="${NLAUNCH:-10}"
REPEAT=3
WITH_SOLVER=0
RESUME=0

usage() {
  cat <<'EOF'
Usage:
  bash 1D_solver/microbenchmark/allrun.sh --out /tmp/weno_all --gpu-cc 80

Options:
  --out DIR          result directory
  --gpu-cc CC        CASE_GPU_CC, e.g. 80
  --nx N             default 4194304
  --nlaunch N        default 10
  --repeat N         default 3
  --solver-orders    also run solver_weno_orders (opt-in; uses report/run_weno_nsys.sh)
  --resume           reuse an existing result directory
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="${2:?missing --out value}"; shift 2;;
    --gpu-cc) GPU_CC="${2:?missing --gpu-cc value}"; shift 2;;
    --nx) NX="${2:?missing --nx value}"; shift 2;;
    --nlaunch) NLAUNCH="${2:?missing --nlaunch value}"; shift 2;;
    --repeat) REPEAT="${2:?missing --repeat value}"; shift 2;;
    --solver-orders) WITH_SOLVER=1; shift;;
    --resume) RESUME=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  OUT="$SCRIPT_DIR/weno_allrun_${stamp}"
fi

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

run_args=(
  --out "$OUT"
  --gpu-cc "$GPU_CC"
  --nx "$NX"
  --nlaunch "$NLAUNCH"
  --repeat "$REPEAT"
)
[ "$WITH_SOLVER" -eq 1 ] && run_args+=(--solver-orders)
[ "$RESUME" -eq 1 ] && run_args+=(--resume)

bash "$SCRIPT_DIR/run_weno_nsys.sh" "${run_args[@]}"

python3 "$SCRIPT_DIR/analyze_weno_static.py" \
  --exe "$SCRIPT_DIR/build/weno_micro" \
  --build-log "$OUT/build.log" \
  --output "$OUT/static_instruction_summary.csv" \
  > "$OUT/static_instruction_summary.stdout.csv"

python3 "$SCRIPT_DIR/summarize_weno_fixed_cost.py" \
  "$OUT/summary.csv" \
  --output "$OUT/fixed_cost_summary.md" \
  > "$OUT/fixed_cost_summary.stdout.md"

cat > "$OUT/README.allrun.txt" <<EOF
allrun output directory: $OUT

primary files:
  summary.csv
  build.log
  static_instruction_summary.csv
  fixed_cost_summary.md

optional:
  solver_weno_orders/summary.csv   only when --solver-orders is used
EOF

echo
echo "allrun complete:"
echo "  $OUT/summary.csv"
echo "  $OUT/static_instruction_summary.csv"
echo "  $OUT/fixed_cost_summary.md"
