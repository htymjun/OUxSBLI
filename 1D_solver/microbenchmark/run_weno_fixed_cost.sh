#!/usr/bin/env bash
# Focused driver for the fixed-cost cuts discussed in weno_micro_a100_report.tex.
#
# It runs only the modes needed to answer:
#   * how much seq -> serial costs
#   * how much serial -> warp recovers
#   * how much barrier placement matters
#   * how much shared-memory footprint matters
# and then emits a markdown summary next to summary.csv.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT=""
GPU_CC="${CASE_GPU_CC:-native}"
NX="${NX:-4194304}"
NLAUNCH="${NLAUNCH:-10}"
REPEAT=3

usage() {
  cat <<'EOF'
Usage:
  bash 1D_solver/microbenchmark/run_weno_fixed_cost.sh --out /tmp/weno_fixed

Options:
  --out DIR       output directory
  --gpu-cc CC     CASE_GPU_CC, e.g. 80
  --nx N          default 4194304
  --nlaunch N     default 10
  --repeat N      interleaved rounds, default 3
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="${2:?missing --out value}"; shift 2;;
    --gpu-cc) GPU_CC="${2:?missing --gpu-cc value}"; shift 2;;
    --nx) NX="${2:?missing --nx value}"; shift 2;;
    --nlaunch) NLAUNCH="${2:?missing --nlaunch value}"; shift 2;;
    --repeat) REPEAT="${2:?missing --repeat value}"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  OUT="$SCRIPT_DIR/weno_fixed_cost_${stamp}"
fi

MODES=(
  weight_poly32_seq
  weight_poly32_serial
  weight_poly32_warp
  weight_poly32_serial_oncebar
  weight_poly32_warp_oncebar
  weight_poly32_wsmem_serial
  weight_poly32_wsmem_warp
  weight_poly32_wsmem_tile2_serial
  weight_poly32_wsmem_tile2_warp
  weight_poly32_halfwarp_serial
  weight_poly32_halfwarp_shfl
  weight_poly32_seq7
  weight_poly32_serial7
  weight_poly32_warp7
  weight_poly32_seq9
  weight_poly32_serial9
  weight_poly32_warp9
  w32_poly64_seq
  w32_poly64_serial
  w32_poly64_warp
  w32_poly64_seq7
  w32_poly64_serial7
  w32_poly64_warp7
  w32_poly64_seq9
  w32_poly64_serial9
  w32_poly64_warp9
  w64_only_seq
  poly32_only_seq
  w64_only_seq7
  poly32_only_seq7
  w64_only_seq9
  poly32_only_seq9
)

MODE_STR="${MODES[*]}"

bash "$SCRIPT_DIR/run_weno_nsys.sh" \
  --out "$OUT" \
  --gpu-cc "$GPU_CC" \
  --nx "$NX" \
  --nlaunch "$NLAUNCH" \
  --repeat "$REPEAT" \
  --no-solver-orders \
  --modes "$MODE_STR"

python3 "$SCRIPT_DIR/summarize_weno_fixed_cost.py" \
  "$OUT/summary.csv" \
  --output "$OUT/fixed_cost_summary.md"

echo
echo "wrote:"
echo "  $OUT/summary.csv"
echo "  $OUT/fixed_cost_summary.md"
