#!/usr/bin/env bash
# Thin dispatcher for the 1D microbenchmark Nsight Compute sweeps.
#
# Prefer running the focused scripts directly:
#   bash report/run_keep_ncu.sh --out report/ncu_runs/keep_a100 --gpu-cc 80 --repeat 3
#   bash report/run_slau_ncu.sh --out report/ncu_runs/slau_a100 --gpu-cc 80 --repeat 3
#   bash report/run_weno_ncu.sh --out report/ncu_runs/weno_a100 --gpu-cc 80 --repeat 3
#
# This wrapper can run either or both:
#   bash report/run_summary_ncu.sh --suite all --out report/ncu_runs/a100 --gpu-cc 80 --repeat 3
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE=all
OUT=""
PASS_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  bash report/run_summary_ncu.sh [--suite keep|slau|weno|all] --out DIR [common options]

The old monolithic sweep has been split:
  keep: report/run_keep_ncu.sh
  slau: report/run_slau_ncu.sh
  weno: report/run_weno_ncu.sh

When --suite all is used, --out DIR becomes:
  DIR/keep
  DIR/slau
  DIR/weno

Use --resume with the same --out directory to continue an interrupted run.
All unknown options are passed through to the selected focused script(s).
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --suite)
      SUITE="${2:?missing --suite value}"
      shift 2;;
    --out)
      OUT="${2:?missing --out value}"
      shift 2;;
    -h|--help)
      usage
      exit 0;;
    *)
      PASS_ARGS+=("$1")
      shift;;
  esac
done

case "$SUITE" in
  keep|slau|weno|all)
    ;;
  *)
    echo "bad --suite: $SUITE (use keep, slau, weno, or all)" >&2
    exit 2;;
esac

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  OUT="$SCRIPT_DIR/ncu_runs/summary_${stamp}"
fi

run_keep() {
  bash "$SCRIPT_DIR/run_keep_ncu.sh" --out "$1" "${PASS_ARGS[@]}"
}

run_slau() {
  bash "$SCRIPT_DIR/run_slau_ncu.sh" --out "$1" "${PASS_ARGS[@]}"
}

run_weno() {
  bash "$SCRIPT_DIR/run_weno_ncu.sh" --out "$1" "${PASS_ARGS[@]}"
}

case "$SUITE" in
  keep)
    run_keep "$OUT";;
  slau)
    run_slau "$OUT";;
  weno)
    run_weno "$OUT";;
  all)
    run_keep "$OUT/keep"
    run_slau "$OUT/slau"
    run_weno "$OUT/weno";;
esac
