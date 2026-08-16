#!/usr/bin/env bash
# Thin dispatcher for the 1D microbenchmark Nsight Systems sweeps.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUITE=all
OUT=""
PASS_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  bash report/run_summary_nsys.sh [--suite keep|slau|all] --out DIR [common options]

Focused scripts:
  keep: report/run_keep_nsys.sh
  slau: report/run_slau_nsys.sh

When --suite all is used, --out DIR becomes:
  DIR/keep
  DIR/slau

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
  keep|slau|all)
    ;;
  *)
    echo "bad --suite: $SUITE (use keep, slau, or all)" >&2
    exit 2;;
esac

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  OUT="$SCRIPT_DIR/nsys_runs/summary_${stamp}"
fi

run_keep() {
  bash "$SCRIPT_DIR/run_keep_nsys.sh" --out "$1" "${PASS_ARGS[@]}"
}

run_slau() {
  bash "$SCRIPT_DIR/run_slau_nsys.sh" --out "$1" "${PASS_ARGS[@]}"
}

case "$SUITE" in
  keep)
    run_keep "$OUT";;
  slau)
    run_slau "$OUT";;
  all)
    run_keep "$OUT/keep"
    run_slau "$OUT/slau";;
esac
