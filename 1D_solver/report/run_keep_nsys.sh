#!/usr/bin/env bash
# KEEP-only Nsight Systems sweep for environments where ncu replay is unstable.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/_nsys_common.sh"

MEASURE_FUSED=0
CHECK_FUSED_SASS=1

usage() {
  cat <<'EOF'
Usage:
  bash report/run_keep_nsys.sh --out report/nsys_runs/keep_a100 --gpu-cc 80 --repeat 3

KEEP sweep:
  ORDER=2,4,6
  KEEP_PREC=fp64,fp32,df
  VISC_PREC=fp32 for all KEEP_PREC values
  modes=seq,warp
  fused is skipped when its normalized SASS matches seq.

KEEP-specific options:
  --measure-fused       always profile fused with nsys
  --no-sass-skip        do not compare seq/fused SASS; same as --measure-fused
EOF
  nsys_common_usage_tail
}

nsys_common_defaults
while [ $# -gt 0 ]; do
  case "$1" in
    --measure-fused)
      MEASURE_FUSED=1
      CHECK_FUSED_SASS=0
      shift;;
    --no-sass-skip)
      MEASURE_FUSED=1
      CHECK_FUSED_SASS=0
      shift;;
    -h|--help)
      usage
      exit 0;;
    *)
      nsys_parse_common_arg "$@" || rc=$?
      rc=${rc:-0}
      if [ "$rc" -eq 0 ]; then
        echo "unknown argument: $1" >&2
        usage >&2
        exit 2
      fi
      shift "$rc"
      unset rc;;
  esac
done

nsys_common_init keep

run_keep_case() {
  local order="$1" keep_prec="$2" visc_prec case_prefix
  visc_prec="$(ncu_visc_for_keep_precision "$keep_prec")"
  case_prefix="keep_${keep_prec}_o${order}"

  nsys_run_case keep "${case_prefix}_seq" seq "$order" "$order" "$keep_prec" "$visc_prec" fp64 fp64 fp64 fp64

  if [ "$MEASURE_FUSED" -eq 1 ]; then
    nsys_run_case keep "${case_prefix}_fused" fused "$order" "$order" "$keep_prec" "$visc_prec" fp64 fp64 fp64 fp64
  elif [ "$CHECK_FUSED_SASS" -eq 1 ]; then
    if ncu_sass_same_seq_fused "$order" "$order" "$keep_prec" "$visc_prec" fp64; then
      echo "[$(date +%H:%M:%S)] skip fused ${case_prefix}: SASS matches seq" | tee -a "$OUT/progress.log"
      ncu_note_skipped_fused keep "${case_prefix}_fused" "$order" "$order" "$keep_prec" "$visc_prec" fp64
    else
      echo "[$(date +%H:%M:%S)] profile fused ${case_prefix}: SASS differs from seq" | tee -a "$OUT/progress.log"
      nsys_run_case keep "${case_prefix}_fused" fused "$order" "$order" "$keep_prec" "$visc_prec" fp64 fp64 fp64 fp64
    fi
  fi

  nsys_run_case keep "${case_prefix}_warp" warp "$order" "$order" "$keep_prec" "$visc_prec" fp64 fp64 fp64 fp64
}

for order in 2 4 6; do
  for keep_prec in fp64 fp32 df; do
    run_keep_case "$order" "$keep_prec"
  done
done

nsys_finish
