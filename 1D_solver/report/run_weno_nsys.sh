#!/usr/bin/env bash
# WENO-Z Nsight Systems sweep across the WENO_ORDER dimension.
#
# This is the A100-friendly companion to run_weno_ncu.sh. It keeps the same
# case set but avoids ncu replay, and writes .nsys-rep plus CSV/SQLite exports
# that can be parsed by the report scripts and by AI agents.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/_nsys_common.sh"

usage() {
  cat <<'EOF'
Usage:
  bash report/run_weno_nsys.sh --out report/nsys_runs/weno_a100 --gpu-cc 80 --repeat 3

Sweep (ORDER=6 throughout; WENO_ORDER selects the WENO-Z width):

  WENO_ORDER=5   split / smem / warp FP64 baselines, warp all-DF, warp per-variable
  WENO_ORDER=7   same
  WENO_ORDER=9   same
EOF
  nsys_common_usage_tail
}

nsys_common_defaults
while [ $# -gt 0 ]; do
  case "$1" in
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

nsys_common_init weno

for WO in 5 7 9; do
  export BENCH_WENO_ORDER="$WO"
  nsys_run_case "slau_weno${WO}" "split_fp64_o6"   slau_weno      6 6 fp64 fp64 fp64 fp64 fp64 fp64
  nsys_run_case "slau_weno${WO}" "smem_fp64_o6"    slau_weno_smem 6 6 fp64 fp64 fp64 fp64 fp64 fp64
  nsys_run_case "slau_weno${WO}" "warp_fp64_o6"    slau_weno_warp 6 6 fp64 fp64 fp64 fp64 fp64 fp64
  nsys_run_case "slau_weno${WO}" "warp_all_df_o6"  slau_weno_warp 6 6 fp64 fp64 fp64 df   df   df
  nsys_run_case "slau_weno${WO}" "warp_rho64_o6"   slau_weno_warp 6 6 fp64 fp64 fp64 fp64 df   df
  nsys_run_case "slau_weno${WO}" "warp_u64_o6"     slau_weno_warp 6 6 fp64 fp64 fp64 df   fp64 df
  nsys_run_case "slau_weno${WO}" "warp_p64_o6"     slau_weno_warp 6 6 fp64 fp64 fp64 df   df   fp64
done
unset BENCH_WENO_ORDER

nsys_finish
