#!/usr/bin/env bash
# SLAU-only Nsight Systems sweep for environments where ncu replay is unstable.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/_nsys_common.sh"

usage() {
  cat <<'EOF'
Usage:
  bash report/run_slau_nsys.sh --out report/nsys_runs/slau_a100 --gpu-cc 80 --repeat 3

SLAU sweep:
  MUSCL: ORDER=4,6 split/smem FP64 baselines and warp FP64/DF variants
  WENO : ORDER=6 split/smem FP64 baselines and warp FP64/DF variants
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

nsys_common_init slau

# SLAU + MUSCL reconstruction.
nsys_run_case slau_muscl split_fp64_o4       slau      4 4 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_muscl smem_fp64_o4        slau_smem 4 4 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_muscl warp_fp64_o4        slau_warp 4 4 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_muscl warp_all_df_o4      slau_warp 4 4 fp64 fp64 fp64 df   df   df
nsys_run_case slau_muscl warp_rho64_o4       slau_warp 4 4 fp64 fp64 fp64 fp64 df   df
nsys_run_case slau_muscl split_fp64_o6       slau      6 6 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_muscl smem_fp64_o6        slau_smem 6 6 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_muscl warp_fp64_o6        slau_warp 6 6 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_muscl warp_all_df_o6      slau_warp 6 6 fp64 fp64 fp64 df   df   df
nsys_run_case slau_muscl warp_rho64_o6       slau_warp 6 6 fp64 fp64 fp64 fp64 df   df
nsys_run_case slau_muscl warp_u64_o6         slau_warp 6 6 fp64 fp64 fp64 df   fp64 df
nsys_run_case slau_muscl warp_p64_o6         slau_warp 6 6 fp64 fp64 fp64 df   df   fp64

# SLAU + WENO5-Z reconstruction.
nsys_run_case slau_weno split_fp64_o6        slau_weno      6 6 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_weno smem_fp64_o6         slau_weno_smem 6 6 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_weno warp_fp64_o6         slau_weno_warp 6 6 fp64 fp64 fp64 fp64 fp64 fp64
nsys_run_case slau_weno warp_all_df_o6       slau_weno_warp 6 6 fp64 fp64 fp64 df   df   df
nsys_run_case slau_weno warp_rho64_o6        slau_weno_warp 6 6 fp64 fp64 fp64 fp64 df   df
nsys_run_case slau_weno warp_u64_o6          slau_weno_warp 6 6 fp64 fp64 fp64 df   fp64 df
nsys_run_case slau_weno warp_p64_o6          slau_weno_warp 6 6 fp64 fp64 fp64 df   df   fp64

nsys_finish
