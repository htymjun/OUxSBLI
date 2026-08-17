#!/usr/bin/env bash
# WENO-Z Nsight Compute sweep across the WENO_ORDER dimension.
#
# Companion to run_slau_ncu.sh, which covers MUSCL and WENO5-Z only. This one
# sweeps WENO_ORDER in {5, 7, 9} -- the knob added so the FP32/FP64 co-issue
# experiment has two halves large enough to overlap: at a flat 4 divisions per
# call the reconstruction carries 131 / 238 / 366 FP64-pipe SASS instructions
# per face (see report/rtx4060_weno_order_sweep.md).
#
# The case label carries the stencil width (e.g. slau_weno7), because bench.sh's
# CSV has no weno_order column -- WENO_ORDER is an env override
# (BENCH_WENO_ORDER), matching how BENCH_SCHEME/TVD/RECON already work.
#
# WENO9-Z now also has a double-float twin, so all three stencil widths run the
# same FP64/DF case set.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/_ncu_common.sh"

usage() {
  cat <<'EOF'
Usage:
  bash report/run_weno_ncu.sh --out report/ncu_runs/weno_a100 --gpu-cc 80 --repeat 3

Sweep (ORDER=6 throughout; ORDER selects the KEEP stencil and must stay 6 for
any WENO build, while WENO_ORDER selects the WENO-Z width):

  WENO_ORDER=5   split / smem / warp FP64 baselines, warp all-DF, warp per-variable
  WENO_ORDER=7   same
  WENO_ORDER=9   same
EOF
  ncu_common_usage_tail
}

ncu_common_defaults
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help)
      usage
      exit 0;;
    *)
      ncu_parse_common_arg "$@" || rc=$?
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

ncu_common_init weno

# ---------------------------------------------------------------------------
# FP64 baselines and the warp FP64/DF split, per stencil width.
# ---------------------------------------------------------------------------
for WO in 5 7 9; do
  export BENCH_WENO_ORDER="$WO"
  ncu_run_case "slau_weno${WO}" "split_fp64_o6"   slau_weno      6 6 fp64 fp64 fp64 fp64 fp64 fp64
  ncu_run_case "slau_weno${WO}" "smem_fp64_o6"    slau_weno_smem 6 6 fp64 fp64 fp64 fp64 fp64 fp64
  ncu_run_case "slau_weno${WO}" "warp_fp64_o6"    slau_weno_warp 6 6 fp64 fp64 fp64 fp64 fp64 fp64
  ncu_run_case "slau_weno${WO}" "warp_all_df_o6"  slau_weno_warp 6 6 fp64 fp64 fp64 df   df   df
  ncu_run_case "slau_weno${WO}" "warp_rho64_o6"   slau_weno_warp 6 6 fp64 fp64 fp64 fp64 df   df
  ncu_run_case "slau_weno${WO}" "warp_u64_o6"     slau_weno_warp 6 6 fp64 fp64 fp64 df   fp64 df
  ncu_run_case "slau_weno${WO}" "warp_p64_o6"     slau_weno_warp 6 6 fp64 fp64 fp64 df   df   fp64
done
unset BENCH_WENO_ORDER

ncu_finish
