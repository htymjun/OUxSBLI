#!/usr/bin/env bash
# Per-kernel benchmark for the 1D flux variants.
#
# Why this exists: the earlier rtx4060_*.md sweeps all ran at the stock
# nx=4096, where every configuration sits 4-35x above its own compute floor and
# the grid is only 1.3-2.7 waves/SM. At that size the numbers are dominated by
# launch overhead and by block-count rounding across the SMs, not by the kernel.
# Any performance claim needs a grid large enough that those vanish.
#
# Default nx=4194304 gives ~1365 waves/SM and ~470 MB of device arrays
# (112 B/cell), comfortable on an 8 GB card.
#
#   report/bench.sh <mode> <ORDER> <VISC_ORDER> [KEEP_PREC] [VISC_PREC] [NX] [NT] [PRESS_PREC]
#
# Extra modes beyond the five KERNEL_MODEs: `split_visc` times the split
# path's second kernel (calc_Ev; total split time = split + split_visc rows),
# and `fill` times calc_quantities_shadow32, the fused primitives+FP32-shadow
# pass (built as seq; works for any KEEP_PREC/VISC_PREC/PRESS_PREC now that
# seq accepts PRESS_PREC='fp32' directly).
#
# Emits one CSV row on stdout; run it in a loop and prepend the header from
# --header. Correctness is NOT checked here -- use the check_*.py scripts at
# nx=4096 for that.
set -u
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MPIBIN=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/mpi/bin
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/12.5/hpcx/hpcx-2.19/ompi

if [ "${1:-}" = "--header" ]; then
  echo "mode,order,visc_order,keep_prec,visc_prec,press_prec,nx,time_us,time_min_us,spread_us,fp64_pipe_inst,fp64_pipe_pct,dram_pct,sm_pct,occupancy_pct,waves_sm,regs,shared_b"
  exit 0
fi

MODE=$1; O=$2; VO=$3; KP=${4:-fp64}; VP=${5:-fp32}; NX=${6:-4194304}; NT=${7:-20}; PP=${8:-fp64}
case $MODE in
  seq) K=calc_seq_x;; fused) K=calc_fused_x;; warp) K=calc_warp_x;; warp_fused) K=calc_warp_fused_x;; split) K=calc_keep_x;;
  split_visc) K=calc_ev;; fill) K=calc_quantities_shadow32;;
  *) echo "bad mode: $MODE" >&2; exit 1;;
esac
# split_visc / fill are measurement aliases, not KERNEL_MODEs.
CMODE=$MODE
[ "$MODE" = split_visc ] && CMODE=split
# calc_quantities_shadow32 lives in calc_flux_base.f90.fypp, shared by all
# four non-split modes; seq now accepts PRESS_PREC='fp32' directly, so no
# mode substitution is needed there. With KP=fp64/VP=fp64/PP=fp64 (no shadow
# needed at all) this mode correctly reports NCU_FAIL -- there is no fused
# kernel to measure, which is the intended "killed the fill pass" outcome.
if [ "$MODE" = fill ]; then
  CMODE=seq
fi

python3 - "$CMODE" "$O" "$VO" "$KP" "$VP" "$NX" "$NT" "$R" "$PP" <<'EOF'
import re, sys, pathlib
mode,o,vo,kp,vp,nx,nt,root,pp = sys.argv[1:10]
R = pathlib.Path(root)
p = R/'ST/config.fypp'; s = p.read_text()
s = re.sub(r"(?m)^#:set ORDER\s*=.*$",       f"#:set ORDER        = {o}", s)
s = re.sub(r"(?m)^#:set VISC_ORDER\s*=.*$",  f"#:set VISC_ORDER   = {vo}", s)
s = re.sub(r"(?m)^#:set KEEP_PREC\s*=.*$",   f"#:set KEEP_PREC    = '{kp}'", s)
s = re.sub(r"(?m)^#:set VISC_PREC\s*=.*$",   f"#:set VISC_PREC    = '{vp}'", s)
s = re.sub(r"(?m)^#:set KERNEL_MODE\s*=.*$", f"#:set KERNEL_MODE  = '{mode}'", s)
s = re.sub(r"(?m)^#:set PRESS_PREC\s*=.*$",  f"#:set PRESS_PREC   = '{pp}'", s)
p.write_text(s)
g = R/'ST/mod_globals.f90'; s = g.read_text()
s = re.sub(r"(?m)^(\s*integer, parameter :: nx = )\d+", rf"\g<1>{nx}", s)
s = re.sub(r"(?m)^(\s*integer, parameter :: nt = )\d+", rf"\g<1>{nt}", s)
g.write_text(s)
EOF

cd "$R/ST" || exit 1
cmake --build build -j > /tmp/ouxsbli_bench_build.log 2>&1 || {
  echo "$MODE,$O,$VO,$KP,$VP,$PP,$NX,BUILD_FAIL,,,,,,,,,," ; exit 1; }

RES=$(cuobjdump -res-usage build/a.out 2>/dev/null | grep -A1 -i "${K}_" | tr '\n' ' ')
REG=$(echo "$RES" | grep -oP 'REG:\K[0-9]+'    | head -1)
SHM=$(echo "$RES" | grep -oP 'SHARED:\K[0-9]+' | head -1)

# --launch-skip 20 skips the first RK stages (warm-up); 10 profiled launches.
CSVF=$(mktemp)
# ncu writes its own ==PROF== chatter to stdout; strip it so the CSV parses.
(cd build && ncu --target-processes all --kernel-name "regex:$K" \
      --launch-skip 20 --launch-count 10 \
      --metrics gpu__time_duration.sum,launch__waves_per_multiprocessor,\
sm__inst_executed_pipe_fp64.sum,\
sm__pipe_fp64_cycles_active.avg.pct_of_peak_sustained_active,\
gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed,\
sm__throughput.avg.pct_of_peak_sustained_elapsed,\
sm__warps_active.avg.pct_of_peak_sustained_active \
      --csv --page raw $MPIBIN/mpirun -n 1 ./a.out 2>/dev/null) \
  | grep -v '^==' > "$CSVF"

python3 "$R/report/_bench_parse.py" "$MODE" "$O" "$VO" "$KP" "$VP" "$PP" "$NX" "${REG:-}" "${SHM:-}" "$CSVF"
rm -f "$CSVF"
