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
#   report/bench.sh <mode> <ORDER> <VISC_ORDER> [KEEP_PREC] [VISC_PREC] [NX] [NT] [PRESS_PREC] [SLAU_RHO_PREC] [SLAU_U_PREC] [SLAU_P_PREC]
#
# Extra modes beyond the five KERNEL_MODEs: `split_visc` times the split
# path's second kernel (calc_Ev; total split time = split + split_visc rows),
# `quant` times calc_quantities_T_1D, and `fill` times
# calc_quantities_shadow32, the fused primitives+FP32-shadow pass (built as
# seq; works for any KEEP_PREC/VISC_PREC/PRESS_PREC now that seq accepts
# PRESS_PREC='fp32' directly). SLAU modes: `slau`/`slau_warp`
# use MUSCL reconstruction, `slau_weno`/`slau_weno_warp` use WENO5-Z
# reconstruction. The SLAU Riemann flux itself remains FP64.
#
# Emits one CSV row on stdout; run it in a loop and prepend the header from
# --header. Correctness is NOT checked here -- use the check_*.py scripts at
# nx=4096 for that.
# MPI is taken from PATH after `module load`: `MPIEXEC` defaults to `mpirun`.
# Override with e.g. `MPIEXEC=mpiexec` or `MPIEXEC_NP_FLAG=-np` if needed.
set -u
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MPIEXEC="${MPIEXEC:-${MPIRUN:-mpirun}}"
MPIEXEC_NP_FLAG="${MPIEXEC_NP_FLAG:--n}"
NCU="${NCU:-ncu}"
NCU_LAUNCH_SKIP="${NCU_LAUNCH_SKIP:-20}"
NCU_LAUNCH_COUNT="${NCU_LAUNCH_COUNT:-10}"
NCU_LAUNCH_MODE="${NCU_LAUNCH_MODE:-auto}"
NCU_RANKS="${NCU_RANKS:-1}"

if [ "${1:-}" = "--header" ]; then
  echo "mode,order,visc_order,keep_prec,visc_prec,press_prec,slau_rho_prec,slau_u_prec,slau_p_prec,nx,time_us,time_min_us,spread_us,fp64_pipe_inst,fp64_pipe_pct,dram_pct,sm_pct,occupancy_pct,waves_sm,regs,shared_b"
  exit 0
fi

MODE=$1; O=$2; VO=$3; KP=${4:-fp64}; VP=${5:-fp32}; NX=${6:-4194304}; NT=${7:-20}; PP=${8:-fp64}
SR=${9:-fp64}; SU=${10:-fp64}; SPREC=${11:-fp64}
REPORT_MODE="${BENCH_REPORT_MODE:-$MODE}"
case $MODE in
  seq) K=calc_seq_x;; fused) K=calc_fused_x;; warp) K=calc_warp_x;; warp_fused) K=calc_warp_fused_x;; split) K=calc_keep_x;;
  split_visc) K=calc_ev;; fill) K=calc_quantities_shadow32;; quant) K=calc_quantities_t_1d;;
  slau) K=calc_slau_x;;
  slau_warp) K=calc_slau_warp_x;;
  slau_weno) K=calc_slau_x;;
  slau_weno_warp) K=calc_slau_warp_x;;
  *) echo "bad mode: $MODE" >&2; exit 1;;
esac
K="${BENCH_KERNEL:-$K}"
# split_visc / fill are measurement aliases, not KERNEL_MODEs.
CMODE=$MODE
[ "$MODE" = split_visc ] && CMODE=split
[ "$MODE" = quant ] && CMODE=seq
# calc_quantities_shadow32 lives in calc_flux_base.f90.fypp, shared by all
# four non-split modes; seq now accepts PRESS_PREC='fp32' directly, so no
# mode substitution is needed there. With KP=fp64/VP=fp64/PP=fp64 (no shadow
# needed at all) this mode correctly reports NCU_FAIL -- there is no fused
# kernel to measure, which is the intended "killed the fill pass" outcome.
if [ "$MODE" = fill ]; then
  CMODE=seq
fi
SCHEME="${BENCH_SCHEME:-KEEP}"
TVD="${BENCH_TVD:-none}"
RECON="${BENCH_RECON:-MUSCL}"
if [ "$MODE" = slau ]; then
  SCHEME=SLAU
  TVD=tvd
  CMODE=split
fi
if [ "$MODE" = slau_warp ]; then
  SCHEME=SLAU
  TVD=tvd
  CMODE=warp
fi
if [ "$MODE" = slau_weno ]; then
  SCHEME=SLAU
  TVD=tvd
  RECON=WENO
  CMODE=split
fi
if [ "$MODE" = slau_weno_warp ]; then
  SCHEME=SLAU
  TVD=tvd
  RECON=WENO
  CMODE=warp
fi

python3 - "$CMODE" "$O" "$VO" "$KP" "$VP" "$NX" "$NT" "$R" "$PP" "$SCHEME" "$TVD" "$RECON" "$SR" "$SU" "$SPREC" <<'EOF'
import re, sys, pathlib
mode,o,vo,kp,vp,nx,nt,root,pp,scheme,tvd,recon,sr,su,sprec = sys.argv[1:16]
R = pathlib.Path(root)
p = R/'ST/config.fypp'; s = p.read_text()
s = re.sub(r"(?m)^#:set SCHEME\s*=.*$",      f"#:set SCHEME       = '{scheme}'", s)
s = re.sub(r"(?m)^#:set TVD\s*=.*$",         f"#:set TVD          = '{tvd}'", s)
s = re.sub(r"(?m)^#:set SLAU_RECON\s*=.*$",  f"#:set SLAU_RECON   = '{recon}'", s)
s = re.sub(r"(?m)^#:set ORDER\s*=.*$",       f"#:set ORDER        = {o}", s)
s = re.sub(r"(?m)^#:set VISC_ORDER\s*=.*$",  f"#:set VISC_ORDER   = {vo}", s)
s = re.sub(r"(?m)^#:set KEEP_PREC\s*=.*$",   f"#:set KEEP_PREC    = '{kp}'", s)
s = re.sub(r"(?m)^#:set VISC_PREC\s*=.*$",   f"#:set VISC_PREC    = '{vp}'", s)
s = re.sub(r"(?m)^#:set KERNEL_MODE\s*=.*$", f"#:set KERNEL_MODE  = '{mode}'", s)
s = re.sub(r"(?m)^#:set PRESS_PREC\s*=.*$",  f"#:set PRESS_PREC   = '{pp}'", s)
s = re.sub(r"(?m)^#:set SLAU_MUSCL_RHO_PREC\s*=.*$", f"#:set SLAU_MUSCL_RHO_PREC = '{sr}'", s)
s = re.sub(r"(?m)^#:set SLAU_MUSCL_U_PREC\s*=.*$",   f"#:set SLAU_MUSCL_U_PREC   = '{su}'", s)
s = re.sub(r"(?m)^#:set SLAU_MUSCL_P_PREC\s*=.*$",   f"#:set SLAU_MUSCL_P_PREC   = '{sprec}'", s)
p.write_text(s)
g = R/'ST/mod_globals.f90'; s = g.read_text()
s = re.sub(r"(?m)^(\s*integer, parameter :: nx = )\d+", rf"\g<1>{nx}", s)
s = re.sub(r"(?m)^(\s*integer, parameter :: nt = )\d+", rf"\g<1>{nt}", s)
g.write_text(s)
EOF

cd "$R/ST" || exit 1
if ! cmake --build build -j > /tmp/ouxsbli_bench_build.log 2>&1; then
  cmake --build build --target clean >> /tmp/ouxsbli_bench_build.log 2>&1
  cmake --build build -j >> /tmp/ouxsbli_bench_build.log 2>&1 || {
    echo "$REPORT_MODE,$O,$VO,$KP,$VP,$PP,$SR,$SU,$SPREC,$NX,BUILD_FAIL,,,,,,,,,," ; exit 1; }
fi

RES=$(cuobjdump -res-usage build/a.out 2>/dev/null | grep -A1 -i "${K}_" | tr '\n' ' ')
REG=$(echo "$RES" | grep -oP 'REG:\K[0-9]+'    | head -1)
SHM=$(echo "$RES" | grep -oP 'SHARED:\K[0-9]+' | head -1)

bench_launch_attempts() {
  case "$NCU_LAUNCH_MODE" in
    auto)
      if [ "$NCU_RANKS" = 1 ]; then
        printf '%s\n' direct mpi
      else
        printf '%s\n' mpi
      fi
      ;;
    direct)
      printf '%s\n' direct
      ;;
    mpi)
      printf '%s\n' mpi
      ;;
    *)
      echo "bad NCU_LAUNCH_MODE: $NCU_LAUNCH_MODE (use auto, direct, or mpi)" >&2
      exit 2
      ;;
  esac
}

bench_setup_launch() {
  local launch_mode="$1"
  if [ "$launch_mode" = direct ]; then
    BENCH_TARGET_PROCESSES=application-only
    BENCH_LAUNCH_DESC=direct
    BENCH_LAUNCH_CMD=(./a.out)
    return 0
  fi

  BENCH_TARGET_PROCESSES=all
  BENCH_LAUNCH_DESC="${MPIEXEC} ${MPIEXEC_NP_FLAG} ${NCU_RANKS} ./a.out"
  BENCH_LAUNCH_CMD=("$MPIEXEC" "$MPIEXEC_NP_FLAG" "$NCU_RANKS" ./a.out)
}

bench_line_has_valid_time() {
  local line="$1"
  local c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 time_us time_min_us spread_us rest

  [ -n "$line" ] || return 1
  case "$line" in
    *NCU_FAIL*|*RUN_FAIL*|*BUILD_FAIL*)
      return 1
      ;;
  esac

  IFS=, read -r c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 time_us time_min_us spread_us rest <<< "$line"
  [ -n "${time_us:-}" ] || return 1
  [ -n "${time_min_us:-}" ] || return 1
  case "$time_us,$time_min_us,$spread_us" in
    *nan*|*NaN*|*inf*|*Inf*)
      return 1
      ;;
  esac
  return 0
}

# --launch-skip 20 skips the first RK stages (warm-up); 10 profiled launches.
CSVF=$(mktemp)
RAWCSV="${NCU_RAW_CSV:-}"
if [ -z "$RAWCSV" ]; then
  RAWCSV=$(mktemp)
  CLEAN_RAWCSV=1
else
  mkdir -p "$(dirname "$RAWCSV")"
  CLEAN_RAWCSV=0
fi

NCU_ARGS=(--kernel-name "regex:$K"
      --launch-skip "$NCU_LAUNCH_SKIP" --launch-count "$NCU_LAUNCH_COUNT"
      --metrics gpu__time_duration.sum,launch__waves_per_multiprocessor,\
sm__inst_executed_pipe_fp64.sum,\
sm__pipe_fp64_cycles_active.avg.pct_of_peak_sustained_active,\
gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed,\
sm__throughput.avg.pct_of_peak_sustained_elapsed,\
sm__warps_active.avg.pct_of_peak_sustained_active \
      --csv --page raw)
if [ -n "${NCU_EXPORT:-}" ]; then
  mkdir -p "$(dirname "$NCU_EXPORT")"
  NCU_ARGS+=(--export "$NCU_EXPORT" --force-overwrite)
fi

LINE=""
for LAUNCH_MODE in $(bench_launch_attempts); do
  bench_setup_launch "$LAUNCH_MODE"
  : > "$RAWCSV"
  (cd build && "$NCU" \
    --target-processes "$BENCH_TARGET_PROCESSES" \
    "${NCU_ARGS[@]}" \
    "${BENCH_LAUNCH_CMD[@]}") > "$RAWCSV" 2>&1 || true

  # ncu writes its own ==PROF== chatter to stdout; strip it so the CSV parses.
  grep -v '^==' "$RAWCSV" > "$CSVF"
  LINE="$(python3 "$R/report/_bench_parse.py" "$REPORT_MODE" "$O" "$VO" "$KP" "$VP" "$PP" "$SR" "$SU" "$SPREC" "$NX" "${REG:-}" "${SHM:-}" "$CSVF")"
  if bench_line_has_valid_time "$LINE"; then
    break
  fi

  if [ "$NCU_LAUNCH_MODE" = auto ]; then
    echo "bench retry: launch_mode=$LAUNCH_MODE produced no finite kernel time for $REPORT_MODE/$O/$KP" >&2
  fi
done

printf '%s\n' "$LINE"
rm -f "$CSVF"
if [ "$CLEAN_RAWCSV" -eq 1 ]; then
  rm -f "$RAWCSV"
fi
