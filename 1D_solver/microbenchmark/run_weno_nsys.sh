#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLVER_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

NSYS="${NSYS:-nsys}"
FC="${FC:-nvfortran}"
GPU_CC="${CASE_GPU_CC:-native}"
NX="${NX:-4194304}"
NREPEAT="${NREPEAT:-1}"
NREPEAT_LIST="${NREPEAT_LIST:-1 2 4 8}"
REPEAT=1
OUT=""
RESUME=0
KEEP_GOING=0
NSYS_LAUNCH_SKIP="${NSYS_LAUNCH_SKIP:-0}"
NSYS_LAUNCH_COUNT="${NSYS_LAUNCH_COUNT:-0}"
NSYS_TRACE="${NSYS_TRACE:-cuda,nvtx,osrt}"

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
)

usage() {
  cat <<'EOF'
Usage:
  bash 1D_solver/microbenchmark/run_weno_nsys.sh --out 1D_solver/report/weno_micro_a100 --gpu-cc 80

Options:
  --out DIR             output directory
  --gpu-cc CC           CASE_GPU_CC for nvfortran, e.g. 80 or native
  --nx N                number of input cells, default 4194304
  --nrepeat N           inner kernel repeat count when --nrepeat-list is empty
  --nrepeat-list "..."  sweep inner repeat counts, default "1 2 4 8"
  --repeat N            process-level repeats, default 1
  --modes "A B ..."     override benchmark mode list. By default this script
                        runs the WENO weight/poly modes needed for the
                        FP32/FP64 overlap and shared-memory-cost study.
  --resume              skip completed rows in an existing --out
  --keep-going          continue after a failed mode
  --trace LIST          nsys trace list, default cuda,nvtx,osrt

The selected modes compare identical arithmetic in:
  weight_poly32_*       FP64 WENO-Z weights + FP32 candidate polynomials.

For each pair, *_seq executes both pieces sequentially in one thread and
*_warp lets the lower warp role compute the FP64 part while the upper warp role
computes the FP32 part for the same face.
weight_poly32_serial uses the same 256-thread/shared-memory layout as
weight_poly32_warp, but runs both FP64 and FP32 pieces on the lower warp role;
serial -> warp is the cleanest isolation of simultaneous FP64/FP32 execution.
weight_poly32_*_oncebar moves the barrier outside the inner repeat loop.
weight_poly32_wsmem_* stores only FP64 weights in shared memory, keeping FP32
polynomials in registers. weight_poly32_wsmem_tile2_* processes two faces per
thread role to amortize the single barrier. weight_poly32_halfwarp_* is a
shared-memory-free warp-shuffle diagnostic; because it splits lanes inside one
warp, it should not be interpreted as the main simultaneous-warp strategy.

The older var_* physical-variable split modes are still available via --modes,
but they are not part of the default run because the current microbenchmark
focus is the cleaner WENO weight/poly decomposition.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="${2:?missing --out value}"; shift 2;;
    --gpu-cc) GPU_CC="${2:?missing --gpu-cc value}"; shift 2;;
    --nx) NX="${2:?missing --nx value}"; shift 2;;
    --nrepeat) NREPEAT="${2:?missing --nrepeat value}"; shift 2;;
    --nrepeat-list) NREPEAT_LIST="${2:?missing --nrepeat-list value}"; shift 2;;
    --repeat) REPEAT="${2:?missing --repeat value}"; shift 2;;
    --modes) read -r -a MODES <<< "${2:?missing --modes value}"; shift 2;;
    --resume) RESUME=1; shift;;
    --keep-going) KEEP_GOING=1; shift;;
    --trace) NSYS_TRACE="${2:?missing --trace value}"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  OUT="$SOLVER_ROOT/report/weno_micro_${stamp}"
fi
mkdir -p "$OUT/raw" "$OUT/reports" "$OUT/stats" "$OUT/sqlite"
OUT="$(cd "$OUT" && pwd)"
if [ -z "${NREPEAT_LIST// }" ]; then
  NREPEAT_LIST="$NREPEAT"
fi

BUILD_DIR="$SCRIPT_DIR/build"
EXE="$BUILD_DIR/weno_micro"
if [ "$RESUME" -eq 0 ]; then
  : > "$OUT/progress.log"
  printf 'mode,repeat,nx,nrepeat,kernel,total_instances,used_instances,time_us,time_min_us,time_max_us,spread_us,total_us,trace_csv,kern_sum_csv,report,sqlite,status\n' > "$OUT/summary.csv"
  printf 'mode\trepeat\tnx\tnrepeat\tcommand\n' > "$OUT/commands.tsv"
fi

{
  echo "date=$(date -Is)"
  echo "out=$OUT"
  echo "gpu_cc=$GPU_CC"
  echo "nx=$NX"
  echo "nrepeat=$NREPEAT"
  echo "nrepeat_list=$NREPEAT_LIST"
  echo "repeat=$REPEAT"
  echo "nsys_trace=$NSYS_TRACE"
  echo
  echo "[gpu]"
  nvidia-smi --query-gpu=name,driver_version,compute_cap,pci.bus_id --format=csv 2>/dev/null || true
  echo
  echo "[nsys]"
  "$NSYS" --version 2>&1 || true
  echo
  echo "[compiler]"
  "$FC" --version 2>&1 | head -20 || true
  echo
  echo "[git]"
  git -C "$SOLVER_ROOT/.." rev-parse --short HEAD 2>/dev/null || true
  git -C "$SOLVER_ROOT/.." status --short 2>/dev/null || true
} > "$OUT/env.txt"

cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_Fortran_COMPILER="$FC" -DCASE_GPU_CC="$GPU_CC" > "$OUT/cmake_configure.log" 2>&1
cmake --build "$BUILD_DIR" -j > "$OUT/build.log" 2>&1

complete_row() {
  local mode="$1" rep="$2" nrep="$3"
  [ "$RESUME" -eq 1 ] || return 1
  awk -F, -v m="$mode" -v r="$rep" -v nr="$nrep" \
    'NR > 1 && $1 == m && $2 == r && $4 == nr && $NF !~ /FAIL/ { found = 1 } END { exit found ? 0 : 1 }' "$OUT/summary.csv"
}

for nrep in $NREPEAT_LIST; do
  for mode in "${MODES[@]}"; do
    for rep in $(seq 1 "$REPEAT"); do
      if complete_row "$mode" "$rep" "$nrep"; then
        echo "skip completed $mode nrepeat $nrep repeat $rep" | tee -a "$OUT/progress.log"
        continue
      fi
      base="${mode}__nr${nrep}__r${rep}"
      report_base="$OUT/reports/$base"
      report="$report_base.nsys-rep"
      sqlite="$OUT/sqlite/$base.sqlite"
      trace_csv="$OUT/stats/${base}__trace_cuda_gpu_trace_base.csv"
      kern_csv="$OUT/stats/${base}__kern_sum_cuda_gpu_kern_sum_base.csv"
      raw_log="$OUT/raw/$base.log"
      echo "[$(date +%H:%M:%S)] nsys $mode nrepeat $nrep repeat $rep/$REPEAT" | tee -a "$OUT/progress.log"
      printf '%s\t%s\t%s\t%s\t%s %s %s %s\n' "$mode" "$rep" "$NX" "$nrep" "$EXE" "$mode" "$NX" "$nrep" >> "$OUT/commands.tsv"
      set +e
      "$NSYS" profile --trace="$NSYS_TRACE" --sample=none --cpuctxsw=none --backtrace=none \
        --force-overwrite=true --output="$report_base" "$EXE" "$mode" "$NX" "$nrep" > "$raw_log" 2>&1
      rc=$?
      set -e
      if [ "$rc" -ne 0 ]; then
        printf '%s,%s,%s,%s,weno_%s,0,0,NSYS_FAIL,,,,,,,,PROFILE_FAIL\n' "$mode" "$rep" "$NX" "$nrep" "$mode" >> "$OUT/summary.csv"
        [ "$KEEP_GOING" -eq 1 ] || exit "$rc"
        continue
      fi
      "$NSYS" export --type sqlite --force-overwrite=true --output "$sqlite" "$report" > "$OUT/raw/${base}__export_sqlite.log" 2>&1 || true
      "$NSYS" stats --report cuda_gpu_trace:base --format csv --output "$OUT/stats/${base}__trace" "$report" > "$OUT/raw/${base}__stats_trace.log" 2>&1 || true
      "$NSYS" stats --report cuda_gpu_kern_sum:base --format csv --output "$OUT/stats/${base}__kern_sum" "$report" > "$OUT/raw/${base}__stats_kern_sum.log" 2>&1 || true
      python3 "$SCRIPT_DIR/parse_weno_nsys.py" \
        --mode "$mode" --repeat "$rep" --nx "$NX" --nrepeat "$nrep" \
        --kernel "weno_.*${mode%%_*}" \
        --trace-csv "$trace_csv" --kern-sum-csv "$kern_csv" \
        --report "$report" --sqlite "$sqlite" \
        --skip "$NSYS_LAUNCH_SKIP" --count "$NSYS_LAUNCH_COUNT" >> "$OUT/summary.csv"
    done
  done
done

echo "summary: $OUT/summary.csv"
