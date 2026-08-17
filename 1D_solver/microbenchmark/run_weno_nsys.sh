#!/usr/bin/env bash
# WENO microbenchmark driver.
#
# Primary measurement is weno_micro.f90's OWN CUDA-event timing: the binary
# takes `nlaunch` timed launches and prints time_min_us / time_avg_us. No
# profiler is required, which matters because ncu and nsys have both hung
# repeatedly on this hardware while the unprofiled binary ran fine
# (report/rtx4060_weno_division_reduction.md). Pass --nsys to additionally
# capture a timeline. By default this script also runs the solver-level
# WENO_ORDER={5,7,9} Nsight Systems sweep in a sibling output directory, so one
# command gives both the WENO5-Z isolated microbenchmark and the WENO7/9 solver
# data.
#
# Two harness rules are enforced here rather than left to the caller:
#
#   * nrepeat is pinned to 1. At nrepeat>1 the comparison silently breaks --
#     nvfortran hoists the inner loop of any mode that writes only private
#     registers (weight_poly32_seq, halfwarp_serial: ~2-4 us per extra repeat)
#     while modes that touch shared memory or a barrier keep paying (~700-1400
#     us). The old --nrepeat-list sweep produced rows that cannot be compared
#     and has been removed.
#   * Rounds are interleaved (round-major, not mode-major) so GPU clock drift
#     hits every mode equally. Clocks cannot be locked without root on some of
#     these boxes; compare time_min_us across rounds, never a single run.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLVER_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

NSYS="${NSYS:-nsys}"
FC="${FC:-nvfortran}"
GPU_CC="${CASE_GPU_CC:-native}"
NX="${NX:-4194304}"
NLAUNCH="${NLAUNCH:-10}"
REPEAT=3
OUT=""
RESUME=0
KEEP_GOING=0
USE_NSYS=0
NSYS_TRACE="${NSYS_TRACE:-cuda,nvtx,osrt}"
RUN_SOLVER_WENO_ORDERS=1
SOLVER_WENO_NT="${SOLVER_WENO_NT:-20}"
SOLVER_WENO_REPEAT="${SOLVER_WENO_REPEAT:-}"

# Both families, because the point of the study is which half goes to FP32:
#   var_*           demote the WENO-Z WEIGHTS (and thus every division)
#   weight_poly32_* demote the candidate POLYNOMIALS
# The 7/9 suffix selects the WENO-Z stencil width (WENO7-Z / WENO9-Z); unsuffixed
# is WENO5-Z. Only the five core co-issue modes have wider-stencil versions --
# the shared-memory/barrier layout variants stay WENO5-only.
# The var_* modes were absent from the old default list because they returned
# NaN; that was an nvfortran `value` miscompile, since fixed, and they turned
# out to be the interesting ones (weights are ~96% of the FP64-pipe work).
MODES=(
  var_fp64_seq
  var_fp64_warp
  var_fp32_seq
  var_fp32_warp
  var_rho64_seq
  var_u64_seq
  var_p64_seq
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
  var_fp64_seq7
  var_fp32_seq7
  weight_poly32_seq7
  weight_poly32_serial7
  weight_poly32_warp7
  var_fp64_seq9
  var_fp32_seq9
  weight_poly32_seq9
  weight_poly32_serial9
  weight_poly32_warp9
)

usage() {
  cat <<'EOF'
Usage:
  bash 1D_solver/microbenchmark/run_weno_nsys.sh --out 1D_solver/report/weno_micro_a100 --gpu-cc 80

Options:
  --out DIR             output directory (default report/weno_micro_<stamp>)
  --gpu-cc CC           CASE_GPU_CC for nvfortran, e.g. 80 or native. Set this
                        explicitly: a stale cached value silently builds for the
                        wrong arch and every launch fails with
                        cudaErrorInvalidPtx (218).
  --nx N                number of input cells, default 4194304
  --nlaunch N           timed launches per run, min is reported, default 10
  --repeat N            interleaved rounds, default 3
  --modes "A B ..."     override the mode list
  --nsys                also capture an nsys timeline per row (slow, optional)
  --trace LIST          nsys trace list, default cuda,nvtx,osrt
  --no-solver-orders    do not run the solver-level WENO_ORDER=5/7/9 sweep
  --solver-nt N         time steps for the solver-level sweep, default 20
  --resume              skip completed rows in an existing --out
  --keep-going          continue after a failed mode

Outputs:
  summary.csv                  WENO5-Z isolated microbenchmark
  solver_weno_orders/summary.csv
                               SLAU+WENO5/7/9 solver-level sweep, produced by
                               1D_solver/report/run_weno_nsys.sh

Modes:
  var_*            split by physical VARIABLE -- which of rho/u/p keep FP64
                   weights. var_fp32_* puts all three sets of WENO-Z weights
                   (and every division) in FP32; var_rho64/u64/p64 keep one.
  weight_poly32_*  split by ALGEBRAIC PIECE -- FP64 weights + FP32 candidate
                   polynomials. *_seq is one thread doing both; *_warp gives the
                   FP64 part to the lower warp role and the FP32 part to the
                   upper role for the same face; *_serial is the identical
                   256-thread/shared-memory/barrier layout with both pieces on
                   the lower role, so serial -> warp is the clean isolation of
                   simultaneous FP64/FP32 execution.
                   *_oncebar hoists the barrier; *_wsmem_* keeps FP32
                   polynomials in registers and only weights in shared memory;
                   *_wsmem_tile2_* does two faces per role to amortise the
                   barrier; *_halfwarp_* is a shared-memory-free shuffle
                   diagnostic and splits lanes WITHIN a warp, so it is not the
                   main simultaneous-warp strategy.

  The binary also accepts var_rho64_warp, var_u64_warp and var_p64_warp; they
  are left out of the default list because var_fp64_warp and var_fp32_warp
  already bracket the seq-vs-warp comparison. Add them with --modes if you want
  the per-variable warp split too.

Correctness: every row records checksum_all and nonfinite. A nonzero nonfinite,
or a checksum that moves between modes by more than rounding, invalidates the
timing -- weno_micro error-stops on non-finite output rather than reporting it.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="${2:?missing --out value}"; shift 2;;
    --gpu-cc) GPU_CC="${2:?missing --gpu-cc value}"; shift 2;;
    --nx) NX="${2:?missing --nx value}"; shift 2;;
    --nlaunch) NLAUNCH="${2:?missing --nlaunch value}"; shift 2;;
    --repeat) REPEAT="${2:?missing --repeat value}"; shift 2;;
    --modes) read -r -a MODES <<< "${2:?missing --modes value}"; shift 2;;
    --nsys) USE_NSYS=1; shift;;
    --trace) NSYS_TRACE="${2:?missing --trace value}"; shift 2;;
    --no-solver-orders) RUN_SOLVER_WENO_ORDERS=0; shift;;
    --solver-nt) SOLVER_WENO_NT="${2:?missing --solver-nt value}"; shift 2;;
    --resume) RESUME=1; shift;;
    --keep-going) KEEP_GOING=1; shift;;
    --nrepeat|--nrepeat-list)
      echo "error: $1 was removed. nrepeat is pinned to 1 -- at nrepeat>1" >&2
      echo "       nvfortran hoists the inner loop for register-only modes but" >&2
      echo "       not for shared-memory ones, so the rows are not comparable." >&2
      exit 2;;
    -h|--help) usage; exit 0;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done

NREPEAT=1   # see the header; not configurable on purpose
[ -n "$SOLVER_WENO_REPEAT" ] || SOLVER_WENO_REPEAT="$REPEAT"

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  OUT="$SOLVER_ROOT/report/weno_micro_${stamp}"
fi
mkdir -p "$OUT/raw"
[ "$USE_NSYS" -eq 1 ] && mkdir -p "$OUT/reports" "$OUT/stats" "$OUT/sqlite"
OUT="$(cd "$OUT" && pwd)"

BUILD_DIR="$SCRIPT_DIR/build"
EXE="$BUILD_DIR/weno_micro"
if [ "$RESUME" -eq 0 ]; then
  : > "$OUT/progress.log"
  printf 'mode,round,nx,nrepeat,nlaunch,time_min_us,time_avg_us,checksum_all,nonfinite,status\n' > "$OUT/summary.csv"
  printf 'mode\tround\tcommand\n' > "$OUT/commands.tsv"
fi

{
  echo "date=$(date -Is)"
  echo "out=$OUT"
  echo "gpu_cc=$GPU_CC"
  echo "nx=$NX"
  echo "nrepeat=$NREPEAT (pinned)"
  echo "nlaunch=$NLAUNCH"
  echo "repeat=$REPEAT"
  echo "use_nsys=$USE_NSYS"
  echo "run_solver_weno_orders=$RUN_SOLVER_WENO_ORDERS"
  echo "solver_weno_nt=$SOLVER_WENO_NT"
  echo "solver_weno_repeat=$SOLVER_WENO_REPEAT"
  echo
  echo "[gpu]"
  nvidia-smi --query-gpu=name,driver_version,compute_cap,pci.bus_id --format=csv 2>/dev/null || true
  echo
  echo "[compiler]"
  "$FC" --version 2>&1 | head -20 || true
  [ "$USE_NSYS" -eq 1 ] && { echo; echo "[nsys]"; "$NSYS" --version 2>&1 || true; }
  echo
  echo "[git]"
  git -C "$SOLVER_ROOT/.." rev-parse --short HEAD 2>/dev/null || true
  git -C "$SOLVER_ROOT/.." status --short 2>/dev/null || true
} > "$OUT/env.txt"

cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_Fortran_COMPILER="$FC" -DCASE_GPU_CC="$GPU_CC" > "$OUT/cmake_configure.log" 2>&1
cmake --build "$BUILD_DIR" -j > "$OUT/build.log" 2>&1
# The target architecture is echoed by microbenchmark/CMakeLists.txt; surface it
# so a wrong-arch build is visible here and not only as a launch failure.
grep -i "target architecture" "$OUT/cmake_configure.log" | tee -a "$OUT/progress.log" || true

complete_row() {
  local mode="$1" rnd="$2"
  [ "$RESUME" -eq 1 ] || return 1
  awk -F, -v m="$mode" -v r="$rnd" \
    'NR > 1 && $1 == m && $2 == r && $NF == "OK" { found = 1 } END { exit found ? 0 : 1 }' "$OUT/summary.csv"
}

field() { grep -oP "$1=\\s*\\K[-0-9.eE+]+" "$2" | tail -1; }

# Round-major so clock drift affects every mode equally.
for rnd in $(seq 1 "$REPEAT"); do
  for mode in "${MODES[@]}"; do
    if complete_row "$mode" "$rnd"; then
      echo "skip completed $mode round $rnd" | tee -a "$OUT/progress.log"
      continue
    fi
    base="${mode}__r${rnd}"
    raw_log="$OUT/raw/$base.log"
    echo "[$(date +%H:%M:%S)] $mode round $rnd/$REPEAT" | tee -a "$OUT/progress.log"
    printf '%s\t%s\t%s %s %s %s %s\n' "$mode" "$rnd" "$EXE" "$mode" "$NX" "$NREPEAT" "$NLAUNCH" >> "$OUT/commands.tsv"

    set +e
    "$EXE" "$mode" "$NX" "$NREPEAT" "$NLAUNCH" > "$raw_log" 2>&1
    rc=$?
    set -e
    if [ "$rc" -ne 0 ]; then
      printf '%s,%s,%s,%s,%s,RUN_FAIL,,,,FAIL\n' "$mode" "$rnd" "$NX" "$NREPEAT" "$NLAUNCH" >> "$OUT/summary.csv"
      echo "  FAILED (rc=$rc), see $raw_log" | tee -a "$OUT/progress.log"
      [ "$KEEP_GOING" -eq 1 ] || exit "$rc"
      continue
    fi
    tmin="$(field time_min_us "$raw_log")"
    tavg="$(field time_avg_us "$raw_log")"
    csum="$(field checksum_all "$raw_log")"
    nbad="$(field nonfinite "$raw_log")"
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
      "$mode" "$rnd" "$NX" "$NREPEAT" "$NLAUNCH" \
      "${tmin:-NA}" "${tavg:-NA}" "${csum:-NA}" "${nbad:-NA}" \
      "$([ "${nbad:-1}" = "0" ] && echo OK || echo NONFINITE)" >> "$OUT/summary.csv"
    echo "  time_min_us=${tmin:-NA}  nonfinite=${nbad:-NA}" | tee -a "$OUT/progress.log"

    if [ "$USE_NSYS" -eq 1 ]; then
      report_base="$OUT/reports/$base"
      set +e
      "$NSYS" profile --trace="$NSYS_TRACE" --sample=none --cpuctxsw=none --backtrace=none \
        --force-overwrite=true --output="$report_base" \
        "$EXE" "$mode" "$NX" "$NREPEAT" "$NLAUNCH" > "$OUT/raw/${base}__nsys.log" 2>&1
      nrc=$?
      set -e
      if [ "$nrc" -eq 0 ]; then
        "$NSYS" export --type sqlite --force-overwrite=true --output "$OUT/sqlite/$base.sqlite" \
          "$report_base.nsys-rep" > "$OUT/raw/${base}__export.log" 2>&1 || true
        "$NSYS" stats --report cuda_gpu_kern_sum:base --format csv \
          --output "$OUT/stats/${base}__kern_sum" "$report_base.nsys-rep" \
          > "$OUT/raw/${base}__stats.log" 2>&1 || true
      else
        echo "  nsys failed (rc=$nrc); CUDA-event timing above is unaffected" | tee -a "$OUT/progress.log"
      fi
    fi
  done
done

if [ "$RUN_SOLVER_WENO_ORDERS" -eq 1 ]; then
  solver_out="$OUT/solver_weno_orders"
  echo
  echo "[$(date +%H:%M:%S)] solver-level WENO_ORDER=5/7/9 sweep -> $solver_out" | tee -a "$OUT/progress.log"
  solver_args=(
    --out "$solver_out"
    --gpu-cc "$GPU_CC"
    --nx "$NX"
    --nt "$SOLVER_WENO_NT"
    --repeat "$SOLVER_WENO_REPEAT"
  )
  [ "$RESUME" -eq 1 ] && solver_args+=(--resume)
  [ "$KEEP_GOING" -eq 1 ] && solver_args+=(--keep-going)
  NSYS="$NSYS" NSYS_TRACE="$NSYS_TRACE" \
    bash "$SOLVER_ROOT/report/run_weno_nsys.sh" "${solver_args[@]}" \
    2>&1 | tee "$OUT/solver_weno_orders.log"
fi

echo
echo "summary: $OUT/summary.csv"
python3 - "$OUT/summary.csv" <<'EOF'
import csv, sys, collections, statistics as st
rows = [r for r in csv.DictReader(open(sys.argv[1])) if r["status"] == "OK"]
if not rows:
    sys.exit(0)
d = collections.defaultdict(list)
for r in rows:
    try:
        d[r["mode"]].append(float(r["time_min_us"]))
    except ValueError:
        pass
if not d:
    sys.exit(0)
best = min(min(v) for v in d.values())
print(f"\n{'mode':34s} {'min us':>11} {'rounds':>7} {'vs fastest':>11}")
for m, v in sorted(d.items(), key=lambda kv: min(kv[1])):
    print(f"{m:34s} {min(v):11.1f} {len(v):7d} {min(v)/best:10.2f}x")
bad = [r for r in csv.DictReader(open(sys.argv[1])) if r["status"] != "OK"]
if bad:
    print(f"\n{len(bad)} row(s) NOT ok: " + ", ".join(sorted({r['mode'] for r in bad})))
EOF

if [ -s "$OUT/solver_weno_orders/summary.csv" ]; then
  echo
  echo "solver WENO_ORDER summary: $OUT/solver_weno_orders/summary.csv"
  python3 - "$OUT/solver_weno_orders/summary.csv" <<'EOF'
import csv, math, sys
rows = list(csv.DictReader(open(sys.argv[1])))
groups = {}
for r in rows:
    try:
        t = float(r.get("time_us", "nan"))
    except ValueError:
        continue
    if not math.isfinite(t):
        continue
    groups.setdefault((r.get("table_id",""), r.get("case_id","")), []).append(t)
if not groups:
    sys.exit(0)
print(f"\n{'table':12s} {'case':20s} {'min us':>11} {'runs':>5}")
for (table, case), vals in sorted(groups.items()):
    print(f"{table:12s} {case:20s} {min(vals):11.1f} {len(vals):5d}")
EOF
fi
