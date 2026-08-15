#!/usr/bin/env bash
# Collect every Nsight Compute row needed to refresh
# paper/fp32fp64_1d_summary.tex on a new GPU.
#
# Typical use:
#
#   # A100
#   bash report/run_summary_ncu.sh --out report/ncu_runs/a100 --gpu-cc 80 --repeat 3
#
#   # RTX 6000 Blackwell / other Blackwell systems
#   bash report/run_summary_ncu.sh --out report/ncu_runs/rtx6000_blackwell --gpu-cc native --repeat 3
#
# Full Nsight Compute reports are collected once per case by default
# (`--full-repeat first`). Use `--full-repeat all` to collect full reports for
# every repeat, or `--no-full-ncu` to collect only timing CSV/reports.
#
# MPI is taken from PATH after `module load`: `MPIEXEC` defaults to `mpirun`.
# Override with e.g. `MPIEXEC=mpiexec` or `MPIEXEC_NP_FLAG=-np` if needed.
#
# The output directory is self-contained:
#   summary.csv       all parsed rows, with table_id/case_id/repeat columns
#   commands.tsv      exact bench.sh argument list for each row
#   env.txt           GPU/toolchain metadata
#   raw/*.log         stdout/stderr capture from bench.sh
#   csv/*.csv         raw ncu CSV from the timing pass
#   timing_reports/*  .ncu-rep reports for the lightweight timing metrics
#   full_reports/*    .ncu-rep reports from a separate --set full pass
#   full_logs/*.log   console output from the --set full pass
set -euo pipefail

SOLVER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$SOLVER_ROOT/.." && pwd)"
BENCH="$SOLVER_ROOT/report/bench.sh"
MPIEXEC="${MPIEXEC:-${MPIRUN:-mpirun}}"
MPIEXEC_NP_FLAG="${MPIEXEC_NP_FLAG:--n}"

OUT=""
GPU_CC="${CASE_GPU_CC:-native}"
NX="${NX:-4194304}"
NT="${NT:-20}"
REPEAT=1
KEEP_GOING=0
COLLECT_FULL=1
FULL_REPEAT="${FULL_REPEAT:-first}"
FULL_NCU_SET="${FULL_NCU_SET:-full}"
FULL_NCU_LAUNCH_SKIP="${FULL_NCU_LAUNCH_SKIP:-20}"
FULL_NCU_LAUNCH_COUNT="${FULL_NCU_LAUNCH_COUNT:-1}"

usage() {
  sed -n '2,/^set -euo pipefail/p' "$0" | sed '$d;s/^# \{0,1\}//'
}

while [ $# -gt 0 ]; do
  case "$1" in
    --out)
      OUT="${2:?missing --out value}"; shift 2;;
    --gpu-cc)
      GPU_CC="${2:?missing --gpu-cc value}"; shift 2;;
    --nx)
      NX="${2:?missing --nx value}"; shift 2;;
    --nt)
      NT="${2:?missing --nt value}"; shift 2;;
    --repeat)
      REPEAT="${2:?missing --repeat value}"; shift 2;;
    --keep-going)
      KEEP_GOING=1; shift;;
    --no-full-ncu)
      COLLECT_FULL=0; shift;;
    --full-repeat)
      FULL_REPEAT="${2:?missing --full-repeat value}"; shift 2;;
    --full-set)
      FULL_NCU_SET="${2:?missing --full-set value}"; shift 2;;
    --full-launch-skip)
      FULL_NCU_LAUNCH_SKIP="${2:?missing --full-launch-skip value}"; shift 2;;
    --full-launch-count)
      FULL_NCU_LAUNCH_COUNT="${2:?missing --full-launch-count value}"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 2;;
  esac
done

if [ -z "$OUT" ]; then
  stamp="$(date +%Y%m%d_%H%M%S)"
  gpu_slug="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 | tr '[:upper:] /' '[:lower:]__' | tr -cd 'a-z0-9_-' || true)"
  OUT="$SOLVER_ROOT/report/ncu_runs/${gpu_slug:-gpu}_${stamp}"
fi

mkdir -p "$OUT/raw" "$OUT/csv" "$OUT/timing_reports" "$OUT/full_reports" "$OUT/full_logs"
OUT="$(cd "$OUT" && pwd)"

CFG="$SOLVER_ROOT/ST/config.fypp"
GLOB="$SOLVER_ROOT/ST/mod_globals.f90"
CFG_BAK="$(mktemp)"
GLOB_BAK="$(mktemp)"
cp "$CFG" "$CFG_BAK"
cp "$GLOB" "$GLOB_BAK"
restore_inputs() {
  cp "$CFG_BAK" "$CFG"
  cp "$GLOB_BAK" "$GLOB"
  rm -f "$CFG_BAK" "$GLOB_BAK"
}
trap restore_inputs EXIT

{
  echo "date=$(date -Is)"
  echo "repo=$REPO_ROOT"
  echo "solver_root=$SOLVER_ROOT"
  echo "out=$OUT"
  echo "gpu_cc=$GPU_CC"
  echo "nx=$NX"
  echo "nt=$NT"
  echo "repeat=$REPEAT"
  echo "collect_full=$COLLECT_FULL"
  echo "full_repeat=$FULL_REPEAT"
  echo "full_ncu_set=$FULL_NCU_SET"
  echo "full_ncu_launch_skip=$FULL_NCU_LAUNCH_SKIP"
  echo "full_ncu_launch_count=$FULL_NCU_LAUNCH_COUNT"
  echo
  echo "[gpu]"
  nvidia-smi --query-gpu=name,driver_version,compute_cap,pci.bus_id --format=csv 2>/dev/null || true
  echo
  echo "[ncu]"
  "${NCU:-ncu}" --version 2>&1 || true
  echo
  echo "[compiler]"
  "${FC:-mpif90}" --version 2>&1 | head -20 || true
  echo
  echo "[mpi]"
  echo "MPIEXEC=$MPIEXEC"
  echo "MPIEXEC_NP_FLAG=$MPIEXEC_NP_FLAG"
  command -v "$MPIEXEC" 2>/dev/null || true
  "$MPIEXEC" --version 2>&1 | head -20 || true
  echo
  echo "[git]"
  git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || true
  git -C "$REPO_ROOT" status --short 2>/dev/null || true
} > "$OUT/env.txt"

cmake_args=(-S "$SOLVER_ROOT/ST" -B "$SOLVER_ROOT/ST/build" -DCASE_GPU_CC="$GPU_CC")
if [ -n "${FC:-}" ]; then
  cmake_args+=("-DCMAKE_Fortran_COMPILER=$FC")
fi
cmake "${cmake_args[@]}" > "$OUT/cmake_configure.log" 2>&1

header="$("$BENCH" --header)"
echo "table_id,case_id,repeat,$header" > "$OUT/summary.csv"
echo -e "table_id\tcase_id\trepeat\tmode\torder\tvisc_order\tkeep_prec\tvisc_prec\tnx\tnt\tpress_prec\tslau_rho_prec\tslau_u_prec\tslau_p_prec" > "$OUT/commands.tsv"
echo -e "table_id\tcase_id\trepeat\tstatus\treport_base\tlog" > "$OUT/full_reports.tsv"

sanitize() {
  printf '%s' "$1" | tr ' /,=:' '______' | tr -cd 'A-Za-z0-9_.-'
}

kernel_for_mode() {
  case "$1" in
    seq) echo calc_seq_x;;
    fused) echo calc_fused_x;;
    warp) echo calc_warp_x;;
    warp_fused) echo calc_warp_fused_x;;
    split) echo calc_keep_x;;
    split_visc) echo calc_ev;;
    fill) echo calc_quantities_shadow32;;
    quant) echo calc_quantities_t_1d;;
    slau|slau_weno) echo calc_slau_x;;
    slau_warp|slau_weno_warp) echo calc_slau_warp_x;;
    *) echo "bad mode: $1" >&2; return 1;;
  esac
}

should_collect_full() {
  local rep="$1"
  [ "$COLLECT_FULL" -eq 1 ] || return 1
  case "$FULL_REPEAT" in
    first) [ "$rep" -eq 1 ];;
    all) return 0;;
    none) return 1;;
    *)
      echo "bad --full-repeat: $FULL_REPEAT (use first, all, or none)" >&2
      return 2;;
  esac
}

collect_full_report() {
  local table_id="$1" case_id="$2" mode="$3" rep="$4" stem="$5"
  local kernel report_base log status
  kernel="$(kernel_for_mode "$mode")"
  report_base="$OUT/full_reports/${stem}"
  log="$OUT/full_logs/${stem}.log"

  echo "[$(date +%H:%M:%S)] full ncu ${table_id}/${case_id} repeat ${rep}" | tee -a "$OUT/progress.log"
  set +e
  (cd "$SOLVER_ROOT/ST/build" && "${NCU:-ncu}" \
      --target-processes all --kernel-name "regex:$kernel" \
      --launch-skip "$FULL_NCU_LAUNCH_SKIP" --launch-count "$FULL_NCU_LAUNCH_COUNT" \
      --set "$FULL_NCU_SET" --import-source yes \
      --export "$report_base" --force-overwrite --log-file "$log" \
      "$MPIEXEC" "$MPIEXEC_NP_FLAG" 1 ./a.out)
  status=$?
  set -e

  if [ "$status" -eq 0 ]; then
    echo -e "${table_id}\t${case_id}\t${rep}\tOK\t${report_base}\t${log}" >> "$OUT/full_reports.tsv"
  else
    echo -e "${table_id}\t${case_id}\t${rep}\tFAIL\t${report_base}\t${log}" >> "$OUT/full_reports.tsv"
    if [ "$KEEP_GOING" -eq 0 ]; then
      echo "full ncu failed: ${table_id}/${case_id}; see $log" >&2
      exit "$status"
    fi
  fi
}

run_case() {
  local table_id="$1" case_id="$2" mode="$3" order="$4" visc_order="$5" keep_prec="$6" visc_prec="$7" press_prec="$8" sr="$9" su="${10}" sprec="${11}"
  local rep raw stem line status

  for rep in $(seq 1 "$REPEAT"); do
    stem="$(sanitize "${table_id}__${case_id}__r${rep}")"
    raw="$OUT/raw/${stem}.log"
    timing_report="$OUT/timing_reports/${stem}"
    raw_csv="$OUT/csv/${stem}.csv"
    echo -e "${table_id}\t${case_id}\t${rep}\t${mode}\t${order}\t${visc_order}\t${keep_prec}\t${visc_prec}\t${NX}\t${NT}\t${press_prec}\t${sr}\t${su}\t${sprec}" >> "$OUT/commands.tsv"
    echo "[$(date +%H:%M:%S)] ${table_id}/${case_id} repeat ${rep}/${REPEAT}" | tee -a "$OUT/progress.log"

    set +e
    line="$(NCU_EXPORT="$timing_report" NCU_RAW_CSV="$raw_csv" "$BENCH" \
      "$mode" "$order" "$visc_order" "$keep_prec" "$visc_prec" "$NX" "$NT" \
      "$press_prec" "$sr" "$su" "$sprec" 2>&1 | tee "$raw" | tail -n 1)"
    status=$?
    set -e

    if [ "$status" -ne 0 ]; then
      echo "${table_id},${case_id},${rep},${mode},${order},${visc_order},${keep_prec},${visc_prec},${press_prec},${sr},${su},${sprec},${NX},RUN_FAIL,,,,,,,,,," >> "$OUT/summary.csv"
      if [ "$KEEP_GOING" -eq 0 ]; then
        echo "failed: ${table_id}/${case_id}; see $raw" >&2
        exit "$status"
      fi
    else
      echo "${table_id},${case_id},${rep},${line}" >> "$OUT/summary.csv"
      if should_collect_full "$rep"; then
        collect_full_report "$table_id" "$case_id" "$mode" "$rep" "$stem"
      fi
    fi
  done
}

# Table: kernel mode comparison, KEEP=fp64, VISC=fp32.
for order in 2 4 6; do
  for mode in seq fused warp; do
    run_case "kernel_modes" "${mode}_o${order}" "$mode" "$order" "$order" fp64 fp32 fp64 fp64 fp64 fp64
  done
done

# Table: ORDER=6 precision settings. `seq` is the neutral same-thread baseline.
run_case "precision_o6" "keep64_visc64" seq 6 6 fp64 fp64 fp64 fp64 fp64 fp64
run_case "precision_o6" "keep64_visc32" seq 6 6 fp64 fp32 fp64 fp64 fp64 fp64
run_case "precision_o6" "term_visc32"   seq 6 6 term fp32 fp64 fp64 fp64 fp64
run_case "precision_o6" "keep32_visc32" seq 6 6 fp32 fp32 fp64 fp64 fp64 fp64

# Table: PRESS_PREC='fp32' comparison.
for order in 2 4 6; do
  run_case "press32" "seq_pp64_o${order}"        seq        "$order" "$order" fp64 fp32 fp64 fp64 fp64 fp64
  run_case "press32" "seq_pp32_o${order}"        seq        "$order" "$order" fp64 fp32 fp32 fp64 fp64 fp64
  run_case "press32" "warp_fused_pp32_o${order}" warp_fused "$order" "$order" fp64 fp32 fp32 fp64 fp64 fp64
done

# Table: calc_quantities_shadow32. `quant` measures the ordinary
# calc_quantities_T_1D decode pass; `fill` measures calc_quantities_shadow32.
run_case "shadow" "quant_baseline_o6" quant 6 6 fp64 fp64 fp64 fp64 fp64 fp64
run_case "shadow" "shadow3_u_t_mu_o6" fill  6 6 fp64 fp32 fp64 fp64 fp64 fp64
run_case "shadow" "shadow4_plus_p_o6" fill  6 6 fp64 fp32 fp32 fp64 fp64 fp64
run_case "shadow" "shadow5_term_plus_rho_o6" fill 6 6 term fp32 fp32 fp64 fp64 fp64

# Table: SLAU + MUSCL reconstruction.
run_case "slau_muscl" "split_fp64_o4"       slau      4 4 fp64 fp64 fp64 fp64 fp64 fp64
run_case "slau_muscl" "warp_fp64_o4"        slau_warp 4 4 fp64 fp64 fp64 fp64 fp64 fp64
run_case "slau_muscl" "warp_all_df_o4"      slau_warp 4 4 fp64 fp64 fp64 df   df   df
run_case "slau_muscl" "warp_rho64_o4"       slau_warp 4 4 fp64 fp64 fp64 fp64 df   df
run_case "slau_muscl" "split_fp64_o6"       slau      6 6 fp64 fp64 fp64 fp64 fp64 fp64
run_case "slau_muscl" "warp_fp64_o6"        slau_warp 6 6 fp64 fp64 fp64 fp64 fp64 fp64
run_case "slau_muscl" "warp_all_df_o6"      slau_warp 6 6 fp64 fp64 fp64 df   df   df
run_case "slau_muscl" "warp_rho64_o6"       slau_warp 6 6 fp64 fp64 fp64 fp64 df   df
run_case "slau_muscl" "warp_u64_o6"         slau_warp 6 6 fp64 fp64 fp64 df   fp64 df
run_case "slau_muscl" "warp_p64_o6"         slau_warp 6 6 fp64 fp64 fp64 df   df   fp64

# Table: SLAU + WENO5-Z reconstruction.
run_case "slau_weno" "split_fp64_o6"        slau_weno      6 6 fp64 fp64 fp64 fp64 fp64 fp64
run_case "slau_weno" "warp_fp64_o6"         slau_weno_warp 6 6 fp64 fp64 fp64 fp64 fp64 fp64
run_case "slau_weno" "warp_all_df_o6"       slau_weno_warp 6 6 fp64 fp64 fp64 df   df   df
run_case "slau_weno" "warp_rho64_o6"        slau_weno_warp 6 6 fp64 fp64 fp64 fp64 df   df
run_case "slau_weno" "warp_u64_o6"          slau_weno_warp 6 6 fp64 fp64 fp64 df   fp64 df
run_case "slau_weno" "warp_p64_o6"          slau_weno_warp 6 6 fp64 fp64 fp64 df   df   fp64

cp "$CFG" "$OUT/final_config.fypp"
cp "$GLOB" "$OUT/final_mod_globals.f90"

echo
echo "Wrote:"
echo "  $OUT/summary.csv"
echo "  $OUT/commands.tsv"
echo "  $OUT/env.txt"
echo "  $OUT/raw/"
echo "  $OUT/csv/"
echo "  $OUT/timing_reports/"
echo "  $OUT/full_reports/"
echo "  $OUT/full_logs/"
