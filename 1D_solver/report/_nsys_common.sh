#!/usr/bin/env bash
# Shared helpers for 1D_solver/report/run_*_nsys.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Reuse the existing fypp configuration/build helpers and kernel-name mapping.
# shellcheck source=/dev/null
source "$SCRIPT_DIR/_ncu_common.sh"

nsys_common_defaults() {
  SOLVER_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
  REPO_ROOT="$(cd "$SOLVER_ROOT/.." && pwd)"
  MPIEXEC="${MPIEXEC:-${MPIRUN:-mpirun}}"
  MPIEXEC_NP_FLAG="${MPIEXEC_NP_FLAG:--n}"

  OUT=""
  GPU_CC="${CASE_GPU_CC:-native}"
  NX="${NX:-4194304}"
  NT="${NT:-20}"
  REPEAT=1
  KEEP_GOING=0
  RESUME=0
  NSYS="${NSYS:-nsys}"
  NSYS_TRACE="${NSYS_TRACE:-cuda,nvtx,openacc,osrt}"
  NSYS_LAUNCH_MODE="${NSYS_LAUNCH_MODE:-auto}"
  NSYS_RANKS="${NSYS_RANKS:-1}"
  NSYS_LAUNCH_SKIP="${NSYS_LAUNCH_SKIP:-20}"
  NSYS_LAUNCH_COUNT="${NSYS_LAUNCH_COUNT:-10}"
  NSYS_EXPORT_SQLITE="${NSYS_EXPORT_SQLITE:-1}"
}

nsys_parse_common_arg() {
  case "$1" in
    --out)
      OUT="${2:?missing --out value}"; return 2;;
    --gpu-cc)
      GPU_CC="${2:?missing --gpu-cc value}"; return 2;;
    --nx)
      NX="${2:?missing --nx value}"; return 2;;
    --nt)
      NT="${2:?missing --nt value}"; return 2;;
    --repeat)
      REPEAT="${2:?missing --repeat value}"; return 2;;
    --resume)
      RESUME=1; return 1;;
    --keep-going)
      KEEP_GOING=1; return 1;;
    --trace)
      NSYS_TRACE="${2:?missing --trace value}"; return 2;;
    --launch-mode)
      NSYS_LAUNCH_MODE="${2:?missing --launch-mode value}"; return 2;;
    --ranks)
      NSYS_RANKS="${2:?missing --ranks value}"; return 2;;
    --launch-skip)
      NSYS_LAUNCH_SKIP="${2:?missing --launch-skip value}"; return 2;;
    --launch-count)
      NSYS_LAUNCH_COUNT="${2:?missing --launch-count value}"; return 2;;
    --no-sqlite)
      NSYS_EXPORT_SQLITE=0; return 1;;
  esac
  return 0
}

nsys_common_usage_tail() {
  cat <<'EOF'

Common options:
  --out DIR                 output directory
  --gpu-cc CC              CMake CASE_GPU_CC value, e.g. 80 or native
  --nx N                   grid size, default 4194304
  --nt N                   time steps, default 20
  --repeat N               profile repeats, default 1
  --resume                 resume an interrupted run in the same --out directory
  --keep-going             continue after failures
  --trace LIST             nsys trace list, default cuda,nvtx,openacc,osrt
  --launch-mode MODE       auto, direct, or mpi; default auto
  --ranks N                MPI ranks when launch-mode=mpi; default 1
  --launch-skip N          target-kernel instances to skip in CSV reduction, default 20
  --launch-count N         target-kernel instances to reduce after skip, default 10
  --no-sqlite              skip nsys export --type sqlite

Output includes:
  reports/*.nsys-rep       original Nsight Systems reports
  sqlite/*.sqlite          SQLite exports for agent/script inspection
  stats/*.csv              nsys stats CSV reports
  summary.csv              reduced per-case target-kernel timings
EOF
}

nsys_common_init() {
  local suite_name="$1"

  if [ "$RESUME" -eq 1 ] && [ -z "$OUT" ]; then
    echo "--resume requires the original --out directory" >&2
    exit 2
  fi

  if [ -z "$OUT" ]; then
    local stamp gpu_slug
    stamp="$(date +%Y%m%d_%H%M%S)"
    gpu_slug="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 | tr '[:upper:] /' '[:lower:]__' | tr -cd 'a-z0-9_-' || true)"
    OUT="$SOLVER_ROOT/report/nsys_runs/${suite_name}_${gpu_slug:-gpu}_${stamp}"
  fi

  mkdir -p "$OUT/raw" "$OUT/reports" "$OUT/sqlite" "$OUT/stats" "$OUT/csv" "$OUT/sass"
  OUT="$(cd "$OUT" && pwd)"
  if [ "$RESUME" -eq 0 ]; then
    : > "$OUT/progress.log"
  fi

  CFG="$SOLVER_ROOT/ST/config.fypp"
  GLOB="$SOLVER_ROOT/ST/mod_globals.f90"
  CFG_BAK="$(mktemp)"
  GLOB_BAK="$(mktemp)"
  cp "$CFG" "$CFG_BAK"
  cp "$GLOB" "$GLOB_BAK"
  trap nsys_restore_inputs EXIT

  nsys_write_env "$suite_name"
  nsys_configure_cmake
  nsys_init_output_files
}

nsys_restore_inputs() {
  cp "$CFG_BAK" "$CFG"
  cp "$GLOB_BAK" "$GLOB"
  rm -f "$CFG_BAK" "$GLOB_BAK"
}

nsys_write_env() {
  local suite_name="$1"
  {
    echo "date=$(date -Is)"
    echo "suite=$suite_name"
    echo "repo=$REPO_ROOT"
    echo "solver_root=$SOLVER_ROOT"
    echo "out=$OUT"
    echo "gpu_cc=$GPU_CC"
    echo "nx=$NX"
    echo "nt=$NT"
    echo "repeat=$REPEAT"
    echo "resume=$RESUME"
    echo "nsys_trace=$NSYS_TRACE"
    echo "nsys_launch_mode=$NSYS_LAUNCH_MODE"
    echo "nsys_ranks=$NSYS_RANKS"
    echo "nsys_launch_skip=$NSYS_LAUNCH_SKIP"
    echo "nsys_launch_count=$NSYS_LAUNCH_COUNT"
    echo "nsys_export_sqlite=$NSYS_EXPORT_SQLITE"
    echo
    echo "[gpu]"
    nvidia-smi --query-gpu=name,driver_version,compute_cap,pci.bus_id --format=csv 2>/dev/null || true
    echo
    echo "[nsys]"
    "$NSYS" --version 2>&1 || true
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
}

nsys_configure_cmake() {
  local cmake_args=(-S "$SOLVER_ROOT/ST" -B "$SOLVER_ROOT/ST/build" -DCASE_GPU_CC="$GPU_CC")
  if [ -n "${FC:-}" ]; then
    cmake_args+=("-DCMAKE_Fortran_COMPILER=$FC")
  fi
  if [ "$RESUME" -eq 1 ]; then
    cmake "${cmake_args[@]}" >> "$OUT/cmake_configure.log" 2>&1
  else
    cmake "${cmake_args[@]}" > "$OUT/cmake_configure.log" 2>&1
  fi
}

nsys_init_output_files() {
  nsys_init_output_file "$OUT/summary.csv" \
    "table_id,case_id,repeat,mode,order,visc_order,keep_prec,visc_prec,press_prec,slau_rho_prec,slau_u_prec,slau_p_prec,nx,kernel,total_instances,used_instances,time_us,time_min_us,time_max_us,spread_us,total_us,trace_csv,kern_sum_csv,report,sqlite,status"
  nsys_init_output_file "$OUT/commands.tsv" \
    $'table_id\tcase_id\trepeat\tmode\torder\tvisc_order\tkeep_prec\tvisc_prec\tnx\tnt\tpress_prec\tslau_rho_prec\tslau_u_prec\tslau_p_prec'
}

nsys_init_output_file() {
  local path="$1" header_line="$2"
  if [ "$RESUME" -eq 1 ] && [ -s "$path" ]; then
    return 0
  fi
  printf '%s\n' "$header_line" > "$path"
}

nsys_timing_complete() {
  local table_id="$1" case_id="$2" rep="$3"
  [ "$RESUME" -eq 1 ] || return 1
  [ -s "$OUT/summary.csv" ] || return 1

  awk -F, -v t="$table_id" -v c="$case_id" -v r="$rep" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "time_us") time_col = i
        if ($i == "status") status_col = i
      }
      next
    }
    NR > 1 && $1 == t && $2 == c && $3 == r &&
    time_col > 0 && status_col > 0 &&
    $time_col != "" && $time_col !~ /FAIL|[Nn][Aa][Nn]|[Ii][Nn][Ff]/ &&
    $status_col !~ /FAIL/ {
      found = 1
    }
    END { exit found ? 0 : 1 }
  ' "$OUT/summary.csv"
}

nsys_prune_rows() {
  local table_id="$1" case_id="$2" rep="$3" tmp
  [ "$RESUME" -eq 1 ] || return 0
  [ -s "$OUT/summary.csv" ] || return 0

  tmp="$(mktemp)"
  awk -F, -v t="$table_id" -v c="$case_id" -v r="$rep" \
    'NR == 1 || !($1 == t && $2 == c && $3 == r)' \
    "$OUT/summary.csv" > "$tmp"
  mv "$tmp" "$OUT/summary.csv"

  if [ -s "$OUT/commands.tsv" ]; then
    tmp="$(mktemp)"
    awk -F '\t' -v t="$table_id" -v c="$case_id" -v r="$rep" \
      'NR == 1 || !($1 == t && $2 == c && $3 == r)' \
      "$OUT/commands.tsv" > "$tmp"
    mv "$tmp" "$OUT/commands.tsv"
  fi
}

nsys_launch_attempts() {
  case "$NSYS_LAUNCH_MODE" in
    auto)
      if [ "$NSYS_RANKS" = 1 ]; then
        printf '%s\n' direct mpi
      else
        printf '%s\n' mpi
      fi
      ;;
    direct)
      printf '%s\n' direct;;
    mpi)
      printf '%s\n' mpi;;
    *)
      echo "bad NSYS_LAUNCH_MODE: $NSYS_LAUNCH_MODE (use auto, direct, or mpi)" >&2
      return 2;;
  esac
}

nsys_setup_launch() {
  local launch_mode="$1"
  if [ "$launch_mode" = direct ]; then
    NSYS_LAUNCH_CMD=(./a.out)
    return 0
  fi
  NSYS_LAUNCH_CMD=("$MPIEXEC" "$MPIEXEC_NP_FLAG" "$NSYS_RANKS" ./a.out)
}

nsys_find_stats_csv() {
  local base="$1"
  find "$OUT/stats" -maxdepth 1 -type f -name "$(basename "$base")*.csv" | sort | head -1
}

nsys_make_stats() {
  local report="$1" stem="$2" trace_base kern_base
  trace_base="$OUT/stats/${stem}__trace"
  kern_base="$OUT/stats/${stem}__kern_sum"

  "$NSYS" stats -q --force-export true --force-overwrite true \
    --report cuda_gpu_trace:base --format csv --output "$trace_base" "$report" > "$OUT/raw/${stem}__stats_trace.log" 2>&1 || true
  "$NSYS" stats -q --force-export true --force-overwrite true \
    --report cuda_gpu_kern_sum:base --format csv --output "$kern_base" "$report" > "$OUT/raw/${stem}__stats_kern_sum.log" 2>&1 || true

  NSYS_TRACE_CSV="$(nsys_find_stats_csv "$trace_base")"
  NSYS_KERN_SUM_CSV="$(nsys_find_stats_csv "$kern_base")"
}

nsys_run_case() {
  local table_id="$1" case_id="$2" mode="$3" order="$4" visc_order="$5" keep_prec="$6" visc_prec="$7" press_prec="$8" sr="$9" su="${10}" sprec="${11}"
  local rep stem report_base report raw_log build_log sqlite kernel line launch_mode status

  kernel="$(ncu_kernel_for_mode "$mode")"

  for rep in $(seq 1 "$REPEAT"); do
    stem="$(ncu_sanitize "${table_id}__${case_id}__r${rep}")"
    report_base="$OUT/reports/${stem}"
    report="${report_base}.nsys-rep"
    raw_log="$OUT/raw/${stem}.log"
    build_log="$OUT/raw/${stem}__build.log"
    sqlite="$OUT/sqlite/${stem}.sqlite"

    if nsys_timing_complete "$table_id" "$case_id" "$rep"; then
      echo "[$(date +%H:%M:%S)] resume skip complete ${table_id}/${case_id} repeat ${rep}/${REPEAT}" | tee -a "$OUT/progress.log"
      continue
    fi

    nsys_prune_rows "$table_id" "$case_id" "$rep"
    echo -e "${table_id}\t${case_id}\t${rep}\t${mode}\t${order}\t${visc_order}\t${keep_prec}\t${visc_prec}\t${NX}\t${NT}\t${press_prec}\t${sr}\t${su}\t${sprec}" >> "$OUT/commands.tsv"
    echo "[$(date +%H:%M:%S)] nsys ${table_id}/${case_id} repeat ${rep}/${REPEAT}" | tee -a "$OUT/progress.log"

    ncu_configure_case "$mode" "$order" "$visc_order" "$keep_prec" "$visc_prec" "$press_prec" "$sr" "$su" "$sprec"
    if ! ncu_build_case "$build_log"; then
      echo "${table_id},${case_id},${rep},${mode},${order},${visc_order},${keep_prec},${visc_prec},${press_prec},${sr},${su},${sprec},${NX},${kernel},0,0,NSYS_FAIL,,,,,,,,BUILD_FAIL" >> "$OUT/summary.csv"
      if [ "$KEEP_GOING" -eq 0 ]; then
        echo "build failed: ${table_id}/${case_id}; see $build_log" >&2
        exit 1
      fi
      continue
    fi

    status=1
    while IFS= read -r launch_mode; do
      nsys_setup_launch "$launch_mode"
      set +e
      (cd "$SOLVER_ROOT/ST/build" && "$NSYS" profile \
        --trace="$NSYS_TRACE" --sample=none --cpuctxsw=none --backtrace=none \
        --force-overwrite=true --output="$report_base" \
        "${NSYS_LAUNCH_CMD[@]}") > "$raw_log" 2>&1
      status=$?
      set -e
      if [ "$status" -eq 0 ] && [ -s "$report" ]; then
        break
      fi
      if [ "$NSYS_LAUNCH_MODE" = auto ]; then
        echo "[$(date +%H:%M:%S)] retry nsys ${table_id}/${case_id} launch=${launch_mode} failed (status=${status})" | tee -a "$OUT/progress.log" >&2
      fi
    done < <(nsys_launch_attempts)

    if [ "$status" -ne 0 ] || [ ! -s "$report" ]; then
      echo "${table_id},${case_id},${rep},${mode},${order},${visc_order},${keep_prec},${visc_prec},${press_prec},${sr},${su},${sprec},${NX},${kernel},0,0,NSYS_FAIL,,,,,,,,PROFILE_FAIL" >> "$OUT/summary.csv"
      if [ "$KEEP_GOING" -eq 0 ]; then
        echo "nsys profile failed: ${table_id}/${case_id}; see $raw_log" >&2
        exit "$status"
      fi
      continue
    fi

    if [ "$NSYS_EXPORT_SQLITE" -eq 1 ]; then
      "$NSYS" export --type=sqlite --force-overwrite=true --output="$sqlite" "$report" > "$OUT/raw/${stem}__export_sqlite.log" 2>&1 || true
    else
      sqlite=""
    fi

    NSYS_TRACE_CSV=""
    NSYS_KERN_SUM_CSV=""
    nsys_make_stats "$report" "$stem"

    line="$("$SOLVER_ROOT/report/_nsys_parse.py" \
      --mode "$mode" --order "$order" --visc-order "$visc_order" \
      --keep-prec "$keep_prec" --visc-prec "$visc_prec" --press-prec "$press_prec" \
      --slau-rho-prec "$sr" --slau-u-prec "$su" --slau-p-prec "$sprec" \
      --nx "$NX" --kernel "$kernel" --skip "$NSYS_LAUNCH_SKIP" --count "$NSYS_LAUNCH_COUNT" \
      --trace-csv "$NSYS_TRACE_CSV" --kern-sum-csv "$NSYS_KERN_SUM_CSV" \
      --report "$report" --sqlite "$sqlite")"
    echo "${table_id},${case_id},${rep},${line}" >> "$OUT/summary.csv"

    if [[ "$line" == *NSYS_FAIL* ]] && [ "$KEEP_GOING" -eq 0 ]; then
      echo "nsys stats did not find kernel ${kernel}: ${table_id}/${case_id}; see $OUT/stats/" >&2
      exit 1
    fi
  done
}

nsys_finish() {
  cp "$CFG" "$OUT/final_config.fypp"
  cp "$GLOB" "$OUT/final_mod_globals.f90"

  echo
  echo "Wrote:"
  echo "  $OUT/summary.csv"
  echo "  $OUT/commands.tsv"
  echo "  $OUT/env.txt"
  echo "  $OUT/raw/"
  echo "  $OUT/reports/"
  echo "  $OUT/sqlite/"
  echo "  $OUT/stats/"
  echo "  $OUT/csv/"
  echo "  $OUT/sass/"
}
