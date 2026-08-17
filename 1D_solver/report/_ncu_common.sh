#!/usr/bin/env bash
# Shared helpers for 1D_solver/report/run_*_ncu.sh.

ncu_common_defaults() {
  SOLVER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  REPO_ROOT="$(cd "$SOLVER_ROOT/.." && pwd)"
  BENCH="$SOLVER_ROOT/report/bench.sh"
  MPIEXEC="${MPIEXEC:-${MPIRUN:-mpirun}}"
  MPIEXEC_NP_FLAG="${MPIEXEC_NP_FLAG:--n}"
  NCU_LAUNCH_MODE="${NCU_LAUNCH_MODE:-auto}"
  NCU_RANKS="${NCU_RANKS:-1}"

  OUT=""
  GPU_CC="${CASE_GPU_CC:-native}"
  NX="${NX:-4194304}"
  NT="${NT:-20}"
  REPEAT=1
  KEEP_GOING=0
  RESUME=0
  COLLECT_FULL=1
  FULL_REPEAT="${FULL_REPEAT:-first}"
  FULL_NCU_MODE="${FULL_NCU_MODE:-sections}"
  FULL_NCU_SET="${FULL_NCU_SET:-full}"
  FULL_NCU_SECTIONS="${FULL_NCU_SECTIONS:-LaunchStats Occupancy SpeedOfLight WorkloadDistribution ComputeWorkloadAnalysis MemoryWorkloadAnalysis InstructionStats SchedulerStats WarpStateStats}"
  FULL_NCU_IMPORT_SOURCE="${FULL_NCU_IMPORT_SOURCE:-no}"
  FULL_NCU_REPLAY_MODE="${FULL_NCU_REPLAY_MODE:-kernel}"
  FULL_NCU_FALLBACK_SECTIONS="${FULL_NCU_FALLBACK_SECTIONS:-1}"
  FULL_NCU_REQUIRED="${FULL_NCU_REQUIRED:-0}"
  FULL_NCU_LAUNCH_SKIP="${FULL_NCU_LAUNCH_SKIP:-20}"
  FULL_NCU_LAUNCH_COUNT="${FULL_NCU_LAUNCH_COUNT:-1}"
}

ncu_parse_common_arg() {
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
    --keep-going)
      KEEP_GOING=1; return 1;;
    --resume)
      RESUME=1; return 1;;
    --no-full-ncu)
      COLLECT_FULL=0; return 1;;
    --full-repeat)
      FULL_REPEAT="${2:?missing --full-repeat value}"; return 2;;
    --full-mode)
      FULL_NCU_MODE="${2:?missing --full-mode value}"; return 2;;
    --full-set)
      FULL_NCU_SET="${2:?missing --full-set value}"; return 2;;
    --full-sections)
      FULL_NCU_SECTIONS="${2:?missing --full-sections value}"; return 2;;
    --full-import-source)
      FULL_NCU_IMPORT_SOURCE="${2:?missing --full-import-source value}"; return 2;;
    --full-replay-mode)
      FULL_NCU_REPLAY_MODE="${2:?missing --full-replay-mode value}"; return 2;;
    --no-full-fallback)
      FULL_NCU_FALLBACK_SECTIONS=0; return 1;;
    --require-full-ncu)
      FULL_NCU_REQUIRED=1; return 1;;
    --full-launch-skip)
      FULL_NCU_LAUNCH_SKIP="${2:?missing --full-launch-skip value}"; return 2;;
    --full-launch-count)
      FULL_NCU_LAUNCH_COUNT="${2:?missing --full-launch-count value}"; return 2;;
  esac
  return 0
}

ncu_common_usage_tail() {
  cat <<'EOF'

Common options:
  --out DIR                 output directory
  --gpu-cc CC              CMake CASE_GPU_CC value, e.g. 80 or native
  --nx N                   grid size, default 4194304
  --nt N                   time steps, default 20
  --repeat N               timing repeats, default 1
  --resume                 resume an interrupted run in the same --out directory
  --keep-going             continue after timing failures
  --no-full-ncu            collect timing CSV/reports only
  --full-repeat MODE       first, all, or none; default first
  --full-mode MODE         sections or set; default sections
  --full-set SET           ncu section set for --full-mode set; default full
  --full-sections "..."    section identifiers for --full-mode sections
  --require-full-ncu       abort if a full report fails

MPI is taken from PATH after module load. Override with MPIEXEC=mpiexec or
MPIEXEC_NP_FLAG=-np if needed. Under ncu, NCU_LAUNCH_MODE=auto prefers direct
launch for 1-rank jobs and falls back to mpirun only if needed.
EOF
}

ncu_common_init() {
  local suite_name="$1"

  if [ "$RESUME" -eq 1 ] && [ -z "$OUT" ]; then
    echo "--resume requires the original --out directory" >&2
    exit 2
  fi

  if [ -z "$OUT" ]; then
    local stamp gpu_slug
    stamp="$(date +%Y%m%d_%H%M%S)"
    gpu_slug="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 | tr '[:upper:] /' '[:lower:]__' | tr -cd 'a-z0-9_-' || true)"
    OUT="$SOLVER_ROOT/report/ncu_runs/${suite_name}_${gpu_slug:-gpu}_${stamp}"
  fi

  mkdir -p "$OUT/raw" "$OUT/csv" "$OUT/timing_reports" "$OUT/full_reports" "$OUT/full_logs" "$OUT/sass"
  OUT="$(cd "$OUT" && pwd)"
  if [ "$RESUME" -eq 0 ]; then
    : > "$OUT/progress.log"
  fi

  FULL_NCU_SECTIONS_REQUESTED="$FULL_NCU_SECTIONS"
  FULL_NCU_SECTIONS_SKIPPED=""
  if [ "$COLLECT_FULL" -eq 1 ] &&
     { [ "$FULL_NCU_MODE" = sections ] ||
       { [ "$FULL_NCU_MODE" = set ] && [ "$FULL_NCU_SET" = full ] && [ "$FULL_NCU_FALLBACK_SECTIONS" -eq 1 ]; }; }; then
    ncu_select_full_sections "$FULL_NCU_SECTIONS"
  fi

  CFG="$SOLVER_ROOT/ST/config.fypp"
  GLOB="$SOLVER_ROOT/ST/mod_globals.f90"
  CFG_BAK="$(mktemp)"
  GLOB_BAK="$(mktemp)"
  cp "$CFG" "$CFG_BAK"
  cp "$GLOB" "$GLOB_BAK"
  trap ncu_restore_inputs EXIT

  ncu_write_env "$suite_name"
  ncu_configure_cmake
  ncu_init_output_files
}

ncu_restore_inputs() {
  cp "$CFG_BAK" "$CFG"
  cp "$GLOB_BAK" "$GLOB"
  rm -f "$CFG_BAK" "$GLOB_BAK"
}

ncu_select_full_sections() {
  local requested="$1" selected="" skipped="" section
  for section in $requested; do
    if "${NCU:-ncu}" --section "$section" --list-metrics > /dev/null 2>&1; then
      selected="${selected:+$selected }$section"
    else
      skipped="${skipped:+$skipped }$section"
    fi
  done

  FULL_NCU_SECTIONS="$selected"
  FULL_NCU_SECTIONS_SKIPPED="$skipped"
  if [ -n "$FULL_NCU_SECTIONS_SKIPPED" ]; then
    echo "Skipping unavailable ncu sections: $FULL_NCU_SECTIONS_SKIPPED" | tee -a "$OUT/progress.log" >&2
  fi
  if [ "$COLLECT_FULL" -eq 1 ] && [ "$FULL_NCU_MODE" = sections ] && [ -z "$FULL_NCU_SECTIONS" ]; then
    echo "No requested ncu sections are available; falling back to --set basic" | tee -a "$OUT/progress.log" >&2
    FULL_NCU_MODE=set
    FULL_NCU_SET=basic
  fi
}

ncu_write_env() {
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
    echo "bench_weno_order=${BENCH_WENO_ORDER:-5}"
    echo "repeat=$REPEAT"
    echo "resume=$RESUME"
    echo "collect_full=$COLLECT_FULL"
    echo "full_repeat=$FULL_REPEAT"
    echo "full_ncu_mode=$FULL_NCU_MODE"
    echo "full_ncu_set=$FULL_NCU_SET"
    echo "full_ncu_sections_requested=$FULL_NCU_SECTIONS_REQUESTED"
    echo "full_ncu_sections=$FULL_NCU_SECTIONS"
    echo "full_ncu_sections_skipped=$FULL_NCU_SECTIONS_SKIPPED"
    echo "full_ncu_import_source=$FULL_NCU_IMPORT_SOURCE"
    echo "full_ncu_replay_mode=$FULL_NCU_REPLAY_MODE"
    echo "full_ncu_fallback_sections=$FULL_NCU_FALLBACK_SECTIONS"
    echo "full_ncu_required=$FULL_NCU_REQUIRED"
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
    echo "NCU_LAUNCH_MODE=$NCU_LAUNCH_MODE"
    echo "NCU_RANKS=$NCU_RANKS"
    command -v "$MPIEXEC" 2>/dev/null || true
    "$MPIEXEC" --version 2>&1 | head -20 || true
    echo
    echo "[git]"
    git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || true
    git -C "$REPO_ROOT" status --short 2>/dev/null || true
  } > "$OUT/env.txt"
}

ncu_configure_cmake() {
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

ncu_init_output_files() {
  local header
  header="$("$BENCH" --header)"
  ncu_init_output_file "$OUT/summary.csv" "table_id,case_id,repeat,$header"
  ncu_init_output_file "$OUT/commands.tsv" $'table_id\tcase_id\trepeat\tmode\torder\tvisc_order\tkeep_prec\tvisc_prec\tnx\tnt\tpress_prec\tslau_rho_prec\tslau_u_prec\tslau_p_prec'
  ncu_init_output_file "$OUT/full_reports.tsv" $'table_id\tcase_id\trepeat\tprofile\tstatus\treport_base\tlog'
}

ncu_init_output_file() {
  local path="$1" header_line="$2"
  if [ "$RESUME" -eq 1 ] && [ -s "$path" ]; then
    return 0
  fi
  printf '%s\n' "$header_line" > "$path"
}

ncu_sanitize() {
  printf '%s' "$1" | tr ' /,=:' '______' | tr -cd 'A-Za-z0-9_.-'
}

ncu_kernel_for_mode() {
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
    slau_smem|slau_weno_smem) echo calc_slau_smem_x;;
    slau_warp|slau_weno_warp) echo calc_slau_warp_x;;
    *) echo "bad mode: $1" >&2; return 1;;
  esac
}

ncu_timing_complete() {
  local table_id="$1" case_id="$2" rep="$3"
  [ "$RESUME" -eq 1 ] || return 1
  [ -s "$OUT/summary.csv" ] || return 1

  awk -F, -v t="$table_id" -v c="$case_id" -v r="$rep" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "time_us") {
          time_col = i
          break
        }
      }
      next
    }
    time_col == 0 {
      exit 1
    }
    NR > 1 && $1 == t && $2 == c && $3 == r &&
    $time_col != "" &&
    $time_col !~ /(^|_)FAIL$/ &&
    $time_col !~ /^[Nn][Aa][Nn]$/ &&
    $time_col !~ /^[Ii][Nn][Ff]$/ &&
    $time_col !~ /^-[Ii][Nn][Ff]$/ {
      found = 1
    }
    END { exit found ? 0 : 1 }
  ' "$OUT/summary.csv"
}

ncu_full_profile_complete() {
  local table_id="$1" case_id="$2" rep="$3" profile="$4" report_base="$5"
  [ "$RESUME" -eq 1 ] || return 1
  [ -s "$OUT/full_reports.tsv" ] || return 1
  [ -s "${report_base}.ncu-rep" ] || return 1

  awk -F '\t' -v t="$table_id" -v c="$case_id" -v r="$rep" -v p="$profile" '
    NR > 1 && $1 == t && $2 == c && $3 == r && $4 == p && $5 == "OK" {
      found = 1
    }
    END { exit found ? 0 : 1 }
  ' "$OUT/full_reports.tsv"
}

ncu_prune_timing_rows() {
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

ncu_prune_full_profile_rows() {
  local table_id="$1" case_id="$2" rep="$3" profile="$4" tmp
  [ "$RESUME" -eq 1 ] || return 0
  [ -s "$OUT/full_reports.tsv" ] || return 0

  tmp="$(mktemp)"
  awk -F '\t' -v t="$table_id" -v c="$case_id" -v r="$rep" -v p="$profile" \
    'NR == 1 || !($1 == t && $2 == c && $3 == r && $4 == p)' \
    "$OUT/full_reports.tsv" > "$tmp"
  mv "$tmp" "$OUT/full_reports.tsv"
}

ncu_should_collect_full() {
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

ncu_full_collection_complete() {
  local table_id="$1" case_id="$2" rep="$3" stem="$4"
  local section section_slug report_base profile

  [ "$COLLECT_FULL" -eq 1 ] || return 0
  case "$FULL_REPEAT" in
    first)
      [ "$rep" -eq 1 ] || return 0
      ;;
    all)
      ;;
    none)
      return 0
      ;;
  esac

  case "$FULL_NCU_MODE" in
    sections)
      for section in $FULL_NCU_SECTIONS; do
        section_slug="$(ncu_sanitize "$section")"
        report_base="$OUT/full_reports/${stem}__section_${section_slug}"
        profile="section:${section}"
        ncu_full_profile_complete "$table_id" "$case_id" "$rep" "$profile" "$report_base" || return 1
      done
      return 0
      ;;
    set)
      report_base="$OUT/full_reports/${stem}__set_${FULL_NCU_SET}"
      profile="set:${FULL_NCU_SET}"
      ncu_full_profile_complete "$table_id" "$case_id" "$rep" "$profile" "$report_base" && return 0
      if [ "$FULL_NCU_SET" = full ] && [ "$FULL_NCU_FALLBACK_SECTIONS" -eq 1 ] && [ -n "$FULL_NCU_SECTIONS" ]; then
        for section in $FULL_NCU_SECTIONS; do
          section_slug="$(ncu_sanitize "$section")"
          report_base="$OUT/full_reports/${stem}__section_${section_slug}"
          profile="section:${section}"
          ncu_full_profile_complete "$table_id" "$case_id" "$rep" "$profile" "$report_base" || return 1
        done
        return 0
      fi
      return 1
      ;;
  esac
  return 1
}

ncu_configure_case() {
  local mode="$1" order="$2" visc_order="$3" keep_prec="$4" visc_prec="$5" press_prec="$6" sr="$7" su="$8" sprec="$9"
  local cmode="$mode" scheme=KEEP tvd=none recon=MUSCL
  local weno_order="${BENCH_WENO_ORDER:-5}"

  [ "$mode" = split_visc ] && cmode=split
  [ "$mode" = quant ] && cmode=seq
  [ "$mode" = fill ] && cmode=seq
  case "$mode" in
    slau)
      scheme=SLAU; tvd=tvd; cmode=split;;
    slau_smem)
      scheme=SLAU; tvd=tvd; cmode=smem;;
    slau_warp)
      scheme=SLAU; tvd=tvd; cmode=warp;;
    slau_weno)
      scheme=SLAU; tvd=tvd; recon=WENO; cmode=split;;
    slau_weno_smem)
      scheme=SLAU; tvd=tvd; recon=WENO; cmode=smem;;
    slau_weno_warp)
      scheme=SLAU; tvd=tvd; recon=WENO; cmode=warp;;
  esac

  python3 - "$cmode" "$order" "$visc_order" "$keep_prec" "$visc_prec" "$NX" "$NT" "$SOLVER_ROOT" "$press_prec" "$scheme" "$tvd" "$recon" "$sr" "$su" "$sprec" "$weno_order" <<'EOF'
import pathlib
import re
import sys

mode,o,vo,kp,vp,nx,nt,root,pp,scheme,tvd,recon,sr,su,sprec,weno_order = sys.argv[1:17]
R = pathlib.Path(root)
p = R/'ST/config.fypp'
s = p.read_text()
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
s = re.sub(r"(?m)^#:set WENO_ORDER\s*=.*$", f"#:set WENO_ORDER   = {weno_order}", s)
p.write_text(s)
g = R/'ST/mod_globals.f90'
s = g.read_text()
s = re.sub(r"(?m)^(\s*integer, parameter :: nx = )\d+", rf"\g<1>{nx}", s)
s = re.sub(r"(?m)^(\s*integer, parameter :: nt = )\d+", rf"\g<1>{nt}", s)
g.write_text(s)
EOF
}

ncu_build_case() {
  local build_log="$1"
  if ! cmake --build "$SOLVER_ROOT/ST/build" -j > "$build_log" 2>&1; then
    cmake --build "$SOLVER_ROOT/ST/build" --target clean >> "$build_log" 2>&1
    cmake --build "$SOLVER_ROOT/ST/build" -j >> "$build_log" 2>&1 || return 1
  fi
}

ncu_prepare_case_for_full() {
  local stem="$1" mode="$2" order="$3" visc_order="$4" keep_prec="$5" visc_prec="$6" press_prec="$7" sr="$8" su="$9" sprec="${10}"
  ncu_configure_case "$mode" "$order" "$visc_order" "$keep_prec" "$visc_prec" "$press_prec" "$sr" "$su" "$sprec"
  ncu_build_case "$OUT/raw/${stem}__resume_build.log"
}

ncu_collect_full_report() {
  local table_id="$1" case_id="$2" mode="$3" rep="$4" stem="$5"
  local kernel status
  kernel="$(ncu_kernel_for_mode "$mode")"

  ncu_launch_attempts() {
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
        return 2
        ;;
    esac
  }

  ncu_setup_launch() {
    local launch_mode="$1"
    if [ "$launch_mode" = direct ]; then
      NCU_TARGET_PROCESSES=application-only
      NCU_LAUNCH_DESC=direct
      NCU_LAUNCH_CMD=(./a.out)
      return 0
    fi

    NCU_TARGET_PROCESSES=all
    NCU_LAUNCH_DESC="${MPIEXEC} ${MPIEXEC_NP_FLAG} ${NCU_RANKS} ./a.out"
    NCU_LAUNCH_CMD=("$MPIEXEC" "$MPIEXEC_NP_FLAG" "$NCU_RANKS" ./a.out)
  }

  ncu_run_full_ncu() {
    local profile="$1" report_base="$2" log="$3"
    local launch_mode attempt_status
    shift 3

    if ncu_full_profile_complete "$table_id" "$case_id" "$rep" "$profile" "$report_base"; then
      echo "[$(date +%H:%M:%S)] resume skip full ncu ${profile} ${table_id}/${case_id} repeat ${rep}" | tee -a "$OUT/progress.log"
      return 0
    fi
    ncu_prune_full_profile_rows "$table_id" "$case_id" "$rep" "$profile"

    echo "[$(date +%H:%M:%S)] full ncu ${profile} ${table_id}/${case_id} repeat ${rep}" | tee -a "$OUT/progress.log"
    status=1
    while IFS= read -r launch_mode; do
      ncu_setup_launch "$launch_mode"
      set +e
      (cd "$SOLVER_ROOT/ST/build" && "${NCU:-ncu}" \
        --target-processes "$NCU_TARGET_PROCESSES" --kernel-name "regex:$kernel" \
        --launch-skip "$FULL_NCU_LAUNCH_SKIP" --launch-count "$FULL_NCU_LAUNCH_COUNT" \
        --replay-mode "$FULL_NCU_REPLAY_MODE" \
        "$@" --import-source "$FULL_NCU_IMPORT_SOURCE" \
        --export "$report_base" --force-overwrite --log-file "$log" \
        "${NCU_LAUNCH_CMD[@]}")
      attempt_status=$?
      set -e

      if [ "$attempt_status" -eq 0 ]; then
        status=0
        break
      fi

      if [ "$NCU_LAUNCH_MODE" = auto ]; then
        echo "[$(date +%H:%M:%S)] retry full ncu ${profile} ${table_id}/${case_id} with launch=${launch_mode} failed (status=${attempt_status})" | tee -a "$OUT/progress.log" >&2
      fi
      status="$attempt_status"
    done < <(ncu_launch_attempts)

    if [ "$status" -eq 0 ]; then
      echo -e "${table_id}\t${case_id}\t${rep}\t${profile}\tOK\t${report_base}\t${log}" >> "$OUT/full_reports.tsv"
    elif [ -s "${report_base}.ncu-rep" ]; then
      echo -e "${table_id}\t${case_id}\t${rep}\t${profile}\tWARN_EXIT_${status}\t${report_base}\t${log}" >> "$OUT/full_reports.tsv"
      status=0
    else
      echo -e "${table_id}\t${case_id}\t${rep}\t${profile}\tFAIL\t${report_base}\t${log}" >> "$OUT/full_reports.tsv"
    fi
    return "$status"
  }

  ncu_collect_full_sections() {
    local section section_slug report_base log section_status overall_status=0
    for section in $FULL_NCU_SECTIONS; do
      section_slug="$(ncu_sanitize "$section")"
      report_base="$OUT/full_reports/${stem}__section_${section_slug}"
      log="$OUT/full_logs/${stem}__section_${section_slug}.log"
      set +e
      ncu_run_full_ncu "section:${section}" "$report_base" "$log" --section "$section"
      section_status=$?
      set -e
      if [ "$section_status" -ne 0 ]; then
        overall_status="$section_status"
        echo "full ncu section failed: ${table_id}/${case_id} ${section}; see $log" >&2
        if [ "$FULL_NCU_REQUIRED" -eq 1 ] && [ "$KEEP_GOING" -eq 0 ]; then
          return "$overall_status"
        fi
      fi
    done
    return "$overall_status"
  }

  case "$FULL_NCU_MODE" in
    sections)
      set +e
      ncu_collect_full_sections
      status=$?
      set -e
      ;;
    set)
      set +e
      ncu_run_full_ncu "set:${FULL_NCU_SET}" "$OUT/full_reports/${stem}__set_${FULL_NCU_SET}" "$OUT/full_logs/${stem}__set_${FULL_NCU_SET}.log" --set "$FULL_NCU_SET"
      status=$?
      set -e
      if [ "$status" -ne 0 ] && [ "$FULL_NCU_SET" = full ] && [ "$FULL_NCU_FALLBACK_SECTIONS" -eq 1 ]; then
        echo "full ncu --set full failed for ${table_id}/${case_id}; retrying split sections" | tee -a "$OUT/progress.log" >&2
        set +e
        ncu_collect_full_sections
        status=$?
        set -e
      fi
      ;;
    *)
      echo "bad --full-mode: $FULL_NCU_MODE (use sections or set)" >&2
      exit 2
      ;;
  esac

  if [ "$status" -ne 0 ] && [ "$FULL_NCU_REQUIRED" -eq 1 ] && [ "$KEEP_GOING" -eq 0 ]; then
    echo "full ncu failed: ${table_id}/${case_id}; see $OUT/full_logs/" >&2
    exit "$status"
  fi
}

ncu_run_case() {
  local table_id="$1" case_id="$2" mode="$3" order="$4" visc_order="$5" keep_prec="$6" visc_prec="$7" press_prec="$8" sr="$9" su="${10}" sprec="${11}"
  local rep raw stem line status timing_report raw_csv bench_mode bench_scheme bench_tvd bench_recon bench_kernel

  bench_mode="$mode"
  bench_scheme=KEEP
  bench_tvd=none
  bench_recon=MUSCL
  bench_kernel="$(ncu_kernel_for_mode "$mode")"
  case "$mode" in
    slau)
      bench_mode=split
      bench_scheme=SLAU
      bench_tvd=tvd
      bench_recon=MUSCL;;
    slau_smem)
      bench_mode=smem
      bench_scheme=SLAU
      bench_tvd=tvd
      bench_recon=MUSCL;;
    slau_warp)
      bench_mode=warp
      bench_scheme=SLAU
      bench_tvd=tvd
      bench_recon=MUSCL;;
    slau_weno)
      bench_mode=split
      bench_scheme=SLAU
      bench_tvd=tvd
      bench_recon=WENO;;
    slau_weno_smem)
      bench_mode=smem
      bench_scheme=SLAU
      bench_tvd=tvd
      bench_recon=WENO;;
    slau_weno_warp)
      bench_mode=warp
      bench_scheme=SLAU
      bench_tvd=tvd
      bench_recon=WENO;;
  esac

  for rep in $(seq 1 "$REPEAT"); do
    stem="$(ncu_sanitize "${table_id}__${case_id}__r${rep}")"
    raw="$OUT/raw/${stem}.log"
    timing_report="$OUT/timing_reports/${stem}"
    raw_csv="$OUT/csv/${stem}.csv"

    if ncu_timing_complete "$table_id" "$case_id" "$rep"; then
      if ncu_full_collection_complete "$table_id" "$case_id" "$rep" "$stem"; then
        echo "[$(date +%H:%M:%S)] resume skip complete ${table_id}/${case_id} repeat ${rep}/${REPEAT}" | tee -a "$OUT/progress.log"
        continue
      fi

      echo "[$(date +%H:%M:%S)] resume skip timing; complete missing full ncu ${table_id}/${case_id} repeat ${rep}/${REPEAT}" | tee -a "$OUT/progress.log"
      if ncu_prepare_case_for_full "$stem" "$mode" "$order" "$visc_order" "$keep_prec" "$visc_prec" "$press_prec" "$sr" "$su" "$sprec"; then
        ncu_collect_full_report "$table_id" "$case_id" "$mode" "$rep" "$stem"
      else
        echo "resume rebuild failed: ${table_id}/${case_id}; see $OUT/raw/${stem}__resume_build.log" >&2
        if [ "$FULL_NCU_REQUIRED" -eq 1 ] && [ "$KEEP_GOING" -eq 0 ]; then
          exit 1
        fi
      fi
      continue
    fi

    ncu_prune_timing_rows "$table_id" "$case_id" "$rep"
    echo -e "${table_id}\t${case_id}\t${rep}\t${mode}\t${order}\t${visc_order}\t${keep_prec}\t${visc_prec}\t${NX}\t${NT}\t${press_prec}\t${sr}\t${su}\t${sprec}" >> "$OUT/commands.tsv"
    echo "[$(date +%H:%M:%S)] ${table_id}/${case_id} repeat ${rep}/${REPEAT}" | tee -a "$OUT/progress.log"

    set +e
    line="$(BENCH_REPORT_MODE="$mode" BENCH_KERNEL="$bench_kernel" \
      BENCH_SCHEME="$bench_scheme" BENCH_TVD="$bench_tvd" BENCH_RECON="$bench_recon" \
      BENCH_WENO_ORDER="${BENCH_WENO_ORDER:-5}" \
      NCU_EXPORT="$timing_report" NCU_RAW_CSV="$raw_csv" "$BENCH" \
      "$bench_mode" "$order" "$visc_order" "$keep_prec" "$visc_prec" "$NX" "$NT" \
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
      if ncu_should_collect_full "$rep"; then
        ncu_collect_full_report "$table_id" "$case_id" "$mode" "$rep" "$stem"
      fi
    fi
  done
}

ncu_visc_for_keep_precision() {
  # KEEP microbenchmarks focus on whether the convective work can coexist with
  # FP32-side work. Keep the viscous side on FP32 for fp64/fp32/df KEEP cases.
  echo fp32
}

ncu_dump_kernel_sass() {
  local mode="$1" out_prefix="$2" kernel raw norm
  kernel="$(ncu_kernel_for_mode "$mode")"
  raw="${out_prefix}.sass"
  norm="${out_prefix}.norm.sass"

  cuobjdump -sass "$SOLVER_ROOT/ST/build/a.out" > "$raw" 2>/dev/null || return 1
  awk -v k="$kernel" '
    /^        Function :/ {
      on = index($0, k) > 0
    }
    on {
      print
    }
  ' "$raw" |
    sed -E \
      -e '/Function :/d' \
      -e '/\.headerflags/d' \
      -e '/code for/d' \
      -e 's@/\*[0-9a-fA-Fx ]+\*/@@g' \
      -e 's/[[:space:]]+/ /g' \
      -e 's/^ //' \
      -e '/^$/d' > "$norm"
  [ -s "$norm" ]
}

ncu_sass_same_seq_fused() {
  local order="$1" visc_order="$2" keep_prec="$3" visc_prec="$4" press_prec="$5"
  local sr=fp64 su=fp64 sprec=fp64 stem seq_prefix fused_prefix
  stem="$(ncu_sanitize "sass_o${order}_${keep_prec}_${visc_prec}_${press_prec}")"
  seq_prefix="$OUT/sass/${stem}__seq"
  fused_prefix="$OUT/sass/${stem}__fused"

  ncu_configure_case seq "$order" "$visc_order" "$keep_prec" "$visc_prec" "$press_prec" "$sr" "$su" "$sprec"
  ncu_build_case "$OUT/sass/${stem}__seq_build.log" || return 1
  ncu_dump_kernel_sass seq "$seq_prefix" || return 1

  ncu_configure_case fused "$order" "$visc_order" "$keep_prec" "$visc_prec" "$press_prec" "$sr" "$su" "$sprec"
  ncu_build_case "$OUT/sass/${stem}__fused_build.log" || return 1
  ncu_dump_kernel_sass fused "$fused_prefix" || return 1

  if cmp -s "${seq_prefix}.norm.sass" "${fused_prefix}.norm.sass"; then
    echo -e "order\tvisc_order\tkeep_prec\tvisc_prec\tpress_prec\tstatus\tseq_sass\tfused_sass" > "$OUT/sass/${stem}.tsv"
    echo -e "${order}\t${visc_order}\t${keep_prec}\t${visc_prec}\t${press_prec}\tSAME\t${seq_prefix}.norm.sass\t${fused_prefix}.norm.sass" >> "$OUT/sass/${stem}.tsv"
    return 0
  fi

  echo -e "order\tvisc_order\tkeep_prec\tvisc_prec\tpress_prec\tstatus\tseq_sass\tfused_sass" > "$OUT/sass/${stem}.tsv"
  echo -e "${order}\t${visc_order}\t${keep_prec}\t${visc_prec}\t${press_prec}\tDIFFER\t${seq_prefix}.norm.sass\t${fused_prefix}.norm.sass" >> "$OUT/sass/${stem}.tsv"
  return 1
}

ncu_note_skipped_fused() {
  local table_id="$1" case_id="$2" order="$3" visc_order="$4" keep_prec="$5" visc_prec="$6" press_prec="$7"
  if [ ! -s "$OUT/fused_skips.tsv" ]; then
    echo -e "table_id\tcase_id\tmode\torder\tvisc_order\tkeep_prec\tvisc_prec\tpress_prec\treason" > "$OUT/fused_skips.tsv"
  fi
  echo -e "${table_id}\t${case_id}\tfused\t${order}\t${visc_order}\t${keep_prec}\t${visc_prec}\t${press_prec}\tSASS_SAME_AS_SEQ" >> "$OUT/fused_skips.tsv"
}

ncu_finish() {
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
  echo "  $OUT/sass/"
  if [ -s "$OUT/fused_skips.tsv" ]; then
    echo "  $OUT/fused_skips.tsv"
  fi
}
