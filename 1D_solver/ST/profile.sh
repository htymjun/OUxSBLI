#!/usr/bin/env bash
set -euo pipefail

timestamp="$(date +"%Y%m%d_%H%M%S")"
MPIEXEC="${MPIEXEC:-${MPIRUN:-mpirun}}"
MPIEXEC_NP_FLAG="${MPIEXEC_NP_FLAG:--n}"
PROFILE_RANKS="${PROFILE_RANKS:-1}"
PROFILE_LAUNCH_MODE="${PROFILE_LAUNCH_MODE:-auto}"

build_launch_cmd() {
  local mode="$1"
  if [ "$mode" = direct ]; then
    PROFILE_TARGET_PROCESSES=application-only
    PROFILE_CMD=(./build/a.out)
    return 0
  fi

  PROFILE_TARGET_PROCESSES=all
  PROFILE_CMD=("$MPIEXEC" "$MPIEXEC_NP_FLAG" "$PROFILE_RANKS" ./build/a.out)
}

case "$PROFILE_LAUNCH_MODE" in
  auto)
    if [ "$PROFILE_RANKS" = 1 ]; then
      build_launch_cmd direct
    else
      build_launch_cmd mpi
    fi
    ;;
  direct|mpi)
    build_launch_cmd "$PROFILE_LAUNCH_MODE"
    ;;
  *)
    echo "bad PROFILE_LAUNCH_MODE: $PROFILE_LAUNCH_MODE (use auto, direct, or mpi)" >&2
    exit 2
    ;;
esac

nsys profile \
    -t cuda,nvtx,openacc,osrt \
    -f true \
    -o "my_report_${timestamp}" \
    "${PROFILE_CMD[@]}" &

ncu \
    --set full \
    --target-processes "$PROFILE_TARGET_PROCESSES" \
    --import-source yes \
    -o "my_report_ncu_${timestamp}" \
    "${PROFILE_CMD[@]}" &
