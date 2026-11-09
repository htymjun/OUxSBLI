#!/bin/bash
# Switch the TBL case between its two run stages and rebuild.
#
#   ./stage.sh A     spin-up   : fresh start, 60k steps  = 12.1 flow-throughs, 11 snapshots
#   ./stage.sh B     sampling  : restart,    120k steps  = 24.2 flow-throughs, 101 snapshots
#
# Stage A brings a laminar-plus-noise initial condition through transition and
# leaves a statistically stationary layer in recal/. Stage B restarts from it
# and writes the dense snapshot series the statistics are computed from; only
# stage B output should be averaged.
#
# The two stages MUST use an identical grid -- recal/Q00001.dat is a raw dump
# with no metadata, so changing blt, Lx/Ly/Lz, nx/ny/nz or the stretching
# constant s between them silently corrupts the restart. Only endT, np,
# step_offset and RESTART may differ.
set -euo pipefail
cd "$(dirname "$0")"

STAGE="${1:-}"
case "$STAGE" in
  A) ENDT='3.6d-4'; NP=10;  OFFSET=0;  RESTART=False ;;
  B) ENDT='7.2d-4'; NP=100; OFFSET=20; RESTART=True  ;;
  *) echo "usage: $0 {A|B}" >&2; exit 2 ;;
esac

sed -i -E "s/^( *real\(8\), parameter :: endT  = ).*/\1${ENDT}/"            mod_globals.f90
sed -i -E "s/^( *integer, parameter :: np    = ).*/\1${NP}/"                mod_globals.f90
sed -i -E "s/^( *integer, parameter :: step_offset   = ).*/\1${OFFSET}/"    mod_globals.f90
sed -i -E "s/^(#:set RESTART = ).*/\1${RESTART}/"                           config.fypp

if [ "$STAGE" = "B" ] && [ ! -f recal/Q00001.dat ]; then
  echo "error: stage B needs recal/Q00001.dat and recal/Qm.dat from stage A" >&2
  exit 1
fi

echo "stage $STAGE: endT=$ENDT np=$NP step_offset=$OFFSET RESTART=$RESTART"
grep -E 'endT|np    =|step_offset' mod_globals.f90 | grep parameter
grep RESTART config.fypp

# CMake's Fortran search does not know about nvfortran, so without this it
# silently picks gfortran and then fails at find_package(MPI). The case
# CMakeLists sets CMAKE_Fortran_COMPILER only if it is not already defined, and
# project() in this directory has already run by then.
export FC="${FC:-nvfortran}"

cmake -B build >/dev/null && cmake --build build -j
echo "built build/a.out for stage $STAGE"
