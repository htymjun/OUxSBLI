#!/bin/bash
# Switch SWTBLI between its run stages and rebuild. delta/u0 = 1.9678e-6 s, dt = 1.08e-8 s.
#
# Grid generation 2 (current mod_globals.f90: dx0 = 0.0403 delta, 18 delta high,
# ximp_in = 31, 4 delta span). Initial field: regrid.py applied to the end of B0.
#   ./stage.sh T2         smoke test, 100 steps from recal_gen2_init, theta = 9.5
#   ./stage.sh G0         theta = 9.5, 50 delta/U: flush the interpolation transient, VTK every 5 delta/U
#   ./stage.sh C1         theta = 11,  150 delta/U, VTK every 5 delta/U
#   ./stage.sh C2 N       theta = 12,  segment N (>=1) of 300 delta/U, VTK every 10 delta/U,
#                         until x_sh and x_0 pass the stationarity test (track.json)
#   ./stage.sh P N        theta = 12,  statistics segment N of 75 delta/U, VTK every 2.5 delta/U
#
# Grid generation 1 (A1, A1b, A2a, A2b, A2c N, B0; dx0 = 0.0604 delta, 12.2 delta
# high, ximp_in = 25) cannot be rerun with the current mod_globals.f90: its
# parameters are archived in runs/<stage>/mod_globals.f90.
#
# endT is chosen slightly above nt*np*dt so that int(endT/(np*dt)) cannot round
# down. Only endT, np, step_offset, theta_deg, the span and RESTART change between
# stages; everything else (grid, blt, dt) must stay fixed or the raw restart files
# in recal/ no longer match.
set -euo pipefail
cd "$(dirname "$0")"

STAGE="${1:-}"
SEG="${2:-1}"
NZ=114; LZ='4.d0'; RESTART=True
case "$STAGE" in
  T2) ENDT='1.0802d-6'; NP=2;  OFF=0;    TH='9.5d0' ;;
  G0) ENDT='9.839d-5';  NP=10; OFF=6000; TH='9.5d0' ;;
  C1) ENDT='2.9517d-4'; NP=30; OFF=6100; TH='11.d0' ;;
  C2) ENDT='5.9034d-4'; NP=30; OFF=$((6200 + 100 * (SEG - 1))); TH='12.d0' ;;
  P)  ENDT='1.4743d-4'; NP=30; OFF=$((8000 + 100 * (SEG - 1))); TH='12.d0' ;;
  A1|A1b|A2a|A2b|A2c|B0|T)
    echo "$STAGE is a grid-generation-1 stage; its parameters are in runs/$STAGE*/mod_globals.f90" >&2; exit 2 ;;
  *) echo "usage: $0 {T2|G0|C1|C2 N|P N}" >&2; exit 2 ;;
esac

sed -i -E "s/^( *real\(8\), parameter :: endT  = ).*/\1${ENDT}/"            mod_globals.f90
sed -i -E "s/^( *integer, parameter :: np    = ).*/\1${NP}/"                mod_globals.f90
sed -i -E "s/^( *integer, parameter :: step_offset   = ).*/\1${OFF}/"       mod_globals.f90
sed -i -E "s/^( *real\(8\), parameter :: theta_deg = ).*/\1${TH}/"          mod_globals.f90
sed -i -E "s/^( *real\(8\), parameter :: Lz1 = ).*/\1${LZ} * blt/"          mod_globals.f90
sed -i -E "s/^( *integer, parameter :: nz1 = ).*/\1${NZ}/"                  mod_globals.f90
sed -i -E "s/^(#:set RESTART = ).*/\1${RESTART}/"                           config.fypp

echo "stage $STAGE (segment $SEG): endT=$ENDT np=$NP step_offset=$OFF theta=$TH nz1=$NZ Lz1=$LZ*blt RESTART=$RESTART"
grep -E "parameter :: (endT|np |step_offset|theta_deg|Lz1|nz1)" mod_globals.f90

export FC="${FC:-nvfortran}"
cmake -B build >/dev/null && cmake --build build -j 16 >/dev/null
echo "built build/a.out for stage $STAGE"
