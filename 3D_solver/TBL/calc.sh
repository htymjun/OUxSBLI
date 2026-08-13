#!/bin/bash
# Local (interactive) run of one TBL stage. For a batch system use
# job_miyabi.sh or job_tsubame.sh instead.
#
#   ./calc.sh A     spin-up,  60k steps -- fresh start
#   ./calc.sh B     sampling, 120k steps -- restarts from recal/
#
# Needs the NVIDIA HPC SDK on PATH, including its own mpif90: the system
# /usr/bin/mpif90 is a gfortran wrapper and find_package(MPI) will reject it.
#   export NVROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7
#   export PATH=$NVROOT/comm_libs/openmpi/openmpi-3.1.5/bin:$NVROOT/compilers/bin:$PATH
set -euo pipefail
cd "$(dirname "$0")"

STAGE="${1:-A}"
bash ./stage.sh "$STAGE"

mkdir -p data recal
cp mod_globals.f90 set.f90 config.fypp ./data/
nohup mpiexec -n 2 ./build/a.out > "stage${STAGE}.log" 2>&1 &
echo "stage $STAGE running, pid $!  --  tail -f stage${STAGE}.log"
