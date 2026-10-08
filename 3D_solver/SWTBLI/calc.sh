#!/bin/bash
mkdir -p data
cp mod_globals.f90 ./data
cp set.f90 ./data
# run_retry.sh restarts the run if nvfortran's -Mchkptr aborts it with its spurious
# "Null pointer for" message (see tools/run_retry.sh); output still goes to nohup.out
nohup bash "$(dirname "$0")/../../tools/run_retry.sh" mpiexec 4 ./build/a.out &
