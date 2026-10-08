#!/bin/bash
mkdir -p data
# snapshot the exact inputs alongside the run output (provenance), before launch
cp ./set.f90 ./mod_globals.f90 ./config.fypp ./data/
# run_retry.sh restarts the run if nvfortran's -Mchkptr aborts it with its spurious
# "Null pointer for" message (see tools/run_retry.sh); output still goes to nohup.out
nohup bash "$(dirname "$0")/../../tools/run_retry.sh" mpirun 2 ./build/a.out &
