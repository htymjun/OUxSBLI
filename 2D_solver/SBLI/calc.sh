#!/bin/bash
mkdir -p data
# snapshot the exact inputs alongside the run output (provenance), before launch
cp ./set.f90 ./mod_globals.f90 ./config.fypp ./data/
nohup mpirun -n 2 ./build/a.out &
