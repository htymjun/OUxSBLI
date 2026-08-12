#!/bin/bash
mkdir -p data
nohup mpiexec -n 2 ./build/a.out &
cp mod_globals.f90 ./data
cp set.f90 ./data
