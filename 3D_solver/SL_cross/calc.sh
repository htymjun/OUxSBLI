#!/bin/bash
set -eu
cd -- "$(dirname -- "$0")"
mkdir -p data recal
cp mod_globals.f90 set.f90 config.fypp data/
nohup mpiexec -n 2 ./build/a.out > nohup.out 2>&1 &
