#!/bin/bash
set -eu
cd -- "$(dirname -- "$0")"
mkdir -p data recal
cp mod_globals.f90 set.f90 config.fypp data/
# x方向MPI分割: mod_globals.f90 の npx(xスラブ数)に合わせて 2*npx ランクで起動する
npx=$(sed -n 's/^ *integer, parameter :: npx *= *\([0-9][0-9]*\).*/\1/p' mod_globals.f90)
nohup mpiexec -n $((2*npx)) ./build/a.out > nohup.out 2>&1 &
