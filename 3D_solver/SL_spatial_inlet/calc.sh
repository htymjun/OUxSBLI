#!/bin/bash
mkdir -p data
# x方向MPI分割: mod_globals.f90 の npx(xスラブ数)に合わせて 2*npx ランクで起動する
npx=$(sed -n 's/^ *integer, parameter :: npx *= *\([0-9][0-9]*\).*/\1/p' mod_globals.f90)
nohup mpiexec -n $((2*npx)) ./build/a.out &
cp mod_globals.f90 ./data
cp set.f90 ./data
