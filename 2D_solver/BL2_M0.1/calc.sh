#!/bin/bash
mkdir -p data
nohup mpirun -n 2 ./build/a.out &
cp ./set.f90 ./data
cp ./mod_globals.f90 ./data
