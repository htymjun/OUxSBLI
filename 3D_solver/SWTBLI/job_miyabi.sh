#!/bin/bash
#PBS -q small-g
#PBS -l select=1:mpiprocs=4
#PBS -l walltime=24:00:00
#PBS -W group_list=gv82
#PBS -j oe

module purge
module load nvidia/25.9 nv-hpcx/25.9

cd ${PBS_O_WORKDIR}
mpiexec ./a.out

exit 0
