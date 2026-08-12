#!/bin/bash
#PBS -q regular-g
#PBS -l select=1:mpiprocs=2
#PBS -l walltime=00:10:00
#PBS -W group_list=gv82
#PBS -j oe

module purge
module load nvidia/25.9 nv-hpcx/25.9

cd ${PBS_O_WORKDIR}
ncu --set full --import-source yes --target-processes all -o my_report_ncu mpiexec ./build/a.out

exit 0
