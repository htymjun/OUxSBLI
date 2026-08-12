#!/bin/bash
#$ -cwd
#$ -l node_h=1
#$ -l h_rt=0:05:00
#$ -p -5

module load nvhpc/25.1_cuda12.6 openmpi/5.0.7-nvhpc

ncu --set basic --replay-mode kernel --target-processes all \
-o my_report_tsubame -f \
mpiexec -npernode 4 -n 4 -x LD_LIBRARY_PATH ./a.out

