#!/bin/bash
#$ -cwd
#$ -l node_h=1
#$ -l h_rt=0:10:00
#$ -p -5

module load nvhpc/25.1_cuda12.6 openmpi/5.0.7-nvhpc

compute-sanitizer --tool memcheck --report-api-error no mpiexec -npernode 4 -n 4 -x LD_LIBRARY_PATH ./a.out

