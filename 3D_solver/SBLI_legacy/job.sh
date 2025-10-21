#!/bin/bash
#PJM -L rscgrp=debug-a
#PJM -L gpu=8
#PJM --mpi proc=16
#PJM -L elapse=00:30:00
#PJM -g gi81

module load nvidia/22.7 cuda/11.4 ompi-cuda

mpiexec -machinefile $PJM_O_NODEINF -n \
$PJM_MPI_PROC -npernode 16 ./a.out

