#!/bin/bash
#PJM -L rscgrp=tutorial1-a
#PJM -L node=1
#PJM --mpi proc=16
#PJM -L elapse=00:10:00
#PJM -g gt01

module load nvidia/22.7 cuda/11.4 ompi-cuda

export UCX_IB_GPU_DIRECT_RDMA=y
export UCX_PROTO_ENABLE=y

mpiexec -machinefile $PJM_O_NODEINF -n \
$PJM_MPI_PROC -npernode 16 ./a.out

