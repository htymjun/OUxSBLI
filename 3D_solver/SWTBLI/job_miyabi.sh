#!/bin/bash
#PBS -q regular-g
#PBS -l select=1:mpiprocs=4
#PBS -l walltime=24:00:00
#PBS -W group_list=gv82
#PBS -j oe
#
# SWTBLI on Miyabi (one GH200 per node): 4 ranks, both domains share the GPU.
#   qsub -v STAGE=G0 job_miyabi.sh            (T2, G0, C1)
#   qsub -v STAGE=C2,SEG=2 job_miyabi.sh      (C2 N, P N)
module purge
module load nvidia/25.9 nv-hpcx/25.9
cd ${PBS_O_WORKDIR}
./run_stage.sh "${STAGE:-G0}" "${SEG:-1}"
exit 0
