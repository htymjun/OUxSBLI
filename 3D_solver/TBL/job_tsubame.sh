#!/bin/bash
#$ -cwd
#$ -l node_h=1
#$ -l h_rt=24:00:00
#$ -N tbl
#$ -j y
#
# TBL (M=2.5 supersonic turbulent boundary layer) on TSUBAME.
#
# Two MPI ranks, ONE GPU: main.f90 pairs the ranks, so the even rank computes
# on the GPU and the odd rank runs the recycling-rescaling on the host and
# writes the VTK snapshots. Do not raise this to 4 the way SWTBLI does --
# rescaling only supports the single rerank pair.
#
# Submit as:
#   qsub -v STAGE=A job_tsubame.sh     # spin-up,  60k steps
#   qsub -v STAGE=B job_tsubame.sh     # sampling, 120k steps  (needs recal/ from A)
#
# On the RTX 4060 Laptop this was 4.8 steps/s, i.e. 3.5 h for stage A and 7 h
# for stage B. Scale h_rt by your GPU; the grid is 6.4 M points and needs
# ~2 GB of device memory.
set -euo pipefail

module load nvhpc/25.1_cuda12.6 openmpi/5.0.7-nvhpc

STAGE="${STAGE:-A}"
bash ./stage.sh "${STAGE}"

mkdir -p data recal
mpiexec -npernode 2 -n 2 -x LD_LIBRARY_PATH ./build/a.out

# Reduce 13 GB of snapshots to a ~13 MB moment file here, before copying
# anything back. Needs only numpy -- see ouxsbli/analysis/vtr_raw.py.
if [ "${STAGE}" = "B" ]; then
  python3 -m ouxsbli.analysis.tbl_stats accumulate data \
      --out tbl_acc.npz --halves || echo "post-processing skipped (no numpy?)"
fi
