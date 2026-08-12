#!/bin/bash
#PBS -q regular-g
#PBS -l select=1:mpiprocs=2
#PBS -l walltime=00:30:00
#PBS -W group_list=gv82
#PBS -j oe

module purge
module load nvidia/25.9 nv-hpcx/25.9

cd ${PBS_O_WORKDIR}

# Name the capture after the dispatch variant actually compiled in, so the
# three OVERLAP_MODE runs don't overwrite each other.
MODE=$(sed -n "s/^#:set OVERLAP_MODE *= *'\([a-z]*\)'.*/\1/p" config.fypp)
OUT="my_report_nsys_${MODE:-unknown}_$(date +%Y%m%d)"

nsys profile -t cuda,nvtx -f true -o "${OUT}" mpiexec ./build/a.out

# The .nsys-rep is written by the cluster's nsys (2025.x) and cannot be opened
# by older local installs. Export sqlite alongside it -- this is what makes the
# per-stream busy/idle timeline analysis possible off-cluster, and it needs no
# re-run to produce.
nsys export --type sqlite --force-overwrite true -o "${OUT}.sqlite" "${OUT}.nsys-rep"

exit 0
