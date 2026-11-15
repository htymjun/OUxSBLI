#!/bin/bash
# Local run of one SWTBLI stage (4 MPI ranks: TBL pair on GPU0, SBLI pair on GPU1
# when two GPUs are visible). See stage.sh for the stages and run_stage.sh for
# what is archived under runs/.
#   ./calc.sh T2 | G0 | C1 | C2 N | P N
set -euo pipefail
cd "$(dirname "$0")"
[ $# -ge 1 ] || { echo "usage: $0 {T2|G0|C1|C2 N|P N}" >&2; exit 2; }
nohup ./run_stage.sh "$@" > "run_${1}${2:+_$2}.log" 2>&1 &
echo "stage $* running in the background, log: run_${1}${2:+_$2}.log"
