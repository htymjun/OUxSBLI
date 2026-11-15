#!/bin/bash
# Build, run, analyse and archive one SWTBLI stage on the local workstation (4 ranks;
# with two visible GPUs the TBL pair runs on GPU0 and the SBLI pair on GPU1).
#   ./run_stage.sh T2 | G0 | C1 | C2 N | P N
# Output: runs/<stage>[_N]/{data,recal,analysis,gif,run.log}; recal/ stays in place
# for the next stage. T2 (smoke test) and G0 start from recal_gen2_init/ (regrid.py).
# After a production stage post_stage.sh runs the analysis and thins data/.
set -uo pipefail
cd "$(dirname "$0")"
STAGE="$1"; SEG="${2:-1}"
TAG="$STAGE"; case "$STAGE" in C2|P) TAG="${STAGE}_$SEG" ;; esac
bash ./stage.sh "$STAGE" "$SEG" || { echo "STAGE BUILD FAILED $TAG"; exit 1; }
if [ "$STAGE" = T2 ] || [ "$STAGE" = G0 ]; then
  [ -f recal_gen2_init/Q00002.dat ] || { echo "recal_gen2_init/ missing: run regrid.py runs/B0 mod_globals.f90 recal_gen2_init"; exit 1; }
  rm -rf recal; cp -r recal_gen2_init recal; rm -f recal/regrid_source_mod_globals.f90
fi
rm -rf data; mkdir -p data/1 data/3 recal "runs/$TAG"
cp mod_globals.f90 set.f90 config.fypp "runs/$TAG/"
echo "RUN START $TAG $(date '+%F %T')"
mpiexec -n 4 ./build/a.out > "runs/$TAG/run.log" 2>&1
rc=$?
if [ $rc -ne 0 ] || grep -qiE "nan|abort|error" "runs/$TAG/run.log"; then
  echo "RUN FAILED $TAG rc=$rc $(date '+%F %T')"; exit 1
fi
rm -rf "runs/$TAG/data" "runs/$TAG/recal"
mv data "runs/$TAG/data"
if [ "$STAGE" = T2 ]; then
  echo "RUN DONE $TAG (smoke test, not archived as a stage) $(date '+%F %T')"; exit 0
fi
cp -r recal "runs/$TAG/recal"
echo "RUN DONE $TAG $(date '+%F %T')"
bash ./post_stage.sh "$TAG"
