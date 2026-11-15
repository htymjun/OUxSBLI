#!/bin/bash
# Analyse one finished grid-generation-2 SWTBLI stage and reclaim disk space.
#   ./post_stage.sh TAG          (called by run_stage.sh; TAG = G0, C1, C2_N, P_N)
#
# 1. swtbli_stats: turb, track (continuing the previous stage's track.npz), report,
#    wall accumulate (wall_acc.npz: merge the P_N files with `wall report`), gif.
# 2. data/: only the last KEEP snapshots of each domain stay (an SBLI snapshot is
#    ~1.5 GB on this grid); everything else needed later is in analysis/ and gif/.
# 3. recal/: kept for this stage and the previous one only (3 GB each).
# The stage order is recorded in runs/STAGES_GEN2 (the stage before G0 is B0).
set -uo pipefail
cd "$(dirname "$0")"
TAG="$1"
KEEP="${KEEP:-3}"
RD="runs/$TAG"
ROOT="$(cd ../.. && pwd)"
REF="runs/A1b/analysis/ref_rows.npz"
LIST="runs/STAGES_GEN2"
PREV="B0"; [ -s "$LIST" ] && PREV="$(tail -n 1 "$LIST")"
[ "$PREV" = "$TAG" ] && PREV="$( (grep -vx "$TAG" "$LIST" || true) | tail -n 1)"
[ -n "$PREV" ] || PREV="B0"
st() { (cd "$ROOT" && python3 -m ouxsbli.analysis.swtbli_stats "$@"); }
A="$(pwd)"

echo "POST $TAG (previous stage $PREV) $(date '+%F %T')"
st turb   "$A/$RD"                                   > "$RD/turb.log"   2>&1 || echo "  turb failed (see $RD/turb.log)"
st track  "$A/$RD" --ref "$A/$REF" --prev "$A/runs/$PREV/analysis/track.npz" \
                                                     > "$RD/track.log"  2>&1 || echo "  track failed (see $RD/track.log)"
st report "$A/$RD" --ref "$A/$REF"                   > "$RD/report.log" 2>&1 || echo "  report failed (see $RD/report.log)"
st wall accumulate "$A/$RD"                          > "$RD/wall.log"   2>&1 || echo "  wall accumulate failed (see $RD/wall.log)"
st gif    "$A/$RD" --track "$A/$RD/analysis/track.npz" > "$RD/gif.log"  2>&1 || echo "  gif failed (see $RD/gif.log)"

if [ -f "$RD/analysis/report.json" ] && [ -f "$RD/analysis/track.npz" ] && [ -f "$RD/analysis/wall_acc.npz" ]; then
  for d in "$RD/data/1" "$RD/data/3"; do
    ls "$d"/Q*.vtr 2>/dev/null | sort | head -n -"$KEEP" | xargs -r rm -f
  done
  echo "  data/ thinned to the last $KEEP snapshots"
else
  echo "  analysis incomplete: data/ left untouched"
fi

grep -qx "$TAG" "$LIST" 2>/dev/null || echo "$TAG" >> "$LIST"
# recal/ of every generation-2 stage except the last two
head -n -2 "$LIST" | while read -r old; do
  [ -d "runs/$old/recal" ] && rm -rf "runs/$old/recal" && echo "  removed runs/$old/recal"
done
du -sh "$RD" runs | sed 's/^/  /'
echo "POST DONE $TAG $(date '+%F %T')"
