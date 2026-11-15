#!/bin/bash
# Drive grid generation 2 to the statistics: after any running chain.sh, continue
# C2 N (theta = 12, 300 delta/U each) until x_sh and x_0 both pass the stationarity
# test of track.json (at most C2_MAX segments), then P 1..P_N (75 delta/U each) and
# the final wall report (P segments + the last, stationary C2 segment, against B0).
#   CHAIN_PID=<pid of chain.sh> nohup ./auto_gen2.sh > auto_gen2.log 2>&1 &
cd "$(dirname "$0")"
C2_MAX="${C2_MAX:-4}"
P_N="${P_N:-4}"
ROOT="$(cd ../.. && pwd)"
A="$(pwd)"
# CHAIN_PID: the chain.sh process to wait for (pgrep -f would also match shells
# whose command line merely mentions chain.sh)
if [ -n "${CHAIN_PID:-}" ]; then while kill -0 "$CHAIN_PID" 2>/dev/null; do sleep 60; done; fi
grep -q "CHAIN DONE" chain.log || { echo "chain.sh did not finish cleanly; stopping"; exit 1; }

stationary() {
  python3 -c "import json,sys; t=json.load(open('runs/$1/analysis/track.json')); \
sys.exit(0 if t['x_sh_stationarity']['stationary'] and t['x_0_stationarity']['stationary'] else 1)"
}
N=$(grep -c '^C2_' runs/STAGES_GEN2)
last="C2_$N"
until stationary "$last"; do
  N=$((N + 1))
  [ "$N" -gt "$C2_MAX" ] && { echo "not stationary after $C2_MAX C2 segments; stopping before P $(date '+%F %T')"; exit 1; }
  ./run_stage.sh C2 "$N" || { echo "C2 $N failed"; exit 1; }
  last="C2_$N"
done
echo "stationary at $last $(date '+%F %T')"
for k in $(seq 1 "$P_N"); do
  ./run_stage.sh P "$k" || { echo "P $k failed"; exit 1; }
done
npz="$A/runs/$last/analysis/wall_acc.npz"
for k in $(seq 1 "$P_N"); do npz="$npz $A/runs/P_$k/analysis/wall_acc.npz"; done
mkdir -p runs/final
(cd "$ROOT" && python3 -m ouxsbli.analysis.swtbli_stats wall report $npz \
   --compare "$A/runs/B0/analysis/wall_acc.npz" --out "$A/runs/final") > runs/final/wall.log 2>&1
echo "AUTO DONE $(date '+%F %T')"
