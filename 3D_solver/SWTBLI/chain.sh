#!/bin/bash
# Run several SWTBLI stages back to back, stopping at the first failure.
#   nohup ./chain.sh G0 C1 "C2 1" "C2 2" > chain.log 2>&1 &
cd "$(dirname "$0")"
for s in "$@"; do
  ./run_stage.sh $s || { echo "CHAIN STOPPED at $s $(date '+%F %T')"; exit 1; }
done
echo "CHAIN DONE $(date '+%F %T')"
