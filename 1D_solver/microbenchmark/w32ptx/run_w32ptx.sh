#!/usr/bin/env bash
# Sweep driver for the w32ptx PTX ladder. Modelled on ../run_weno_nsys.sh and
# keeping its measurement protocol, because that protocol encodes hard-won
# lessons:
#
#  * ROUND-MAJOR interleaving, not mode-major. GPU clocks cannot be locked
#    without root on these boxes, so a mode-major sweep charges the last modes
#    for the thermal drift accumulated by the first. Compare time_min_us across
#    interleaved rounds.
#  * --nrepeat is REJECTED. It does not amplify: measured here, NVVM hoists the
#    loop-invariant face body out of the repeat loop and leaves only the ~6 DADD
#    accumulator, so the slope is ~3.5 us/iteration for every mode regardless of
#    its arithmetic. The Fortran harness pins nrepeat for the same reason.
#  * --gpu-cc is MANDATORY with no env fallback. A wrong architecture builds
#    cleanly and fails at every launch with cudaErrorInvalidPtx (218).
#  * A non-finite output invalidates the row, and the output buffer is poisoned
#    before the run so a kernel that never launched cannot pass as zeros.
set -uo pipefail

OUT=""; GPU_CC=""; NX=4194304; NLAUNCH=10; ROUNDS=3; MODES_ARG=""; GROUP="all"

usage() {
  cat <<'EOF'
usage: run_w32ptx.sh --out DIR --gpu-cc CC [options]
  --out DIR         output directory (required)
  --gpu-cc CC       compute capability, e.g. 80 or 89 (required)
  --nx N            problem size (default 4194304)
  --nlaunch N       timed launches per row (default 10)
  --repeat N        interleaved rounds (default 3)
  --group G         ref | rcp | abl | all   (default all)
  --modes "A B ..." explicit mode list, overrides --group
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --out) OUT="$2"; shift 2;;
    --gpu-cc) GPU_CC="$2"; shift 2;;
    --nx) NX="$2"; shift 2;;
    --nlaunch) NLAUNCH="$2"; shift 2;;
    --repeat) ROUNDS="$2"; shift 2;;
    --group) GROUP="$2"; shift 2;;
    --modes) MODES_ARG="$2"; shift 2;;
    --nrepeat|--nrepeat-list)
      echo "error: --nrepeat is rejected. It does not amplify arithmetic here --" >&2
      echo "       NVVM hoists the loop-invariant body and leaves only the acc" >&2
      echo "       chain (~3.5 us/iter for every mode), so rows stop being" >&2
      echo "       comparable while looking fine. See the header comment." >&2
      exit 2;;
    -h|--help) usage; exit 0;;
    *) echo "error: unknown option $1" >&2; usage; exit 2;;
  esac
done

[ -n "$OUT" ]    || { echo "error: --out is required" >&2; exit 2; }
[ -n "$GPU_CC" ] || { echo "error: --gpu-cc is required (no default on purpose)" >&2; exit 2; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT/raw"
SUMMARY="$OUT/summary.csv"
PROGRESS="$OUT/progress.log"

log() { echo "$*" | tee -a "$PROGRESS"; }

log "=== build for sm_${GPU_CC} ==="
make -C "$HERE" ARCH="$GPU_CC" > "$OUT/build.log" 2>&1 || { echo "build FAILED, see $OUT/build.log" >&2; exit 1; }
grep -E '[1-9][0-9]* bytes spill' "$OUT/build.log" && { echo "error: register spill -- timings would be invalid" >&2; exit 1; }
EXE="$HERE/w32ptx"

{
  echo "date=$(date -Is)"
  echo "gpu_cc_requested=$GPU_CC"
  echo "nx=$NX nlaunch=$NLAUNCH rounds=$ROUNDS"
  echo "git_head=$(git -C "$HERE" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "--- nvcc ---";      nvcc --version 2>&1
  echo "--- nvidia-smi ---"; nvidia-smi --query-gpu=name,compute_cap,driver_version --format=csv 2>&1
} > "$OUT/env.txt"

ALL_MODES=$("$EXE" --list-kernels | cut -f1)
if [ -n "$MODES_ARG" ]; then
  MODES=$MODES_ARG
else
  case "$GROUP" in
    ref) MODES=$(echo "$ALL_MODES" | grep -E '^w32_poly64_.*_ref$');;
    rcp) MODES=$(echo "$ALL_MODES" | grep -E '^w32_poly64_.*_(ref|rcp|divfull|rcpn)$');;
    abl) MODES=$(echo "$ALL_MODES" | grep -E '^(w32_only|weno64)');;
    all) MODES=$ALL_MODES;;
    *) echo "error: unknown --group $GROUP" >&2; exit 2;;
  esac
fi
log "modes: $(echo $MODES | wc -w)"

echo "mode,round,nx,nrepeat,nlaunch,time_min_us,time_avg_us,checksum_all,checksum_bits,nonfinite,status" > "$SUMMARY"

field() { sed -n "s/^$2=[[:space:]]*//p" "$1" | head -1; }

for rnd in $(seq 1 "$ROUNDS"); do
  for m in $MODES; do
    raw="$OUT/raw/${m}__r${rnd}.log"
    "$EXE" "$m" "$NX" 1 "$NLAUNCH" > "$raw" 2>&1
    rc=$?
    if [ $rc -ne 0 ]; then
      echo "$m,$rnd,$NX,1,$NLAUNCH,,,,,,RUN_FAIL_rc$rc" >> "$SUMMARY"
      log "round $rnd  $m  RUN_FAIL rc=$rc"
      continue
    fi
    tmin=$(field "$raw" time_min_us); tavg=$(field "$raw" time_avg_us)
    csum=$(field "$raw" checksum_all); bits=$(field "$raw" checksum_bits)
    nf=$(field "$raw" nonfinite)
    st=OK; [ "${nf:-1}" = "0" ] || st=NONFINITE
    echo "$m,$rnd,$NX,1,$NLAUNCH,$tmin,$tavg,$csum,$bits,$nf,$st" >> "$SUMMARY"
    log "round $rnd  $(printf '%-30s' "$m") ${tmin} us  $st"
  done
done

log ""
log "=== ranked (best time_min_us across rounds) ==="
python3 - "$SUMMARY" <<'PY' | tee -a "$PROGRESS"
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1])))
best, bits = {}, {}
for r in rows:
    if r["status"] != "OK":
        continue
    t = float(r["time_min_us"])
    m = r["mode"]
    if m not in best or t < best[m]:
        best[m] = t
    bits.setdefault(m, set()).add(r["checksum_bits"])

for m, s in bits.items():
    if len(s) > 1:
        print(f"WARNING: {m} produced differing checksum_bits across rounds: {s}")

print(f"{'mode':<32}{'time_min_us':>12}{'vs _ref':>10}")
for m, t in sorted(best.items(), key=lambda kv: kv[1]):
    ref = m.rsplit("_", 1)[0] + "_ref"
    sp = f"{best[ref]/t:.3f}x" if ref in best and ref != m else "-"
    print(f"{m:<32}{t:>12.3f}{sp:>10}")

# Marginal cost of the FP64 half: t(w32_poly64) - t(w32_only) at the same rung.
print("\n=== marginal cost of the FP64 polynomial half ===")
print(f"{'width':>5} {'rung':<8}{'w32_poly64':>12}{'w32_only':>10}{'marginal':>10}")
for w in (5, 7, 9):
    for rung in ("ref", "rcp"):
        a, b = f"w32_poly64_seq{w}_{rung}", f"w32_only_seq{w}_{rung}"
        if a in best and b in best:
            print(f"{w:>5} {rung:<8}{best[a]:>12.3f}{best[b]:>10.3f}{best[a]-best[b]:>10.3f}")
PY

log ""
log "wrote $SUMMARY"
