#!/bin/bash
# Run the solver under MPI and restart it if it aborts with the spurious
# "Null pointer for <array>" message of nvfortran's -Mchkptr.
#
#   run_retry.sh <launcher> <nranks> <executable>      (run from the case directory)
#   e.g.  nohup bash ../../tools/run_retry.sh mpiexec 2 ./build/a.out &
#
# Why: the shared CMake flags build with -Mchkptr. Its runtime check (pgf90_ptrchk)
# tests only the low 32 bits of a 64-bit pointer, so a device array that happens
# to be allocated at an address that is a multiple of 4 GiB is reported as a null
# pointer although it is valid. Which array hits this, and whether any does, changes
# from run to run, and it can only happen at the first use of each array, i.e.
# before the time loop. Re-running therefore gives the same result as a clean run.
#
# A retry happens only for that message (any other failure is returned as is), at
# most MAXTRY times (default 10). Output the failed attempt appended to data/ or
# recal/ (the entropy / kinetic-energy logs are opened with position="append") is
# undone first: files that existed before the attempt are truncated back to their
# previous size, files it created are removed. Files outside data/ and recal/ are
# not touched.
#
# Standard output carries the solver output of every attempt (use nohup as usual);
# retry.log records each retry.
MAXTRY=${MAXTRY:-10}
if [ $# -ne 3 ]; then
  echo "usage: $0 <launcher> <nranks> <executable>" >&2
  exit 2
fi
LAUNCHER=$1; NP=$2; EXE=$3

mkdir -p data recal
TMP=.run_retry.$$
trap 'rm -f "$TMP".before "$TMP".log' EXIT

snapshot() {   # "path size" for every file that exists now
  find data recal -type f -printf '%p %s\n' 2>/dev/null | sort
}
restore() {    # undo what the failed attempt wrote, using the snapshot taken before it
  local f sz
  while read -r f sz; do
    [ -f "$f" ] && [ "$(stat -c %s "$f")" -gt "$sz" ] && truncate -s "$sz" "$f"
  done < "$TMP".before
  find data recal -type f 2>/dev/null | sort | while read -r f; do
    grep -q "^$f " "$TMP".before || rm -f "$f"   # did not exist before the attempt
  done
}

for try in $(seq 1 "$MAXTRY"); do
  snapshot > "$TMP".before
  "$LAUNCHER" -n "$NP" "$EXE" 2>&1 | tee "$TMP".log
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ] && grep -q "Null pointer for" "$TMP".log; then
    echo "try $try: spurious -Mchkptr abort, restoring data/ and retrying" | tee -a retry.log
    restore
    continue
  fi
  exit "$rc"
done
echo "giving up after $MAXTRY spurious -Mchkptr aborts" | tee -a retry.log
exit 1
