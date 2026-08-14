#!/usr/bin/env python3
"""CSV reducer for report/bench.sh -- not useful on its own.

time_min_us is the primary number to compare. This is a laptop GPU whose clocks
cannot be locked (nvidia-smi -lgc is not permitted for a non-root user), so a
long sweep drifts as the card heats: mean times can move by more than 2x while
the pipe-utilisation ratios stay put. The per-launch minimum is the closest
available estimate of unthrottled time; always cross-check it against the
spread column and against a repeated, interleaved run.

`ncu --csv --page raw` emits a WIDE table: row 0 = metric names, row 1 = units,
rows 2+ = one profiled launch each. Average the numeric columns over the
launches and also report min/max so the spread is visible -- the launch-to-launch
spread on this case is ~1%, which is larger than several of the deltas earlier
sweeps reported as signal.
"""
import csv
import statistics
import sys

MODE, O, VO, KP, VP, PP, NX, REG, SHM, PATH = sys.argv[1:11]

with open(PATH) as fh:
    rows = [r for r in csv.reader(fh) if r]
# The solver's own stdout is interleaved with ncu's CSV, so locate the header
# row by content rather than assuming it is row 0.
hdr = next((i for i, r in enumerate(rows)
            if any(c.startswith("gpu__time_duration") for c in r)), None)
if hdr is None or len(rows) < hdr + 3:
    print(f"{MODE},{O},{VO},{KP},{VP},{PP},{NX},NCU_FAIL,,,,,,,,,{REG},{SHM}")
    sys.exit(0)

names, units, data = rows[hdr], rows[hdr + 1], rows[hdr + 2:]


def col(prefix):
    """Mean over profiled launches of the first column whose name starts with prefix."""
    for j, n in enumerate(names):
        if n.startswith(prefix):
            vals = []
            for r in data:
                try:
                    vals.append(float(r[j].replace(",", "")))
                except (ValueError, IndexError):
                    pass
            if vals:
                return statistics.mean(vals), min(vals), max(vals), units[j]
    return None, None, None, None


t, tmin, tmax, unit = col("gpu__time_duration.sum")
if t is None:
    print(f"{MODE},{O},{VO},{KP},{VP},{PP},{NX},NCU_FAIL,,,,,,,,,{REG},{SHM}")
    sys.exit(0)

scale = {"ns": 1e-3, "nsecond": 1e-3, "us": 1.0, "usecond": 1.0,
         "ms": 1e3, "msecond": 1e3}.get((unit or "ns").strip(), 1e-3)
t_us, tmin_us, spread = t * scale, tmin * scale, (tmax - tmin) * scale


def num(prefix, fmt="{:.1f}"):
    v = col(prefix)[0]
    return "" if v is None else fmt.format(v)


print(",".join([
    MODE, O, VO, KP, VP, PP, NX,
    f"{t_us:.2f}", f"{tmin_us:.2f}", f"{spread:.2f}",
    num("sm__inst_executed_pipe_fp64.sum", "{:.0f}"),
    num("sm__pipe_fp64_cycles_active"),
    num("gpu__dram_throughput"),
    num("sm__throughput"),
    num("sm__warps_active"),
    num("launch__waves_per_multiprocessor", "{:.1f}"),
    REG, SHM,
]))
