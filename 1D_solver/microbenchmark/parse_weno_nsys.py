#!/usr/bin/env python3
import argparse
import csv
import math
import re
import statistics
from pathlib import Path


def read_rows(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(newline="") as fh:
        return [r for r in csv.reader(fh) if r]


def find_header(rows, required):
    for i, row in enumerate(rows):
        low = [c.strip().lower() for c in row]
        if all(any(req in c for c in low) for req in required):
            return i
    return None


def col(header, names):
    low = [c.strip().lower() for c in header]
    for name in names:
        for i, c in enumerate(low):
            if name in c:
                return i
    return None


def number(x):
    try:
        v = float(str(x).replace(",", "").strip())
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def scale_us(name):
    low = name.lower()
    if "(ns)" in low or low.endswith(" ns"):
        return 1.0e-3
    if "(us)" in low or low.endswith(" us"):
        return 1.0
    if "(ms)" in low or low.endswith(" ms"):
        return 1.0e3
    if "(s)" in low or low.endswith(" s"):
        return 1.0e6
    return 1.0e-3


def trace_values(path, pattern):
    rows = read_rows(path)
    hdr = find_header(rows, ["duration", "name"])
    if hdr is None:
        return []
    header = rows[hdr]
    name_col = col(header, ["name", "kernel"])
    dur_col = col(header, ["duration"])
    if name_col is None or dur_col is None:
        return []
    rx = re.compile(pattern, re.I)
    scale = scale_us(header[dur_col])
    vals = []
    for row in rows[hdr + 1:]:
        if len(row) <= max(name_col, dur_col):
            continue
        if not rx.search(row[name_col]):
            continue
        v = number(row[dur_col])
        if v is not None:
            vals.append(v * scale)
    return vals


def fail(args, status):
    print(",".join([
        args.mode, args.repeat, args.nx, args.nrepeat, args.kernel,
        "0", "0", "NSYS_FAIL", "", "", "", "",
        args.trace_csv, args.kern_sum_csv, args.report, args.sqlite, status,
    ]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True)
    ap.add_argument("--repeat", required=True)
    ap.add_argument("--nx", required=True)
    ap.add_argument("--nrepeat", required=True)
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--trace-csv", default="")
    ap.add_argument("--kern-sum-csv", default="")
    ap.add_argument("--report", default="")
    ap.add_argument("--sqlite", default="")
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--count", type=int, default=0)
    args = ap.parse_args()

    vals = trace_values(args.trace_csv, args.kernel)
    total = len(vals)
    used = vals[args.skip:]
    if args.count > 0:
        used = used[:args.count]
    if not used:
        fail(args, "NO_KERNEL")
        return
    avg = statistics.mean(used)
    mn = min(used)
    mx = max(used)
    print(",".join([
        args.mode, args.repeat, args.nx, args.nrepeat, args.kernel,
        str(total), str(len(used)), f"{avg:.3f}", f"{mn:.3f}",
        f"{mx:.3f}", f"{mx-mn:.3f}", f"{sum(used):.3f}",
        args.trace_csv, args.kern_sum_csv, args.report, args.sqlite, "trace",
    ]))


if __name__ == "__main__":
    main()
