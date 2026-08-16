#!/usr/bin/env python3
"""Reduce Nsight Systems CSV exports to one benchmark row.

The nsys-rep file is useful for GUI inspection, but hard to review from an
agent. This reducer consumes nsys stats CSV output, filters the target CUDA
kernel, applies the same warm-up skip/count convention as the ncu scripts, and
prints one CSV row.
"""
import argparse
import csv
import math
import re
import statistics
from pathlib import Path


def read_csv_rows(path):
    if not path:
        return []
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


def number(text):
    try:
        value = float(str(text).replace(",", "").strip())
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def unit_scale_us(name):
    low = name.lower()
    if "(ns)" in low or low.endswith(" ns") or "nsec" in low:
        return 1.0e-3
    if "(us)" in low or low.endswith(" us") or "usec" in low:
        return 1.0
    if "(ms)" in low or low.endswith(" ms") or "msec" in low:
        return 1.0e3
    if "(s)" in low or low.endswith(" s") or "sec" in low:
        return 1.0e6
    # Nsight Systems trace reports normally use ns when the unit is omitted.
    return 1.0e-3


def column_index(header, candidates):
    low = [c.strip().lower() for c in header]
    for candidate in candidates:
        for i, name in enumerate(low):
            if candidate in name:
                return i
    return None


def trace_durations_us(path, kernel_pattern):
    rows = read_csv_rows(path)
    hdr = find_header(rows, ["duration", "name"])
    if hdr is None:
        return []

    header = rows[hdr]
    name_col = column_index(header, ["name", "kernel"])
    dur_col = column_index(header, ["duration"])
    if name_col is None or dur_col is None:
        return []

    scale = unit_scale_us(header[dur_col])
    rx = re.compile(kernel_pattern)
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


def kern_sum_values_us(path, kernel_pattern):
    rows = read_csv_rows(path)
    hdr = find_header(rows, ["name"])
    if hdr is None:
        return None

    header = rows[hdr]
    name_col = column_index(header, ["name", "kernel"])
    total_col = column_index(header, ["total time"])
    avg_col = column_index(header, ["avg"])
    min_col = column_index(header, ["min"])
    max_col = column_index(header, ["max"])
    inst_col = column_index(header, ["instances", "calls", "count"])
    if name_col is None:
        return None

    rx = re.compile(kernel_pattern)
    for row in rows[hdr + 1:]:
        if len(row) <= name_col or not rx.search(row[name_col]):
            continue
        def get(col):
            if col is None or len(row) <= col:
                return None
            v = number(row[col])
            if v is None:
                return None
            return v * unit_scale_us(header[col])

        inst = None
        if inst_col is not None and len(row) > inst_col:
            inst = number(row[inst_col])
        return {
            "instances": int(inst) if inst is not None else 0,
            "total": get(total_col),
            "avg": get(avg_col),
            "min": get(min_col),
            "max": get(max_col),
        }
    return None


def emit_fail(args, status):
    print(",".join([
        args.mode, args.order, args.visc_order, args.keep_prec, args.visc_prec,
        args.press_prec, args.slau_rho_prec, args.slau_u_prec, args.slau_p_prec,
        args.nx, args.kernel, "0", "0",
        "NSYS_FAIL", "", "", "", "", args.trace_csv or "", args.kern_sum_csv or "",
        args.report or "", args.sqlite or "", status,
    ]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True)
    ap.add_argument("--order", required=True)
    ap.add_argument("--visc-order", required=True)
    ap.add_argument("--keep-prec", required=True)
    ap.add_argument("--visc-prec", required=True)
    ap.add_argument("--press-prec", required=True)
    ap.add_argument("--slau-rho-prec", required=True)
    ap.add_argument("--slau-u-prec", required=True)
    ap.add_argument("--slau-p-prec", required=True)
    ap.add_argument("--nx", required=True)
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--skip", type=int, default=20)
    ap.add_argument("--count", type=int, default=10)
    ap.add_argument("--trace-csv", default="")
    ap.add_argument("--kern-sum-csv", default="")
    ap.add_argument("--report", default="")
    ap.add_argument("--sqlite", default="")
    args = ap.parse_args()

    vals = trace_durations_us(args.trace_csv, args.kernel)
    source = "trace"
    total_instances = len(vals)
    used = vals[args.skip:]
    if args.count > 0:
        used = used[:args.count]

    if not used:
        summary = kern_sum_values_us(args.kern_sum_csv, args.kernel)
        if not summary or summary["avg"] is None:
            emit_fail(args, "NO_KERNEL")
            return
        source = "kern_sum"
        total_instances = summary["instances"]
        used_instances = summary["instances"]
        avg = summary["avg"]
        mn = summary["min"] if summary["min"] is not None else avg
        mx = summary["max"] if summary["max"] is not None else avg
        total = summary["total"] if summary["total"] is not None else avg * max(used_instances, 1)
        spread = mx - mn
    else:
        used_instances = len(used)
        avg = statistics.mean(used)
        mn = min(used)
        mx = max(used)
        total = sum(used)
        spread = mx - mn

    print(",".join([
        args.mode, args.order, args.visc_order, args.keep_prec, args.visc_prec,
        args.press_prec, args.slau_rho_prec, args.slau_u_prec, args.slau_p_prec,
        args.nx, args.kernel, str(total_instances), str(used_instances),
        f"{avg:.3f}", f"{mn:.3f}", f"{mx:.3f}", f"{spread:.3f}", f"{total:.3f}",
        args.trace_csv or "", args.kern_sum_csv or "", args.report or "",
        args.sqlite or "", source,
    ]))


if __name__ == "__main__":
    main()
