#!/usr/bin/env python3
"""SASS stream analysis for FP32/FP64 mixing.

This complements `static_w32ptx.py` by looking at the *order* of the emitted
instructions, not just their totals. It reports:

- per-mnemonic counts for the relevant FP32/FP64 ops
- FP32/FP64 adjacency statistics after filtering to arithmetic ops
- short traces that make it easy to inspect whether the two pipes are
  alternating or clustering

Important caveat: SASS is a static schedule, not a runtime trace. The metrics
here therefore show whether ptxas *laid out* FP32 and FP64 ops next to each
other, which is the strongest evidence we can extract from the binary alone.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

try:
    from analyze_weno_static import FP32, FP64, CONV  # noqa: E402
except Exception as exc:  # pragma: no cover
    print(f"error: cannot import mnemonic sets from ../analyze_weno_static.py: {exc}",
          file=sys.stderr)
    raise SystemExit(2)

FUNC = re.compile(r"^\s*Function : (\S+)\s*$")
INSTR = re.compile(r"^\s+/\*[0-9a-f]{4}\*/\s+(.*?);", re.I)

CORE_MNEMONICS = ("FFMA", "FMUL", "FADD", "DFMA", "DMUL", "DADD", "MUFU")
DETAIL_MNEMONICS = (
    "FFMA", "FMUL", "FADD", "FSET", "FMNMX", "MUFU",
    "DFMA", "DMUL", "DADD", "F2F", "FCHK", "CALL", "BSSY", "BSYNC",
)


def opcode(text: str) -> str:
    """Leading mnemonic, stripped of predicate and of .MODIFIER suffixes."""
    t = text.strip()
    if t.startswith("@"):
        t = t.split(None, 1)[1] if " " in t else ""
    t = t.strip()
    if not t:
        return ""
    return t.split()[0].split(".")[0]


def parse_sass(path: Path) -> dict[str, list[str]]:
    per: dict[str, list[str]] = {}
    cur = None
    for line in path.read_text(errors="replace").splitlines():
        m = FUNC.match(line)
        if m:
            cur = m.group(1)
            per.setdefault(cur, [])
            continue
        if cur is None:
            continue
        m = INSTR.match(line)
        if m:
            op = opcode(m.group(1))
            if op:
                per[cur].append(op)
    return per


def cls(op: str) -> str | None:
    if op in FP64:
        return "fp64"
    if op in FP32:
        return "fp32"
    return None


def summarize_ops(ops: list[str]) -> dict[str, int]:
    c = Counter(ops)
    out = {m: c.get(m, 0) for m in DETAIL_MNEMONICS}
    out["total"] = len(ops)
    out["fp32_total"] = sum(c.get(m, 0) for m in FP32)
    out["fp64_total"] = sum(c.get(m, 0) for m in FP64)
    out["conv_total"] = sum(c.get(m, 0) for m in CONV)
    return out


def run_stats(ops: list[str]) -> dict[str, float | int]:
    core = [cls(op) for op in ops if cls(op) in {"fp32", "fp64"}]
    if not core:
        return dict(core_ops=0, switches=0, fp32_fp64=0, fp64_fp32=0,
                    fp32_runs=0, fp64_runs=0, max_fp32_run=0, max_fp64_run=0,
                    alternation_rate=0.0, core_shares=0.0)

    switches = 0
    fp32_fp64 = 0
    fp64_fp32 = 0
    max_fp32_run = 0
    max_fp64_run = 0
    runs: list[tuple[str, int]] = []
    cur = core[0]
    n = 1
    for prev, nxt in zip(core, core[1:]):
        if nxt == cur:
            n += 1
            continue
        runs.append((cur, n))
        if cur == "fp32":
            max_fp32_run = max(max_fp32_run, n)
            if nxt == "fp64":
                fp32_fp64 += 1
        else:
            max_fp64_run = max(max_fp64_run, n)
            if nxt == "fp32":
                fp64_fp32 += 1
        switches += 1
        cur = nxt
        n = 1
    runs.append((cur, n))
    if cur == "fp32":
        max_fp32_run = max(max_fp32_run, n)
    else:
        max_fp64_run = max(max_fp64_run, n)

    fp32_runs = sum(1 for kind, _ in runs if kind == "fp32")
    fp64_runs = sum(1 for kind, _ in runs if kind == "fp64")
    alternation_rate = switches / max(len(core) - 1, 1)
    core_shares = len(core) / len(ops)
    return dict(
        core_ops=len(core),
        switches=switches,
        fp32_fp64=fp32_fp64,
        fp64_fp32=fp64_fp32,
        fp32_runs=fp32_runs,
        fp64_runs=fp64_runs,
        max_fp32_run=max_fp32_run,
        max_fp64_run=max_fp64_run,
        alternation_rate=alternation_rate,
        core_shares=core_shares,
    )


def trace_ops(ops: list[str], limit: int) -> str:
    core = [(i, op, cls(op)) for i, op in enumerate(ops) if cls(op) in {"fp32", "fp64"}]
    if not core:
        return "(no fp32/fp64 ops)"
    items = []
    for idx, op, kind in core[:limit]:
        items.append(f"{idx:04d}:{kind}:{op}")
    if len(core) > limit:
        items.append(f"... ({len(core) - limit} more core ops)")
    return " | ".join(items)


def modes_from_list(kernels: dict[str, list[str]], pattern: str | None) -> list[str]:
    modes = sorted(kernels)
    if pattern:
        rx = re.compile(pattern)
        modes = [m for m in modes if rx.search(m)]
    return modes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sass", default=str(HERE / "w32ptx.sass"))
    ap.add_argument("--csv", help="write the summary here as CSV")
    ap.add_argument("--mode", action="append", help="restrict to one or more mode regexes")
    ap.add_argument("--trace", type=int, default=24, help="number of core ops to print per kernel")
    args = ap.parse_args()

    sass = Path(args.sass)
    if not sass.exists():
        print(f"error: {sass} not found -- run `make ARCH=<cc> sass` first", file=sys.stderr)
        return 2

    kernels = parse_sass(sass)
    patterns = args.mode or []
    if patterns:
        selected = []
        for pat in patterns:
            rx = re.compile(pat)
            selected.extend([m for m in kernels if rx.search(m)])
        # Preserve order while deduplicating.
        seen = set()
        modes = [m for m in selected if not (m in seen or seen.add(m))]
    else:
        modes = sorted(kernels)

    rows = []
    for mode in modes:
        ops = kernels[mode]
        counts = summarize_ops(ops)
        stats = run_stats(ops)
        row = {"mode": mode, **counts, **stats}
        rows.append(row)

    header = [
        "mode", "total", "fp64_total", "fp32_total", "conv_total",
        *DETAIL_MNEMONICS,
        "core_ops", "switches", "fp32_fp64", "fp64_fp32",
        "fp32_runs", "fp64_runs", "max_fp32_run", "max_fp64_run",
        "alternation_rate", "core_shares",
    ]

    print(
        f"{'mode':<30}"
        f"{'total':>7}{'fp64':>7}{'fp32':>7}{'core':>7}"
        f"{'sw':>6}{'a_rate':>9}{'max32':>7}{'max64':>7}"
    )
    for row in rows:
        print(
            f"{row['mode']:<30}"
            f"{row['total']:>7}{row['fp64_total']:>7}{row['fp32_total']:>7}"
            f"{row['core_ops']:>7}{row['switches']:>6}"
            f"{row['alternation_rate']:>9.3f}{row['max_fp32_run']:>7}{row['max_fp64_run']:>7}"
        )
    print("\n=== details ===")
    for row in rows:
        mode = row["mode"]
        print(f"\n[{mode}]")
        print("counts: " + ", ".join(f"{m}={row[m]}" for m in DETAIL_MNEMONICS))
        print(
            "mix: "
            f"core_ops={row['core_ops']}, switches={row['switches']}, "
            f"fp32->fp64={row['fp32_fp64']}, fp64->fp32={row['fp64_fp32']}, "
            f"fp32_runs={row['fp32_runs']}, fp64_runs={row['fp64_runs']}, "
            f"max_fp32_run={row['max_fp32_run']}, max_fp64_run={row['max_fp64_run']}, "
            f"alternation_rate={row['alternation_rate']:.3f}, core_share={row['core_shares']:.3f}"
        )
        print("trace: " + trace_ops(kernels[mode], args.trace))

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header)
            w.writeheader()
            for row in rows:
                w.writerow(row)
        print(f"\nwrote {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
