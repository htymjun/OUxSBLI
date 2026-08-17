#!/usr/bin/env python3
"""Static SASS accounting for the w32ptx kernels, and the fidelity gate against
the Fortran reference.

Reuses the mnemonic classification from ../analyze_weno_static.py -- that
partition (which opcodes land on the FP64 pipe, which on FP32, which are
conversions) is the valuable, hard-won part and must not be forked.

Two things this adds over the Fortran-side script:

1. `divctl` = FCHK + CALL + BSSY + BSYNC. The IEEE FP32 division's cost is split
   between an inline FCHK (which falls into `other`) and an out-of-line slow-path
   body reached by CALL.REL.NOINC (whose instructions are attributed to a
   different `Function :` block and are therefore invisible). Without this column
   the div.rn -> rcp.approx change looks far smaller statically than it is in the
   timing.
2. `mufu` split out of the FP32 column. MUFU issues on the SFU, which is a
   separate pipe from the FMA pipe, so folding it into `fp32` overstates the
   saturating stream.

Symbols are matched EXACTLY, never as a substring: `k_..._rcp` is a prefix of
`k_..._rcpn`, and a regex match silently sums the two kernels.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
try:
    from analyze_weno_static import FP64, FP32, CONV  # noqa: E402
except Exception as exc:  # pragma: no cover
    print(f"error: cannot import mnemonic sets from ../analyze_weno_static.py: {exc}",
          file=sys.stderr)
    raise SystemExit(2)

FUNC = re.compile(r"^\s*Function : (\S+)\s*$")
INSTR = re.compile(r"^\s+/\*[0-9a-f]{4}\*/\s+(.*?);", re.I)
DIVCTL = {"FCHK", "CALL", "BSSY", "BSYNC"}

# The Fortran reference, measured from ../build/weno_micro at cc89 and matching
# ../weno_allrun_20260817_170014/static_instruction_summary.csv (cc80) exactly.
#
# WHAT IS ASSERTED EXACTLY, AND WHY THE REST IS NOT
#
# `dfma_dmul` -- the FP64 candidate-polynomial and combine arithmetic -- matches
# instruction for instruction at all three widths (55 / 112 / 169). So do conv,
# ldg, stg, fchk and call. Those are the counts that define the kernel.
#
# Two quantities differ for understood, benign reasons and are reported rather
# than gated:
#
#  * dadd. nvfortran unrolls `do k=1,nrepeat` by 4 with a predicated remainder,
#    so 13/19/19 of its DADDs are unroll copies never executed at nrepeat=1.
#    The C++ side carries `#pragma unroll 1`. Same work, different scaffolding.
#  * fp32. nvcc finds ~39 more FMA contractions than nvfortran in the dense beta
#    forms (WENO9: FFMA +42, FMUL -27, FADD -54 -- FMNMX, FSETP and MUFU all
#    match exactly). That is a compiler-quality difference, and it makes the C++
#    FP32 stream ~3.9% SHORTER than the Fortran's, which slightly understates the
#    saturating stream. Recorded, not chased: pinning ~200 beta expressions with
#    __fmaf_rn to reproduce one compiler's contraction choices would make the
#    numerics unreadable for no experimental gain.
#
# The consequence to keep in mind: because FP32 contraction differs, the C++
# _ref is a faithful implementation of the same ALGORITHM, not a bit-identical
# copy of the same INSTRUCTION SEQUENCE. Its output therefore differs from the
# Fortran's at the FP32 weight level (~1e-7 relative), which is the same order as
# the FP32 split's own deviation from FP64. The correctness basis is accordingly
# (a) the --accuracy order/deviation gate and (b) bit-identity BETWEEN rungs
# inside this binary, where the beta code is compiled identically and the
# reciprocal is the only difference.
#
#   width -> dict of exactly-asserted quantities, plus reference-only values
FORTRAN_REF = {
    5: dict(dfma_dmul=55, conv=36, ldg=18, stg=6, fchk=18, call=24,
            ref_dadd=19, ref_fp32=426),
    7: dict(dfma_dmul=112, conv=48, ldg=24, stg=6, fchk=24, call=30,
            ref_dadd=34, ref_fp32=708),
    9: dict(dfma_dmul=169, conv=60, ldg=30, stg=6, fchk=30, call=36,
            ref_dadd=37, ref_fp32=1050),
}
EXACT_KEYS = ["dfma_dmul", "conv", "ldg", "stg", "fchk", "call"]


def opcode(text: str) -> str:
    """Leading mnemonic, stripped of predicate and of .MODIFIER suffixes."""
    t = text.strip()
    if t.startswith("@"):
        t = t.split(None, 1)[1] if " " in t else ""
    t = t.strip()
    if not t:
        return ""
    return t.split()[0].split(".")[0]


def parse_sass(path: Path) -> dict[str, Counter]:
    per: dict[str, Counter] = {}
    cur = None
    for line in path.read_text(errors="replace").splitlines():
        m = FUNC.match(line)
        if m:
            cur = m.group(1)
            per.setdefault(cur, Counter())
            continue
        if cur is None:
            continue
        m = INSTR.match(line)
        if m:
            op = opcode(m.group(1))
            if op:
                per[cur][op] += 1
    return per


def classify(c: Counter) -> dict:
    total = sum(c.values())
    fp64 = sum(v for k, v in c.items() if k in FP64)
    mufu = c.get("MUFU", 0)
    fp32 = sum(v for k, v in c.items() if k in FP32)
    conv = sum(v for k, v in c.items() if k in CONV)
    divctl = sum(v for k, v in c.items() if k in DIVCTL)
    return dict(
        total=total, fp64=fp64, fp32=fp32, mufu=mufu, fp32_fma=fp32 - mufu,
        conv=conv, divctl=divctl, fchk=c.get("FCHK", 0), call=c.get("CALL", 0),
        ldg=c.get("LDG", 0), stg=c.get("STG", 0),
        spill=c.get("LDL", 0) + c.get("STL", 0),
        other=total - fp64 - fp32 - conv,
        dfma_dmul=c.get("DFMA", 0) + c.get("DMUL", 0),
        dadd=c.get("DADD", 0),
    )


def width_of(name: str) -> int | None:
    m = re.search(r"seq(\d)", name)
    return int(m.group(1)) if m else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sass", default=str(HERE / "w32ptx.sass"))
    ap.add_argument("--bin", default=str(HERE / "w32ptx"))
    ap.add_argument("--csv", help="write the table here as CSV")
    args = ap.parse_args()

    sass = Path(args.sass)
    if not sass.exists():
        print(f"error: {sass} not found -- run `make ARCH=<cc> sass` first", file=sys.stderr)
        return 2

    per = parse_sass(sass)

    order = []
    try:
        out = subprocess.run([args.bin, "--list-kernels"], capture_output=True, text=True,
                             check=True).stdout
        for line in out.splitlines():
            if "\t" in line:
                mode, kern = line.split("\t", 1)
                order.append((mode.strip(), kern.strip()))
    except Exception:
        order = [(k, k) for k in sorted(per)]

    cols = ["total", "fp64", "fp32_fma", "mufu", "conv", "divctl", "fchk", "call",
            "ldg", "stg", "spill", "other"]
    print(f"{'mode':<30}" + "".join(f"{c:>9}" for c in cols))
    rows = []
    for mode, kern in order:
        c = per.get(kern)
        if c is None:
            print(f"{mode:<30}  (not found in SASS)")
            continue
        d = classify(c)
        rows.append((mode, d))
        print(f"{mode:<30}" + "".join(f"{d[c]:>9}" for c in cols))

    # ---- fidelity gate: the _ref rungs must reproduce the Fortran exactly ----
    print("\n=== fidelity gate: w32_poly64_seq{5,7,9}_ref vs Fortran weno_micro ===")
    print(f"{'width':>5} {'quantity':<10} {'fortran':>9} {'cpp':>9}   verdict")
    ok = True
    by_mode = dict(rows)
    for w, ref in sorted(FORTRAN_REF.items()):
        mode = f"w32_poly64_seq{w}_ref"
        d = by_mode.get(mode)
        if d is None:
            print(f"{w:>5} (missing {mode})")
            ok = False
            continue
        for qty in EXACT_KEYS:
            want, got = ref[qty], d[qty]
            good = got == want
            ok = ok and good
            print(f"{w:>5} {qty:<10} {want:>9} {got:>9}   {'PASS' if good else 'FAIL'}")
        # Reported, not gated -- see the FORTRAN_REF comment for why.
        f32 = d["fp32_fma"] + d["mufu"]
        print(f"{w:>5} {'dadd':<10} {ref['ref_dadd']:>9} {d['dadd']:>9}   "
              f"info (nvfortran repeat-loop unroll, delta {d['dadd'] - ref['ref_dadd']:+d})")
        print(f"{w:>5} {'fp32':<10} {ref['ref_fp32']:>9} {f32:>9}   "
              f"info (nvcc FMA contraction, "
              f"{100.0 * (f32 - ref['ref_fp32']) / ref['ref_fp32']:+.1f}%)")

    # ---- the asm rungs must have removed every division branch ----
    print("\n=== asm rungs: division control flow must be gone, no spill ===")
    for mode, d in rows:
        if not re.search(r"_(rcp|rcpn|divfull)$", mode):
            continue
        good = d["fchk"] == 0 and d["call"] == 0 and d["spill"] == 0
        ok = ok and good
        print(f"{mode:<30} fchk={d['fchk']:<3} call={d['call']:<3} spill={d['spill']:<3} "
              f"{'PASS' if good else 'FAIL'}")

    if args.csv:
        import csv as _csv
        with open(args.csv, "w", newline="") as fh:
            w = _csv.writer(fh)
            w.writerow(["mode"] + cols)
            for mode, d in rows:
                w.writerow([mode] + [d[c] for c in cols])
        print(f"\nwrote {args.csv}")

    print(f"\n-> {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
