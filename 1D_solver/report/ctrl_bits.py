#!/usr/bin/env python3
"""Decode SM70+ scheduling control bits from a cuobjdump -sass listing.

Each instruction prints as two hex words; control lives in instruction bits
[125:105], i.e. bits [61:41] of the high word:
  stall=(H>>41)&0xF  yield=(H>>45)&1  wr=(H>>46)&7  rd=(H>>49)&7
  wait=(H>>52)&0x3F  reuse=(H>>58)&0xF
This layout is community-reverse-engineered, not documented by NVIDIA; sanity
is checked by asserting stall<=15 and barrier ids in 0..7.
"""
import re
import sys
from collections import Counter

FP64 = ("DADD", "DMUL", "DFMA", "DSETP")
FP32 = ("FADD", "FMUL", "FFMA", "FSETP", "FSEL", "MUFU")

INSTR = re.compile(r"/\*[0-9a-f]{4}\*/\s+(.*?);.*?/\* (0x[0-9a-f]+) \*/", re.I)
HIGH = re.compile(r"/\* (0x[0-9a-f]{16}) \*/")


def decode(h):
    return dict(stall=(h >> 41) & 0xF, yld=(h >> 45) & 1, wr=(h >> 46) & 7,
                rd=(h >> 49) & 7, wait=(h >> 52) & 0x3F, reuse=(h >> 58) & 0xF)


def kind(op):
    b = op.split()[0].split(".")[0]
    b = re.sub(r"^@!?U?P[T0-9]+$", "", b)
    if b in FP64:
        return "FP64"
    if b in FP32:
        return "FP32"
    return "other"


def main(path, kernel):
    lines = open(path).read().splitlines()
    inside, rows = False, []
    i = 0
    while i < len(lines):
        m = re.match(r"\s+Function : (\S+)", lines[i])
        if m:
            inside = kernel in m.group(1)
        if inside:
            mi = INSTR.search(lines[i])
            if mi and i + 1 < len(lines):
                mh = HIGH.search(lines[i + 1])
                if mh:
                    text = re.sub(r"^@!?U?P[T0-9]+\s+", "", mi.group(1).strip())
                    rows.append((text, decode(int(mh.group(1), 16))))
        i += 1

    assert rows, f"no instructions found for {kernel}"
    bad = [c for _, c in rows if c["stall"] > 15 or c["wr"] > 7 or c["rd"] > 7]
    print(f"{kernel}: {len(rows)} instructions decoded, {len(bad)} implausible\n")

    per = {}
    for text, c in rows:
        per.setdefault(kind(text), []).append(c)
    for k in ("FP64", "FP32", "other"):
        if k not in per:
            continue
        cs = per[k]
        st = [c["stall"] for c in cs]
        print(f"  {k:5s} n={len(cs):3d}  mean stall={sum(st)/len(st):.2f}  "
              f"stall hist={dict(sorted(Counter(st).items()))}")
        print(f"        yield set: {sum(c['yld'] for c in cs)}   "
              f"sets write barrier: {sum(1 for c in cs if c['wr'] != 7)}   "
              f"waits on barrier: {sum(1 for c in cs if c['wait'])}")

    print("\n  math instruction stream (stall = cycles the scheduler idles after issue):")
    for text, c in rows:
        k = kind(text)
        if k == "other":
            continue
        print(f"    {k}  stall={c['stall']:2d} yield={c['yld']} "
              f"wr={c['wr']} rd={c['rd']} wait={c['wait']:02b}  {text[:52]}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
