#!/usr/bin/env python3
"""Split a cuobjdump -sass listing per kernel and classify FP64 / FP32 / convert opcodes."""
import re
import sys
from collections import Counter, OrderedDict

FP64 = {"DADD", "DMUL", "DFMA", "DSETP", "DMMA"}
FP32 = {"FADD", "FMUL", "FFMA", "FSETP", "FSEL", "FMNMX", "MUFU"}
CONV = {"F2F", "F2I", "I2F"}

INSTR = re.compile(r"^\s+/\*[0-9a-f]{4}\*/\s+(.*?);", re.I)


def opcode(text):
    t = text.strip()
    t = re.sub(r"^@!?U?P[T0-9]\s+", "", t)          # strip predicate guard
    return t.split()[0].split(".")[0] if t.split() else ""


def parse(path):
    kernels = OrderedDict()
    cur = None
    for line in open(path):
        m = re.match(r"\s+Function : (\S+)", line)
        if m:
            cur = m.group(1)
            kernels[cur] = []
            continue
        if cur is None:
            continue
        m = INSTR.match(line)
        if m:
            kernels[cur].append(m.group(1).strip())
    return kernels


def klass(op):
    if op in FP64:
        return "FP64"
    if op in FP32:
        return "FP32"
    if op in CONV:
        return "CONV"
    return "other"


def report(name, instrs):
    ops = [opcode(i) for i in instrs]
    counts = Counter(ops)
    cls = Counter(klass(o) for o in ops)
    print(f"\n{'='*72}\n{name}   ({len(ops)} instructions)\n{'='*72}")
    print(f"  FP64 (D*): {cls['FP64']:4d}   FP32 (F*): {cls['FP32']:4d}   "
          f"convert: {cls['CONV']:4d}   other: {cls['other']:4d}")
    for grp, members in (("FP64", FP64), ("FP32", FP32), ("CONV", CONV)):
        detail = {k: v for k, v in counts.items() if k in members}
        if detail:
            print(f"    {grp}: " + ", ".join(f"{k}={v}" for k, v in
                                             sorted(detail.items(), key=lambda x: -x[1])))
    # math-only instruction sequence, to see interleaving vs clustering
    seq = [(i, klass(o)) for i, o in enumerate(ops) if klass(o) in ("FP64", "FP32", "CONV")]
    if seq:
        letters = "".join({"FP64": "D", "FP32": "S", "CONV": "c"}[c] for _, c in seq)
        print(f"  math sequence (D=fp64, S=fp32, c=convert), in program order:")
        for i in range(0, len(letters), 64):
            print(f"    {letters[i:i+64]}")
        # count transitions between D and S runs, ignoring converts
        ds = [c for c in letters if c in "DS"]
        runs = [(ds[0], 1)] if ds else []
        for c in ds[1:]:
            if c == runs[-1][0]:
                runs[-1] = (c, runs[-1][1] + 1)
            else:
                runs.append((c, 1))
        print(f"  D/S runs: " + " ".join(f"{c}x{n}" for c, n in runs))
        print(f"  -> {len(runs)} runs over {len(ds)} fp instructions "
              f"(1 run = fully clustered, {len(ds)} runs = fully alternating)")


if __name__ == "__main__":
    ks = parse(sys.argv[1])
    want = sys.argv[2:] or list(ks)
    for name in want:
        for k in ks:
            if name in k:
                report(k, ks[k])
