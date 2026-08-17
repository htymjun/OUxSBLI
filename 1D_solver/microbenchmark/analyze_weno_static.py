#!/usr/bin/env python3
"""Static SASS/PTXAS summary for weno_micro modes.

This is meant for the two checks the paper now calls out explicitly:
1. static instruction accounting per kernel
2. shared-memory / register fixed-cost accounting from ptxas

The important detail is that several user-facing modes share one kernel body
(`var_fp64_seq` and `var_fp32_seq` both launch `weno_var_seq`, etc.), so the
report always prints both the mode name and the underlying kernel symbol.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, OrderedDict
from pathlib import Path

FP64 = {"DADD", "DMUL", "DFMA", "DSETP", "DMMA"}
FP32 = {"FADD", "FMUL", "FFMA", "FSETP", "FSEL", "FMNMX", "MUFU"}
CONV = {"F2F", "F2I", "I2F"}
INSTR = re.compile(r"^\s+/\*[0-9a-f]{4}\*/\s+(.*?);", re.I)
PTXAS_COMPILE = re.compile(r"Compiling entry function '([^']+)'")
PTXAS_USED = re.compile(
    r"Used\s+(\d+)\s+registers(?:,\s+(\d+)\s+bytes\s+smem)?(?:,\s+(\d+)\s+bytes\s+cmem\[0\])?(?:,\s+(\d+)\s+bytes\s+cmem\[2\])?"
)

MODE_TO_KERNEL = {
    "var_fp64_seq": "weno_micro_kernels_weno_var_seq_",
    "var_fp32_seq": "weno_micro_kernels_weno_var_seq_",
    "var_rho64_seq": "weno_micro_kernels_weno_var_seq_",
    "var_u64_seq": "weno_micro_kernels_weno_var_seq_",
    "var_p64_seq": "weno_micro_kernels_weno_var_seq_",
    "var_fp64_warp": "weno_micro_kernels_weno_var_warp_",
    "var_fp32_warp": "weno_micro_kernels_weno_var_warp_",
    "var_rho64_warp": "weno_micro_kernels_weno_var_warp_",
    "var_u64_warp": "weno_micro_kernels_weno_var_warp_",
    "var_p64_warp": "weno_micro_kernels_weno_var_warp_",
    "weight_poly32_seq": "weno_micro_kernels_weno_weight_poly_seq_",
    "weight_poly32_serial": "weno_micro_kernels_weno_weight_poly_serial_warp_",
    "weight_poly32_warp": "weno_micro_kernels_weno_weight_poly_warp_",
    "weight_poly32_serial_oncebar": "weno_micro_kernels_weno_weight_poly_serial_oncebar_",
    "weight_poly32_warp_oncebar": "weno_micro_kernels_weno_weight_poly_warp_oncebar_",
    "weight_poly32_wsmem_serial": "weno_micro_kernels_weno_weight_poly_wsmem_serial_",
    "weight_poly32_wsmem_warp": "weno_micro_kernels_weno_weight_poly_wsmem_warp_",
    "weight_poly32_wsmem_tile2_serial": "weno_micro_kernels_weno_weight_poly_wsmem_tile2_serial_",
    "weight_poly32_wsmem_tile2_warp": "weno_micro_kernels_weno_weight_poly_wsmem_tile2_warp_",
    "weight_poly32_halfwarp_serial": "weno_micro_kernels_weno_weight_poly_halfwarp_serial_",
    "weight_poly32_halfwarp_shfl": "weno_micro_kernels_weno_weight_poly_halfwarp_shfl_",
    "var_fp64_seq7": "weno_micro_kernels_weno_var_seq7_",
    "var_fp32_seq7": "weno_micro_kernels_weno_var_seq7_",
    "weight_poly32_seq7": "weno_micro_kernels_weno_weight_poly_seq7_",
    "weight_poly32_serial7": "weno_micro_kernels_weno_weight_poly_serial_warp7_",
    "weight_poly32_warp7": "weno_micro_kernels_weno_weight_poly_warp7_",
    "var_fp64_seq9": "weno_micro_kernels_weno_var_seq9_",
    "var_fp32_seq9": "weno_micro_kernels_weno_var_seq9_",
    "weight_poly32_seq9": "weno_micro_kernels_weno_weight_poly_seq9_",
    "weight_poly32_serial9": "weno_micro_kernels_weno_weight_poly_serial_warp9_",
    "weight_poly32_warp9": "weno_micro_kernels_weno_weight_poly_warp9_",
    "w64_only_seq": "weno_micro_kernels_weno_w64_only_seq_",
    "poly32_only_seq": "weno_micro_kernels_weno_poly32_only_seq_",
    "w64_only_seq7": "weno_micro_kernels_weno_w64_only_seq7_",
    "poly32_only_seq7": "weno_micro_kernels_weno_poly32_only_seq7_",
    "w64_only_seq9": "weno_micro_kernels_weno_w64_only_seq9_",
    "poly32_only_seq9": "weno_micro_kernels_weno_poly32_only_seq9_",
    # family G: double-float hybrid. poly64_only_seq9 is the FP64 reference the
    # DF arms replace, so it is the denominator when reading C off the
    # fp32_instructions column.
    "poly64_only_seq9": "weno_micro_kernels_weno_poly64_only_seq9_",
    "polydf_only_seq9": "weno_micro_kernels_weno_polydf_only_seq9_",
    "polydfr_only_seq9": "weno_micro_kernels_weno_polydfr_only_seq9_",
    "w64_polydfr_seq9": "weno_micro_kernels_weno_w64_polydfr_seq9_",
    "w64_polydfr_dfin_seq9": "weno_micro_kernels_weno_w64_polydfr_dfin_seq9_",
    "polydfr_dfin_only_seq9": "weno_micro_kernels_weno_polydfr_dfin_only_seq9_",
    "w32_poly64_seq": "weno_micro_kernels_weno_w32_poly64_seq_",
    "w32_poly64_serial": "weno_micro_kernels_weno_w32_poly64_serial_",
    "w32_poly64_warp": "weno_micro_kernels_weno_w32_poly64_warp_",
    "w32_poly64_seq7": "weno_micro_kernels_weno_w32_poly64_seq7_",
    "w32_poly64_serial7": "weno_micro_kernels_weno_w32_poly64_serial7_",
    "w32_poly64_warp7": "weno_micro_kernels_weno_w32_poly64_warp7_",
    "w32_poly64_seq9": "weno_micro_kernels_weno_w32_poly64_seq9_",
    "w32_poly64_serial9": "weno_micro_kernels_weno_w32_poly64_serial9_",
    "w32_poly64_warp9": "weno_micro_kernels_weno_w32_poly64_warp9_",
    "w32mixh_poly64_seq9": "weno_micro_kernels_weno_w32mix_poly64_seq9_",
    "w32mix1_poly64_seq9": "weno_micro_kernels_weno_w32mix_poly64_seq9_",
    "var3_fp64_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var3_fp32_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var3_k1_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var3_k1_warp9": "weno_micro_kernels_weno_varsplit_warp9_",
    "var4_fp64_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var4_fp32_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var4_k1_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var4_k2_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var4_k1_warp9": "weno_micro_kernels_weno_varsplit_warp9_",
    "var5_fp64_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var5_fp32_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var5_k1_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var5_k2_seq9": "weno_micro_kernels_weno_varsplit_seq9_",
    "var5_k2_warp9": "weno_micro_kernels_weno_varsplit_warp9_",
}


def opcode(text: str) -> str:
    t = text.strip()
    t = re.sub(r"^@!?U?P[T0-9]\s+", "", t)
    return t.split()[0].split(".")[0] if t.split() else ""


def classify(op: str) -> str:
    if op in FP64:
        return "FP64"
    if op in FP32:
        return "FP32"
    if op in CONV:
        return "CONV"
    return "other"


def parse_sass_listing(text: str) -> OrderedDict[str, list[str]]:
    kernels: OrderedDict[str, list[str]] = OrderedDict()
    cur: str | None = None
    for line in text.splitlines():
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


def parse_ptxas_build_log(path: Path) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    cur: str | None = None
    for line in path.read_text().splitlines():
        m = PTXAS_COMPILE.search(line)
        if m:
            cur = m.group(1)
            out.setdefault(cur, {})
            continue
        m = PTXAS_USED.search(line)
        if m and cur is not None:
            out[cur] = {
                "registers": int(m.group(1)),
                "smem_bytes": int(m.group(2) or 0),
                "cmem0_bytes": int(m.group(3) or 0),
                "cmem2_bytes": int(m.group(4) or 0),
            }
            cur = None
    return out


def summarize_kernel(instrs: list[str]) -> dict[str, int]:
    ops = [opcode(i) for i in instrs]
    cls = Counter(classify(op) for op in ops)
    letters = [("D" if classify(op) == "FP64" else "S" if classify(op) == "FP32" else "c")
               for op in ops if classify(op) in ("FP64", "FP32", "CONV")]
    ds = [c for c in letters if c in ("D", "S")]
    runs = 0
    prev = None
    for c in ds:
        if c != prev:
            runs += 1
            prev = c
    return {
        "total_instructions": len(ops),
        "fp64_instructions": cls["FP64"],
        "fp32_instructions": cls["FP32"],
        "convert_instructions": cls["CONV"],
        "other_instructions": cls["other"],
        "ds_runs": runs,
        "ds_length": len(ds),
    }


def load_sass(args: argparse.Namespace) -> OrderedDict[str, list[str]]:
    if args.sass_file:
        return parse_sass_listing(Path(args.sass_file).read_text())
    if not args.exe:
        return OrderedDict()
    cmd = [args.cuobjdump, "-sass", str(args.exe)]
    proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return parse_sass_listing(proc.stdout)


def resolve_modes(requested: list[str]) -> list[str]:
    if not requested:
        return list(MODE_TO_KERNEL)
    bad = [m for m in requested if m not in MODE_TO_KERNEL]
    if bad:
        raise SystemExit(f"unknown modes: {', '.join(bad)}")
    return requested


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exe", type=Path, help="weno_micro executable for cuobjdump -sass")
    ap.add_argument("--build-log", type=Path, help="build.log with ptxas register/smem lines")
    ap.add_argument("--sass-file", type=Path, help="pre-dumped cuobjdump -sass text")
    ap.add_argument("--cuobjdump", default="cuobjdump", help="cuobjdump executable")
    ap.add_argument("--output", type=Path, help="write CSV here")
    ap.add_argument("--modes", nargs="*", default=[], help="subset of mode names")
    args = ap.parse_args()

    modes = resolve_modes(args.modes)
    sass = load_sass(args)
    ptxas = parse_ptxas_build_log(args.build_log) if args.build_log else {}

    rows = []
    for mode in modes:
        kernel = MODE_TO_KERNEL[mode]
        row = {"mode": mode, "kernel": kernel}
        row.update(ptxas.get(kernel, {
            "registers": -1, "smem_bytes": -1, "cmem0_bytes": -1, "cmem2_bytes": -1
        }))
        if kernel in sass:
            row.update(summarize_kernel(sass[kernel]))
        else:
            row.update({
                "total_instructions": -1,
                "fp64_instructions": -1,
                "fp32_instructions": -1,
                "convert_instructions": -1,
                "other_instructions": -1,
                "ds_runs": -1,
                "ds_length": -1,
            })
        rows.append(row)

    fieldnames = [
        "mode", "kernel", "registers", "smem_bytes", "cmem0_bytes", "cmem2_bytes",
        "total_instructions", "fp64_instructions", "fp32_instructions",
        "convert_instructions", "other_instructions", "ds_runs", "ds_length",
    ]
    if args.output:
        with args.output.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

    w = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
