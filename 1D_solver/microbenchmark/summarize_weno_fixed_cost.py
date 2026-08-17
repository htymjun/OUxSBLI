#!/usr/bin/env python3
"""Summarize the fixed-cost sweep from summary.csv into markdown.

The goal is not to invent a new metric, but to put the existing controlled
pairs next to the interpretation used in the paper:
  seq -> serial          structure cost with no warp-level overlap
  serial -> warp         warp-level overlap benefit at matched structure
  serial -> oncebar      barrier-placement effect
  oncebar -> wsmem       shared-memory footprint / handoff effect
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_best(path: Path) -> dict[str, float]:
    best: dict[str, float] = {}
    with path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            mode = row["mode"]
            t = float(row["time_min_us"])
            if mode not in best or t < best[mode]:
                best[mode] = t
    return best


def ratio(a: float, b: float) -> float:
    return a / b


def delta_pct(a: float, b: float) -> float:
    return 100.0 * (a / b - 1.0)


def fmt(x: float) -> str:
    return f"{x:.3f}"


def emit_line(lines: list[str], text: str = "") -> None:
    lines.append(text)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("summary_csv", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    best = load_best(args.summary_csv)
    lines: list[str] = []
    emit_line(lines, "# WENO Fixed-Cost Summary")
    emit_line(lines)
    emit_line(lines, f"source: `{args.summary_csv}`")
    emit_line(lines)

    def have(*modes: str) -> bool:
        return all(m in best for m in modes)

    if have(
        "weight_poly32_seq",
        "w64_only_seq",
        "poly32_only_seq",
        "weight_poly32_serial",
        "weight_poly32_warp",
        "weight_poly32_serial_oncebar",
        "weight_poly32_warp_oncebar",
        "weight_poly32_wsmem_serial",
        "weight_poly32_wsmem_warp",
        "weight_poly32_wsmem_tile2_serial",
        "weight_poly32_wsmem_tile2_warp",
    ):
        emit_line(lines, "## WENO5 weight\\_poly32 family")
        emit_line(lines)
        emit_line(lines, "| comparison | lhs [us] | rhs [us] | lhs/rhs | interpretation |")
        emit_line(lines, "|---|---:|---:|---:|---|")
        emit_line(
            lines,
            f"| `seq -> serial` | {fmt(best['weight_poly32_serial'])} | {fmt(best['weight_poly32_seq'])} | {fmt(ratio(best['weight_poly32_serial'], best['weight_poly32_seq']))} | structure cost of 256-thread/shared-memory/barrier layout without warp-level overlap |",
        )
        emit_line(
            lines,
            f"| `serial -> warp` | {fmt(best['weight_poly32_serial'])} | {fmt(best['weight_poly32_warp'])} | {fmt(ratio(best['weight_poly32_serial'], best['weight_poly32_warp']))} | warp-level overlap benefit at matched structure |",
        )
        emit_line(
            lines,
            f"| `serial -> serial_oncebar` | {fmt(best['weight_poly32_serial'])} | {fmt(best['weight_poly32_serial_oncebar'])} | {fmt(ratio(best['weight_poly32_serial'], best['weight_poly32_serial_oncebar']))} | barrier placement effect in the serial control |",
        )
        emit_line(
            lines,
            f"| `warp -> warp_oncebar` | {fmt(best['weight_poly32_warp'])} | {fmt(best['weight_poly32_warp_oncebar'])} | {fmt(ratio(best['weight_poly32_warp'], best['weight_poly32_warp_oncebar']))} | barrier placement effect in the warp form |",
        )
        emit_line(
            lines,
            f"| `serial_oncebar -> wsmem_serial` | {fmt(best['weight_poly32_serial_oncebar'])} | {fmt(best['weight_poly32_wsmem_serial'])} | {fmt(ratio(best['weight_poly32_serial_oncebar'], best['weight_poly32_wsmem_serial']))} | effect of keeping FP32 polynomials out of shared memory |",
        )
        emit_line(
            lines,
            f"| `warp_oncebar -> wsmem_warp` | {fmt(best['weight_poly32_warp_oncebar'])} | {fmt(best['weight_poly32_wsmem_warp'])} | {fmt(ratio(best['weight_poly32_warp_oncebar'], best['weight_poly32_wsmem_warp']))} | same effect in the warp form |",
        )
        emit_line(
            lines,
            f"| `wsmem_serial -> wsmem_tile2_serial` | {fmt(best['weight_poly32_wsmem_tile2_serial'])} | {fmt(best['weight_poly32_wsmem_serial'])} | {fmt(ratio(best['weight_poly32_wsmem_tile2_serial'], best['weight_poly32_wsmem_serial']))} | two-face tiling tradeoff in the serial control |",
        )
        emit_line(
            lines,
            f"| `wsmem_warp -> wsmem_tile2_warp` | {fmt(best['weight_poly32_wsmem_tile2_warp'])} | {fmt(best['weight_poly32_wsmem_warp'])} | {fmt(ratio(best['weight_poly32_wsmem_tile2_warp'], best['weight_poly32_wsmem_warp']))} | two-face tiling tradeoff in the warp form |",
        )
        emit_line(lines)
        emit_line(
            lines,
            f"`serial -> warp` is a {delta_pct(best['weight_poly32_serial'], best['weight_poly32_warp']):.1f}% gain, while `seq -> serial` is a {delta_pct(best['weight_poly32_serial'], best['weight_poly32_seq']):.1f}% penalty.",
        )
        emit_line(lines)
        hidden_delta = best["weight_poly32_seq"] - best["w64_only_seq"]
        poly_only = best["poly32_only_seq"]
        emit_line(lines, "### single-thread overlap check")
        emit_line(lines)
        emit_line(lines, "| quantity | value [us] | meaning |")
        emit_line(lines, "|---|---:|---|")
        emit_line(lines, f"| `weight_poly32_seq - w64_only_seq` | {fmt(hidden_delta)} | observed extra time after adding the FP32 polynomial half to the FP64-weight stream |")
        emit_line(lines, f"| `poly32_only_seq` | {fmt(poly_only)} | standalone FP32 polynomial-half ablation |")
        emit_line(lines, f"| `poly32_only_seq / (weight_poly32_seq - w64_only_seq)` | {fmt(poly_only / hidden_delta if hidden_delta else float('inf'))} | values >> 1 suggest the polynomial half is being hidden rather than merely being cheap |")
        emit_line(lines)

    core_groups = [
        ("weight_poly32", ["", "7", "9"]),
        ("w32_poly64", ["", "7", "9"]),
    ]
    emit_line(lines, "## Order Sweep")
    emit_line(lines)
    emit_line(lines, "| family | order | seq [us] | serial [us] | warp [us] | serial/warp | warp/seq |")
    emit_line(lines, "|---|---:|---:|---:|---:|---:|---:|")
    for family, suffixes in core_groups:
        for suffix in suffixes:
            order = "5" if suffix == "" else suffix
            seq = f"{family}_seq{suffix}"
            serial = f"{family}_serial{suffix}"
            warp = f"{family}_warp{suffix}"
            if not have(seq, serial, warp):
                continue
            emit_line(
                lines,
                f"| `{family}` | {order} | {fmt(best[seq])} | {fmt(best[serial])} | {fmt(best[warp])} | {fmt(ratio(best[serial], best[warp]))} | {fmt(ratio(best[warp], best[seq]))} |",
            )
    emit_line(lines)

    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
