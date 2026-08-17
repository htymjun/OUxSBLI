# CLAUDE.md — w32ptx

`w32_poly64_seq{5,7,9}` extracted from `../weno_micro.f90` (3280 lines) into a
standalone CUDA C++ benchmark, so the FP32 half can be attacked with **inline
PTX**. Harness-wide rules live in `../CLAUDE.md`; this file covers only what is
different here.

## Why C++ and not Fortran

**nvfortran has no inline-asm statement.** `nvfortran -help` offers only
`-Mkeepasm` (host assembly); `asm volatile` is a CUDA C++ construct. `nvcc` 12.5
ships in the same HPC SDK, so the extraction is the way to reach PTX at all.
Consequence: nothing here ports back into the solver's CUDA Fortran directly —
this is a probe for whether the headroom exists, the same role
`3D_solver_mixed/fp32fp64_report/`'s `.cu` probes played.

## Build and run

```bash
make ARCH=89                     # RTX 4060;  ARCH=80 for A100. NO DEFAULT.
make ARCH=89 regs                # per-kernel FCHK/CALL/MUFU/LDG/STG/spill table
make ARCH=89 ptx                 # NVVM output + a decimal-.f32-literal check
./w32ptx <mode> [nx] [nrepeat] [nlaunch]      # defaults 4194304 / 1 / 10
./w32ptx accuracy                # order + deviation gate, both biases
./w32ptx --list-kernels
python3 static_w32ptx.py         # static counts + the Fortran fidelity gate
bash run_w32ptx.sh --out DIR --gpu-cc 80 --repeat 3
```

`ARCH` is a required make variable with no default, which removes the stale
`CASE_GPU_CC` trap by construction rather than catching its symptom; the binary
additionally compares its compiled arch against the device and refuses to launch
on mismatch.

## Structure

`modes.inc` is an X-macro table included three times (kernel definitions, launch
wrappers, dispatch table), so the mode→kernel map cannot drift — unlike the
Fortran side, where `../analyze_weno_static.py` keeps it by hand.

**Every rung shares one `face_body` and one `kernel_body`.** The reciprocal is a
*policy* (`RecipRef` / `RecipRcp` / `RecipFull` / `RecipRcpN`) and is the only
difference between rungs. Betas, tau, `dd`, eps, the clamp, the polynomials, the
combine and the whole epilogue are written once. That is what makes a timing
difference attributable to the reciprocal rather than to codegen luck.

| rung | reciprocal | role |
|---|---|---|
| `_ref` | IEEE `div.rn.f32` | baseline — what nvfortran `-fast` emits |
| `_rcp` | `rcp.approx.ftz` + mul for ratios, + one Newton step for the normalisation | **primary** |
| `_divfull` | `div.full.ftz.f32` | attribution control: branch-free, more arithmetic |
| `_rcpn` | Newton-refined everywhere | near-IEEE fallback |

Ablations: `w32_only_seq{5,7,9}_{ref,rcp}` (FP32 weight stream alone, combined in
FP32) gives the marginal cost of the FP64 half as
`t(w32_poly64) − t(w32_only)`; both terms sit far above the DRAM floor, which is
what makes it usable where `poly64_only` (202.8 µs against a ~200 µs floor) is
not. `weno64_seq{5,7,9}` is the accuracy reference, never a timing baseline.

## What `_rcp` actually does — the measured result

The FP32 divisions are not merely instructions. `nvfortran -fast` does not relax
FP32 division, so each `min(tau/(b+eps), cap)` is the full IEEE sequence:
`MUFU.RCP` + `FCHK` + 5 `FFMA` + a predicated `BRA` + `CALL.REL.NOINC` into a
~101-instruction out-of-line stub + `BSYNC`. That is `6 × (r+1)` = **24/30/36
CALL sites per face**, i.e. 24/30/36 basic-block boundaries, and **ptxas cannot
co-schedule across a basic-block boundary** — the same mechanism that made
deleting the solver's order-degrading boundary ladder the thing that unblocked
ILP above `ORDER=2`.

`_rcp` removes all of them. Static, WENO9, per face:

| | `_ref` | `_rcp` | Δ |
|---|---:|---:|---:|
| total | 1712 | **1160** | −32% |
| FP32-pipe (incl. MUFU) | 1011 | **828** | −183 |
| FP64-pipe (DADD+DMUL+DFMA) | 187 | 187 | 0 |
| F2F | 60 | 60 | 0 |
| FCHK / CALL / BSSY+BSYNC | 30 / 36 / 78 | **0 / 0 / ~6** | — |
| `other` | 454 | **85** | −369 |

FP64 and CONV are untouched by construction — only the FP32 half changed.
Accuracy is unchanged: `_rcp` sits at 1.02×/0.99× of `_ref`'s deviation from
FP64, and reproduces `_ref`'s L1 order floor exactly.

**This is work removal plus barrier removal, not scheduling.** `_divfull` is the
control that separates the two: branch-free like `_rcp` but retaining most of the
arithmetic (902 vs 792 FP32-pipe).

## Where it can and cannot be measured

On the **RTX 4060 (1:64)** every rung lands within ±1% (WENO9: 3338 vs 3343 µs).
That is expected, not a null result: with FP64 at 1/64 rate the 187 FP64
instructions dominate, so removing 552 instructions of FP32 and control changes
nothing. `w32_only` is no better — at 75 MB / ~200 GB/s it is memory-bound.
**The local box validates correctness and static counts only.**

The headroom is on **A100 (1:2)**, where `w32_poly64_seq9` measured 602 µs
against a 593 µs serial bound and a ~434 µs perfect-overlap floor — about 13% of
the available FP32/FP64 overlap realised. Note the ceiling is modest and already
measured: `RESULTS_df_h100.md` puts FP64 retention under a saturated FFMA stream
at **20.3%** on 1:2-class hardware (vs 94.8% on Blackwell 1:64), so pure
reordering is worth ~2%. The prize here is the −183 FP32 instructions, which the
cost model puts at **~1.15×**, more if the 36 removed barriers also help ρ.

## Fidelity to the Fortran

Asserted **exactly** by `static_w32ptx.py`, all three widths:
`DFMA+DMUL` = 55/112/169, `CONV` = 36/48/60, `LDG` = 18/24/30, `STG` = 6,
`FCHK` = 18/24/30, `CALL` = 24/30/36.

Two quantities differ for understood reasons and are reported, not gated:

- **DADD** (19/34/37 → 6/15/18). nvfortran unrolls `do k=1,nrepeat` by 4 with a
  predicated remainder, so 13/19/19 of its DADDs are unroll copies never executed
  at `nrepeat=1`.
- **FP32** (−1.4% / −3.4% / −3.7%). nvcc finds ~39 more FMA contractions in the
  dense beta forms (WENO9: FFMA +42, FMUL −27, FADD −54; FMNMX, FSETP and MUFU
  match exactly).

So `_ref` is a faithful implementation of the same **algorithm**, not a
bit-identical copy of the same **instruction sequence**, and its output differs
from the Fortran's at the FP32 weight level (~1e-7 relative) — the same order as
the split's own deviation from FP64. **A byte-exact oracle against the Fortran
was therefore not built**: it would fail by construction for a reason already
understood and quantified. The correctness basis is instead (a) the exact static
counts above, (b) `./w32ptx accuracy`, and (c) bit-identity *between rungs*,
where the beta code is compiled identically.

`checksum_all` is also **not** compiler-invariant — measured: local cc89
nvfortran 24.7 gives `2.0342306819608565E+07` at seq9 where the recorded
A100/24.3 run gives `...613412E+07`, agreeing only to ~12 digits, because
Fortran `sum()` at `-fast` is a blocked vector reduction. Use `checksum_bits`
(FNV-1a-64 over the raw bytes) for identity comparisons.

## Traps found here, each of which cost a debugging session

- **`rcp.approx.f32` is 7 instructions; `rcp.approx.ftz.f32` is 1.** The `.ftz`
  *instruction modifier* (not the `-ftz` compile flag) is what removes ptxas's
  inline predicated range-scaling for denormal inputs and overflowing results:
  measured in isolation, `FSETP.GEU + FSETP.GT + FSEL + FSEL + FMUL + MUFU.RCP +
  FMUL` versus a bare `MUFU.RCP`. At 36 reciprocals per face that is ~216
  instructions, and without it the FP32 stream fell by only 29 instead of 183 —
  i.e. the rung's entire reason for existing nearly vanished silently.
  Compiling with `-ftz=true` does **not** fix it (guards remain, total grows);
  `-ftz=false` is correct and also closer to the reference.
  Numerically inert here: every reciprocal argument is `beta+eps >= 1e-20` or
  `sum(alpha) >= 1`, far above the 1.18e-38 min normal.
- **PTX forbids operand negation on `fma`** — `fma.rn.f32 e, -%1, r0, 1;` fails
  with "Operand negation not allowed for instruction 'fma'", even though SASS
  supports it (the reference division sequence contains `FFMA R0, -R63, R64, 1`).
  Write `neg.f32` as its own instruction; ptxas folds it back into the FFMA
  source modifier, so it is free in SASS.
- **Brace-scope every multi-instruction asm body.** Without `{ }` the `.reg`
  names are redeclared on the second instantiation and ptxas fails with
  "Redefinition of variable".
- **Match SASS symbols EXACTLY, never as a substring.** `k_..._rcp` is a prefix
  of `k_..._rcpn`, so an `awk '/Function : k_..._rcp/'` silently concatenates the
  two kernels and every count comes out as their sum — LDG 36 instead of 18, STG
  12 instead of 6. This produced two entirely fictitious "findings" (a
  `div.approx.f32` body duplication and an NVVM loop-peeling bug) before being
  caught. `make regs` and `static_w32ptx.py` both compare `$NF == name`.
- **`nrepeat` does not amplify arithmetic.** NVVM hoists the loop-invariant face
  body out of the repeat loop and leaves only the ~6 DADD accumulator: measured
  slope ~3.5 µs/iteration for *every* mode regardless of its arithmetic. This is
  the nvcc counterpart of the hoisting the Fortran harness documents.
  `run_w32ptx.sh` rejects `--nrepeat`. Amplification needs a runtime-opaque index
  offset (`x[f*n + i + k*koff]`, `koff = 0` at runtime).
- **Never `-use_fast_math`.** It implies `-prec-div=false`, which turns the
  *baseline*'s `div.rn.f32` into an approximate divide — `_ref` and `_rcp` become
  the same kernel and the ladder measures nothing. `static_w32ptx.py` gates on
  `_ref` retaining its 18/24/30 `FCHK`.
- **Do not gate FP32-split accuracy against FP64 with an absolute threshold.** On
  a smooth field the dense beta quadratic forms cancel catastrophically (terms
  ~1e6 producing a beta ~1e-8), so FP32 weights carry ~1e-7 relative error *by
  construction* — `_ref` itself measures 9e-08…3e-07, which is why
  `../../report/check_weno_order.py` carries `W32_FLOOR = 3e-7`. Gate each rung
  against `_ref`, and gate order separately.

## Not done

- **`_allfma`** (force `fma.rn.f32` density) and a **`bar.warp.sync` forced-
  interleave** rung. Both are pure scheduling, and both were measured negative on
  1:2 hardware already (`RESULTS_df_h100.md`: all-FMA 51.0% vs ADD-heavy 54.1%
  retention; FFMA-saturated only 20.3%). Note also that **`asm` is opaque to
  NVVM, not to ptxas** — ptxas strength-reduces `fma(x, 1.0, y)` back to `FADD`
  *inside* an asm block, so all-FMA still needs the 1.0/0.0 constants to arrive
  as kernel arguments.
- **`dfin`-style pre-converted FP32 input.** Only 30 of the 60 F2F are removable
  (ptxas already CSEs the 8 conversions shared between the two biases at a face),
  and F2F fits at ~0.25 µs rather than 0.86, so the prize is ~7 µs against +50%
  read traffic.
- The A100 measurement itself.
