# Precision modes of KEEP: FP64, DF, difference-first mixed

Ports of the diffFP32 study into `3D_solver` (configuration: `PREC`, `KEEP_PREC`,
`WENO_PREC='MIXED'`; implementation notes: "Precision modes" in `3D_solver/CLAUDE.md`).

Machine: NVIDIA RTX PRO 6000 Blackwell Max-Q (sm_120, FP64:FP32 = 1:64), driver 595.58.03,
NVIDIA HPC SDK 26.5, 2 MPI ranks (one compute rank). Branch `fp32fp64-dfmix`.

## Correctness against FP64

DF and the mixed forms against FP64 on the same build. Metric: max |q - q_FP64| / max |q_FP64| of the
last snapshot (`OUTPUT_PRECISION=8`).

| case | steps | DF | MIXED | MIXED_CELL |
|---|---|---|---|---|
| ETGV 66^3 KEEP6 RK4, dt 0.01 | 200 | 1.2e-13 (rho) | 1.1e-8 | 1.5e-8 |
| ETGV 66^3 KEEP2 RK4 | 40 | 2.9e-14 | | |
| KHI 130^3 KEEP6 RK3 | 40 | 3.8e-13 | | |
| EVC KEEP4 RK4 | 40 | 1.0e-12 (velocity) | | |

KEEP6 face flux alone (unit test, 6-point stencil with 1 % variations): MIXED 4e-10 relative,
DF 6e-15. Default FP64 builds: ETGV output bit-identical to d911d5a.

## diffFP32's rules (`ouxsbli/tests/test_prec_difffp32.py`, slow)

3D EVC (Euler vortex, density amplitude ~2e-7, one period), L2(rho) at n = 64/128/256:

| KEEP | mode | n=64 | n=128 | n=256 | order | max / FP64 | rule |
|---|---|---|---|---|---|---|---|
| 4 | FP64 | 3.210e-08 | 2.533e-09 | 1.994e-10 | 3.67 | 1.00 | |
| 4 | DF | 3.210e-08 | 2.533e-09 | 1.993e-10 | 3.67 | 1.00 | pass |
| 4 | MIXED | 3.210e-08 | 2.543e-09 | 3.086e-10 | 3.35 | 1.55 | fail |
| 4 | MIXED_CELL | 3.210e-08 | 2.541e-09 | 3.075e-10 | 3.35 | 1.54 | fail |
| 6 | FP64 | 3.782e-09 | 1.478e-10 | 1.375e-12 | 5.71 | 1.00 | |
| 6 | DF | 3.782e-09 | 1.478e-10 | 1.384e-12 | 5.71 | 1.01 | pass |
| 6 | MIXED | 3.774e-09 | 3.406e-10 | 2.601e-10 | 1.93 | 189 | fail |
| 6 | MIXED_CELL | 3.755e-09 | 2.936e-10 | 3.040e-10 | 1.81 | 221 | fail |

The mixed forms floor at L2(rho) ~ 3e-10 (2.5e-10 relative): FP32 rounding of the deviation
terms, the same mechanism as diffFP32's WENO/TENO floor. diffFP32's KEEP2 mixed matched FP64,
most likely because its 2nd-order errors (>= 2e-10 at n = 256) stayed above such a floor;
KEEP4/6 reach it at n = 128-256. The face form (conservative) and the cell form floor alike.

ETGV 66^3, RK4, dt = 0.01, 20000 steps (t = 200), KEEP2 and KEEP6: no mode leaves
test_etgv's bounds (ke/ke0 <= 1.01, |entropy| < 3e-3); FP64, DF and the mixed forms are all
inside them over the whole run. max |ke/ke0 - FP64| over the run: KEEP2 DF 6.1e-3, MIXED
5.2e-3, MIXED_CELL 5.8e-3; KEEP6 DF 2.2e-3, MIXED 6.2e-3, MIXED_CELL 5.9e-3 (round-off
growth of the inviscid TGV; DF also departs from FP64, as in diffFP32).

Existing suites: `pytest -m "not slow"` 33 passed / 14 skipped (same as d911d5a);
`OUXSBLI_CONFIG=keep_prec=mixed|keep_prec=mixed_cell|prec=df pytest test_etgv.py` pass;
`OUXSBLI_CONFIG=keep_prec=mixed pytest test_nstgv_ke_eps.py` (fused KEEP+NS) passes;
`OUXSBLI_CONFIG=weno_prec=mixed pytest test_weno.py -m "not slow"` passes (EVC WENO+LF order
>= 3.0, OS Rankine-Hugoniot). WENO reconstruction alone (`weno_order.cuf`, O(1) data): MIXED
floors at 2.7e-10 (FP64 5e-13 at N = 1024), W32P64 follows FP64. diffFP32's
`tools/check_identities.py`: all identities hold.

## Cost (`python -m ouxsbli.analysis.bench_prec --n 256 --steps 20`)

ETGV 256^3, RK4, 20 steps; GPU time per step [ms] from nsys. "ops" is the SASS census of
`calc_keep_x_in` (FP64 pipe: D* and F2F/I2F with an F64 operand; FP32: FADD/FMUL/FFMA/FSETP/
FMNMX/FSEL/MUFU). The models scale the census by an FP64:FP32 throughput of 1:64 (this GPU)
and 1:2 (GH200/H100); they ignore memory traffic.

| KEEP | mode | KEEP kernels | all kernels | KEEP / FP64 | FP64 ops | FP32 ops | 1:64 model | 1:2 model |
|---|---|---|---|---|---|---|---|---|
| 2 | FP64 | 12.202 | 29.578 | 1.00 | 27 | 8 | 1.00 | 1.00 |
| 2 | DF | 12.467 | 30.024 | 1.02 | 12 | 447 | 0.70 | 7.60 |
| 2 | MIXED | 17.966 | 35.257 | 1.47 | 57 | 71 | 2.14 | 2.98 |
| 2 | MIXED_CELL | 15.669 | 32.912 | 1.28 | 50 | 125 | 1.92 | 3.63 |
| 6 | FP64 | 30.305 | 47.931 | 1.00 | 113 | 8 | 1.00 | 1.00 |
| 6 | DF | 22.164 | 40.125 | 0.73 | 12 | 2003 | 0.38 | 8.66 |
| 6 | MIXED | 30.552 | 48.175 | 1.01 | 111 | 285 | 1.02 | 2.17 |
| 6 | MIXED_CELL | 27.755 | 45.327 | 0.92 | 98 | 467 | 0.93 | 2.83 |

Reading:

* On this GPU DF is the fastest KEEP6 (0.73x), as in diffFP32 (KEEP2 0.44x there; here KEEP2
  is memory bound at 1.0x).
* The mixed forms do not remove FP64-pipe work for KEEP6: each face takes 6-point FP64
  differences of 7 variables and narrows them (~90 DADD/F2F), as many FP64-pipe instructions
  as the FP64 KEEP6 product sums themselves (113). The deviation arithmetic then costs 2.5-4x
  the FP32 instructions of a plain FP32 KEEP.
* 1:2 model: native FP64 is cheapest; MIXED (2.2x) and MIXED_CELL (2.8x) beat DF (8.7x), so
  the expected order mixed < DF holds there, but neither is expected to beat FP64. To be
  confirmed on Miyabi (GH200) with the same command.
* Next lever for the mixed forms: share the narrowed differences between the faces of a tile
  (one DADD+F2F per cell and variable instead of per face and stencil point), which would
  cut the KEEP6 FP64-pipe count to ~40.
