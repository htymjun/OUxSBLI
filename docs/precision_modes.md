# Precision modes of KEEP: FP64, DF, difference-first mixed

Ports of the diffFP32 study into `3D_solver` (configuration: `PREC`, `KEEP_PREC`,
`WENO_PREC='MIXED'`; implementation notes: "Precision modes" in `3D_solver/CLAUDE.md`).

Machine: NVIDIA RTX PRO 6000 Blackwell Max-Q (sm_120, FP64:FP32 = 1:64), driver 595.58.03,
NVIDIA HPC SDK 26.5, 2 MPI ranks (one compute rank). Branch `fp32fp64-dfmix`.

## Correctness against FP64

DF and the mixed forms against FP64 on the same build. Metric: max |q - q_FP64| / max |q_FP64| of the
last snapshot (`OUTPUT_PRECISION=8`). DF (`PREC='DF'`), 40 steps (SWLBLI 2 steps, analytic
start; NSTGV at 65^3):

| case | features | rho | p | velocity |
|---|---|---|---|---|
| ETGV | Euler KEEP6, RK4, periodic | 5.5e-14 | 2.9e-14 | 4.0e-14 |
| EVC | Euler KEEP6 | 3.8e-14 | 3.9e-14 | 1.4e-12 |
| KHI | Euler KEEP6, RK3 | 1.9e-13 | 3.8e-13 | 9.9e-14 |
| DHIT | NS fused KEEP6 | 8.6e-14 | 6.6e-14 | 4.4e-14 |
| IVST | NS, walls | 8.7e-14 | 1.6e-13 | 5.4e-14 |
| STZ | NS KEEP4, COMMZ (MPI='GPU') | 1.2e-15 | 4.2e-15 | (u = 0) |
| TBL | NS KEEP6, RESCALE, walls | 2.0e-13 | 2.3e-13 | 1.4e-13 |
| NSTGV, 1 step | NS fused SLAU (plain SLAU) | 1.1e-14 | 9.1e-15 | 1.1e-14 |
| BL, 1 step | NS SLAU, walls (plain SLAU) | 7.4e-15 | 8.2e-15 | 2.4e-14 |
| OS, 1 step | Euler SLAU, shock (plain SLAU) | 2.3e-15 | 9.8e-15 | 6.2e-15 |
| SWTBLI, 1 step | two domains, RESCALE (plain SLAU) | 1.1e-14 | 3.8e-14 | 1.4e-14 |
| SWLBLI | NS SLAU, BC_FORCING (plain SLAU) | 2.1e-14 | 2.8e-14 | 1.3e-14 |
| DHIT, 1 step | LES | 9.3e-15 | 1.4e-14 | 1.2e-14 |
| EVC, 1 step | Hybrid / WENO+LF | 4.8e-15 / 5.4e-15 | 6.3e-15 / 7.1e-15 | 1.2e-13 / 7.5e-14 |
| NSTGV, 1 step | VISC_ORDER=4 (separate viscous kernels) | 1.1e-14 | 9.1e-15 | 1.1e-14 |

With the configured HR-SLAU2, NSTGV / BL / OS / SWTBLI differ at 1e-9..1e-3: its wiggle
detector is a 0/1 switch on the sign of products of pressure differences, which are exactly 0
in FP64 at symmetry planes and in uniform regions and +-1e-30 under DF. FP64 default builds:
output bit-identical to d911d5a in every case (SWTBLI's second domain is nondeterministic
already in d911d5a); 2D_solver and 3D_solver_curv SASS identical.

Mixed forms (FP64 state), ETGV 66^3 KEEP6, 200 steps: MIXED 1.1e-8 (rho), MIXED_CELL 1.5e-8.
KEEP6 face flux alone (unit test, 6-point stencil with 1 % variations): MIXED 4e-10 relative.

## 2D_solver: DF and FP32 (`PREC = 'DF' | 'FP32'`)

Machine: RTX 4060 Laptop (sm_89), WSL2, NVIDIA HPC SDK 25.5, 2 MPI ranks. Each configuration
runs 40 steps (`prec2d` runs; nx, ny of the case defaults); metric as above. `slau1` =
`SLAU_VARIANT='SLAU'`. FP64 output of all 21 configurations is bit-identical to the output before the
precision work (`recal/Q00001.dat` md5).

| case | features | DF rho | DF p | DF velocity | FP32 rho | FP32 p | FP32 velocity |
|---|---|---|---|---|---|---|---|
| EVC | Euler KEEP2 / KEEP4 / KEEP6 (RK4) | 2.9e-14 / 3.3e-14 / 3.3e-14 | 2.0e-14 / 3.0e-14 / 3.3e-14 | 2.4e-13 / 1.0e-12 / 1.4e-12 | 8.7e-7 / 9.2e-7 / 8.8e-7 | 9.8e-7 / 8.4e-7 / 8.9e-7 | 1.7e-5 / 2.9e-5 / 2.8e-5 |
| EVC | KEEP6, RK3 | 2.2e-13 | 7.6e-14 | 1.5e-12 | 2.0e-6 | 1.9e-6 | 2.6e-5 |
| EVC | SLAU slau1 / Hybrid / WENO+LF | 2.3e-14 / 2.9e-14 / 2.6e-14 | 1.1e-14 / 3.3e-14 / 1.4e-14 | 1.8e-13 / 1.3e-12 / 1.7e-13 | 6.6e-7 / 8.0e-7 / 3.4e-7 | 2.3e-7 / 8.3e-7 / 3.0e-7 | 1.4e-5 / 2.8e-5 / 6.4e-6 |
| DSL | Euler KEEP6, periodic | 5.2e-14 | 2.1e-14 | 1.0e-13 | 1.3e-6 | 1.3e-6 | 3.5e-6 |
| ST | NS KEEP6, walls, `calc_div` | 6.8e-14 | 7.1e-14 | 2.2e-14 | 1.8e-6 | 2.5e-6 | 1.9e-6 |
| ST | `VISC_ORDER=2` / RK4 | 6.9e-14 / 1.8e-14 | 8.3e-14 / 1.7e-14 | 2.2e-14 / 1.2e-14 | 1.4e-6 / 5.4e-7 | 2.2e-6 / 6.5e-7 | 9.7e-7 / 3.0e-7 |
| OS | Euler SLAU slau1, shock, TVD hybrid | 4.8e-12 | 6.0e-12 | 1.7e-12 | 1.1e-5 | 1.4e-5 | 8.9e-6 |
| SBLI | NS SLAU slau1, shock / BL | 3.2e-12 | 3.3e-12 | 1.7e-12 | 6.0e-6 | 8.5e-6 | 3.3e-6 |
| SBLI | NS KEEP+SLAU (configured) | 2.3e-13 | 2.8e-13 | 2.0e-13 | 5.4e-6 | 7.7e-6 | 4.8e-6 |
| BL | NS SLAU slau1, constant-p far field | 1.2e-13 | 6.7e-14 | 1.8e-13 | 2.2e-6 | 1.4e-6 | 4.9e-6 |
| SWLBLI | NS SLAU slau1, 2700 x 480 | 1.4e-13 | 2.2e-13 | 1.1e-13 | 2.0e-6 | 4.1e-6 | 3.1e-6 |

The two cases above 1e-12 under DF (OS, SBLI slau1) exceed it at 2-4 points (OS near the top
boundary, SBLI at the leading edge); the median difference is 1e-15. With the configured HR-SLAU2,
DF / FP32 differ by 1e-9..1e-6 (OS, BL) and 1e-3 (SWLBLI, after 40 steps from the analytic
start), for the wiggle-detector reason given above; EVC with HR-SLAU2 5e-11 (velocity) under DF.

TVD-RK3 stages 2 and 3 (`calc_step2_3`) are evaluated in increment form,
Q^(n+1) = Q^(*) + (alpha (Q^n - Q^(*)) - gamma R) / (alpha + beta); the former
(alpha Q^n + beta Q^(*) - gamma R) / (alpha + beta) multiplied the whole state by the rounded
1/3. Uniform gas (EVC with beta = 0, 150 steps): in FP32 the old form moved rho by 1 ulp within
the first 50 steps, the increment form keeps it unchanged. FP64: RK4 configurations unchanged,
RK3 configurations differ by 3e-15..9e-14 (HR-SLAU2: BL 2e-6, SWLBLI 8e-4).
Existing suites with `OUXSBLI_CONFIG=prec=df`: test_st, test_os, test_bl, test_evc 9 passed.
With `prec=fp32`, test_evc (KEEP4/6 orders; the vortex's density amplitude 2e-7 is at FP32
rounding) and test_bl's velocity profile fail; their tolerances are left as they are.

FP32 SASS (`cuobjdump -sass`, `bench_prec.FP64_OP`): no FP64 instruction in any kernel of any
of the configurations above. 3D_solver TBL (RESCALE, 65 x 49 x 22, 20 steps; `calc_rescale`'s
`calc_mean` / `copy` and the TBL boundary conditions as hand-written kernels): FP64 output
bit-identical to before; DF vs FP64 rho 5.7e-14, p 3.5e-13, velocity 1.1e-13. With HPC SDK 25.5,
the 3D cases whose `set.f90` keeps state arithmetic in `!$cuf` loops (e.g. SWTBLI, BL, OS, STZ)
still do not build with `PREC='DF'`; the table at the top was measured with HPC SDK 26.5.

## diffFP32's rules (`ouxsbli/tests/test_prec_difffp32.py`, slow)

3D EVC (Euler vortex, density amplitude ~2e-7, one period), L2(rho) at n = 64/128/256:

| KEEP | mode | n=64 | n=128 | n=256 | order | max / FP64 | rule |
|---|---|---|---|---|---|---|---|
| 2 | FP64 | 5.078e-07 | 1.056e-07 | 3.191e-08 | 2.00 | 1.00 | |
| 2 | DF, MIXED, MIXED_CELL | 5.078e-07 | 1.056e-07 | 3.191e-08 | 2.00 | 1.00 | pass |
| 4 | FP64 | 3.210e-08 | 2.533e-09 | 1.994e-10 | 3.67 | 1.00 | |
| 4 | DF | 3.210e-08 | 2.533e-09 | 1.993e-10 | 3.67 | 1.00 | pass |
| 4 | MIXED | 3.210e-08 | 2.543e-09 | 3.086e-10 | 3.35 | 1.55 | fail |
| 4 | MIXED_CELL | 3.210e-08 | 2.541e-09 | 3.075e-10 | 3.35 | 1.54 | fail |
| 6 | FP64 | 3.782e-09 | 1.478e-10 | 1.375e-12 | 5.71 | 1.00 | |
| 6 | DF | 3.782e-09 | 1.478e-10 | 1.384e-12 | 5.71 | 1.01 | pass |
| 6 | MIXED | 3.774e-09 | 3.406e-10 | 2.601e-10 | 1.93 | 189 | fail |
| 6 | MIXED_CELL | 3.755e-09 | 2.936e-10 | 3.040e-10 | 1.81 | 221 | fail |

Decision: `KEEP_PREC='MIXED'` (face form) is kept; the cell form (`MIXED_CELL`, rows above
and below) was removed after this comparison -- same density floor, 5x larger velocity
error, and not conservative.

Where the mixed floor comes from (EVC KEEP6, n = 256):

* The algebra is exact: `type(dev)` with real(8) components reproduces FP64 (L2 1.3754e-12).
  The floor is FP32 rounding of the deviation terms.
* It accumulates over the run: L2(rho_MIXED - rho_FP64) after 1/4, 1/2, 3/4, 1 period is
  0.86, 1.46, 2.07, 2.54 x 1e-10 (MIXED_CELL 0.96, 1.72, 2.40, 2.97 x 1e-10). Each step adds
  ~eps32 x |one direction's flux difference|; in this low-Mach vortex the x/y/z flux
  differences are ~1e-2 relative and nearly cancel (div u ~ 0), so the rounding of each is
  large against the net change. diffFP32's entropy wave has no such cancellation.
* So the limit is the error level, not the order: KEEP2 (3e-8 at n = 256) never sees it,
  KEEP4 reaches it at n = 256, KEEP6 at n = 128.
* Face vs cell form: same density floor; velocity error after one period 1.1e-7 (MIXED)
  vs 5.5e-7 (MIXED_CELL). MIXED is conservative and runs in every KEEP kernel (walls, COMMZ,
  fused NS); MIXED_CELL is Euler/periodic only and 8-13 % cheaper (KEEP6 0.92 vs 1.01 of
  FP64). Recommendation: keep MIXED (face form).

ETGV 66^3, RK4, dt = 0.01, 20000 steps (t = 200), KEEP2 and KEEP6: no mode leaves
test_etgv's bounds (ke/ke0 <= 1.01, |entropy| < 3e-3); FP64, DF and the mixed forms are all
inside them over the whole run. max |ke/ke0 - FP64| over the run: KEEP2 DF 6.1e-3, MIXED
5.2e-3, MIXED_CELL 5.8e-3; KEEP6 DF 2.2e-3, MIXED 6.2e-3, MIXED_CELL 5.9e-3 (round-off
growth of the inviscid TGV; DF also departs from FP64, as in diffFP32).

Existing suites in DF (`OUXSBLI_CONFIG=prec=df`, commit 4ff0ee0): test_etgv 1 passed,
test_evc_3d 2 passed (KEEP4/6 orders), test_os_3d 1 passed (Rankine-Hugoniot), test_bl_3d 6
passed (Blasius), test_nstgv_ke_eps 2 passed (KEEP and SLAU budgets); test_tbl skipped (its
reference data is not in this checkout).
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
