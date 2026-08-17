# CLAUDE.md — 1D_solver

Guidance for the 1D solver. Project-wide concerns live in the repository-root
`CLAUDE.md`; the WENO microbenchmark has its own `microbenchmark/CLAUDE.md`.

## What this is

A 1D compressible Navier-Stokes solver (Sod shock tube) that exists as a
**testbed for simultaneous FP32/FP64 execution**, not as a production solver.
It is the smallest thing that runs the same KEEP convective + Sutherland
viscous physics as `3D_solver`, so precision and instruction-scheduling
experiments can be read out of the SASS instead of inferred from a 500³
profile. It also carries a SLAU + MUSCL/WENO path used for the same study.

Two mechanisms are under study: **instruction-level** (one thread issues FP64
convective and FP32 viscous work — `KERNEL_MODE='fused'`) and **warp-level**
(half the block does KEEP, half viscous, combined through shared memory —
`'warp'`, and `'warp_fused'` which also moves the KEEP pressure terms to the
viscous role). Motivating measurements: `3D_solver_mixed/fp32fp64_report/`
(1:64 on Blackwell/Ada vs 1:2 on GH200 — the conclusions invert). Findings from
this directory: `report/rtx4060_weno_division_reduction.md`,
`rtx4060_weno_order_sweep.md`, `rtx4060_weno_micro_split.md`,
`cpasync_cuda_fortran.md`.

## Layout

| Path | Role |
|---|---|
| `src/main.f90` | driver (MPI single rank, kept for parity with 2D/3D) |
| `src/preprocess.f90.fypp` | pre-run setup |
| `src/calc_time_dev.f90.fypp` | TVD-RK3 stepping (`RK=3` only); owns launch geometry and passes `ng` to `set_bc` |
| `src/calc_flux_base.f90.fypp` | flux dispatcher over `KERNEL_MODE` and `SCHEME`; also `calc_quantities_shadow32` |
| `src/calc_seq_kernel.f90.fypp` | **`seq` — one thread does KEEP then viscous. The baseline to beat** |
| `src/calc_fused_kernel.f90.fypp` | `fused` — same maths hand-inlined and source-interleaved |
| `src/calc_warp_kernel.f90.fypp` | `warp` and `warp_fused` — warp-role split |
| `src/calc_keep_kernel.f90.fypp` | split-path convective kernel |
| `src/calc_keep_1d.f90.fypp` | KEEP2/4/6 device functions (`include` fragment) |
| `src/calc_keep_1d_df.f90.fypp` | double-float KEEP twin |
| `src/calc_slau_kernel.f90.fypp` | SLAU kernels: `calc_slau_x`, `calc_slau_smem_x`, `calc_slau_warp_x`, and the MUSCL/WENO double-float twins |
| `src/calc_slau_1d.f90.fypp` | SLAU1 flux device functions (`include` fragment) |
| `src/calc_visc.f90.fypp` | split-path viscous kernel (always FP64) |
| `src/calc_visc_1d.f90.fypp` | VISC2/4/6 device functions (`include` fragment) |
| `src/calc_steps.f90.fypp` | RK stage kernels; `ng`-dependent index range |
| `src/fltflt*.f90` | double-float (`'df'`) emulation type, operators and interfaces (5 files) |
| `src/print_1d.f90` | writes `Q.dat` |
| `ST/` | the only case: Sod shock tube |
| `report/` | verification scripts, benchmark/profiling drivers, findings |
| `microbenchmark/` | standalone WENO harness — see `microbenchmark/CLAUDE.md` |

`mod_constant.f90`, `calc_physical_quantities.f90`, `calc_muscl.f90.fypp` and
`calc_weno.f90` come from the repo-root `src/`, shared with 2D/3D.

## Build

**`cmake -B build` on its own picks up gfortran and fails.** `ST/CMakeLists.txt`
calls `project()` before including the parent `CMakeLists.txt` that sets the
compiler, so the compiler must be given on the command line:

```bash
cd 1D_solver/ST
cmake -B build -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/mpi/bin/mpif90
cmake --build build -j
cd build && mpirun -n 1 ./a.out     # writes Q.dat
```

Unlike the 2D cases, 1D output needs only **one** rank — there is no separate
writer rank.

If `mpirun` dies in `opal_init`, the system `/usr/bin/mpirun` is shadowing the
one the binary was linked against, and the SDK's OpenMPI has a stale baked-in
prefix. Use its own launcher with `OPAL_PREFIX` set:

```bash
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/12.5/hpcx/hpcx-2.19/ompi
/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/mpi/bin/mpirun -n 1 ./a.out
```

## config.fypp

| Variable | Values | Effect |
|---|---|---|
| `SCHEME` | `'KEEP'`, `'SLAU'` | convective scheme; `'SLAU'` constrains `KERNEL_MODE`/`ORDER`/`TVD`/`SLAU_VARIANT` (compile-time `$:error`) |
| `ORDER` | `2`, `4`, `6` | KEEP stencil width; also gates SLAU reconstruction (`2` = first-order upwind) |
| `VISC_ORDER` | `2`, `4`, `6` | viscous stencil width; defaults to `ORDER`, independent of it |
| `TVD` | `'none'`, `'tvd'` | SLAU requires `'tvd'` (minmod-limited MUSCL) |
| `SLAU_VARIANT` | `'SLAU'` | only SLAU1 is ported to 1D; HRSLAU2 needs a wiggle sensor that is not |
| `SLAU_RECON` | `'MUSCL'`, `'WENO'` | SLAU face-state reconstruction, `ORDER=6` only |
| `WENO_ORDER` | `5`, `7`, `9` | WENO-Z stencil width when `SLAU_RECON='WENO'`; decoupled from `ORDER` |
| `RECON_SPAN` | *derived* | `WENO_ORDER+1` for SLAU+WENO, else 0; feeds every template's `ng` (see below) |
| `KERNEL_MODE` | `'split'`, `'seq'`, `'fused'`, `'smem'`, `'warp'`, `'warp_fused'` | two launches / plain / hand-inlined / SLAU face states via shared memory / warp-specialised / warp-specialised with the pressure split. A typo is a build error, not a silent fall-through |
| `KEEP_PREC` | `'fp64'`, `'fp32'`, `'term'`, `'df'` | convective flux working precision |
| `VISC_PREC` | `'fp64'`, `'fp32'` | viscous precision in fused-style kernels; split stays FP64 reference |
| `PRESS_PREC` | `'fp64'`, `'fp32'` | precision of the split-off KEEP pressure terms; `'fp32'` needs `KEEP_PREC` in `('fp64','term')` and `VISC_PREC='fp32'` |
| `SLAU_MUSCL_{RHO,U,P}_PREC` | `'fp64'`, `'df'` | per-variable reconstruction precision — **applies to WENO too**, so the name is a misnomer. `'df'` needs `KERNEL_MODE='warp'` |
| `RK` | `3` | anything else is a compile-time `$:error` |
| `VISC`, `RESCALE`, `RESTART` | — | reach `mod_constant` as kind tags; **nothing in `1D_solver/src` branches on them** |

`KEEP_PREC='term'` keeps the nearest-neighbour terms in FP64 and drops the
wide-stencil corrections to FP32, reading them from FP32 mirrors of the shared
tile. Needs `ORDER >= 4`. `KEEP_PREC /= VISC_PREC` is what puts work on both
pipes at once.

**FP32 shadow arrays.** When a fused-style kernel needs FP32 copies of its
inputs (`VISC_PREC='fp32'`, `PRESS_PREC='fp32'`, `KEEP_PREC='term'`), the
FP64→FP32 conversions are fused into the per-stage primitives decode by
`calc_quantities_shadow32` (`src/calc_flux_base.f90.fypp`), which reproduces the
shared `calc_quantities_T_1D` bit-for-bit and additionally writes whichever of
`u32`/`t32`/`mu32`/`rho32`/`p32` the config needs. When none are needed the
shared routine is called unmodified — that is what keeps zero-shadow configs
bit-identical. Both `F2F` directions issue on the FP64 pipe, so converting once
here beats converting at every tile load.

**`VISC_ORDER` is a stencil WIDTH knob, not an order-of-accuracy knob.** The
face flux is genuinely 2nd/4th/6th order *at the face*, but `calc_R` turns it
into a divergence with a one-cell difference carrying its own `(h²/24)·F'''`
term, so `d/dx` of the viscous flux is **2nd order at every width** —
`report/check_visc_order.py` measures both. What `VISC_ORDER` is *for* here is
growing the FP32 instruction count.

## Ghost cells, not an order-degrading ladder

`set_bc` fills `ng` ghost cells at each end (zero-gradient), making one uniform
interior stencil valid at every face. There is **no boundary branch in any flux
kernel**.

| | range |
|---|---|
| ghost cells (`set_bc` writes) | `1 … ng` and `nx-ng+1 … nx` |
| physical cells (RK updates) | `ng+1 … nx-ng` |
| faces the flux kernels compute | `ng … nx-ng` |
| cells the stencil reads at face `f` | `f-io … f+io+1` |

**`ng` is derived identically in eight templates** —
`calc_{keep,seq,fused,warp,slau}_kernel`, `calc_steps`, `calc_time_dev`,
`calc_visc` — all as `max(ORDER, VISC_ORDER, RECON_SPAN) // 2`, with the policy
input `RECON_SPAN` in `config.fypp` (0 for every non-SLAU/WENO build, which
keeps those bit-identical). Widening `ng` in one kernel alone would make that
kernel read uninitialised cells, so any change must land in all eight. `io` is
**per-kernel**, not universally `ng-1`: `calc_slau_kernel` uses
`(WENO_ORDER-1)//2` for WENO and `ORDER//2-1` for MUSCL, `calc_keep_kernel`
uses `ORDER//2-1`, and only `seq`/`fused`/`warp` use `ng-1`.

Deleting the old order-degrading boundary ladder is **what makes the ILP
mechanism work above ORDER=2** — ptxas cannot co-schedule across a basic-block
boundary, and the ladder put the FP32 and FP64 halves in different blocks.

Cheapest regression check available: `ORDER=VISC_ORDER=2` must reproduce a
pre-change `Q.dat` **bit-for-bit**, because `ng=1` degenerates to the original
index ranges. Capture a baseline before touching index arithmetic.

## Higher-order KEEP and viscous

`src/calc_keep_1d.f90.fypp` and `src/calc_visc_1d.f90.fypp` generate KEEP2/4/6
and VISC2/4/6 as 1D specialisations of the repo-root
`calc_scheme_math.f90.fypp` and `3D_solver/src/calc_visc_cent.f90.fypp`. In 1D
the transverse velocities and the `Normal` argument collapse away, so the
shared files are not reused directly. Order is selected by the **KIND of an
unused leading argument** (`integer(2)/(4)/(8)`) through a generic interface;
only the configured order is generated, keeping dead bodies out of the SASS
listing.

Two deliberate deviations from the 3D source: pair sums/differences instead of
a serial `fma()` accumulator (a serial chain has no internal ILP, and
independent FP32 instructions are the whole point here), and exact fractions
evaluated at the working kind (`real(125.d0/1920.d0, VK)`) instead of truncated
decimals, which lose digits at `real(4)`.

## WENO-Z

`src/calc_weno.f90` (repo root) holds WENO5-Z (`delta6_weno`), WENO7-Z
(`delta8_weno`) and WENO9-Z (`delta10_weno`), selected by `WENO_ORDER`. The
WENO7/9 bodies and all the double-float twins are **generated, not
hand-written**:

```bash
python3 report/check_weno_order.py             # verify — run after ANY change
python3 report/check_weno_order.py --gen 4     # regenerate weno7z_* bodies
python3 report/check_weno_order.py --gen 5     # regenerate weno9z_* bodies
python3 report/check_weno_order.py --gen-df 4  # regenerate the WENO7-Z DF twin
python3 report/check_weno_order.py --gen-df 5  # regenerate the WENO9-Z DF twin
```

`check_weno_order.py` derives every coefficient from its defining conditions in
**exact rational arithmetic** and then measures the convergence rate on a smooth
solution (observed 5.00 / 7.02 / 8.95). **`check_vs_numpy.py` does not guard
this** — it compares the Fortran against a numpy transcription of the *same*
formulas, so it checks that two implementations agree, not that the scheme
achieves its design order. That gap hid a real bug: `weno5z_right`'s optimal
weights were copied from `weno5z_left` unmirrored, making `v^+` **third order
instead of fifth**, with both implementations carrying it so the check reported
PASS throughout. `v^+` mirrors `d`: near-face candidate 3/10, furthest 1/10.

At a flat 4 divisions per call, widening grows the smoothness-indicator work
quadratically — 131 → 238 → 366 FP64-pipe instructions per face for WENO5/7/9,
`io`/`ng` = 2/3, 3/4, 4/5. That is the intended way to give the co-issue
experiment two halves large enough to overlap:
`report/rtx4060_weno_order_sweep.md`.

## Verification

`report/` holds five checks; run all after touching a flux stencil.

```bash
# 1. coefficients are genuinely the claimed order (pure numpy, no GPU)
python3 report/check_keep_order.py
python3 report/check_visc_order.py
python3 report/check_weno_order.py

# 2. the Fortran on the GPU computes those same coefficients
#    (set nt=1 in ST/mod_globals.f90, rebuild, run, then:)
python3 report/check_vs_numpy.py ST/build/Q.dat <ORDER> [VISC_ORDER]
#    --press fp32                                for PRESS_PREC='fp32' builds
#    --scheme SLAU --recon WENO --weno-order N   for SLAU/WENO builds
#    No KEEP_PREC='fp32'/'term' model exists — use check 3 for those.

# 3. the full run matches the exact Riemann solution (nt=2500)
python3 report/check_sod.py ST/build/Q.dat [reference_Q.dat]
```

Check 1 is the one that catches a mistyped coefficient; **a shock tube cannot**,
because every scheme is first order at a discontinuity. For the same reason
check 3's density L1 barely moves with order (0.185% → 0.182%) — do not read it
as an order test. `check_visc_order.py` prints two tables and the second is
*supposed* to show 2/2/2 (see the `VISC_ORDER` note above).

Benchmark with `report/bench.sh` or the `run_*_ncu.sh` / `run_*_nsys.sh`
drivers, never with the stock `nx=4096` — see Notice. SASS helpers:
`report/sass_analyze.py`, `report/ctrl_bits.py`.

## Rules established by measurement — do not re-derive

- **FP64-pipe instruction count (`DADD+DMUL+DFMA+F2F`) is the static metric that
  tracks time.** The SASS interleaving `runs` count does **not** predict time
  and must never be an optimisation target.
- **Do not add `volatile` scalars or `x - x` zero-dependencies to force
  FP32/FP64 alternation.** Tried: `DADD` and both `F2F` directions issue on the
  FP64 unit, so it added +64% FP64 work, and chaining the halves destroys the
  independence co-issue needs — 9–28% slower.
- **Relocating work is free; removing it is what pays.** Warp specialisation at
  matched precision is worth ~1.00×. `seq`/`fused` with `PRESS_PREC='fp32'`
  match or beat `warp_fused` with the same demotion, so the win is the precision
  demotion, not the role split.
- **FP64 division is never free — `grep -c MUFU.RCP64H` on a `cuobjdump -sass`
  listing.** `div.rn.f64` expands to `MUFU.RCP64H` + ~9 `DFMA` **plus a
  `CALL`/`RET`** (~20 instructions and a basic-block boundary), and `-fast` does
  **not** fold division by a literal: `x/6.0d0` emits `MUFU.RCP64H R9, 6`.
  Removing them is the only lever here that is not GPU-dependent — it pays on
  A100/GH200 as much as on Ada. `report/rtx4060_weno_division_reduction.md`.
- **Batched inversion (`1/(c0*c1*c2)`) is FP64-only.** In `fltflt` the product
  underflows — every `c_i` is `eps=1e-20` on a plateau, `1e-60` flushes to zero
  in FP32 exponent range, giving `1/0 = Inf` then `0*Inf = NaN`. The DF twins
  therefore take only the polynomial-scaling rewrite and are deliberately *not*
  expression-identical to the FP64 path; both files say so.
- **`cp.async` lives in the `wmma` module** (`pipelineMemcpyAsync` /
  `pipelineCommit` / `pipelineWaitPrior`), not `cudadevice` — whose
  `__pgi_memcpy_async*` are the *host* `cudaMemcpyAsync` specifics and emit no
  `LDGSTS`. It is **load-only**: `cp.async` is hardware-restricted to
  global→shared, so a register or `shared` source compiles fine and then fails
  at launch with **CUDA 717**. `report/cpasync_cuda_fortran.md`.

Four `report/*.md` files cited by earlier revisions of this document
(`SASS_ilp.md`, `rtx4060_variant_optimization.md`,
`rtx4060_warp_fused_optimization.md`, `rtx4060_precision_vs_specialization.md`)
**were never committed**. The conclusions above survive from those sweeps; the
backing CSVs do not. Do not go looking for them.

## Notice

* `nx=4096` is far too small to measure throughput: ~0.1 waves/SM, ~11%
  achieved occupancy, launch-latency bound. Use `report/bench.sh` (default
  `nx=4194304`) for any performance claim; keep `nx=4096` for SASS and
  correctness only.
* **GPU clocks cannot be locked without root**, so long sweeps throttle: mean
  kernel times drift >2× while pipe-utilisation ratios stay put. Compare
  `time_min_us` from interleaved, repeated rounds.
* **Always pass `-DCASE_GPU_CC=<cc>`.** It defaults to `ccnative`, so the target
  architecture otherwise follows the build machine — and that matters here,
  because FP64:FP32 is 1:64 on consumer Ada/Blackwell and 1:2 on GH200, with the
  mixed-precision conclusions inverting between them. A stale cached value also
  produces a clean build that then fails at every launch with
  `cudaErrorInvalidPtx`.
* `maxregcount:96` (raised from 64), set in `1D_solver/CMakeLists.txt`. Check
  `ptxas` output for spills after any change that adds live values.
  `calc_slau_x` currently uses 60 registers with zero spill at every
  `WENO_ORDER`, so the cap is not presently binding.
* **`ncu` and `nsys` hang on the solver binary on some machines** while the
  unprofiled binary runs fine. If that happens, fall back to end-to-end wall
  clock, and kill leftover profiler processes **by PID** — `pkill -f Nsight`
  also matches the shell running the `pkill` and will kill your own script.
* A crashed run leaves the previous `Q.dat` in place. Always `rm -f Q.dat`
  before a run in a sweep script, or a crash is silently reported as a pass.
* `nt` in `ST/mod_globals.f90` is switched between `1` (SASS/profiling,
  `check_vs_numpy.py`) and `2500` (physics, `check_sod.py`). Check which is set
  before interpreting `Q.dat`.
