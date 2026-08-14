# CLAUDE.md — 1D_solver

Guidance for the 1D solver. Project-wide concerns live in the repository-root
`CLAUDE.md`.

## What this is

A 1D compressible Navier-Stokes solver (Sod shock tube) that exists as a
**testbed for simultaneous FP32/FP64 execution**, not as a production solver.
It is the smallest thing that runs the same KEEP convective + Sutherland
viscous physics as `3D_solver`, so precision and instruction-scheduling
experiments can be read directly out of the SASS instead of being inferred from
a 500³ profile.

Two mechanisms are under study:

- **Instruction-level** — one thread issues FP64 convective and FP32 viscous
  work; both pipes are busy at once. `src/calc_fused_kernel.f90.fypp`
  (`KERNEL_MODE='fused'`) is hand-inlined so KEEP and viscous arithmetic can be
  source-ordered for SASS experiments.
- **Warp-level** — half of each block computes KEEP and the other half computes
  the viscous flux, then the two halves combine through shared memory. This is
  `src/calc_warp_kernel.f90.fypp` (`KERNEL_MODE='warp'`). Its extension
  `'warp_fused'` also moves the explicit KEEP pressure terms (`KEEPP`) to the
  viscous-role warps; with `PRESS_PREC='fp32'` those terms run on the FP32 pipe,
  which is the only way the role split can *reduce* FP64-pipe work rather than
  merely relocate it.

Background measurements that motivate the whole exercise are in
`3D_solver_mixed/fp32fp64_report/` (Blackwell/Ada at FP64:FP32 = 1:64 vs GH200
at 1:2 — the conclusions are opposite between the two).

## Layout

| Path | Role |
|---|---|
| `src/main.f90` | driver (MPI single rank, kept for structural parity with 2D/3D) |
| `src/calc_time_dev.f90.fypp` | TVD-RK3 time stepping (`RK=3` only) |
| `src/calc_flux_base.f90.fypp` | flux dispatcher: `split` (two launches), `seq`, `fused`, `warp` |
| `src/calc_seq_kernel.f90.fypp` | **`seq` — one thread does KEEP then viscous. The baseline, and the one to beat** |
| `src/calc_fused_kernel.f90.fypp` | `fused` — same maths hand-inlined and source-interleaved |
| `src/calc_warp_kernel.f90.fypp` | `warp` — warps 0–3 do KEEP, warps 4–7 do viscous for the same faces |
| `src/calc_warp_kernel.f90.fypp` | `warp_fused` — warp mode with the explicit KEEP pressure terms (`KEEPNP`/`KEEPP` split) on the viscous role; `PRESS_PREC` selects their precision |
| `src/calc_keep_kernel.f90.fypp` | split-path convective kernel |
| `src/calc_keep_1d.f90.fypp` | KEEP2/KEEP4/KEEP6 device functions (`include` fragment) |
| `src/calc_visc.f90.fypp` | split-path viscous kernel (always FP64) |
| `src/calc_visc_1d.f90.fypp` | VISC2/VISC4/VISC6 device functions (`include` fragment) |
| `src/calc_steps.f90.fypp` | RK stage kernels; order-independent arithmetic, `ng`-dependent index range |
| `src/print_1d.f90` | writes `Q.dat` |
| `ST/` | the only case: Sod shock tube, `nx=4096` |
| `report/` | SASS analysis, verification scripts, findings |

`mod_constant.f90` and `calc_physical_quantities.f90` come from the repo-root
`src/`, shared with 2D/3D.

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

Only these actually change generated code in 1D. `VISC`, `SCHEME`, `TVD`,
`SLAU_VARIANT`, `RESCALE`, `RESTART` reach `mod_constant.f90.fypp` as kind tags
but **nothing in `1D_solver/src` branches on them** — there is no Euler/LES path
and no SLAU kernel here.

| Variable | Values | Effect |
|---|---|---|
| `ORDER` | `2`, `4`, `6` | KEEP stencil width |
| `VISC_ORDER` | `2`, `4`, `6` | viscous stencil width; defaults to `ORDER`, independent of it |
| `KERNEL_MODE` | `'split'`, `'seq'`, `'fused'`, `'warp'`, `'warp_fused'` | two launches / one plain kernel / hand-inlined / warp-specialised / warp-specialised with the pressure split. Validated — a typo is now a build error, not a silent fall-through to `split` |
| `KEEP_PREC` | `'fp64'`, `'fp32'`, `'term'` | working precision of the convective flux |
| `VISC_PREC` | `'fp64'`, `'fp32'` | working precision of the viscous flux in `fused`/`warp`; split stays FP64 reference |
| `PRESS_PREC` | `'fp64'`, `'fp32'` | precision of the split-off KEEP pressure terms (`KEEPNP`/`KEEPP`), available on `seq`/`fused`/`warp_fused`. `'fp32'` requires `KEEP_PREC` in `('fp64','term')` and `VISC_PREC='fp32'`; any other combination (incl. `KERNEL_MODE` in `('split','warp')`) is a compile-time `$:error` |
| `RK` | `3` | anything else is a compile-time `$:error` |

`KEEP_PREC='term'` keeps the leading (nearest-neighbour, l=1) terms in FP64 and
drops the wide-stencil corrections (l≥2) to FP32, reading them from FP32 mirrors
of the shared tile. It needs `ORDER >= 4` — KEEP2 has no corrections to demote,
and asking for it is a compile-time `$:error`.

`KEEP_PREC /= VISC_PREC` is what puts work on both hardware pipes at once. The
split path's `calc_Ev` is always FP64 and serves as the reference.

**FP32 shadow arrays (`calc_quantities_shadow32`).** Whenever a fused-style
kernel needs FP32 copies of its inputs (`VISC_PREC='fp32'`, `PRESS_PREC='fp32'`,
or `KEEP_PREC='term'`), the FP64→FP32 conversions are fused directly into the
per-stage primitives-decode kernel: `calc_quantities_shadow32`
(`src/calc_flux_base.f90.fypp`) reproduces the repo-root, shared-with-2D/3D
`calc_quantities_T_1D`'s exact arithmetic bit-for-bit (that shared file is
never modified) and additionally writes whichever of `u32`/`t32`/`mu32`/
`rho32`/`p32` the active config needs, in one launch. When none are needed,
`calc_quantities_T_1D` is still called directly and unmodified — this is what
keeps every zero-shadow config bit-for-bit identical. Both F2F directions
issue on the FP64 pipe — the bound resource — so converting once here instead
of at every flux-kernel tile load is a net win (`report/rtx4060_precision_vs_specialization.md`:
this single kernel costs 3518–3653 µs/stage depending on which arrays it
writes, **less** than the 3307 µs `calc_quantities_T_1D` cost on its own
*plus* the old separate `fill_shadow32` pass it replaced — that base cost was
never measured before this fusion and had been silently missing from every
earlier end-to-end claim in this directory). The flux kernels load the
shadows with plain 4-byte loads; the FP32 bits are identical either way.

**`VISC_ORDER` is a stencil WIDTH knob, not an order-of-accuracy knob.** The face
flux it builds is genuinely 2nd/4th/6th-order accurate *at the face*, but `calc_R`
turns it into a divergence with a one-cell difference, which carries its own
irreducible `(h²/24)·F'''` term — so `d/dx` of the viscous flux is **2nd order at
every width**. `report/check_visc_order.py` prints both tables so the fact is
measured rather than assumed. `3D_solver`'s `calc_visc_high` has exactly the same
property; 1D mirrors it deliberately. What `VISC_ORDER` is *for* here is that it
grows the FP32 instruction count, which is what feeds the other hardware pipe.

The ghost-cell depth `ng = max(ORDER, VISC_ORDER)//2` and the shared-memory halo
`io = ng - 1` are **derived inside the templates**, not set in `config.fypp`, so
they cannot drift out of sync with the stencil widths. This deviates from the 3D
convention (`ORDER_IO`) on purpose.

## Ghost cells, not an order-degrading ladder

`set_bc` fills `ng` ghost cells at each end (zero-gradient/flat extension), which
makes one uniform interior stencil valid at every face the RK update consumes.
There is **no boundary branch in any flux kernel**.

| | range |
|---|---|
| ghost cells (`set_bc` writes) | `1 … ng` and `nx-ng+1 … nx` |
| physical cells (RK updates) | `ng+1 … nx-ng` |
| faces the flux kernels compute | `ng … nx-ng` |
| cells the stencil reads at face `f` | `f-io … f+io+1` |

The arithmetic is tight: over all computed faces the stencil reach is exactly
`1 … nx` at every order. `ng=1` reproduces the pre-ghost-cell ranges (faces
`1..nx-1`, cells `2..nx-1`) exactly, so `ORDER=VISC_ORDER=2` is bit-identical to
the old order-degrading code — which is the cheapest regression check available
here, and worth running after any change to the index arithmetic.

This is not only tidier than degrading the order near the edge; **it is what makes
the ILP mechanism work at all above ORDER=2.** ptxas cannot co-schedule across a
basic-block boundary, and the old `if/elseif/else` ladder put the FP32 viscous
half and the FP64 convective half in different blocks: the emitted SASS was
`S×12` then `D×151`. Interleaving survived only at ORDER=2, i.e. exactly the case
where the ladder collapsed to one unconditional statement. Deleting it also
removed the dead boundary arms from the instruction stream (ORDER=6 FP64: 151 →
95). See `report/SASS_ilp.md` §7.

`ST/set.f90` takes `ng` as a runtime argument rather than deriving it, so that
file stays out of the fypp pipeline; `calc_time_dev` passes the compile-time
constant and also owns the `blocksE`/`blocksEv`/`blocks` launch geometry, which
`main.f90` used to compute.

## Higher-order KEEP

`src/calc_keep_1d.f90.fypp` generates KEEP2/KEEP4/KEEP6 as 1D specialisations of
`src/calc_scheme_math.f90.fypp` (repo root). In 1D the transverse velocities
vanish and the face normal is always +1, so the `v`/`w` terms and the `Normal`
argument collapse away. The shared file is **not** reused: its
`VELS = [...] if DIM == 3 else [...]` has no `DIM=1` branch and its `' + '.join()`
expressions would emit empty strings.

The shared versions hand-write `fma()` chains to pin the association order; the
1D versions use the equivalent plain expressions instead, because libm's `fma`
is real(8)-only and the functions must also compile at `real(4)`. `-Mfma`
contracts them to D/FFMA regardless.

Order is selected by the **KIND of the unused leading `id_acc` argument**
(`integer(2)/(4)/(8)`) through a generic `KEEP` interface — the same
compile-time dispatch trick `3D_solver` uses. Only the configured order is
generated: with no boundary ladder nothing calls the lower orders any more, and
leaving them out keeps dead function bodies out of the SASS listing the study
reads.

## Higher-order viscous

`src/calc_visc_1d.f90.fypp` is the viscous sibling of `calc_keep_1d.f90.fypp`,
generating `VISC2`/`VISC4`/`VISC6` with the same `id_vacc`-KIND dispatch. It is a
1D specialisation of `3D_solver/src/calc_visc_cent.f90.fypp` — with no transverse
velocities there are no cross-derivative terms, leaving only a face interpolation
of `mu` and `u` and a face derivative of `u` and `T`:

| p | `interp` | `diff` (× 1/dx) |
|---|---|---|
| 2 | `[1,1]/2` | `[-1,1]` |
| 4 | `[-1,9,9,-1]/16` | `[1,-27,27,-1]/24` |
| 6 | `[3,-25,150,150,-25,3]/256` | `[-9,125,-2250,2250,-125,9]/1920` |

Two deliberate deviations from the 3D source. It writes **pair sums/differences**
rather than 3D's serial `fma()` accumulator: 3D chains six dependent FMAs into
one register to minimise live values, but that is a serial chain with no internal
ILP, and here the whole point is to leave independent FP32 instructions for the
scheduler to slot into the FP64 stream. And coefficients are written as exact
fractions evaluated at the working kind (`real(125.d0/1920.d0, VK)`) rather than
3D's truncated decimals like `0.065104167d0`, which lose digits at `real(4)`.

Because the split path must stay FP64 while the fused path follows `VISC_PREC`,
the fragment emits precision-suffixed specifics (`VISC6_d`, `VISC6_s`) from one
fypp loop, and each including module names the one it wants in its own generic
`VISC` interface.

## Which variant to use, and what the ILP experiment actually showed

Measured at `nx=4194304` (`report/rtx4060_variant_optimization.md`): **all three
variants are within 3% of each other**, and `seq` — which does nothing clever —
is the reference. `fused` matches it to 0.1%; `warp` is 0.4–3.0% faster and does
strictly more work per face (a shared-memory round trip plus a second barrier)
using twice the threads.

Two dead ends are recorded so they are not retried:

- **Do not add `volatile` scalars or artificial `x - x` zero-dependencies to force
  the FP32/FP64 alternation.** It was tried. `DADD` and *both* `F2F` directions
  issue on the FP64 unit, so it added +64% FP64-pipe work, and chaining the two
  halves destroys the independence co-issue needs. The SASS `runs` count rose
  while FP64 issue stall more than doubled and the kernel got 9–28% slower.
- **`runs` (the SASS interleaving count) does not predict time** and must not be
  used as an optimization target. The static metric that *does* track time is
  FP64-pipe instruction count = `DADD + DMUL + DFMA + F2F`.

The kernel is **FP64-issue-bound at ~86% pipe utilization** whenever
`KEEP_PREC='fp64'`. The only lever that matters is reducing FP64 work:
`KEEP_PREC='term'` is worth 2.9× at ORDER=6, where scheduling tricks are worth
1.00×. Type conversions alone are 39% of FP64-pipe work at ORDER=2.

**Update (2026-08-14, `report/rtx4060_warp_fused_optimization.md`):** hoisting
those conversions into a shadow pass made every fused-style kernel 4–17%
faster (bit-identical), and `warp_fused` + `PRESS_PREC='fp32'` — FP32 pressure
terms on the viscous-role warps — is now the fastest flux kernel at ORDER≥4:
**0.819×/0.791× seq at ORDER 4/6**.

**Update 2 (2026-08-14, `report/rtx4060_precision_vs_specialization.md`):**
`PRESS_PREC='fp32'` is no longer `warp_fused`-only — `seq`/`fused` support it
too, as a control for exactly this question: does the win come from the
precision demotion or the warp specialization? **Entirely the demotion.**
`seq(PP=fp32)`/`fused(PP=fp32)` match or *beat* `warp_fused(PP=fp32)` at every
order (0.948–0.994× its time) — the warp role split contributes no additional
speedup once the FP64-pipe instructions are actually removed rather than
relocated; it is pure overhead on top of the same demotion. Prefer `seq`/`fused`
with `PRESS_PREC='fp32'` over `warp_fused` for this precision combination.
Separately, the old `fill_shadow32` pass was fused directly into the
quantities-decode kernel (see above), which turned out to be a real end-to-end
win (3.3–9.0% per RK stage) once `calc_quantities_T_1D`'s own previously
unmeasured 3307 µs/stage was accounted for. `KEEP_PREC='term'` composes with
`PRESS_PREC='fp32'` too (~0.82–0.92×, smaller than the non-term ratio since the
two demotions overlap rather than stack).

Benchmark with `report/bench.sh`, never with the stock `nx=4096` — see the
performance note under Notice.

## Verification

`report/` holds four checks; run all of them after touching a flux stencil.

```bash
# 1. the coefficients are genuinely 2nd/4th/6th order (pure numpy, no GPU)
python3 report/check_keep_order.py
python3 report/check_visc_order.py

# 2. the Fortran on the GPU computes those same coefficients
#    (set nt=1 in ST/mod_globals.f90, rebuild, run, then:)
python3 report/check_vs_numpy.py ST/build/Q.dat <ORDER> [VISC_ORDER]
#    for seq/fused/warp_fused with PRESS_PREC='fp32', append `--press fp32` so
#    the numpy reference also evaluates the KEEPP pressure terms in float32.
#    This script has no KEEP_PREC='fp32'/'term' model at all -- only valid
#    against KEEP_PREC='fp64' builds. For 'term' builds (incl. term x
#    PRESS_PREC='fp32'), use check 3 instead, against exact and against a
#    same-config PRESS_PREC='fp64' reference Q.dat.

# 3. the full run matches the exact Riemann solution (nt=2500)
python3 report/check_sod.py ST/build/Q.dat [reference_Q.dat]
```

Check 1 is the one that catches a mistyped coefficient; a shock tube cannot,
because every scheme is first-order at a discontinuity. `check_visc_order.py`
prints two tables — read the docstring before reacting to the second one, which
is *supposed* to show 2/2/2 (see the `VISC_ORDER` note above). Check 3's density
L1 against the exact solution barely moves with order (0.185% → 0.182% → 0.182%)
for the same discontinuity reason — do not read it as an order test.

The cheapest regression check of all: `ORDER=VISC_ORDER=2` must reproduce a
pre-change `Q.dat` **bit-for-bit**, because `ng=1` degenerates to the original
index ranges. Capture a baseline before touching the index arithmetic.

SASS analysis helpers (`report/sass_analyze.py`, `report/ctrl_bits.py`) and the
findings are described in `report/SASS_ilp.md`.

## Notice

* `nx=4096` is far too small to measure throughput: 0.1 waves/SM, ~11% achieved
  occupancy, launch-latency bound — every configuration sits 4–35× above its own
  compute floor there. Use `report/bench.sh` (default `nx=4194304`, ~114–341
  waves/SM, 0.03% launch spread) for any performance claim, and keep `nx=4096`
  for SASS and correctness only. Four earlier reports drew variant rankings from
  `nx=4096`; two of those rankings did not survive being re-measured.
* **GPU clocks cannot be locked on this machine** (`nvidia-smi -lgc` needs root),
  so a long sweep throttles: mean kernel times drift >2× while pipe-utilization
  ratios stay put. Compare `time_min_us` from interleaved, repeated rounds, and
  distrust any `time_us` whose `spread_us` approaches it.
* `CMakeLists.txt` hard-codes `ccnative` (every other solver uses
  `cc${CASE_GPU_CC}`), so the target architecture silently follows the build
  machine. This matters a lot here: FP64:FP32 is 1:64 on consumer Ada/Blackwell
  and 1:2 on GH200, and the mixed-precision conclusions invert between them.
* `maxregcount:96` (raised from 64). `ORDER=6` now uses 68–76 registers — the
  uniform 6-point stencil plus a 6-point viscous half no longer fits under 64,
  and spilling would put LDL/STL traffic into the very instruction stream the
  SASS study measures. Occupancy is not a concern at `nx=4096` (see above), but
  it would be on a real grid. Check `ptxas` output for spills after any change
  that adds live values.
* A crashed run leaves the previous `Q.dat` in place. Always `rm -f Q.dat`
  before a run in a sweep script, or a crash will be silently reported as a pass
  with stale numbers.
* `nt` in `ST/mod_globals.f90` is switched between `1` (SASS/profiling and
  `check_vs_numpy.py`) and `2500` (physics and `check_sod.py`). Check which is
  set before interpreting `Q.dat`.
