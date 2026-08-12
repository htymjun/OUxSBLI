# Accelerate the overlapped dispatch in `3D_solver_mixed` (GH200)

## Context

`3D_solver_mixed` (FP64 convection / FP32 viscous) selects between `calc_conv_visc_seq` and the
Green-Context overlapped `calc_conv_visc_paired_mixed` via `OVERLAP_GC`. The overlapped path is
slower than sequential on GH200 and the goal is to make it faster.

The new captures in `3D_solver_mixed/nsys_ncu_SoA/` (post-SoA, post-`calc_div_mp`, 513³) — now with
sqlite exports, so the **timeline** is readable, not just per-kernel counters — settle what is
actually wrong. One steady-state RK stage:

| | GC=False | GC=True |
|---|---|---|
| **stage wall** | **63.83 ms** | **84.31 ms** (+32%) |
| serial prologue (quantities 5.72, div 2.83, zero-fill) | 8.55 | 11.62 |
| overlap region | — (serial 42.64) | 11.63 → 71.64 = **60.01** |
| serial epilogue (`calc_step2_3` 11.99, bc 0.63) | 12.63 | 12.63 |

Inside the GC overlap region:

| stream | busy | idle waiting on events |
|---|---|---|
| conv (56 SM) | **60.00 ms** | **0.01 ms** |
| visc (76 SM) | **39.78 ms** | **0.01 ms** — then idle 51.43 → 71.64 (**20.2 ms**) |

## What the measurement actually says

**1. The overlap machinery is already perfect.** Neither stream stalls on an event (0.01 ms).
`wall = prologue + max(conv, visc) + epilogue` holds to within 0.1 ms. There is nothing to fix in the
event/phase structure, and the `cudaStreamWaitEvent(stream_conv, event_Fv/Gv)` calls cost nothing —
contrary to what I suspected before seeing the timeline.

**2. The loss is entirely overhead plus a badly imbalanced split.** The 20.5 ms deficit decomposes as:
- **zero-fill 3.06 ms/stage** (`__pgi_dev_cumemset_16n` ×2 — `F=0.d0`/`G=0.d0`, GC-path only)
- **`_acc` read-modify-write ≈10.5 ms/stage** — scaling `slau_y`/`slau_z` at the rate stock
  `slau_x` actually scales (2.317×) predicts 13.90/17.82 ms; they measure 22.35/19.87 ms
- **imbalance ≈20.2 ms** — conv needs 60.0 ms on 56 SMs while visc needs 39.8 ms on 76 SMs, so the
  viscous partition sits idle for a third of the region. `conv_share=0.4` starves the side that
  needs the SMs.

**3. Work conservation holds, so Green Contexts can only tie sequential — the original conclusion in
`mod_streams_mixed.f90` was right.** Measured against perfectly-linear SM scaling: conv (minus the
`_acc` excess) lands at 101.2% of linear, visc at 104.7%. Rebalancing to the work-optimal 64/68 split
and deleting both overheads gives `8.55 + max(2740/64, 2888/68) + 12.63 = 8.55 + 42.8 + 12.63
= **63.98 ms** vs sequential 63.83 ms` — a dead heat.

> Correction to my earlier reading of the `.ncu-rep` files: there the viscous kernels *appeared* to
> scale sublinearly (1.40–1.65× on 58% of the SMs), which would have meant spare capacity and a real
> win. That was an **ncu replay artifact** — ncu serializes kernels, so a kernel confined to 76 SMs
> still gets the whole memory system to itself. Under genuine concurrency in the nsys timeline both
> families scale ~linearly. Per-kernel ncu counters cannot answer this question; only the timeline can.

**4. So the only mechanism with real headroom is SM *sharing*, not SM partitioning** — and your own
GH200 microbenchmarks in `fp32fp64_report/` already characterise it: FP64:FP32 = 1:1.65, co-residency
is trivially achieved (**all 132 SMs, 100% time overlap in every configuration**), and concurrent
FP32+FP64 beat serial in every configuration by **4.0–27.9%**, best when neither side saturates issue
slots. The instruction probe is the relevant detail: an issue-saturated FFMA partner strangles FP64 to
20.3% retention, while a **low-issue-rate partner leaves 67.9%**. Our viscous kernels sit at 43–59%
issue occupancy (ncu) — the favourable regime.

## Design

Add an `OVERLAP_MODE = 'none' | 'gc' | 'streams'` selector, and for `'streams'` use **plain
non-blocking streams with a per-direction decomposition**:

```
prologue (default stream) : calc_quantities -> calc_div_{ux,vy,wz} -> calc_Ducros
stream 1 : calc_slau_x_in -> calc_Ev6_in_mp     (owns E)
stream 2 : calc_slau_y_in -> calc_Fv6_in_mp     (owns F)
stream 3 : calc_slau_z_in -> calc_Gv6_in_mp     (owns G)
epilogue (default stream) : calc_step* , set_bc
```

Why this shape rather than the current conv-stream/visc-stream split:

- **No cross-stream dependencies at all** — each stream owns one flux array end to end, so there are
  no events, no `cudaStreamWaitEvent`, and no ordering subtleties. Within a stream the order is
  conv-then-visc, exactly as `calc_conv_visc_seq` does per array.
- **Both overheads disappear**: conv writes fresh (stock `intent(out)` kernels, so no `_acc`) and
  nothing reads the array before it is written (so no zero-fill). That is the full 13.6 ms/stage.
- **It is bit-identical to sequential by construction**, since per array the write order is unchanged.
  The 3-phase design inverted that order, which is exactly why it needed the zero-fills, the `_acc`
  rewrite, and two rounds of bug fixes (NaN blowup, then fresh-write-clobbers-viscous).
- **It mixes the pipes naturally.** The three streams have different kernel durations, so at any
  instant one stream is typically in its FP32 viscous kernel while another is in its FP64 convective
  kernel — the co-execution the microbenchmarks measured, with no scheduling effort.
- `calc_div` moves into the prologue where it is already unavoidable serial work, but **on the
  measured timeline it is 2.83 ms/stage that currently overlaps with nothing** — see the optional
  refinement below.

**Primary risk to state plainly:** the microbenchmarks forced co-residency with 1-block/SM launches.
Our kernels launch millions of blocks, and the hardware work distributor may simply drain one grid
before admitting the next, in which case `'streams'` degenerates to sequential (no loss, no gain).
This is cheap to find out empirically and is the main thing the first `'streams'` run answers.

## Implementation

1. **`NSTGV/config.fypp`** — `#:set OVERLAP_MODE = 'streams'`, then
   `#:set OVERLAP_GC = (OVERLAP_MODE != 'none')` so every existing `#:if OVERLAP_GC` site keeps
   working untouched.

2. **`src/mod_streams_mixed.f90` → `.f90.fypp`** (move `_BASE` → `_SHARED_FYPP` in
   `3D_solver_mixed/CMakeLists.txt` so it can read `OVERLAP_MODE`). For `'streams'`: create three
   ordinary streams with `cudaStreamCreateWithFlags(..., cudaStreamNonBlocking)` and skip the
   Green-Context path. Keep `green_ctx_bindings_mixed.f90` compiled unconditionally (see the existing
   `_BASE` comment). Create events with `cudaEventDisableTiming`. **Check every `cires`/`istat`** —
   all return codes are currently discarded, so a failed create silently yields a null handle.

3. **`src/calc_flux_base_mixed.f90.fypp`** — add `calc_conv_visc_streamed` implementing the three-way
   split above: no zero-fill, stock `calc_slau_y_in`/`calc_slau_z_in`, viscous kernels unchanged
   (they already accumulate). Leave `calc_conv_visc_paired_mixed` in place for `'gc'`.

4. **`conv_share` 0.4 → 0.48** (`'gc'` mode only; `round_to_multiple8(int(0.48*132))` = 64, the
   work-optimal 64/68). Rewrite the comment: keep its work-conservation conclusion — now confirmed
   from the timeline rather than assumed — but correct the arithmetic, which derives 53/79 while the
   code actually produces 56/76, and record that the loss was imbalance plus overhead, not the
   partitioning mechanism.

5. **Also fix the `'gc'` path's overhead** (shell-only zero-fill, or reuse the streamed dispatch's
   ordering) so `'gc'` is measured at its true ceiling rather than 20 ms below it.

6. **Instrumentation** — add `-lnvhpcwrapnvtx` and NVTX ranges (this fork has none, unlike
   `3D_solver_fltflt`), plus `NSTGV/job_miyabi_nsys.sh` next to the existing `job_miyabi_ncu.sh`; no
   nsys job script exists in the repo. Keep exporting `.sqlite` alongside `.nsys-rep` — that export
   is what made this analysis possible.

`src/calc_slau_kernel_internal_acc.f90.fypp` becomes unused by `'streams'`; leave the file and its
CMake entry so `'gc'` still builds.

## Expected outcome

Stage-level, against the measured 63.83 ms sequential baseline:

| variant | predicted | note |
|---|---|---|
| `'gc'` today | 84.31 ms | measured |
| `'gc'` rebalanced + de-overheaded | ~64.0 ms | ties sequential — work conservation |
| `'streams'` @ 4% co-exec gain | ~62.1 ms | microbenchmark worst case |
| `'streams'` @ 15% | ~57.4 ms | ~10% faster |
| `'streams'` @ 28% | ~51.9 ms | microbenchmark best case, ~19% faster |

So: fixing `'gc'` is worth ~24% *to that path* but only reaches parity; `'streams'` is the only route
to an actual win, with a realistic ~10% and an optimistic ~19%.

## Verification

- Build all three modes: `cd 3D_solver_mixed/NSTGV && cmake -B build -DCMAKE_Fortran_COMPILER=mpif90 && cmake --build build -j`
  (NVHPC's own mpif90 on PATH: `/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/comm_libs/openmpi4/bin`).
- **Correctness gate:** `entropy.d` / `kinetic_energy.d` **bit-identical** across `'none'`, `'gc'`,
  `'streams'` at the local 129³ / np=nt=10 config. For `'streams'` this should be exact by
  construction, since per-array write order matches sequential.
- Local wall-clock is not a proxy (24 SMs, split collapses to 16/8, 129³ is launch-overhead dominated:
  seq 6 s vs gc 18 s this session). Local runs are correctness-only.
- **GH200 verdict run:** all three modes at 513³, nsys **with `--type sqlite` export** plus ncu, into
  `nsys_ncu_SoA/`. The decisive number is the per-stream busy/idle breakdown of the overlap region,
  exactly as tabulated above — for `'streams'`, check whether the three streams' kernels actually
  overlap in time or serialize (the primary risk).

## Adjacent, and probably the bigger prize

The serial prologue+epilogue is **21.2 ms = 33% of the stage** and overlaps with nothing:

- `calc_step2_3` **11.99 ms at ~91% DRAM throughput** — bandwidth-saturated, immune to overlap, and
  the largest single non-flux cost. `CLAUDE.md` already flags `calc_steps.f90` as a main bottleneck.
- `calc_div_*_mp` **2.83 ms** — a pure dependency of the viscous kernels only, so it could run
  concurrently with `calc_Ducros`/the convective kernels instead of blocking them. Worth ~4% of the
  stage on its own and independent of which overlap mode wins.
- `threadsEv = threadsE = dim3(32,1,1)` — single-warp blocks, 4.18 M of them at 513³; Hopper's
  32-blocks-per-SM cap pins `calc_Ev6_in_mp` near 50% occupancy. This interacts directly with the
  streams work, since the microbenchmarks found co-execution is best at moderate occupancy on both
  sides.
- `CLAUDE.md:347` points at `3D_solver/nsys_ncu/*.ncu-rep`, which no longer exists.
