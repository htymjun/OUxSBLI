# Warp-level + ILP optimization of the fused KEEP+viscous kernel — Phase A status

## Context

A prior session fused the KEEP convective flux and NS viscous flux into one
kernel per sweep direction — `3D_solver/src/calc_keep_visc_kernel.f90.fypp`
(`calc_keepv_x_in/y_in/z_in`), still **uncommitted** in this worktree along with
its dispatch wiring in `3D_solver/src/calc_flux_base.f90.fypp` (search
`fused_keepv`) and the registration in `3D_solver/CMakeLists.txt`'s
`_SHARED_FYPP` list. It's gated by an fypp guard (`SCHEME=='KEEP' and
VISC=='NS' and ORDER==6 and VISC_ORDER==6 and not BC_X/Y/Z and not COMMZ`), a
no-op everywhere except **DHIT** by default. It's pure FP64 — despite the
kernel originally being framed as "FP32/FP64 concurrent execution
optimization," no FP32 path or stream dispatch exists, and that framing was
explicitly dropped this session (see below).

Full plan and reasoning: `/home/jhatayama/.claude/plans/i-fused-convection-kernels-fancy-curry.md`
(approved plan file from this session — still the authoritative design doc).

**Two things ruled out the FP32/FP64 / stream-concurrency direction**, investigated
via a sibling checkout `~/ouxsbli/mixed/OUxSBLI` (branch `mixed`):
- CUDA streams / Green Contexts to overlap conv/visc (or FP32/FP64) kernels:
  measured **1.00x concurrency factor** on GH200 (zero overlap) — a single flux
  kernel at production grid sizes already floods the GPU across hundreds of
  scheduling waves, so there's no idle SM capacity for a second stream to fill.
  Matches this repo's own `CLAUDE.md` rule against separate streams in
  `calc_flux_base.f90`.
- Standalone FP32/FP64 concurrent-issue microbenchmarks (never wired into the
  real solver): weak/fragile on GH200 specifically (~28% best case, vs. strong
  on a degenerate-FP64 consumer GPU) — not worth pursuing here.

**A more aggressive prior fusion attempt already failed on this exact
x-direction kernel** — an untracked benchmark (`benchmark_dr_vs_efg.md`, found
in sibling clones `~/ouxsbli/SoA/OUxSBLI` and `~/ouxsbli/fusion/OUxSBLI`) used
`__shfl_down` to share flux faces across adjacent threads in a warp, writing
straight to the flux-divergence array. Measured on an RTX4060: **~2x slower**,
at *identical* occupancy to the legacy split kernel — registers rose 68→80 and
local-memory (spill) frame tripled 40→120 bytes. Not bit-reproducible either
(atomicAdd in y/z). **This plan deliberately avoids that mechanism.** The lever
used instead: more independent, communication-free warps per block (no
cross-warp data sharing at all) — exactly how `calc_keepv_y_in`/`_z_in` already
work.

GH200 profiling (`ncu/my_report_ncu_20260426.json`, git-tracked, captured on an
actual GH200 120GB) shows the concrete opportunity: the x-direction kernel
launches **1 warp/block** (`threadsE=dim3(32,1,1)`), which caps occupancy at
**50%** on GH200 (32 max resident blocks/SM × 1 warp ÷ 64 max warps/SM)
*regardless of registers* — measured ~42-49%. The y/z kernels already use 4
warps/block and reach 54-60%, limited by registers/shared-mem instead.

## Completed items (do not redo)

| Item | Status |
|---|---|
| Hand-traced `calc_keepv_x_in`'s shared-mem sizing, offset arithmetic, halo-load loop, and `blocksE`'s grid formula (`src/set_coordinate.f90:35`) — confirmed all already generalize to `threadsE%y>1` with **zero kernel-code changes needed** | ✅ Done |
| **Phase A implemented**: `3D_solver/DHIT_warpilp/` created (copy of `DHIT/`), only `threadsE = dim3(32,2,1)` changed (was `dim3(32,1,1)`) in `mod_globals.f90` | ✅ Done |
| Built + ran both `DHIT` (baseline) and `DHIT_warpilp` at `nx=ny=nz=70`, `np=200` on local RTX4060 | ✅ Done — 9.0s vs 8.3s (not a GH200 occupancy signal, just confirms it runs) |
| Correctness gate on DHIT pair: `kinetic_energy.d`, `entropy.d` byte-identical; all 201 `Q*.vtr` snapshots `cmp`-identical | ✅ **Passed exactly** |
| **Second, independent verification** on a larger/different grid: `3D_solver/NSTGV_keep128` created (copy of `NSTGV`, `SCHEME` switched `'SLAU'→'KEEP'`, `M0` lowered `1.25→0.4` subsonic — matches ETGV's own KEEP+TGV `M0`, since KEEP needs shock-free flow — grid reduced `513³→128³` for local 8GiB RTX4060 feasibility, matching `benchmark_dr_vs_efg.md`'s own precedent) | ✅ Done |
| `3D_solver/NSTGV_keep128_warpilp` created (same, `threadsE=dim3(32,2,1)`) | ✅ Done |
| Built + ran both NSTGV-derived cases (~12 min each on RTX4060, nearly identical: 12m2s vs 12m4s) | ✅ Done |
| Correctness gate on NSTGV_keep128 pair: same triple check | ✅ **Passed exactly** — 0 mismatches across 101 `Q*.vtr` + both `.d` files |

Production `3D_solver/DHIT/` and `3D_solver/NSTGV/` were **not modified** — all
new code lives in untracked sibling directories (`DHIT_warpilp`,
`NSTGV_keep128`, `NSTGV_keep128_warpilp`). `DHIT/data/` had a stale `np=100`
archive from an earlier run; preserved as `DHIT/data_archive_np100_prefusion_check/`
before regenerating a fresh `np=200` baseline for the comparison.

## Current status

| Item | Status |
|---|---|
| Phase A (`threadsE=dim3(32,2,1)` for `calc_keepv_x_in`) | ✅ Implemented, bit-reproducibility verified on 2 independent cases |
| GH200 occupancy/performance verification | 🔲 **Next — blocked on GH200 access** (this dev box only has an RTX4060) |
| Escalate to `dim3(32,4,1)` if N=2 leaves headroom | ⏳ Decide from GH200 data |
| Phase B (ILP) | ⏳ Not started — deliberately deferred until Phase A's real profiling data is in hand |

## Next step: GH200 profiling handoff

Run on the GH200 machine (mirrors the repo's own `job_miyabi_ncu.sh` /
`ncu/export_ai_context.py` convention):

```bash
# Baseline, unmodified (threadsE=dim3(32,1,1))
cd 3D_solver/DHIT/build   # or 3D_solver/NSTGV_keep128/build for the larger-grid case
ncu --kernel-name regex:calc_keepv_x_in --launch-count 30 --set full --import-source yes \
    --target-processes all -o keepv_x_baseline_N1 mpiexec -np 2 ./a.out

# Phase A (threadsE=dim3(32,2,1))
cd 3D_solver/DHIT_warpilp/build   # or 3D_solver/NSTGV_keep128_warpilp/build
ncu --kernel-name regex:calc_keepv_x_in --launch-count 30 --set full --import-source yes \
    --target-processes all -o keepv_x_phaseA_N2 mpiexec -np 2 ./a.out
```

Export: `ncu --import <file>.ncu-rep --csv --page raw > <file>.csv`, then
`ncu/export_ai_context.py` for the JSON summary. Share both baseline and
Phase-A captures back.

**What decides success**: achieved occupancy rises well past the old 50%
ceiling **and** `gpu__time_duration.sum` for `calc_keepv_x_in` drops or holds
flat. If duration regresses despite higher occupancy, check Registers Per
Thread / Static Shared Memory Per Block first (the `__shfl_down` attempt's
failure signature was register/local-memory growth, not occupancy) before
concluding the warp-count lever was wrong.

Warp Cycles/Issued-Instruction and Scheduler "No Eligible %" from this same
capture decide which Phase B (ILP) item — if any — to pursue next (menu in the
approved plan file: breaking the `utau_sum` sequential-accumulation chain,
per-thread work coarsening, loop unrolling, or a speculative `syncwarp()` swap
— each gated on a specific profiling signal, not to be picked blind).

## Critical files

| File | Role |
|---|---|
| `3D_solver/src/calc_keep_visc_kernel.f90.fypp` | The fused kernel itself — unmodified in Phase A |
| `3D_solver/DHIT_warpilp/mod_globals.f90` | Phase A change #1 — `threadsE=dim3(32,2,1)` |
| `3D_solver/NSTGV_keep128/mod_globals.f90`, `config.fypp` | New KEEP+subsonic+128³ baseline case (activates the fusion guard on NSTGV's grid/setup) |
| `3D_solver/NSTGV_keep128_warpilp/mod_globals.f90` | Phase A change #2, same edit on the NSTGV-derived case |
| `src/set_coordinate.f90:35` | `blocksE`'s grid formula — confirmed already generic in `threadsE%y/%z`, no edit needed |
| `/home/jhatayama/.claude/plans/i-fused-convection-kernels-fancy-curry.md` | Full approved plan with complete reasoning, occupancy math, and Phase B menu |

## Verification already done vs. still needed

- ✅ Correctness (bit-reproducibility): done, twice, on RTX4060.
- 🔲 Performance/occupancy on GH200: **not done** — needs the user's remote
  machine. This is the only remaining gate before deciding whether to escalate
  to `dim3(32,4,1)` or move on to Phase B.

## Update — 2026-08-05: GH200 fuse-vs-nofuse capture + FP32/FP64 re-verdict

Two fresh Nsight Compute captures landed in `3D_solver/ncu/`
(`my_report_ncu_fuse_20260805.ncu-rep`, `my_report_ncu_nofuse_20260805.ncu-rep`),
analyzed directly via `ncu --import` (real GH200 120GB, CC9.0 — confirmed from
`--page session` metadata, not assumed). This is **not** the N=1-vs-N=2
`threadsE` occupancy test requested above — it's a different, also-useful
comparison, described precisely below.

**Experimental setup** (identified from ncu-rep process-path metadata, not
guessed): both captures profile the **NSTGV case at full production grid
(513³)**, `VISC='NS'`, `SCHEME` temporarily set to `'KEEP'` on both sides (to
exercise the KEEP convective path at all — production `NSTGV` normally ships
`SCHEME='SLAU'`). The *fuse* capture ran a remote `.../SoA/OUxSBLI` checkout
(this branch, with the still-uncommitted `calc_keep_visc_kernel.f90.fypp`
fusion active — `fused_keepv` guard true). The *nofuse* capture ran a separate
remote `.../stable/OUxSBLI` checkout (the `stable` branch, which has no fusion
code at all, so it naturally runs the legacy split `calc_keep_kernel` +
`calc_visc_high` kernels). Both x-direction kernels launch with identical
`threadsE=(32,1,1)` — 1 warp/block on both sides, confirmed via
`launch__block_dim_*` metrics — so this isolates **just the fusion effect**,
cleanly, at full production grid size.

(Local note: this dev box's `3D_solver/NSTGV/config.fypp`/`mod_globals.f90` had
picked up unstaged leftovers from setting up this experiment — `SCHEME`
SLAU→KEEP and grid 513³→129³ — unrelated to the remote checkouts above, which
have their own independent `nx=513` state. Reverted back to production
defaults on 2026-08-05.)

**Kernel time breakdown** (sum over captured invocations):

| Report | Total captured GPU time | KEEP+visc combined |
|---|---|---|
| fuse | 145.5 ms | 86.2 ms (25.3+28.6+32.4, x/y/z) |
| nofuse | 195.9 ms | 136.5 ms (KEEP 48.3 + visc 88.2) |

Fusing KEEP+viscous cuts their combined time **~37%** (136.5ms → 86.2ms) at
full 513³ production grid, on real GH200 hardware. Other shared kernels
(`calc_steps`, `calc_physical_quantities`, `calc_div_*`) are ~identical between
runs, confirming the difference is isolated to KEEP+viscous as expected.

**Pipe utilization / occupancy** (Hopper: FP32 goes through the FMA pipe, no
separate fp32 counter):

| Kernel | Mem Throughput % | FMA(FP32) %peak | FP64 %peak | Occupancy % |
|---|---|---|---|---|
| calc_physical_quantities | 75.7 | 5.8 | 11.1 | 59.4 |
| calc_div_ux | 68.7 | 18.6 | 10.2 | 78.4 |
| calc_keepv_x (fuse) | 62.4 | 11.6 | 28.7 | 29.8 |
| calc_keep_x (nofuse) | 64.8 | 10.3 | 25.4 | 34.6 |
| calc_visc_ev6 (nofuse) | 52.8 | 10.3 | 14.8 | 49.2 |
| calc_steps_step1 | 91.4 | 7.4 | 4.4 | 68.7 |

Every kernel is dominated by `long_scoreboard` (memory-latency) stalls — 55-93%
of stall cycles across the board (e.g. `calc_steps` 92.9%, `calc_physical_quantities`
82.6%, `calc_keepv_x` 55.7%) — combined with moderate occupancy, partly capped
by the 1-warp/block launch shape for the KEEP/visc kernels specifically.

**FP32/FP64 verdict — re-confirmed, not reversed.** FMA(FP32) utilization is
6–26% of peak and FP64 is 4–29% of peak, **simultaneously**, across every
kernel in both captures — neither pipe is ever saturated while the other
idles; both are underused together because every kernel is memory/latency
bound, not compute-pipe bound. Concurrent FP32/FP64 execution is a lever for
filling an idle compute pipe sitting next to a saturated one — that scenario
doesn't exist anywhere in this data. It would not address the actual
bottleneck (memory latency / occupancy), and this independently reproduces the
prior session's conclusion (§ Context above: 1.00x stream-concurrency factor,
weak GH200 FP32/64 microbenchmarks) with a second, production-grid dataset.
**FP32/FP64 concurrent execution remains ruled out as an optimization
direction for these kernels.**

**Open item, not acted on in this update:** the fusion's ~37% combined-time
reduction is a validated result on real GH200 hardware, but
`calc_keep_visc_kernel.f90.fypp` and its dispatch wiring remain uncommitted —
whether/when to commit that work is left as a decision for later, not resolved
here.

**Phase A's original ask is still outstanding.** Today's data doesn't answer
it — both sides above use `threadsE=(32,1,1)` (N=1). The N=1-vs-N=2 occupancy
comparison via `DHIT_warpilp` / `NSTGV_keep128_warpilp` described in the "Next
step" section above still needs to be captured on GH200.

## Update — 2026-08-05 (continued): SLAU+NS fusion generalized, Roe removed, NS pytest added

Follow-up session, same day: generalized the KEEP+NS fusion recipe to SLAU,
removed the unused Roe scheme, and added the first pytest coverage for
Navier-Stokes physics in this repo (previously zero — every existing
integration test forces `visc="euler"`). Full reasoning and scope decisions:
`/home/jhatayama/.claude/plans/i-got-ncu-rep-files-crystalline-crystal.md`.

### Roe removal

Confirmed unused (no case anywhere sets `SCHEME='Roe'`) — moved rather than
deleted: `3D_solver/src/calc_roe_kernel.f90.fypp`,
`calc_roe_kernel_internal.f90.fypp`, `calc_roe_3d.f90` → `3D_solver/src/roe/`.
Stripped from `calc_flux_base.f90.fypp`'s `use` statements and the three
`('SLAU','Roe','Hybrid')` sensor-guard tuples, from `mod_constant.f90.fypp`'s
`id_scheme` kind-dispatch, and from `CMakeLists.txt`'s `_SHARED_FYPP` list.
Updated CLAUDE.md/README.md/docs/*.md/tutorials/installation.md/`case.py`'s
`_VALUE_NORMALIZE` — no other stale mentions remain outside historical
`docs/plans/*.md` records and the retired `3D_solver/src/roe/` files
themselves. Build-check across all 8 `3D_solver` cases: clean.

### SLAU+NS fusion (`calc_slau_visc_kernel.f90.fypp`)

Generalized `calc_flux_base.f90.fypp`'s guard (renamed `fused_keepv` →
`fused_conv_visc`) from `SCHEME=='KEEP'` to `SCHEME in ('KEEP','SLAU')`,
dispatching to the matching fused kernel. This activates directly on
production `3D_solver/NSTGV` (SLAU, NS, `ORDER=6`, periodic by default) — a
real case, not a demo one. Full design writeup is now in CLAUDE.md's new
"Fused Convective+Viscous Kernels" section; in short, two adaptations were
needed beyond copying KEEP's recipe:

1. SLAU's own reconstruction would, in the unfused kernel, overwrite the raw
   primitive tile with the left state before calling `SLAU()`. The fused
   kernel computes viscous stress/heat-flux from the **raw** tile first and
   keeps the reconstructed left state in local scalars instead of writing it
   back — so the overwrite that motivated the original ordering never
   happens.
2. SLAU never tiles `T` (KEEP does, for free reuse); the fused kernel reads
   `T` straight from DRAM per-point rather than adding an 11th/12th shared
   array. Also had to locally rename SLAU's own `interp2/4/6` reconstruction
   routines to `recon2/4/6` — they collide by name with
   `calc_visc_cent.f90.fypp`'s own `interp6` function once both are
   `#:include`d into one module (a KEEP never hit, having no reconstruction
   step of its own).

**Correctness — bit-reproducible, verified exactly.** Snapshot-before/edit/
rebuild-after methodology (mirroring the KEEP fusion's own
`DHIT`/`DHIT_warpilp` gate): built & ran `3D_solver/NSTGV` at `nx=ny=nz=129`
(21 output blocks) both before and after the fusion edit. `kinetic_energy.d`,
`entropy.d`, and all 21 `Q*.vtr` snapshots are **byte-identical** (`cmp`) —
no atomicAdd/shuffle/reduction in this fusion, so bit-identical output is the
correct bar, not "close enough", exactly as for KEEP's fusion. Also passed:
build-check (all 8 cases), `test_etgv.py`, and the SLAU regression smoke test
on NSTGV @ nx=128 (both before this edit, for Roe removal, and after, for the
fusion itself) — `ke/ke0` final value identical (`0.0192`) across both smoke
runs, consistent with the exact-diff result.

**Performance — mixed, real signal, needs GH200 to resolve.** This dev box
only has an RTX4060, so nothing below is the real verdict — but it's not
nothing either:

- Direct `ncu` comparison at `nx=128` (registers/thread, static shared
  mem/block, `sm__warps_active` occupancy, `gpu__time_duration.sum`, 3-launch
  average):

  | Kernel | Registers | Shared mem | Occupancy limiter (blocks) | Time |
  |---|---|---|---|---|
  | `calc_slau_x_in` (unfused conv) | 68 | 2.792 KB | regs=28, mem=26, **hw=24** | 6.69 ms |
  | `calc_ev6_in` (unfused visc) | 48 | 1.512 KB | regs=40, mem=25, **hw=24** | 3.83 ms |
  | unfused combined | — | — | — | **10.52 ms** |
  | `calc_slauv_x_in` (fused) | 76 | 3.4 KB | regs=24, **mem=22**, hw=24 | 11.44 ms |

  The fused kernel is measurably heavier — added shared-memory tiles (the
  cross-derivative arrays needed for viscous stress, on top of SLAU's own
  10-array reconstruction tiles) push the **binding occupancy constraint**
  from the hardware block limit (24) down to shared memory (22), exactly the
  register/shared-mem-pressure risk flagged before writing this kernel. Net
  effect at the kernel level: **~9% slower** (11.44ms vs. 10.52ms combined) on
  this GPU — a real but modest regression, not a wash.
- However, the *full* SLAU-smoke-test wall-clock (50,000 steps @ nx=128) told
  a much larger story: ~20-25 min unfused vs. ~90+ min fused, a ~4x gap. That
  gap is far bigger than the ~9% kernel-level number above accounts for, and
  this session ran many long (20-90+ min), sustained 100%-GPU-utilization
  jobs back-to-back on a *laptop* RTX4060 — thermal throttling across the
  session is a much more likely explanation for a wall-clock-only 4x gap than
  a 9%-at-the-kernel-level change somehow costing 4x in practice. Flagging
  both numbers rather than picking one: the short, controlled `ncu` capture
  is the more trustworthy of the two, but neither is a substitute for a real
  GH200 measurement.

**Bottom line:** unlike KEEP's fusion (which GH200 has never actually
profiled either, per Phase A's own still-open ask above), this SLAU
generalization has at least one piece of *quantified, controlled* evidence
pointing to a real (if modest) regression from added shared-memory pressure,
on top of the inherent limits of testing on non-target hardware. **This
should be verified on GH200 before deciding whether to keep it** — same
ask as Phase A's, plus a `calc_slauv_x_in`/`calc_slau_x_in`+`calc_ev6_in`
`ncu` comparison at production grid size specifically watching Static
Shared Memory Per Block and the occupancy limiter, not just kernel duration.
Per this repo's own optimization convention (`ouxsbli-optimize` skill), a
confirmed regression on GH200 should be reverted with a written failure
record; a confirmed win should be committed with the before/after numbers.
**Neither has happened yet — this change is intentionally left uncommitted**
pending that GH200 data, exact-diff correctness notwithstanding.

### New Navier-Stokes pytest coverage

`ouxsbli/tests/test_nstgv_ke_eps.py` (two cases: KEEP+NS subsonic `M0=0.4`,
SLAU+NS supersonic `M0=1.25` — matching NSTGV's own production Mach number)
is the first pytest in this repo to actually exercise viscous physics, and
regression-guards both fused kernels directly since it dispatches on
NSTGV's own case. Built on `ouxsbli/analysis/tgv_ke_eps.py`, refactored from
a hardcoded standalone script into an importable `compute()` (mirroring
`dhit_decay_report.py`'s split) parameterized on physical constants rather
than hardcoded ones (`M0` differs between the KEEP/SLAU cases here).

The check is the *full* compressible kinetic-energy budget on a periodic
domain: `-dEk/dtau == (epsE+epsD) - pdiv`, where `epsE`/`epsD` are the
enstrophy/dilatational viscous-dissipation terms (already in the original
script) and `pdiv = <p*div(u)>` is the pressure-dilatation term — a
*reversible* exchange with internal energy, newly added to
`calc_total_dissipation`, needed because it's not negligible for a
compressible (let alone `M0=1.25` supersonic) TGV. An earlier version of this
test checked `-dEk/dtau` against `epsE+epsD` alone; that only holds within
~1% near `t=0` and grows past 100% error within the first several intervals
for the supersonic case once dilatation develops. Adding `pdiv` closes the
budget to within ~0.4% (KEEP) / ~2.5% (SLAU) across the **entire** run for
both schemes, including through the SLAU case's `dEk/dtau` sign change
partway through — confirmed by direct recomputation on the actual captured
data, not assumed. Test tolerance: 5%, comfortable margin over both measured
maxima.

### Verification summary

- Build-check (8/8 cases): clean, both after Roe removal and after the SLAU
  fusion edit.
- Exact-diff (NSTGV @ nx=129, `kinetic_energy.d`/`entropy.d`/21×`Q*.vtr`):
  byte-identical, fused vs. unfused.
- `test_etgv.py`: passed (both edits).
- NSTGV SLAU smoke (`nx=128`, 50k steps): passed, identical `ke/ke0` (both
  edits).
- Full suite (`pytest ouxsbli/tests/`, 24 tests including the 2 new ones):
  **24/24 passed**, ~21.6 min.

### Outstanding

- **GH200 performance verification for this SLAU fusion** — see above;
  this is the actual blocker before committing.
- Phase A's original ask (KEEP `threadsE` N=1-vs-N=2 on GH200) — still
  outstanding, unrelated to and not resolved by this update.
- Hybrid+NS fusion — deferred; no periodic Hybrid+NS case exists today to
  benefit from or validate it, and it would need correct viscous-stress
  computation under both of Hybrid's per-point KEEP/SLAU sensor branches.
- Boundary-aware fusion (BC_X/Y/Z=True, where SBLI/TBL run), 2D solver
  fusion (no 6th-order 2D viscous kernel exists), and curvilinear fusion
  (mismatched conv/visc thread-block shapes) — all explicitly out of scope
  for this pass, each needing its own follow-up plan.

## Update — 2026-08-05 (continued again): GH200 confirms the SLAU fusion win; committed; Phase B register-pressure fix applied to both kernels

New GH200 captures landed (`my_report_ncu_{fuse,nofuse}_20260805_slau.ncu-rep`,
full production grid, 513³, real GH200 120GB — confirmed via session
metadata, same `SoA`-vs-`stable` branch-checkout comparison methodology as
before). **This resolves the open question from the previous update.**

### SLAU fusion: confirmed win, not a regression

| Direction | Fused | Unfused (conv+visc) | Reduction |
|---|---|---|---|
| x | 29.6 ms | 42.3 ms | 30.0% |
| y | 33.3 ms | 48.4 ms | 31.3% |
| z | 36.7 ms | 57.3 ms | 36.1% |
| **Total** | **99.5 ms** | **148.1 ms** | **32.8%** |

This directly contradicts the earlier RTX4060 signal (~9% slower). Occupancy
data explains why: on RTX4060 the fused kernel's added shared-memory tiles
pushed the binding occupancy constraint below the hardware block cap; on
GH200 (far more registers/shared-mem per SM), that never binds — the SLAU
fusion behaves like KEEP's fusion, a clean win from eliminating the
intermediate global-memory round-trip. **RTX4060 was confirmed not to be a
valid proxy for this decision.** Committed as two commits: `127dd00` (Roe
removal) and `24bc57f` (SLAU fusion + NS pytest coverage), both including the
before/after numbers above in the commit message.

### Next lever, identified directly from the GH200 data

All three fused kernels (x, y, z) land at the **same occupancy ceiling**: 24
resident warps/SM out of 64 max (37.5%), matching measured occupancy almost
exactly (36.4/36.6/36.7%) — x is register-bound (74 regs/thread → 24 blocks),
y/z are register- *and* shared-mem-bound simultaneously (6 blocks × 4
warps). Every kernel's dominant stall reason is `long_scoreboard` (memory
latency, 3.9–4.3 cycles/instruction) by a wide margin over any other stall
category — meaning more resident warps would directly help hide that
latency. This matches the exact trigger condition for Phase B item #1 from
the original plan ("Block Limit Registers is the binding occupancy
constraint") — not a blind guess.

**Applied to both fused kernels** (`calc_keep_visc_kernel.f90.fypp` and
`calc_slau_visc_kernel.f90.fypp`, all three directions each): the `utau_sum`
sequential-accumulation chain — previously threaded live from the first
`calc_tau_straight` block through the final heat-flux block before a single
deferred `flux(5) = flux(5) - (utau_sum + kTx)` — now commits each viscous
contribution to `flux(5)` immediately inside its own block
(`flux(5) = flux(5) - utau`), removing `utau_sum`'s cross-block live range.

**Correctness bar shifts here, deliberately.** Unlike Phase A's launch-config
change (byte-identical by construction — same math, different thread
mapping) and the SLAU generalization above (byte-identical — same
computation, just fused into one launch), this changes floating-point
*evaluation order*: `a - b - c - d` is not bit-identical to `a - (b+c+d)` in
IEEE 754 in general. Verified accordingly, not by relaxing scrutiny but by
using the right tool for what actually changed:
- `kinetic_energy.d` (both DHIT/KEEP @ nx=70 and NSTGV/SLAU @ nx=129):
  **exactly identical** (domain-averaged scalar, insensitive at this
  precision).
- Raw `Q*.vtr` field snapshots: **not** byte-identical (as expected), but
  quantified directly — max relative difference ~5.5×10⁻¹⁰ (v-velocity,
  DHIT) and ~9.4×10⁻¹⁰ (w-velocity, NSTGV), i.e. rounding-level, not a
  physics-affecting change. (Some fields matched exactly even at the raw
  level — the reordering doesn't always change the rounded result, and the
  single-precision VTK output quantizes most of what's left.)
- Build-check (8/8), `test_etgv.py`, and the full pytest suite (24/24,
  including both `test_nstgv_ke_eps.py` energy-budget checks) all pass.

**Not yet committed** — this needs the same GH200 `ncu` re-capture as
before (same kernel names, same metrics: registers/thread, occupancy
limiter, `gpu__time_duration.sum`) to confirm the expected occupancy/latency
benefit actually materializes before deciding to keep it. Expect it might be
a *small* win at best: eliminating one live scalar register is unlikely to
single-handedly move the register-limit needle much (74→~73 registers,
roughly) — this change is more of a "measure and see" than a slam dunk, per
its own original framing in the Phase B menu. If GH200 shows no
improvement, this is cheap to revert (a mechanical, well-isolated diff).

### Outstanding (updated)

- **GH200 re-profiling of the register-pressure fix** — the actual blocker
  before committing this specific change. Same `ncu` command pattern as
  before (`--kernel-name regex:calc_keepv_x_in` / `calc_slauv_x_in` etc.,
  `--metrics launch__registers_per_thread,launch__shared_mem_per_block_static,sm__warps_active.avg.pct_of_peak_sustained_active,gpu__time_duration.sum,launch__occupancy_limit_registers,launch__occupancy_limit_shared_mem,launch__occupancy_limit_blocks,launch__occupancy_limit_warps`).
- Phase A's original ask (KEEP `threadsE` N=1-vs-N=2 on GH200) — still
  outstanding, independent of everything above.
- Hybrid+NS fusion, boundary-aware fusion, 2D solver fusion, curvilinear
  fusion — all still out of scope/deferred, unchanged from the previous
  update.
