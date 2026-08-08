# Self-review + real-profiling-backed optimization report + CLAUDE.md update

## Context

The `calc_div`/DRAM-precompute work (see prior sessions) is functionally
done and runtime is acceptable. The user now wants: (1) a bug self-review,
fixing anything found; (2) further optimization ideas for the **high-order**
path specifically (not just minor order-2 memory tweaks); (3) `CLAUDE.md`
updated.

Two corrections from the user's feedback on my first draft of this plan:
- `calc_div_2` is **not** dead code to remove — it's intentionally kept
  because `calc_hybrid`'s `calc_Ducros` shock sensor (used for
  SLAU/Roe/Hybrid schemes) also needs a divergence field, and `ux/vy/wz`
  could be reused there for 2nd-order builds too. Do not suggest removing it.
- Real profiling data already exists at `3D_solver/nsys_ncu/*.ncu-rep`
  (Nsight Compute captures, `CC 9.0` — an H100, not this session's RTX 4060)
  — I should mine that instead of speculating about performance.

### What the profiling data already shows (from `my_report_ncu0727_2.ncu-rep`,
the most recent capture, NSTGV-shaped grid `513³`, via
`ncu --import ... --print-summary per-kernel --section Occupancy`)

| Kernel | Registers | Theoretical/Achieved Occ. | Binding limiter |
|---|---|---|---|
| `calc_div_6_in` | 60 | 50% / 43.25% | **Registers** (block limit = 4, vs. 32 for SM/shared-mem, 8 for warps) |
| `calc_Ev6_in` | 48 | 50% / 49.2% | **Block-count-per-SM hardware cap** (block limit SM=32 blocks, but registers/shared-mem allow 40 — `threadsEv` is only `(32,1,1)`, i.e. 1 warp/block, so the 32-blocks/SM ceiling caps total resident warps at 32 instead of the 64 the registers/smem would otherwise allow) |
| `calc_Fv6_in`/`calc_Gv6_in` | 54/52 | 56.25% / ~53.5% | **Registers** (block limit = 9, vs. 10 for shared-mem) |

`calc_div_6_in` itself costs **~3.83ms/call** (×3 per RK3 step ≈ 11.5ms/step)
— comparable to `calc_quantities_T_3D` (4.41ms) — a real, measurable new
cost from this feature, not hidden by the smem/register wins already made.

This gives 3 concrete, data-backed optimization candidates (see below) —
**report only, do not implement**, since they involve kernel splits/launch
config changes I haven't validated and the user only asked me to *find and
report* optimization points (as opposed to bugs, which the user did ask me
to fix directly).

## Plan

### 1. Self-review pass + fix bugs directly
Re-read the current on-disk state (not memory) of `calc_div.f90`,
`calc_visc_cent.f90.fypp`, `calc_visc_high.f90.fypp`,
`calc_visc_high_internal.f90.fypp`, `load_smem_visc_cent.f90.fypp`,
`calc_flux_base.f90.fypp`, `preprocess.f90.fypp`, `calc_time_dev.f90.fypp`,
`CMakeLists.txt`. Focus on:
- Argument-order consistency across all 4 `calc_tau_straight_s` call sites
  (Fv/Gv × 2 files) and between the order-4/order-6 branches.
- Margin/index consistency between `calc_div`'s per-component gates and
  each consumer's own `interior` condition (bit-exact-verified for order 6
  already; re-derive for order 4 where `io_v=1` shifts every bound).
- `_koff`/`is_les` paths (unreachable for all 8 current cases) — confirm
  they still compile and are internally self-consistent.
- The device-array plumbing added to `CMakeLists.txt`/`preprocess.f90.fypp`/
  `calc_time_dev.f90.fypp`.
Fix anything concretely wrong; note what changed and why.

### 2. Close the order-4 testing gap (the highest-value thing to actually
run, not just read — no current case exercises `VISC_ORDER==4` for NS/LES)
Temporarily add `#:set VISC_ORDER = 4` to `3D_solver/DHIT/config.fypp`
(`VISC_ORDER` is documented in `CLAUDE.md` as an independent override of
`ORDER`), rebuild+run a short DHIT integration (temporary `nt`/`np`
reduction as in prior sessions), confirm finite/sane output. Revert both
temporary edits afterward regardless of outcome. If it fails to compile or
produces garbage, fix it in place and re-verify.

### 3. Write the report (in my final chat response, not a new file)
- Bugs found/fixed in step 1 (or "none found" if the review turns up
  nothing beyond what's already validated).
- Optimization ideas, clearly marked as reported-not-implemented:
  a. **`calc_div` register pressure**: 60 registers/thread, occupancy
     capped by registers at 50% theoretical. Candidate: split into 3
     independent single-component kernels (`calc_div_ux`/`_vy`/`_wz`), each
     gated only by its own margin with zero branching (today's single
     kernel carries 3 independent conditional branches in one register
     footprint) — likely reduces registers substantially at the cost of 2
     extra kernel launches.
  b. **`calc_Ev6_in` block shape**: occupancy is capped by the SM's
     max-blocks-per-SM hardware limit, not registers/shared-mem (both have
     large headroom) — `threadsEv` is `(32,1,1)` (1 warp/block) while
     `threadsFv`/`threadsGv` use 4-warp blocks. Widening `threadsEv`'s y/z
     extent (mod_globals.f90, per-case) could raise achieved occupancy
     without touching kernel code at all.
  c. **`calc_Fv6_in`/`calc_Gv6_in` register pressure**: still register-bound
     (54/52 regs, block limit 9 vs. 10 for shared mem) — shaving a few more
     registers (e.g. via the same kind of restructuring as (a), or
     revisiting whether all `H(:)` LES-path locals are needed when
     `is_les` is false at compile time) could unlock one more resident
     block/SM.
  d. **`calc_Ducros` reuse of `ux/vy/wz`**: `calc_Ducros` (only relevant for
     SBLI among the 8 cases — the one case using both high-order viscous
     and SLAU) currently recomputes its own 9-term velocity-gradient tensor
     from scratch via always-2nd-order central differences to get both
     `div(u)` and `curl(u)` (2.02ms/call). `calc_div` already computes 3 of
     those 9 terms (`ux,vy,wz`, i.e. exactly what's needed for `div(u)`) —
     but at **6th order**, not 2nd, and only in the interior margin (not the
     3-cell-wide edge band `calc_Ducros` also covers). Reusing them would
     change the Ducros sensor's numerical values (mixed-order dilatation vs.
     vorticity) and requires a fallback for the edge band — a real,
     targeted idea, but a physics-affecting change that needs explicit
     buy-in before implementing, not a free win.

### 4. Update `CLAUDE.md`
- Add `calc_div.f90` to the "Notice" section's bottleneck list (runs every
  RK stage right after `calc_quantities_T_3D`).
- Add a short paragraph documenting the `ux,vy,wz` precompute-once
  architecture: what it replaced, which files are involved
  (`calc_div.f90`, `load_smem_visc_cent.f90.fypp`'s reduced outputs,
  `calc_tau_straight_s`/`interp2_${N}$_s` in `calc_visc_cent.f90.fypp`), and
  its one documented limitation (COMMZ+NS/LES unwired — no current case
  needs it).
- Note that `calc_div_2` is intentionally unwired today but retained for a
  possible future `calc_Ducros` reuse path at 2nd order (per user
  correction — don't describe it as dead code).
- Mention `3D_solver/nsys_ncu/` as where Nsight Compute profiling captures
  live, since it's how the optimization ideas above were substantiated.

## Verification
- Rebuild DHIT, NSTGV, SBLI, TBL, ETGV after any fix from step 1/2.
- Re-run the DHIT bit-exact regression check if any fix touches numerics.
- Leave the repo clean afterward (no stray build/data dirs, no leftover
  `config.fypp`/`mod_globals.f90` test edits).
