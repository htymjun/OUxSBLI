# 2026-08-06 ncu analysis → occupancy + smem-diet + ILP optimizations

## Context

The two new GH200 captures (`3D_solver/ncu/my_report_ncu_fuse_{keep,slau}_20260806.ncu-rep`,
NSTGV 513³, CC 9.0, 132 SMs) are the re-profiling that commit `6f0c722` ("Break utau_sum
accumulation chain") said it needed before being judged. Verdict against the 0805 fuse captures:

| Kernel | 0805 (pre-fix) | 0806 (post-fix) | regs pre→post |
|---|---|---|---|
| keepv x/y/z | 8.44 / 9.52 / 10.78 ms | 8.83 / 9.77 / 11.13 ms (**+3–5%**) | 95/94/95 → 96/95/96 |
| slauv x/y/z | 9.87 / 11.08 / 12.22 ms | 9.98 / 11.11 / 12.23 ms (flat) | 74/78/80 → 76/78/80 |

**The utau_sum change regressed KEEP, did nothing for SLAU, and raised registers.** Jun confirmed:
revert it (Step 1).

**Where a KEEP RK stage goes** (SLAU similar, plus calc_ducros 2.03 ms):

| Kernel | Time | DRAM% | Occ ach/theo | Limiter |
|---|---|---|---|---|
| keepv x+y+z (fused conv+visc) | 29.7 ms (~61%) | 46–60 | 30 / 31.25% | **regs (96)**; y/z smem 2nd |
| calc_step1 / calc_step2_3 | 9.0 / 12.0 ms (~20%) | **91** / **79.7** | 69/75 · 58/62.5% | step2_3: regs 42→48 alloc |
| calc_quantities_T | 4.4 ms | 75 | 59/75% | latency (long_scoreboard) |
| calc_div ux/vy/wz | 0.84/1.37/1.49 ms | 64/38/36 | 78 · 51 · 40% | vy/wz: regs (48, 60) |
| calc_ducros (SLAU runs) | 2.0 ms | 53 | 55/62.5% | **L2 88% busy** |

Every hot kernel is long-scoreboard-stalled with DRAM unsaturated → resident warps are the lever;
registers (and for SLAU y/z, shared memory) cap them everywhere.

**Phase A closure** (docs/plans/2026-08-04-…-phase-a-status.md): at 96 regs, `threadsE=(32,2,1)`
gives 10 blocks × 2 warps = the same 20 warps — registers bind before the 32-block HW cap, so the
threadsE lever is moot until regs < ~64. Answered from production data; no dedicated capture needed.

Scope per Jun: revert; quick register caps; fused-kernel launch_bounds sweep; **SLAU y/z smem diet
(privatize right-state tiles — "implement the idea")**; **calc_ducros tiling ("try it once")**;
**calc_quantities_T ILP ("I like ILP")**.

## Key occupancy math (H100/GH200: 64K regs/SM, 64 warps/SM, smem config 135168 B, 1024 B/block driver)

- CUDA Fortran: `attributes(global) launch_bounds(<threads>, <minBlocks/SM>) subroutine …`. Prefer
  parameter expressions (`threadsF%x*threadsF%y*threadsF%z`) — nvfortran already accepts `threadsE%x`
  as a constant in the smem bounds; fall back to fypp literals if rejected.
- Tile unit in fused y/z kernels: 32×9×8 B = 2304 B/array (sy = threads%y + 2·io + 1, io=2).
- SLAU y/z today: 7 halo arrays (rho,u,v,w,p + 2 cross-derivs) = 16128 B **+ 5 right-state arrays
  (rhor…pr, 32×4×8 = 1024 B each) = 5120 B** → 21.25 KB → 6 blocks; regs 78–80 → 6 blocks. Dual-capped, 24 warps.
- **Right-state tiles are thread-private**: written at `idx_r` and read back at the same `idx_r`
  only ([calc_slau_visc_kernel.f90.fypp:210,217](3D_solver/src/calc_slau_visc_kernel.f90.fypp#L210)).
  → local scalars. The `syncthreads()` between recon and SLAU (lines 212/324/436) protects nothing
  after that (already vestigial — inherited from the unfused kernel's tile-overwrite design) → remove.
  Post-diet smem 16128+1024 → **7 blocks = 28 warps (43.75%)** with `launch_bounds(128,7)` (≤72 regs).
  8 blocks is out of reach by 256 B (needs ≤15872 B static). Bit-identical change.
- KEEP y/z: 8 arrays = 18432 B → 6 blocks; regs 96 → 5 blocks binding. `launch_bounds(128,6)` (≤85→80
  regs) → 24 warps (37.5%, from 20). 7 blocks would need dropping the T tile (146 B over) — optional
  risky rung only.
- x-kernels (32 thr, smem trivial): KEEP 96 regs → 20 warps; SLAU 76 → 24. Rungs: `launch_bounds(32,24)`
  ≤85 regs, `(32,32)` ≤64 regs → 32 warps (50%).
- step2_3 (128 thr): 42→48 alloc → 40 warps; `launch_bounds(128,12)` caps 40 regs → 48 warps (75%),
  matching step1 (91% DRAM). calc_div vy `(256,6)`→40 regs (75%); wz `(256,5)`→48 regs first, escalate.
- Risk: forced caps spill to local. DRAM headroom 40–50% absorbs moderate spill, but it's empirical —
  ptxinfo (already on) reports spill bytes at build; GH200 decides each rung. (Smem carveout ↑ to 228 KB
  is an alternative for SLAU y/z but steals L1 the strided mu/T/ux/vy reads rely on — not this pass.)

## Implementation steps (each independently gated; all changes bit-identical by construction)

### Step 1 — Revert 6f0c722 + records
- `git revert 6f0c722`; write `docs/2026-08-06-utau-sum-chain-break-failed-optimization.md` with the
  table above + hypothesis (early commits to flux(5) serialize FADDs on the accumulator, removed
  scheduling freedom; live range didn't shrink — regs rose). Append verdict + Phase A closure to the
  2026-08-04 status doc. Re-run gates once to re-baseline.

### Step 2 — Quick register caps
- [3D_solver/src/calc_steps.f90](3D_solver/src/calc_steps.f90): `launch_bounds(128,12)` on
  `calc_step2_3` only (step1 already 40; RK4 kernels untouched — unprofiled).
- [3D_solver/src/calc_div.f90](3D_solver/src/calc_div.f90): confirm launch shape in-file (ncu: (32,4,2)=256),
  then `launch_bounds(256,6)` on `calc_div_vy_6_in`, `(256,5)` on `calc_div_wz_6_in`.
  Leave `*_4_in`/`*_2` variants alone (untested config branches per CLAUDE.md).

### Step 3 — SLAU right-state privatization + vestigial-barrier removal
[calc_slau_visc_kernel.f90.fypp](3D_solver/src/calc_slau_visc_kernel.f90.fypp), all three directions:
- Replace `rhor/ur/vr/wr/pr` shared tiles with local scalars (mirror of existing `rhol…pl`); delete
  `idx_r/offset_*r/sxr/syr/szr`; drop the second `syncthreads()` and merge the recon+SLAU blocks back
  into one guarded region.
- Registers will rise (~5 live doubles); pair with `launch_bounds(128,7)` on y/z (72-reg cap → 28 warps)
  and the Step-4 knob on x. Verify generated `.f90` in NSTGV's build dir.

### Step 4 — launch_bounds sweep knobs on both fused kernels
- Knobs in `config.fypp` with template defaults via `#:if not defined(...)`:
  `FUSED_LB_MIN_BLOCKS_X` (KEEP default 24, SLAU 32), `FUSED_LB_MIN_BLOCKS_YZ` (KEEP 6, SLAU 7).
- [calc_keep_visc_kernel.f90.fypp](3D_solver/src/calc_keep_visc_kernel.f90.fypp): x `(32, X)`, y/z `(128, YZ)`.
- [calc_slau_visc_kernel.f90.fypp](3D_solver/src/calc_slau_visc_kernel.f90.fypp): x `(32, X)`; y/z from Step 3.
- Sweep rungs on GH200 via per-case clones varying config.fypp (established `NSTGV_keep128_warpilp`
  pattern): KEEP x {24, 32}; SLAU x {32}; KEEP yz {6}; SLAU yz {7}. Check DHIT's `mod_globals.f90`
  thread shapes match NSTGV's before reusing numbers.

### Step 5 — calc_ducros smem tiling (one-shot experiment)
[calc_hybrid.f90](3D_solver/src/calc_hybrid.f90) + call sites in
[calc_flux_base.f90.fypp](3D_solver/src/calc_flux_base.f90.fypp) (`<<<blocks,threads>>>` = (32,4,1), one
k-plane/block → zero z-reuse today):
- Reshape to (32,4,2) (local dim3 at the call site; promote to `mod_globals` convention if it wins),
  cooperative-load u,v,w tiles with ±1 halo: 34×6×4×8 B×3 = 19.58 KB → 6 blocks; `launch_bounds(256,6)`
  (≤40 regs) → 48 warps (75% theo, up from 62.5) and ~2–3× fewer L1/L2 requests on an 88%-L2-bound kernel.
- Est. 2.03 → ~1.4 ms. Affects SLAU/Hybrid dispatch paths (NSTGV, SBLI, TBL — sensor values unchanged,
  bit-identical). If GH200 says flat/worse → failure doc + revert, per "try it once".

### Step 6 — calc_quantities_T ILP
`calc_quantities_T_3D`'s kernel in [src/calc_physical_quantities.f90](src/calc_physical_quantities.f90)
(~line 100–140): 2 points per thread along z (halve grid z), keeping per-point math identical;
`launch_bounds(256,6)` to hold 40 regs. Streaming kernel at 75% DRAM / 34.5 warp-cycles-per-inst —
ILP=2 doubles independent memory chains. Est. 4.44 → ~3.9 ms.

### Gates after every step (ouxsbli-optimize skill)
1. `bash scripts/build_all_cases.sh <repo>` (8/8; record ptxinfo spill bytes of edited kernels).
2. `pytest ouxsbli/tests/test_etgv.py -v`.
3. `python scripts/nstgv_slau_check.py <repo> --nt 500`.
4. Exact-diff (all steps here change no FP math → byte-identical required): short DHIT @ nx=70 and
   NSTGV @ nx=129 runs vs post-revert baseline (`cmp` on `Q*.vtr`, `kinetic_energy.d`, `entropy.d`).

### GH200 handoff (perf verdict is never local — RTX4060 proven invalid as proxy)
Append per-rung capture instructions to the status doc:
```
ncu --kernel-name "regex:calc_keepv_|calc_slauv_|calc_step|calc_div_|calc_ducros|calc_quantities" \
    --launch-count 3 --set full -o <tag>_$(date +%Y%m%d) mpiexec -np 2 ./a.out
```
Decision metrics: `gpu__time_duration.sum`, achieved occupancy, regs/thread, **Local Memory Spilling
Requests** (currently 0 — watch), DRAM%. Baseline = post-revert (≈ 0805 fuse numbers). Per rung:
commit via `scripts/record_optimization_outcome.sh` (win) or failure-doc + revert (loss).

## Expected effect (per RK stage, if rungs win; honest ranges)
step2_3 −1.5–2 ms (stages 2/3) · div −0.5 ms · fused kernels −2–5 ms (main uncertainty: spill cost)
· ducros −0.5 ms (SLAU) · qT −0.5 ms ⇒ plausible ~8–15% stage-time reduction; each piece
independently revertible.

## Future (not this session)
- SLAU/KEEP y/z beyond 7/6 blocks: needs dropping a halo tile (KEEP: T tile → 7 global reads/pt) — only
  after the sweep settles registers.
- E/F/G → in-kernel divergence into R: traffic shift from 91%-DRAM steps kernels into latency-bound fused
  kernels; a shfl/atomics variant already failed 2× (benchmark_dr_vs_efg.md) — any retry must use
  block-overlap recompute.
- set_bc_cyclic6_247 launch shape (121 µs @ 6% occ) — negligible.
