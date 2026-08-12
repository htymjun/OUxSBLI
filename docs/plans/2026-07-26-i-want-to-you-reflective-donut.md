# Reduce register usage of calc_visc_high(.f90.fypp)/calc_visc_high_internal(.f90.fypp) to ≤40 (NS mode, NSTGV 513³) — round 2

## Context

Continuing the same register-reduction task. Round 1 (already implemented and kept in the working tree)
folded in the untracked `_int` prototypes' single-`uvw`-shared-buffer technique and restructured the
`interp6_s`/`diff6_s` scalar-stencil call sites to pass pre-summed pairs instead of 6 raw loads (bit-exact).
A third idea (uncaching `dxi`/`dyj`/`dzk`, recomputing `dx(i)`/`dy(j)`/`dz(k)` at each use site) was tried,
measured, found to be a net regression (+6 aggregate registers across the 12 kernel variants, including a
+6 hit on `calc_Ev6`), and reverted. None of this reached ≤40, but it is a real (if partial) improvement
for the kernels NSTGV actually launches. You've asked to keep pushing — you believe ≤40 is achievable — so
this round tries two more concrete levers before considering anything bigger.

**Current state (round 1 result, kept), for the kernels NSTGV actually launches at runtime** (confirmed via
`3D_solver/NSTGV/build/calc_flux_base.f90:87-92`: only `calc_Ev6_in`/`calc_Fv6_in`/`calc_Gv6_in`, no `_koff`
since `COMMZ=False`):

| Kernel | Before round 1 | After round 1 (current) |
|---|---|---|
| `calc_Ev6_in` | 56 | 56 |
| `calc_Fv6_in` | 52 | 48 |
| `calc_Gv6_in` | 48 | 48 |

(Full 12-kernel table incl. `_koff` and the boundary-aware `calc_visc_high` module is in the prior
conversation; `_koff`/boundary kernels are dead code for NSTGV specifically but still compiled and were
measured too — some regressed slightly, some improved; net aggregate across all 12 was roughly flat.)
No spilling anywhere (`LOCAL:0`, `STACK:0` on every kernel) — must stay that way.

Target: ≤40 registers, NS mode only, restructuring only (confirmed again: no `-Mcuda=maxregcount` fallback).

## Round 2 approach (per your answers)

Apply both, rebuild NSTGV, measure via `cuobjdump -res-usage`, keep only what nets an improvement (same
measure-then-decide discipline as round 1 — a plausible-sounding change here already regressed once):

1. **Recompute `mui` fresh in stage 2 and stage 3 instead of caching it across all 3 stages.** `mui` is
   currently declared at subroutine scope and computed once in stage 1 via `interp${VISC_ORDER}$_s` on `mu`'s
   6-point neighbor stencil, then reused unchanged by stage 2/3's `calc_tau_cross` calls. Since it depends
   only on `mu`'s spatial neighbors (not on stage), recompute the identical `interp${VISC_ORDER}$_s(...)` call
   locally inside stage 2's and stage 3's own `block`, using a stage-local variable instead of the persistent
   one. Bit-exact (same formula, same inputs) — costs a few redundant `mu` reads (likely still L1/L2-hot from
   stage 1) in exchange for not holding `mui` live across 2 `syncthreads()` barriers.

2. **Stop holding `kTx`/`utau_sum` as persistent registers across all 3 stages.** Currently: `kTx` computed
   once in stage 1 and held until the final combination in stage 3; `utau_sum` accumulated stage-by-stage in
   a register (`utau_sum = utau`, then `utau_sum = utau_sum + utau` ×2) and combined with `kTx` in one final
   `__stcg`/`__ldlu` read-modify-write to component 5 of E/F/G. Instead, write each stage's own contribution
   directly to component 5 via its own `__stcg`/`__ldlu` RMW call (mirroring the existing per-stage RMW
   pattern already used for components 2/3/4), immediately after that stage computes it — stage 1 does two
   RMWs (kTx, then utau1) or one combined, stage 2 does one (utau2), stage 3 does one (utau3). This is where
   the ULP-level reassociation you approved comes in: the same four quantities (kTx + 3 partial utaus) end up
   subtracted from E/F/G(...,5) as separate floating-point operations instead of one grouped sum — same math,
   negligibly different rounding, and no longer needs `kTx`/`utau_sum` to survive past the stage that produces
   them.

Apply to `calc_visc_high_internal.f90.fypp` first (the file that actually matters for NSTGV), measure, then
apply the same techniques to `calc_visc_high.f90.fypp` (including its boundary-fallback `else` branches, which
also use `mui`/`kTx`/`utau_sum` the same way) for consistency, per original task scope.

**If ≤40 is reached this way:** stop, report final numbers, done.

**If not:** per your answer, do NOT proceed to splitting each kernel into 3 per-stage launches (which would
touch `calc_flux_base.f90.fypp`'s launch sites, triple kernel-launch count, and redundantly recompute
`mui`/`Tx` per stage-kernel) without checking in again first. Instead stop and report the best achieved
number, with this bigger option described as the next escalation if you want to authorize it.

## Verification

- `cmake --build 3D_solver/NSTGV/build -j` after each of the two changes (or together, then bisect if the
  combined result is mixed — same approach as round 1, where one candidate change had to be measured and
  reverted after the fact).
- `cuobjdump -res-usage 3D_solver/NSTGV/build/a.out`, filtered to `calc_Ev6_in`/`calc_Fv6_in`/`calc_Gv6_in`
  (and the `calc_visc_high` / `_koff` variants for completeness) — confirm `LOCAL:0`/`STACK:0` throughout
  (no spilling introduced).
- No simulation run (per your original instruction) — compile + register-count evidence only. Numerical
  verification (confirming the ULP-level reassociation doesn't affect test tolerances) remains a follow-up
  for a later session via `pytest ouxsbli/tests/` or a short NSTGV run.

## Files touched

- `3D_solver/src/calc_visc_high.f90.fypp`
- `3D_solver/src/calc_visc_high_internal.f90.fypp`

(`load_smem_visc_cent.f90.fypp` already updated in round 1, not touched further this round unless a
follow-on issue turns up.) No changes to `config.fypp`, `mod_globals.f90`, `calc_flux_base.f90.fypp`, or
CMake wiring in this round.
