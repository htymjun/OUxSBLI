# Verify default viscous stencil is Central (not ME4-Base) for 4th/6th order

## Context

The user wants the default viscous discretization for `VISC_ORDER = 4` and `VISC_ORDER = 6`
to be the basic central-difference stencil ("Central"/Cent), not the ME4-Base scheme, unless a
case explicitly opts into ME4-Base.

This request lands on top of an already-staged (uncommitted) refactor — documented in
`docs/plans/2026-06-26-now-this-code-has-stateless-cascade.md` — that:
- unified `calc_visc_cent4.f90` / `calc_visc_cent6.f90` into `3D_solver/src/calc_visc_cent.f90.fypp`
- unified `load_smem_visc_cent4.f90` / `load_smem_visc_cent6.f90` into `3D_solver/src/load_smem_visc_cent.f90.fypp`
- renamed `calc_visc4[_internal].f90.fypp` → `calc_visc_high[_internal].f90.fypp`, which now dispatch
  between ME4-Base and Central via a new `VISC_STENCIL` config.fypp variable
- updated `3D_solver/CMakeLists.txt` accordingly and deleted the old `calc_visc4*`/`*_cent4`/`*_cent6` files

**Finding from investigation:** the dispatch logic in the currently staged
[`calc_visc_high.f90.fypp`](3D_solver/src/calc_visc_high.f90.fypp) and
[`calc_visc_high_internal.f90.fypp`](3D_solver/src/calc_visc_high_internal.f90.fypp) already defaults
`VISC_STENCIL` to `'Cent'` (not `'ME4Base'`) when a case's `config.fypp` doesn't set it:

```fortran
#:if not defined('VISC_STENCIL')
  #:set VISC_STENCIL = 'Cent'
#:endif
...
#:if VISC_ORDER == 4
  ! always Cent4 — ME4-Base is not offered at order 4
#:elif VISC_STENCIL == 'ME4Base'
  ! opt-in ME4-Base
#:else
  ! Central (Cent6) — this is what any non-'ME4Base' value hits, including the default
#:endif
```

None of the case `config.fypp` files (checked all of `3D_solver/*/config.fypp`) set
`VISC_STENCIL = 'ME4Base'`, so every case already gets Central by default:
- `VISC_ORDER == 4` → always Cent4 (hard-wired, no override possible)
- `VISC_ORDER == 6` (DHIT, NSTGV, SBLI, TBL are the `VISC in ('NS','LES')` cases exercising this) → Central (Cent6) by default, ME4-Base only via explicit opt-in

This was confirmed by dry-running `fypp` over all 4 affected cases (DHIT, NSTGV, SBLI, TBL): each
generates `use load_smem_visc_cent6` / `module load_smem_visc_cent6`, not `load_smem_visc_me4_base`.

**Naming inconsistency found (user flagged, to fix):** the module/subroutine family is named
`*_cent*` (`calc_visc_cent.f90.fypp`, `load_smem_visc_cent4`/`load_smem_visc_cent6`), but the design
doc (`docs/plans/2026-06-26-...md`) documents the opt-in config value as `VISC_STENCIL = 'Central'`,
while the actual staged default literal is `'Cent'`. The two names aren't used consistently. Since the
dispatch only ever tests `VISC_STENCIL == 'ME4Base'` (anything else — including the default — falls
through to Central), the mismatch is currently harmless functionally, but it's confusing for anyone
reading config.fypp docs vs. the code. **Fix:** standardize on `'Central'` (matches the documented
config-facing name) as the default literal in both:
- `3D_solver/src/calc_visc_high.f90.fypp` — `#:set VISC_STENCIL = 'Cent'` → `'Central'`
- `3D_solver/src/calc_visc_high_internal.f90.fypp` — same change

This is a one-line default-literal change in each file, purely cosmetic/consistency (no dispatch-logic
change, since the `else` branch still catches it) — no case `config.fypp` needs updating since none of
them set `VISC_STENCIL` explicitly.

**Conclusion:** aside from that naming fix, no further source edit is needed for the default-stencil
behavior itself — it's already correct in the current staged tree. What's still outstanding, matching
the user's "check success of compilation" ask, is that this refactor (per its own status note, "ZERO
FILES WRITTEN YET" as of 2026-06-30) has never actually been build-verified with `nvfortran`/CMake.
fypp-level syntax was just confirmed OK, but full compilation has not been tried — and per the user's
follow-up, both the default (Central) *and* the explicit ME4-Base opt-in path need to be build-verified,
not just the default.

## Plan

1. **Unify the `'Cent'`/`'Central'` naming** — change the default literal in both files from `'Cent'`
   to `'Central'`:
   - `3D_solver/src/calc_visc_high.f90.fypp:6`
   - `3D_solver/src/calc_visc_high_internal.f90.fypp:6`

   No other files need updating for this — no `config.fypp` sets `VISC_STENCIL` explicitly, and the
   dispatch logic only special-cases `'ME4Base'`, so this is a pure naming-consistency fix with no
   behavior change.

2. **Build-verify each of the 4 affected cases with the default (Central)** — the only 3D cases with
   `VISC in ('NS','LES')`, i.e. the ones that actually pull in `calc_visc_high`:
   ```bash
   cd 3D_solver/DHIT  && cmake -B build && cmake --build build -j
   cd 3D_solver/NSTGV && cmake -B build && cmake --build build -j   # has a stale build/ from before the rename — reconfigure will pick up new file lists; if link/compile errors reference old calc_visc4 symbols, rm -rf build and reconfigure clean
   cd 3D_solver/SBLI  && cmake -B build && cmake --build build -j
   cd 3D_solver/TBL   && cmake -B build && cmake --build build -j   # has a stale build/ too; same caveat, and this is the LES/Hybrid path exercising the `is_les` branches
   ```
   NSTGV and TBL already have `build/` directories predating this refactor (their generated
   `calc_visc4.f90`/`calc_visc4_internal.f90` are still sitting there) — `cmake --build` should
   regenerate/relink fine since it's driven by `_SHARED_FYPP`/`_BASE` lists in `CMakeLists.txt`, but if
   stale object/module files cause spurious link errors, do a clean `rm -rf build` first.

3. **Build-verify the explicit ME4-Base opt-in path too**, on one NS case and one LES case (LES has
   extra `mut`/`qc2`/`Hsgs` codepaths worth exercising separately from plain NS):
   ```bash
   # temporarily append to config.fypp, build, then revert the temporary line
   echo "#:set VISC_STENCIL = 'ME4Base'" >> 3D_solver/SBLI/config.fypp
   cd 3D_solver/SBLI && rm -rf build && cmake -B build && cmake --build build -j
   # revert config.fypp, repeat for TBL
   ```
   Confirm via `grep` in the generated `build/calc_visc_high.f90` that `load_smem_visc_me4_base` /
   `calc_visc_me4_base` got selected for these two test builds, then revert both `config.fypp` files
   back to not setting `VISC_STENCIL` (so the shipped default stays Central) and re-run their default
   builds from step 2 to leave the tree in the intended end state.

4. **Triage any compile errors** that come up in steps 2 or 3. Likely candidates to check first if
   something fails:
   - `3D_solver/src/calc_visc_cent.f90.fypp` (parametrized `interp${N}$`/`diff${N}$`, `calc_tau_*`
     for N=4 and N=6) — included into both `calc_visc_high.f90.fypp` and `calc_visc_high_internal.f90.fypp`
   - `3D_solver/src/load_smem_visc_cent.f90.fypp` — generates `load_smem_visc_cent4`/`load_smem_visc_cent6`
     modules; the source plan doc calls out specific `inv_dy`/`inv_dx`/`inv_dz` bugs it claims to have
     already fixed relative to the old `load_smem_visc_cent6.f90` — worth double-checking against the
     generated output for DHIT/SBLI/TBL if a NaN/wrong-result issue (not just compile issue) shows up
   - Interior-index conditions in `calc_visc_high.f90.fypp`/`calc_visc_high_internal.f90.fypp` per
     kernel (E/F/G) differ between Cent4/ME4-Base/Cent6 — verify against the table in the design doc
     if a shared-memory bounds/index mismatch appears

5. **Report results**: which cases built cleanly under Central and under ME4-Base, which needed fixes
   (and what was fixed), and confirm via a quick grep of the generated `build/*.f90` that
   `load_smem_visc_cent6`/`calc_tau_straight` (not `load_smem_visc_me4_base`) is what's actually
   compiled in for the default builds of DHIT/NSTGV/SBLI/TBL.

## Out of scope

- `2D_solver/` — its viscous path (`calc_visc4.f90.fypp`) has no `VISC_STENCIL`/Central option at all
  (it's unconditionally ME4-Base at `VISC_ORDER=4`, with no 6th-order viscous path). The user's request
  and the existing design doc both scope this change to the 3D solver only; extending Central-by-default
  to 2D would be new scope, not this fix.
- `3D_solver_curv/` — untouched by this refactor, no `calc_visc4`/`VISC_STENCIL` references found there.
- Not touching `docs/plans/2026-06-26-now-this-code-has-stateless-cascade.md`'s stale "ME4-Base
  (default)" wording (and its `VISC_STENCIL = 'Central'` opt-in example, which is what the code's
  literal will now actually match) — flagging it here for awareness, but leaving the doc as-is unless
  asked.
