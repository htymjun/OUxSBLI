# Sweep ORDER / VISC / VISC_STENCIL for NSTGV compilation

## Context

`3D_solver/src/calc_visc_high.f90.fypp` / `calc_visc_high_internal.f90.fypp` are new files (renamed from
`calc_visc4[.f90.fypp]` / `calc_visc4_internal.f90.fypp`) that implement an in-progress refactor
(`docs/plans/2026-06-26-now-this-code-has-stateless-cascade.md`, approved 2026-06-30) unifying the
viscous-kernel stencils and adding a new `VISC_STENCIL` dispatch knob. `VISC_STENCIL` is not yet exposed
in any case's `config.fypp` — it only exists as an internal fypp variable with a fallback default inside
those two template files. The goal of this task is to validate that NSTGV compiles across the combinations
of `ORDER`, `VISC`, and `VISC_STENCIL` that exercise every distinct branch this refactor introduces, since
the templates have not yet been build-tested for most of these combinations.

Two pre-existing issues were found during investigation and explicitly scoped **out** of this task per user
decision:
- The fallback default for `VISC_STENCIL` (`'Central'`, set at `calc_visc_high.f90.fypp:6` and
  `calc_visc_high_internal.f90.fypp:6`) does not match what the approved refactor plan specifies
  (`'ME4Base'`). Leave as-is; only set `VISC_STENCIL` explicitly in NSTGV's `config.fypp` for the combos
  that need it.
- `ORDER=2` with `VISC=NS`/`LES` produces a hard compile error: `calc_visc_high.f90.fypp` has no guard
  excluding `VISC_ORDER==2`, so it generates `calc_Ev2`/`calc_Fv2`/`calc_Gv2` (+ `_koff` variants) that
  collide with the identically-named subroutines from the legacy `calc_visc2` module — both are
  unconditionally `use`d together in `3D_solver/src/calc_flux_base.f90.fypp:20-22` whenever
  `VISC in ('NS','LES')`. This is unrelated to `VISC_STENCIL` correctness, so `ORDER=2` combos are
  dropped from the sweep entirely (not attempted, not fixed).

## Sweep matrix (6 combinations)

`VISC_STENCIL` only affects behavior at `ORDER=6` (at `ORDER=4` the dispatch always uses Cent4 regardless
of `VISC_STENCIL`, per `calc_visc_high.f90.fypp:13-19`). The matrix below hits every distinct branch
(`is_les` toggle × order-4 Cent4 path × order-6 ME4Base/Central paths) exactly once:

| # | ORDER | VISC | VISC_STENCIL (config.fypp line added?) | Build dir |
|---|-------|------|------------------------------------------|-----------|
| 1 | 4 | 'NS'  | not set (irrelevant at ORDER=4) | `3D_solver/NSTGV/build_o4_ns` |
| 2 | 4 | 'LES' | not set (irrelevant at ORDER=4) | `3D_solver/NSTGV/build_o4_les` |
| 3 | 6 | 'NS'  | `'ME4Base'` | `3D_solver/NSTGV/build_o6_ns_me4base` |
| 4 | 6 | 'NS'  | `'Central'` | `3D_solver/NSTGV/build_o6_ns_central` |
| 5 | 6 | 'LES' | `'ME4Base'` | `3D_solver/NSTGV/build_o6_les_me4base` |
| 6 | 6 | 'LES' | `'Central'` | `3D_solver/NSTGV/build_o6_les_central` |

## Implementation

1. **Snapshot current `3D_solver/NSTGV/config.fypp`** (currently: `VISC='NS'`, `SCHEME='KEEP'`, `ORDER=2`,
   plus other unrelated unstaged edits) so it can be restored exactly after the sweep — these are the
   user's own in-progress edits and must not be lost.

2. **For each of the 6 combinations, in sequence:**
   - Edit `3D_solver/NSTGV/config.fypp`: set `ORDER` and `VISC` to the row's values; add a
     `#:set VISC_STENCIL = '...'` line (placed after the `VISC`/`ORDER` lines) only for the `ORDER=6` rows,
     omit it entirely for `ORDER=4` rows.
   - Run `cmake -B <build_dir>` then `cmake --build <build_dir> -j` from `3D_solver/NSTGV/` (compiler
     confirmed available: `mpif90` wrapping `nvfortran` 24.7 at
     `/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/bin/nvfortran`, `cmake` at `/usr/bin/cmake`).
   - Record configure/build success or failure, capturing the exact compiler error text on failure.
   - For successful builds, spot-check the generated `<build_dir>/calc_visc_high.f90` picked the intended
     branch (e.g. `grep -n "load_smem_visc_me4_base\|load_smem_visc_cent6" <build_dir>/calc_visc_high.f90`)
     so a "compiles" result isn't masking a wrong-branch dispatch.

3. **Restore `3D_solver/NSTGV/config.fypp`** to its exact pre-sweep content from step 1.

4. **Report a results table**: combination → configure result → build result → dispatch-branch
   verification, plus any compiler errors seen, back to the user. Leave the 6 `build_o*` directories in
   place per the user's request (not cleaned up) so generated sources/binaries remain available for
   inspection.

## Verification

- Each `cmake --build <dir> -j` either exits 0 (success, binary `<dir>/a.out` produced) or exits non-zero
  (failure) — this is the primary pass/fail signal for the sweep.
- The `grep` spot-check in step 2 confirms the fypp preprocessor actually selected the requested
  `VISC_STENCIL`/`ORDER` branch for each build, not just that *some* code compiled.
- After restoring `config.fypp`, run `git diff 3D_solver/NSTGV/config.fypp` to confirm it matches the
  pre-sweep unstaged state exactly (no leftover sweep edits).
