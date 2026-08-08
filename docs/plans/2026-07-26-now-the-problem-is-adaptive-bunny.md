# Merge u/v/w loading into each stage, accumulate LES Hsgs incrementally

## Context

Building on the prior rounds of restructuring `calc_visc_high.f90.fypp` /
`calc_visc_high_internal.f90.fypp` / `load_smem_visc_cent.f90.fypp` (staged
`grad_a`/`grad_b` shared buffers, `mu`/`T` read via scalar
`interp/diff${VISC_ORDER}$_s` calls with no register array): `u, v, w` are
still 3 **persistent** shared arrays, loaded once via a dedicated
`load_smem_visc_cent${VISC_ORDER}$_*_uvw` call (1 more `syncthreads()`)
before stage 1 begins, and held resident for the entire kernel body.

Inspecting actual usage shows this is unnecessary: **each of the 3
stress-component stages only ever touches ONE of u/v/w** — both in the
high-order interior branch (`calc_tau_straight`/`calc_tau_cross`'s "self"
velocity argument: `u` in stage 1, `v` in stage 2, `w` in stage 3 for
`calc_Ev`, and the direction-appropriate analogue for `calc_Fv`/`calc_Gv`)
**and** in the 2nd-order boundary fallback (stage 1's `else` branch only
reads `u(idx)`/`u(idx+1)`, stage 2's only `v(idx)`/`v(idx+1)`, stage 3's
only `w(idx)`/`w(idx+1)`). So instead of 3 persistent arrays, a single
reusable shared buffer — loaded fresh in each stage's own load call,
right alongside that stage's `grad_a`/`grad_b` gradient(s), in the same
loop — is all that's needed. This:
- drops the shared-array count from 5 (`u,v,w,grad_a,grad_b`) to 3
  (`uvw,grad_a,grad_b`) — a further ~40% cut on top of the prior rounds'
  reductions,
- removes the standalone `_uvw` load call and its `syncthreads()` entirely
  (folded into stage 1's load), cutting total `syncthreads()` calls per
  kernel from 6 to 5.

The one place all 3 of u,v,w were needed **together** is the (currently
dead-code) LES `Hsgs` term, which needs `u²+v²+w²` at each stencil point.
Per your instruction, this becomes an **incremental accumulation**: a
persistent `H` array (subroutine-scope, sized like the interior stencil, or
`H(2)` for the boundary fallback) is initialized with `Cp*T + qc2 + 0.5*u²`
in stage 1 (when `uvw` holds `u`), has `0.5*v²` added in stage 2 (when
`uvw` holds `v`), has `0.5*w²` added in stage 3 (when `uvw` holds `w`), and
`Hsgs = -muti * diff${VISC_ORDER}$_s(d, H(...)) / Prt` is computed at the
very end of stage 3 (once all 3 components have contributed) instead of
all at once in stage 1 as today.

## Approach

### 1. `3D_solver/src/load_smem_visc_cent.f90.fypp`

For each direction (`_x`, `_y`, `_z`, `_z_koff`):
- Delete the standalone `..._uvw` subroutine.
- Merge the direct (no-stencil) copy of whichever single velocity
  component that direction's stage 1 needs into `..._stage1` — same
  grid-stride loop, guarded by the existing wide range check (`jk_in_range`
  etc., not the narrower `_grad_y`/`_grad_z` guards), writing into a new
  `uvw` output array alongside `grad_a`/`grad_b`. `_stage1` gains the
  needed extra `Q_*` input (whichever of `Q_2/Q_3/Q_4` isn't already an
  argument for that direction's gradients).
- `..._stage2` and `..._stage3` each gain the same treatment: merge in the
  direct copy of *their* stage's velocity component (a different `Q_*`
  each time — e.g. for `calc_Ev`: stage1→`u` via `Q_2`, stage2→`v` via
  `Q_3`, stage3→`w` via `Q_4`), writing into the same reused `uvw` output
  array (matching how `grad_a` is already reloaded/reused stage-to-stage).

The direction→stage→component mapping (already established in the prior
rounds' consumption table) is:

| Direction | Stage 1 | Stage 2 | Stage 3 |
|---|---|---|---|
| Ev (x) | `u` (Q_2) | `v` (Q_3) | `w` (Q_4) |
| Fv (y) | `v` (Q_3) | `u` (Q_2) | `w` (Q_4) |
| Gv (z) | `w` (Q_4) | `u` (Q_2) | `v` (Q_3) |

### 2. `calc_visc_high.f90.fypp` / `calc_visc_high_internal.f90.fypp`

For each of `calc_Ev/Fv/Gv${VISC_ORDER}$[_koff]` (both files):
- Replace the `u, v, w` shared-array declarations with a single `uvw`
  shared array (same shape as today's `u`).
- Delete the `call load_smem_..._uvw(...)` call; update the `..._stage1/2/3`
  calls to match their new signatures (extra `Q_*` in, `uvw` out).
- Replace every reference to `u(...)`, `v(...)`, or `w(...)` in that
  stage's interior/boundary blocks with `uvw(...)` — each stage
  unambiguously means one physical component, per the table above (e.g. in
  `calc_Ev`'s stage 1, `uvw` means `u`; in stage 2, `uvw` means `v`).
- LES `Hsgs` (dead code, `#:if is_les`): change `H`/`H(2)` from a
  stage-1-local block variable to a subroutine-scope persistent array.
  Initialize with `Cp*T(...) + qc2(...) + 0.5*uvw(...)**2` in stage 1, add
  `0.5*uvw(...)**2` in stage 2 and stage 3 (by then `uvw` holds the next
  component), and move the `Hsgs = -muti * diff${VISC_ORDER}$_s(...) / Prt`
  computation to the end of stage 3, right before it's added into the
  final combine block (which already conditionally includes `Hsgs`).

`calc_visc_high_internal.f90.fypp` mirrors the same change, minus the
boundary-fallback branch (interior-only, as established in prior rounds).

## Verification

1. Rebuild `NSTGV`, `DHIT`, `SBLI` — confirm clean compilation.
2. Re-run the `cuobjdump --dump-resource-usage` register/shared-memory
   comparison (final LTO-linked binary) — expect shared memory to drop
   further (3 arrays vs. 5) and report the resulting register counts
   (registers may shift since `uvw` is now reloaded 3× instead of loaded
   once — report whatever the compiler actually produces, don't assume).
3. Re-run the short `mpirun -n 2` smoke test on `NSTGV` (same method as
   the prior round) to confirm merging the loads and removing a
   `syncthreads()` hasn't reintroduced a race or hang — check for finite,
   sane VTK output within a short window.
4. Spot-check the generated `build/calc_visc_high.f90` for one case: no
   `_uvw` call remains, each stage's load call populates `uvw` correctly
   for its direction (per the table above), and boundary-fallback
   references to `u`/`v`/`w` were correctly renamed per-stage.

