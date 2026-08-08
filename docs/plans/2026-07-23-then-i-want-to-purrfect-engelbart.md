# Add __ldcs / __ldlu cache-hint intrinsics to viscous kernels

## Context

The working tree already has WIP (uncommitted) changes adding `__ldcg`/`__stcg`
cache-hint intrinsics to the E/F/G flux read-modify-write in the viscous
kernels (`calc_visc2.f90.fypp` ×2, `calc_visc_high.f90.fypp`,
`calc_visc_high_internal.f90.fypp`). The user wants to extend this with
`__ldcs`, `__ldlu`, and `__stcs` wherever they are a genuinely better fit,
following the same reasoning style already used (and explained in comments)
in the *committed* `calc_steps.f90`, which uses `__stcs` for RK-stage output
writes that are "write-once, not re-read until the next kernel launch."

Investigation findings (see conversation for full detail):

- **`__ldlu` fits well.** In all four target files, the pattern
  `e2v = __ldcg(E(2,i,j-1,k-1)); call __stcg(E(2,i,j-1,k-1), e2v - txx)`
  reads the existing flux value only to immediately overwrite it. Each
  thread owns a unique flux index — no other thread in the kernel launch
  reads that pre-write value — so the load is a true "last use." Swap the
  **read** from `__ldcg` to `__ldlu`. Leave the **store** as `__stcg`
  unchanged, since `calc_R` in `calc_steps.f90` re-reads that exact address
  in the very next kernel launch (same reasoning `calc_steps.f90` already
  uses to justify `__stcg` there).
- **`__ldcs` fits well** on the 2nd-order-style viscous stencil's "corner"
  neighbor loads — each used exactly once per thread — as opposed to the
  central pair (`mu_c`/`mu_e`, `mu_c`/`mu_n`, `mu_c`/`mu_t`), which is reused
  2-4× per thread and should stay default. This pattern appears in full in
  `calc_visc2.f90.fypp` (2D and 3D — order 2 has no separate interior branch,
  it's always this style) and in the boundary-fallback `else` blocks of
  `calc_visc_high.f90.fypp` (used when the high-order stencil can't reach,
  i.e. cells near domain edges).
- **No genuine `__stcs` site exists in these 5 files** — confirmed with the
  user, who agreed to leave `__stcs` out of scope here. Every global store in
  these files is to `E`/`F`/`G`, which are always re-read almost immediately
  by `calc_R`, the opposite of the write-once/not-soon-reread case `__stcs`
  is for. `calc_steps.f90` already uses `__stcs` correctly and needs no
  change.
- `calc_visc_high_internal.f90.fypp` has **no** boundary-fallback `else`
  block (it's interior-only — the high-order stencil always applies), so
  only the `__ldlu` change applies there. The dense high-order stencil reads
  (`mu4(:) = mu(i-1:i+2,j,k)`, `mu6(:)`, `T4/T6`, ME4Base scalars, etc., in
  `calc_visc_high.f90.fypp`/`_internal`) are heavily reused across adjacent
  x/y/z-threads via overlapping stencil windows — leave those on default
  loads; marking them `__ldcs`/`__ldlu` would hurt, not help, cache reuse.

## Changes

### 1. `__ldlu` for the flux read-before-overwrite (all 4 files)

In every `block ... e2v = __ldcg(E(...)) ... end block` (and F/G analogues),
change the read intrinsic from `__ldcg` to `__ldlu`. Add a one-line comment
in the same style as `calc_steps.f90`, e.g.:
```fortran
! no other thread reads this address before the paired store below overwrites it: __ldlu
e2v = __ldlu(E(2,i,j-1,k-1))
call __stcg(E(2,i,j-1,k-1), e2v - txx)
```
Apply to all occurrences:
- `3D_solver/src/calc_visc2.f90.fypp` — `calc_Ev2`/`calc_Fv2`/`calc_Gv2` (both `suf` variants share this code)
- `2D_solver/src/calc_visc2.f90.fypp` — `calc_Ev2`/`calc_Fv2`
- `3D_solver/src/calc_visc_high.f90.fypp` — end of `calc_Ev${VISC_ORDER}$`/`calc_Fv...`/`calc_Gv...` (one site each, outside the order-dependent branches)
- `3D_solver/src/calc_visc_high_internal.f90.fypp` — inside every order/stencil branch of `calc_Ev.../calc_Fv.../calc_Gv...` (multiple sites per subroutine, since here the E/F/G write is duplicated inside each `if` branch rather than factored out after `endif`)

### 2. `__ldcs` for single-use neighbor loads

Change the "corner" neighbor loads (each read exactly once per thread) to
`__ldcs`; leave the central pair and any load reused ≥2× untouched. Add a
short comment analogous to `calc_visc_high.f90.fypp`'s existing
`!mu_c..mu_nt = mu at ...` style comment, noting these are single-use.

- **`3D_solver/src/calc_visc2.f90.fypp`**:
  - `calc_Ev2`: `mu_cjm, mu_cjp, mu_ejm, mu_ejp, mu_ckm, mu_ckp, mu_ekm, mu_ekp`; `mut_cjm, mut_cjp, mut_ejm, mut_ejp, mut_ckm, mut_ckp, mut_ekm, mut_ekp` (under `is_les`); `q2jm, q2ejm, q2jp, q2ejp, q3jm, q3ejm, q3jp, q3ejp` and `q2km, q2ekm, q2kp, q2ekp, q4km, q4ekm, q4kp, q4ekp`. Leave `mu_c, mu_e, T_c, T_e` default.
  - `calc_Fv2`: same shape with `mu_wc, mu_wn, mu_ec, mu_en, mu_ckm, mu_ckp, mu_nkm, mu_nkp` (+ `mut_*` analogues) and the `Q` neighbor reads in the `mx1/mx2`/`mz1/mz2` blocks. Leave `mu_c, mu_n, T_c, T_n` default.
  - `calc_Gv2`: same shape with `mu_wc, mu_wt, mu_ec, mu_et, mu_sc, mu_st, mu_nc, mu_nt` (+ `mut_*`) and its `Q` neighbor reads. Leave `mu_c, mu_t, T_c, T_t` default.

- **`2D_solver/src/calc_visc2.f90.fypp`**:
  - `calc_Ev2`: `mu_i_jm1, mu_ip1_jm1, mu_i_jp1, mu_ip1_jp1` and the four `Q(i[,±1],2/3,j±1)` reads. Leave `mui, muip1` default (each reused: once in `mudx`, once in `my1`/`my2`).
  - `calc_Fv2`: `mu_im1_j, mu_im1_jp1, mu_ip1_jp1, mu_i_jp1` (note: reuse `mu_i_jp1` carefully — check whether it's used once or twice before marking) and the four `Q` neighbor reads. Leave `muij, muip1j` default.

- **`3D_solver/src/calc_visc_high.f90.fypp`** — the three `else` boundary-fallback blocks only (not the order-4/6/ME4Base interior branches):
  - `calc_Ev...` fallback (around the `mu_cjm...mu_ekp` block): same set as `calc_visc2.f90.fypp`'s `calc_Ev2`.
  - `calc_Fv...` fallback: same set as `calc_Fv2`.
  - `calc_Gv...` fallback: same set as `calc_Gv2`, i.e. `mu_wc, mu_wt, mu_ec, mu_et, mu_sc, mu_st, mu_nc, mu_nt` (+ `mut_*`) and its 8 `Q` neighbor reads (`Q(i-1,2,j,k)`, `Q(i,2,j,k)`, etc.). Leave `mu_c, mu_t, T_c, T_t` default.

- **`3D_solver/src/calc_visc_high_internal.f90.fypp`**: no change (no fallback block exists; only the `__ldlu` change from step 1 applies here).

## Verification

Since this only changes cache-hint intrinsics (no change to arithmetic or
control flow), correctness should be unaffected — verify by:
1. For each touched case's build dir, regenerate and diff the fypp output:
   `cmake --build build -j` for at least one case per file touched
   (e.g. `3D_solver/NSTGV` or `3D_solver/SBLI` for NS; `3D_solver/TBL` for LES,
   to exercise the `is_les` branches; `2D_solver/BL` for the 2D file), and
   spot check `build/calc_visc2.f90` / `build/calc_visc_high*.f90` to confirm
   fypp expanded the new intrinsics correctly for both `is_koff` variants and
   both `is_les` branches.
2. Run the existing Python test suite (`pytest ouxsbli/tests/`) for cases
   that exercise NS/LES viscous kernels, to confirm no numerical regression
   (the intrinsics only affect caching behavior, not results, but this is a
   free correctness check since fypp bugs are easy to introduce when editing
   templates by hand across multiple `#:if`/`#:for` branches).
3. Optionally profile with `nsys`/`ncu` (`bash ../profile.sh` in a case's
   `build/` dir) before/after to confirm the change doesn't regress memory
   throughput — cache-hint tuning is inherently something to measure, not
   just reason about.
