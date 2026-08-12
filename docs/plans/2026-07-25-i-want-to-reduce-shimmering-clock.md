# Simplify viscous loaders (drop pipelining, merge copies, hoist guards) + hoist boundary-guard ifs in convective kernels

## Context

The previous register-blocking pass (already completed and verified) deliberately left the shared-memory loader files (`load_smem_visc_cent.f90.fypp`, `load_smem_visc_me4_base.f90`) untouched, flagging that their `pipelineMemcpyAsync`-based staging was a tradeoff worth revisiting later. The user has now asked to revisit exactly that, plus extend the same "hoist loop-invariant work out of the loop" idea to the convective (KEEP/SLAU/Hybrid) kernels' shared-memory load loops. Concretely, five asks:

1. Combine the loaders' separate copy-to-shared-memory loop passes into fewer passes.
2. Hoist the complex boundary-guard `if` conditions inside loader loops out of the loop, cached in a `logical`.
3. Drop `pipelineMemcpyAsync`/`pipelineCommit`/`pipelineWaitPrior` — the compute isn't heavy enough to justify the async-copy staging (this echoes the user's very first piece of guidance this session).
4. The same complex-`if`-inside-a-loop pattern also exists in the convective kernels' inline shared-memory load loops — hoist those into a `logical` too.
5. In the viscous term, `dx`/`dy`/`dz` reciprocal-spacing gmem is read repeatedly inside the loader subroutines. Per the user's explicit choice, cache these in a register **at the call site**, before calling the loader, and change the loader's signature to accept a scalar instead of the whole array.

Read-only investigation (via 2 Explore agents + direct reads) confirmed the exact shape of all five items — details below.

## What's actually there today

### Loaders: `load_smem_visc_cent.f90.fypp` (→ `load_smem_visc_cent4`/`load_smem_visc_cent6`), `load_smem_visc_me4_base.f90`

Each defines `_x`/`_y`/`_z`/`_z_koff`. Each subroutine currently does **3 separate `do ii=...` loop passes**, one per velocity component (u, then v, then w), each wrapped in its own `pipelineMemcpyAsync` + `pipelineCommit()`, with the complementary cross-direction gradient (e.g. `vy`/`wz` while `u` is in flight, `uy` while `v` is in flight, `uz` while `w` is in flight) computed in a second loop right after each commit — a deliberate compute/copy-overlap pattern from the earlier pipelining task. Every one of these loops repeats the same boundary-guard condition, e.g. in `_x` (`load_smem_visc_cent.f90.fypp:55,63,70,82,90,102,110`):
```fortran
if (1 <= i .and. i <= nx .and. j <= ny .and. k <= nz) then          ! copy guard — j,k fixed for this call
if (1 <= i .and. i <= nx .and. D_MIN <= j .and. j <= ny-D_SFX .and. k <= nz) then   ! vy/uy gradient guard
if (1 <= i .and. i <= nx .and. j <= ny .and. D_MIN <= k .and. k <= nz-D_SFX) then   ! wz/uz gradient guard
```
`j`/`k` (or `i`/`k`, `i`/`j` for `_y`/`_z`) are fixed scalar arguments for the whole subroutine call — only the `1<=i<=nx`-style clause actually varies with the loop variable `ii`. The `j<=ny .and. k<=nz` part is recomputed identically on every one of the ~7 conditional checks in a subroutine, and the gradient stencils (e.g. line 65) read `inv_dy(j)`/`inv_dz(k)` from gmem repeatedly across the separate loops (once in the `vy` loop, again in the later `uy` loop, for the exact same `j`).

`load_smem_visc_me4_base.f90` has the identical structure (same 3-stage pipeline, same repeated guards, literal bounds instead of fypp `D_MIN`/`D_SFX`).

`load_smem_visc2.f90` (used only when `VISC_ORDER==2`, i.e. never called from the files being touched by the register-blocking work) already copies u/v/w in one loop pass but still uses a single-stage `pipelineMemcpyAsync`×3 + one `pipelineCommit`/`pipelineWaitPrior`, and repeats a partially loop-invariant guard (`j<=ny .and. k<=nz` etc.) plus a zero-fill `else` branch. It has no gradient computation and does **not** take `inv_dx`/`inv_dy`/`inv_dz` at all (confirmed by reading the file) — so item 5 does not apply here, only items 2/3 (hoist guard, drop pipeline).

### Call sites: `calc_visc_high.f90.fypp`, `calc_visc_high_internal.f90.fypp`

Each of `calc_Ev`/`calc_Fv`/`calc_Gv` (× `_koff`) calls its loader once, before the branch dispatch on `(i,j,k)`, e.g. (`calc_visc_high.f90.fypp:96-108`):
```fortran
j  = (blockIdx%y-1)*blockDim%y + jt + 1
k  = (blockIdx%z-1)*blockDim%z + kt + 1
offset_yz = ...
call load_smem_visc_cent4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q_2, Q_3, Q_4, u, v, w, uy, vy, uz, wz)
i  = (blockIdx%x-1)*blockDim%x + it
```
`j`/`k` are known before this call — so the "other two directions'" spacing can be read into a register right here and passed as a scalar. Each kernel only ever needs the two *cross* directions for its loader call (`calc_Ev` → `dy`,`dz`; `calc_Fv` → `dx`,`dz`; `calc_Gv` → `dx`,`dy`) and its *own* direction later inside the branch bodies (already register-cached as `dxi`/`dyj`/`dzk` from the prior task) — these never overlap, so there's no naming collision between the new subroutine-scope scalars and the existing block-scope ones.

### Convective kernels: `calc_keep_kernel(.f90.fypp/_internal)`, `calc_slau_kernel(.f90.fypp/_internal)`, `calc_hybrid_kernel(.f90.fypp/_internal)`

No separate loader — each `attributes(global)` subroutine loads `Q_1..Q_5`/`T` inline via one `do ii=...` loop before `syncthreads()`, e.g. (`calc_keep_kernel.f90.fypp:51-59`):
```fortran
do ii = it-io, threadsE%x+io+1, blockDim%x
  i = i_base + ii
  if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
    idx = ii + offset_yz
    rho(idx) = Q_1(i,j,k);   u(idx) = Q_2(i,j,k)
      v(idx) = Q_3(i,j,k);   w(idx) = Q_4(i,j,k)
      p(idx) = Q_5(i,j,k); tmp(idx) =   T(i,j,k)
  endif
enddo
```
Same pattern, same fix (hoist the loop-invariant `j`/`k` clauses), appears **27 times**: `calc_keep_kernel.f90.fypp` (x/y/z, 3 sites), `calc_keep_kernel_internal.f90.fypp` (x/y/z/z_koff, 4), `calc_slau_kernel.f90.fypp` (3), `calc_slau_kernel_internal.f90.fypp` (4), `calc_hybrid_kernel.f90.fypp` (3), `calc_hybrid_kernel_internal.f90.fypp` (4), `calc_roe_kernel.f90.fypp` (3), `calc_roe_kernel_internal.f90.fypp` (4). Per `CLAUDE.md` ("Roe scheme is not used... KEEP, SLAU, Hybrid schemes should be optimized"), **Roe is excluded from this pass** — 21 sites remain in scope. The codebase already uses the `logical` idiom for exactly this kind of thing (`calc_hybrid_kernel.f90.fypp:56` declares `logical :: compute_slau`), so the hoisted guard should follow that same style (e.g. `logical :: jk_in_range`).

## Changes

### 1. Loaders — drop pipelining, merge copy passes, hoist guards, accept scalar spacing args

For each of `load_smem_visc_cent.f90.fypp` (`_x`/`_y`/`_z`/`_z_koff`, both `VISC_ORDER` branches share one template) and `load_smem_visc_me4_base.f90` (`_x`/`_y`/`_z`/`_z_koff`):

- **Remove `use wmma`** and all `pipelineMemcpyAsync`/`pipelineCommit`/`pipelineWaitPrior` calls. Replace each `call pipelineMemcpyAsync(dest, src)` with a plain assignment `dest = src`. Keep the final `call syncthreads()` — still required so every thread sees the fully-populated shared tile before use.
- **Merge the 3 separate copy+gradient loop-pairs into one single loop pass** per subroutine (the async-overlap rationale for staggering them no longer applies once the pipeline is gone) — one `do ii=...` loop that assigns `u`,`v`,`w` and computes both cross-direction gradient pairs (e.g. `uy`,`vy`,`uz`,`wz` for `_x`), each still individually guarded since their valid ranges differ (full-domain for the copy, stencil-margin-restricted for the gradients).
- **Hoist the loop-invariant parts of each guard into `logical`s** computed once before the loop, e.g. for `_x`:
  ```fortran
  logical :: jk_in_range, jk_in_range_grad_y, jk_in_range_grad_z
  jk_in_range         = (j <= ny .and. k <= nz)
  jk_in_range_grad_y  = (D_MIN <= j .and. j <= ny-D_SFX .and. k <= nz)
  jk_in_range_grad_z  = (j <= ny .and. D_MIN <= k .and. k <= nz-D_SFX)
  ```
  then inside the merged loop: `if (1 <= i .and. i <= nx .and. jk_in_range) then` (copy), `if (1 <= i .and. i <= nx .and. jk_in_range_grad_y) then` (vy/uy), `if (1 <= i .and. i <= nx .and. jk_in_range_grad_z) then` (uz/wz). Symmetric for `_y` (invariant over `i`,`k`; loop var `j`) and `_z`/`_z_koff` (invariant over `i`,`j`; loop var `k`).
- **Change the spacing-array dummy args to scalars**: e.g. `_x`'s `real(8), intent(in), device, contiguous :: inv_dy(ny-1)` / `inv_dz(nz-1)` become `real(8), intent(in), value :: inv_dy_j`, `inv_dz_k`; use these directly in the stencil formulas instead of `inv_dy(j)`/`inv_dz(k)`. Symmetric for `_y` (→ `inv_dx_i`, `inv_dz_k`) and `_z`/`_z_koff` (→ `inv_dx_i`, `inv_dy_j`).

Concrete sketch of the fused `_x` subroutine body (replaces the current 3 separate copy-loop/commit/gradient-loop stages with one pass):
```fortran
logical :: jk_in_range, jk_in_range_grad_y, jk_in_range_grad_z
jk_in_range        = (j <= ny .and. k <= nz)
jk_in_range_grad_y = (D_MIN <= j .and. j <= ny-D_SFX .and. k <= nz)
jk_in_range_grad_z = (j <= ny .and. D_MIN <= k .and. k <= nz-D_SFX)
do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
  i = i_base + ii
  idx = ii + offset_yz
  if (1 <= i .and. i <= nx .and. jk_in_range) then
    u(idx) = Q_2(i,j,k)
    v(idx) = Q_3(i,j,k)
    w(idx) = Q_4(i,j,k)
  endif
  if (1 <= i .and. i <= nx .and. jk_in_range_grad_y) then
    vy(idx) = (two_third * (-Q_3(i,j-1,k) + Q_3(i,j+1,k)) - one_twelfth * (-Q_3(i,j-2,k) + Q_3(i,j+2,k))) * inv_dy_j
    uy(idx) = (two_third * (-Q_2(i,j-1,k) + Q_2(i,j+1,k)) - one_twelfth * (-Q_2(i,j-2,k) + Q_2(i,j+2,k))) * inv_dy_j
  endif
  if (1 <= i .and. i <= nx .and. jk_in_range_grad_z) then
    wz(idx) = (two_third * (-Q_4(i,j,k-1) + Q_4(i,j,k+1)) - one_twelfth * (-Q_4(i,j,k-2) + Q_4(i,j,k+2))) * inv_dz_k
    uz(idx) = (two_third * (-Q_2(i,j,k-1) + Q_2(i,j,k+1)) - one_twelfth * (-Q_2(i,j,k-2) + Q_2(i,j,k+2))) * inv_dz_k
  endif
enddo
call syncthreads()
```
One loop, one pass over `ii`, all copies and both gradient pairs computed together — no more staggered copy/commit/compute stages, since there's no async latency left to hide. Symmetric restructuring applies to `_y`, `_z`, `_z_koff`, and to `load_smem_visc_me4_base.f90`'s four subroutines (same fusion, ME4-specific stencil constants).

For `load_smem_visc2.f90` (`_x`/`_y`/`_z`/`_z_koff`): drop `use wmma` and the pipeline calls (plain assignment instead), hoist the loop-invariant guard clause into one `logical` before the loop. No spacing-array args exist here, so item 5 doesn't apply to this file. The existing zero-fill `else` branch stays as-is.

### 2. Call sites — read the loader's spacing scalars into registers before calling it

In `calc_visc_high.f90.fypp` and `calc_visc_high_internal.f90.fypp`, for each of `calc_Ev`/`calc_Fv`/`calc_Gv` (× `_koff`), immediately before the existing `call load_smem_visc_cent${VISC_ORDER}$_x(...)` (or `_y`/`_z`/`me4_base` equivalent — whichever the fypp branch selects), add the two register reads for that call's cross directions and pass them instead of the arrays:
```fortran
! calc_Ev: needs dy, dz for the loader (its own direction, dx, is cached later as dxi inside each branch)
dyj = dy(j)
dzk = dz(k)
call load_smem_visc_cent${VISC_ORDER}$_x(it, jt, kt, j, k, nx, ny, nz, dyj, dzk, Q_2, Q_3, Q_4, u, v, w, uy, vy, uz, wz)
```
(`calc_Fv` → `dxi`,`dzk`; `calc_Gv` → `dxi`,`dyj`). Declare these two new scalars at the subroutine-body level (alongside the existing `integer i, j, k, ...` / `real(8) :: kTx, utau_sum` declarations) since they're used before the per-branch `block`s begin. No collision with the existing block-scoped `dxi`/`dyj`/`dzk` from the prior register-blocking task, since a kernel's own-direction register and its loader's cross-direction registers are always disjoint by construction.

### 3. Convective kernels — hoist the repeated boundary-guard `if`

For `calc_keep_kernel.f90.fypp`, `calc_keep_kernel_internal.f90.fypp`, `calc_slau_kernel.f90.fypp`, `calc_slau_kernel_internal.f90.fypp`, `calc_hybrid_kernel.f90.fypp`, `calc_hybrid_kernel_internal.f90.fypp` (Roe excluded per `CLAUDE.md`): in each `_x`/`_y`/`_z`/`_z_koff` subroutine's load loop, add one `logical` declared and computed once right after `j`/`k` (or `i`/`k`, `i`/`j`) are set, e.g. for an x-sweep:
```fortran
logical :: jk_in_range
jk_in_range = (j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz)
...
do ii = it-io, threadsE%x+io+1, blockDim%x
  i = i_base + ii
  if (i >= 1 .and. i <= nx .and. jk_in_range) then
    ...
```
Symmetric for y-sweeps (hoist `i`,`k` → `ik_in_range`) and z-sweeps/`_z_koff` (hoist `i`,`j` → `ij_in_range`). Purely mechanical, same shape at all 21 sites; no formula/behavior change.

## Verification

1. fypp balance (`#:if`/`#:endif`, `block`/`end block`) on every touched `.fypp` file.
2. fypp-preprocess-check the loader + `calc_visc_high(_internal)` files against: NSTGV's real config (non-LES, ORDER=6), a synthetic LES config (ORDER=6 and ORDER=4, as used in the prior task) — confirms the `is_les` branches and both `VISC_ORDER` expansions still compile.
3. Rebuild NSTGV and SBLI (exercise the Cent6 loader + SLAU convective hoist) in fresh build dirs; confirm clean compilation.
4. Rebuild one KEEP-scheme case (check `config.fypp` for ETGV/DHIT/IVST/KHI/STZ to find one) to exercise the KEEP hoist, and check whether any real case uses `SCHEME='Hybrid'` (TBL per docs, but confirm its actual `config.fypp` — if it turns out ORDER=2, Hybrid's shared-memory-load hoist is still exercised since the guard-hoist is convective-only and independent of `VISC_ORDER`).
5. `cuobjdump --dump-resource-usage` register/stack/spill comparison for the touched kernels against the last verified baselines (`build_pipeline_fixed` for NSTGV/SBLI) — expect zero new spills; register counts may shift.
6. Full simulation run + `cmp`/`md5sum` on `recal/Q00001.dat` for NSTGV (and the KEEP case) against their last verified baselines — must be bit-identical (pure restructuring, no formula changes).
