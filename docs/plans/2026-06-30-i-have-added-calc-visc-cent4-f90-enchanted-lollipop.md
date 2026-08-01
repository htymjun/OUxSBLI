# Plan: Unify cent files + rename calc_visc4 → calc_visc_high with 3-way dispatch

## Context

The user added four new plain-Fortran files for central-difference viscous discretization:
- `calc_visc_cent4.f90` / `calc_visc_cent6.f90` — device helper routines (`interp4/6`, `diff4/6`, `calc_tau_*`), plain Fortran include files used inside the kernel module's `contains` section
- `load_smem_visc_cent4.f90` / `load_smem_visc_cent6.f90` — compiled shared-memory loader modules

The cent4 and cent6 variants are structurally identical, differing only in stencil order (4-point vs 6-point) and io_v (cent4: io_v=1, cent6: io_v=2).

**Resolved dispatch design:**
- `VISC_ORDER = ORDER` always (change default in `calc_flux_base.fypp`)
- **ORDER = 4 → always Cent4** (ME4-Base is no longer used at order 4)
- **ORDER = 6 → ME4-Base (default) or Cent6** via new config variable `VISC_STENCIL`
- ME4-Base is preserved as an order-6 option (same stencil coefficients, subroutine now named `calc_Ev6` instead of `calc_Ev4`)

---

## Key Structural Facts (from code)

- `calc_visc4.f90.fypp` line 16: `include 'calc_visc_me4_base.f90'` (Fortran include inside `contains`)
- `calc_visc4.f90.fypp` and `calc_visc4_internal.f90.fypp` both: `use load_smem_visc_me4_base`
- `calc_flux_base.fypp` calls kernels as `calc_${dir}$v${VISC_ORDER}$` — subroutine names stay consistent with VISC_ORDER; only the `use` module name changes
- `CMakeLists.txt` line 76: still references old `"load_smem_visc4.f90"` (bug fix needed → `load_smem_visc_me4_base.f90`)
- ME4-Base interior: `mu3(3)`, `u(6)` slices, calls `flux4`; condition `3 ≤ i ≤ nx-3`
- Cent4 interior: `mu(4)`, `u(4)` slices, calls `calc_tau_straight(4pt)`; condition `2 ≤ i ≤ nx-2` (wider, io_v=1)
- Cent6 interior: `mu(6)`, `u(6)` slices, calls `calc_tau_straight(6pt)`; condition `3 ≤ i ≤ nx-3`
- Fortran include search path: `target_include_directories` adds `3D_solver/src/` and `build/`; fypp `-I` flag currently only adds `<case-dir>`
- `one_24 = 1/24` is a module-level parameter used by both ME4-Base (`flux4`) and Cent4 (`diff4`)

---

## New Config Variable

`VISC_STENCIL` — only relevant when `VISC_ORDER = 6`. Add to `config.fypp` for cases wanting Cent6:

```python
#:set VISC_STENCIL = 'Central'  # enables Cent6 at ORDER=6
```

Default (when unset): `'ME4Base'`. Each `calc_visc_high*.fypp` sets: `#:if not defined('VISC_STENCIL')` → `'ME4Base'`.

---

## Step-by-Step Changes

### 1. Create `3D_solver/src/calc_visc_cent.f90.fypp`

Merge `calc_visc_cent4.f90` and `calc_visc_cent6.f90`. Used only via fypp `#:include` (not a standalone compilation unit, not listed in CMakeLists).

Parameterize with `VISC_ORDER`:

```fortran
#:set N = VISC_ORDER
  pure attributes(device) function interp${N}$(a) result(ai)
    real(8), intent(in), contiguous :: a(${N}$)
    #:if N == 4
    ai = 0.0625d0 * (-a(1) + 9.d0 * (a(2) + a(3)) - a(4))
    #:else  ! N == 6
    ai = 0.00390625d0 * (3.d0*a(1) - 25.d0*a(2) + 150.d0*(a(3)+a(4)) - 25.d0*a(5) + 3.d0*a(6))
    #:endif
  end function interp${N}$

  pure attributes(device) function diff${N}$(d, a) result(da)
    real(8), intent(in)             :: d
    real(8), intent(in), contiguous :: a(${N}$)
    real(8) da
    #:if N == 4
    da = (a(1) - 27.d0 * (a(2) - a(3)) - a(4)) * one_24 * d
    #:else  ! 6th-order, local coef constants
    da = (-0.0046875d0*a(1) + 0.065104167d0*a(2) - 1.171875d0*(a(3)-a(4)) - 0.065104167d0*a(5) + 0.0046875d0*a(6)) * d
    #:endif
  end function diff${N}$

  ! calc_tau_straight, calc_tau_straight_LES, calc_tau_cross, calc_tau_cross_LES
  ! — same structure as existing cent4/cent6.f90, calling interp${N}$ and diff${N}$
  ! — array sizes all (${N}$)
```

**Delete** `calc_visc_cent4.f90` and `calc_visc_cent6.f90` after this.

### 2. Create `3D_solver/src/load_smem_visc_cent.f90.fypp`

Merge `load_smem_visc_cent4.f90` and `load_smem_visc_cent6.f90`. Listed in `_SHARED_FYPP` in CMakeLists.

```fortran
module load_smem_visc_cent${VISC_ORDER}$
  ...
  ! io_v: 1 for VISC_ORDER==4, 2 for VISC_ORDER==6
  ! gradient stencil: 4-point (cent4 coefficients) or 6-point (cent6 coefficients)
  ! subroutine names: load_smem_visc_cent${VISC_ORDER}$_x/y/z and _z_koff
```

**Delete** `load_smem_visc_cent4.f90` and `load_smem_visc_cent6.f90` after this.

### 3. Create `3D_solver/src/calc_visc_high.f90.fypp` (from calc_visc4.f90.fypp)

Module renamed `calc_visc4` → `calc_visc_high`. Public subroutine names `calc_Ev${VISC_ORDER}$` etc. remain consistent with `calc_flux_base`'s call pattern.

**Top-level additions:**

```fortran
#:if not defined('VISC_STENCIL')
  #:set VISC_STENCIL = 'ME4Base'
#:endif

module calc_visc_high
  use mod_globals, only : gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : id_visc, Cp, gamma_1, Cp_over_Pr, one_third, two_third, one_twelfth
  #:if VISC_ORDER == 4
  use load_smem_visc_cent4
  #:elif VISC_STENCIL == 'ME4Base'
  use load_smem_visc_me4_base
  #:else
  use load_smem_visc_cent6
  #:endif
  ...
  real(8), parameter :: one_24 = 1.d0 / 24.d0   ! used by ME4-Base and diff4
contains
  #:if VISC_ORDER == 4
  #:include 'calc_visc_cent.f90.fypp'            ! generates interp4, diff4, calc_tau_*(4pt)
  #:elif VISC_STENCIL == 'ME4Base'
  include 'calc_visc_me4_base.f90'
  #:else
  #:include 'calc_visc_cent.f90.fypp'            ! generates interp6, diff6, calc_tau_*(6pt)
  #:endif
```

**`io_v` inside each kernel subroutine** (currently `integer, parameter :: io_v = 2`):

```fortran
#:if VISC_ORDER == 4
integer, parameter :: io_v = 1
#:else
integer, parameter :: io_v = 2
#:endif
```

**`load_smem` call** (e.g. line 77 in current calc_Ev4):

```fortran
#:if VISC_ORDER == 4
call load_smem_visc_cent4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
#:elif VISC_STENCIL == 'ME4Base'
call load_smem_visc_me4_base_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
#:else
call load_smem_visc_cent6_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
#:endif
```

**Interior path** (replaces the `if (3 <= i ...)` block in calc_Ev4):

```fortran
#:if VISC_ORDER == 4
  ! Cent4: wider interior range (io_v=1)
  if (2 <= i .and. i <= nx-2 .and. 2 <= j .and. j <= ny-2 .and. 2 <= k .and. k <= nz-2) then
    block
      real(8) :: mu4(4)
      mu4(:) = mu(i-1:i+2, j, k)
      kTx    = Cp_over_Pr * interp4(mu4) * diff4(dx(i), T(i-1:i+2, j, k))
      #:if is_les
      call calc_tau_straight_LES(mu4, mut4, u(idx-1:idx+2), vy(idx-1:idx+2), wz(idx-1:idx+2), dx(i), txx, utxx)
      ...
      #:else
      call calc_tau_straight    (mu4,       u(idx-1:idx+2), vy(idx-1:idx+2), wz(idx-1:idx+2), dx(i), txx, utxx)
      call calc_tau_cross       (mu4,       v(idx-1:idx+2), uy(idx-1:idx+2),                  dx(i), txy, vtxy)
      call calc_tau_cross       (mu4,       w(idx-1:idx+2), uz(idx-1:idx+2),                  dx(i), txz, wtxz)
      #:endif
    end block
#:elif VISC_STENCIL == 'ME4Base'
  ! ME4-Base at ORDER=6: identical to current calc_visc4 interior code
  if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
    block
      real(8) :: mu3(3)
      ...  ! UNCHANGED from current calc_Ev4 interior block
    end block
#:else
  ! Cent6: same range as ME4-Base
  if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
    block
      real(8) :: mu6(6)
      mu6(:) = mu(i-2:i+3, j, k)
      kTx    = Cp_over_Pr * interp6(mu6) * diff6(dx(i), T(i-2:i+3, j, k))
      call calc_tau_straight    (mu6, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
      call calc_tau_cross       (mu6, v(idx-2:idx+3), uy(idx-2:idx+3),                  dx(i), txy, vtxy)
      call calc_tau_cross       (mu6, w(idx-2:idx+3), uz(idx-2:idx+3),                  dx(i), txz, wtxz)
    end block
#:endif
```

The `else` (boundary) block stays: the existing 2-point difference code is shared by all three variants.

Apply the same pattern for `calc_Fv${VISC_ORDER}$` and `calc_Gv${VISC_ORDER}$` with appropriate y/z index permutations.

**Delete** `3D_solver/src/calc_visc4.f90.fypp` after creating this file.

### 4. Create `3D_solver/src/calc_visc_high_internal.f90.fypp` (from calc_visc4_internal.f90.fypp)

Same structure as `calc_visc_high.f90.fypp` but:
- Module name: `calc_visc_high_internal`
- Subroutines: `calc_Ev${VISC_ORDER}$_in` etc.
- Interior-only: no `else` boundary block (thread returns if outside interior range)

**Delete** `3D_solver/src/calc_visc4_internal.f90.fypp` after creating this file.

### 5. Update `3D_solver/src/calc_flux_base.f90.fypp`

**VISC_ORDER default** (lines 2-4):

```fortran
#:if not defined('VISC_ORDER')
  #:set VISC_ORDER = ORDER   ← was: 4 if ORDER >= 4 else 2
#:endif
```

**`use` statements** (lines 19-23):

```fortran
#:if VISC in ('NS', 'LES')
use calc_visc2
use calc_visc_high          ← was: calc_visc4
use calc_visc_high_internal ← was: calc_visc4_internal
#:endif
```

No changes to the kernel call lines — `calc_${dir}$v${VISC_ORDER}$` already generates the right names for both order 4 and 6.

### 6. Update `3D_solver/CMakeLists.txt`

In `_SHARED_FYPP` list:
- `"calc_visc4"` → `"calc_visc_high"`
- `"calc_visc4_internal"` → `"calc_visc_high_internal"`
- Add `"load_smem_visc_cent"` (new compiled module)

In `_BASE` list:
- Fix bug: `"load_smem_visc4.f90"` → `"load_smem_visc_me4_base.f90"` (file was already renamed in git)

In `FYPP_FLAGS` (line 29):
- Add `3D_solver/src/` to fypp search path so `#:include 'calc_visc_cent.f90.fypp'` resolves:

```cmake
set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${CMAKE_CURRENT_SOURCE_DIR}/../src")
```

---

## Files Created / Deleted Summary

| Action | File |
|--------|------|
| Create | `3D_solver/src/calc_visc_high.f90.fypp` |
| Create | `3D_solver/src/calc_visc_high_internal.f90.fypp` |
| Create | `3D_solver/src/calc_visc_cent.f90.fypp` (fypp-include helper, not compiled standalone) |
| Create | `3D_solver/src/load_smem_visc_cent.f90.fypp` (compiled module) |
| Delete | `3D_solver/src/calc_visc4.f90.fypp` |
| Delete | `3D_solver/src/calc_visc4_internal.f90.fypp` |
| Delete | `3D_solver/src/calc_visc_cent4.f90` |
| Delete | `3D_solver/src/calc_visc_cent6.f90` |
| Delete | `3D_solver/src/load_smem_visc_cent4.f90` |
| Delete | `3D_solver/src/load_smem_visc_cent6.f90` |
| Modify | `3D_solver/src/calc_flux_base.f90.fypp` |
| Modify | `3D_solver/CMakeLists.txt` |

---

## Behavioral Impact on Existing Cases

| Case | ORDER | Before | After |
|------|-------|--------|-------|
| ORDER=4 cases (DHIT, ETGV, IVST, KHI, SBLI) | 4 | ME4-Base (calc_Ev4) | Cent4 (calc_Ev4, different stencil) |
| ORDER=6 cases (NSTGV, TBL…) | 6 | ME4-Base (calc_Ev4) | ME4-Base (calc_Ev6, same stencil) |
| New Cent6 cases | 6 | — | Cent6 (calc_Ev6) via VISC_STENCIL='Central' |

ORDER=4 cases will see a stencil change (ME4-Base → Cent4). ORDER=6 cases retain ME4-Base by default.

---

## Verification

Build an ORDER=6 case (ME4-Base default, no behavior change in stencil):
```bash
cd 3D_solver/NSTGV
cmake -B build && cmake --build build -j
# Inspect build/calc_visc_high.f90: should contain 'load_smem_visc_me4_base_x'
```

Build an ORDER=4 case (Cent4):
```bash
cd 3D_solver/DHIT   # or any ORDER=4 case
cmake -B build && cmake --build build -j
# Inspect build/calc_visc_high.f90: should contain 'load_smem_visc_cent4_x' and 'interp4'
```

Test Cent6 manually:
1. In `3D_solver/NSTGV/config.fypp` add `#:set VISC_STENCIL = 'Central'` (ORDER=6)
2. Rebuild; `build/calc_visc_high.f90` should contain `load_smem_visc_cent6_x` and `interp6`

Check `build/load_smem_visc_cent.f90` for each case to confirm correct io_v and stencil coefficients.
