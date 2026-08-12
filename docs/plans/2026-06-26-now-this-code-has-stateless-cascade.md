# Implementation Status: Unify cent files + rename calc_visc4 → calc_visc_high

**Status as of 2026-06-30: APPROVED, ZERO FILES WRITTEN YET. Start here.**

---

## What the user asked for

1. Merge `calc_visc_cent4.f90` and `calc_visc_cent6.f90` into one fypp template
2. Merge `load_smem_visc_cent4.f90` and `load_smem_visc_cent6.f90` into one fypp module
3. Rename `calc_visc4.f90.fypp` → `calc_visc_high.f90.fypp`
4. Rename `calc_visc4_internal.f90.fypp` → `calc_visc_high_internal.f90.fypp`
5. Modify the renamed files to dispatch between ME4-Base, Cent4, and Cent6

## Dispatch rules (confirmed by user)

- `VISC_ORDER = ORDER` always (change default in calc_flux_base)
- **ORDER=4** → always Cent4 (ME4-Base retired at order 4)
- **ORDER=6** → ME4-Base (default) or Cent6 via new `VISC_STENCIL = 'Central'` in config.fypp

---

## File 1 to CREATE: `3D_solver/src/calc_visc_cent.f90.fypp`

**fypp-include file — NOT standalone compiled, NOT in CMakeLists, NOT in any list.**
Used only via `#:include 'calc_visc_cent.f90.fypp'` inside the `contains` section of `calc_visc_high.f90.fypp`.

`VISC_ORDER` is in scope when included (set by config.fypp already). `two_third` and `one_24` are in scope from the surrounding module.

```fortran
#:set N = VISC_ORDER
  pure attributes(device) function interp${N}$(a) result(ai)
    real(8), intent(in), contiguous :: a(${N}$)
    real(8) ai
    #:if N == 4
    ai = 0.0625d0 * (-a(1) + 9.d0 * (a(2) + a(3)) - a(4))
    #:else
    ai = 0.00390625d0 * (3.d0*a(1) - 25.d0*a(2) + 150.d0*(a(3)+a(4)) - 25.d0*a(5) + 3.d0*a(6))
    #:endif
  end function interp${N}$

  pure attributes(device) function diff${N}$(d, a) result(da)
    real(8), intent(in)             :: d
    real(8), intent(in), contiguous :: a(${N}$)
    real(8) da
    #:if N == 4
    da = (a(1) - 27.d0 * (a(2) - a(3)) - a(4)) * one_24 * d
    #:else
    da = (-0.0046875d0*a(1) + 0.065104167d0*a(2) - 1.171875d0*(a(3)-a(4)) - 0.065104167d0*a(5) + 0.0046875d0*a(6)) * d
    #:endif
  end function diff${N}$
```

The `calc_tau_straight`, `calc_tau_straight_LES`, `calc_tau_cross`, `calc_tau_cross_LES` subroutines
all have array size `(${N}$)` and call `interp${N}$`/`diff${N}$`. Copy structure verbatim from
`calc_visc_cent4.f90`, replacing size `4` → `${N}$` and `interp4`/`diff4` → `interp${N}$`/`diff${N}$`.

---

## File 2 to CREATE: `3D_solver/src/load_smem_visc_cent.f90.fypp`

**Compiled standalone module. Add `"load_smem_visc_cent"` to `_SHARED_FYPP` in CMakeLists.**

```fortran
#:include 'config.fypp'
#:set io_v_val = 1 if VISC_ORDER == 4 else 2
#:set D_MIN = 3 if VISC_ORDER == 4 else 4
#:set D_SFX = 2 if VISC_ORDER == 4 else 3

module load_smem_visc_cent${VISC_ORDER}$
  use wmma
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  #:if VISC_ORDER == 4
  use mod_constant, only : two_third, one_twelfth
  #:else
  use mod_constant, only : one_sixty
  #:endif
  implicit none
  private
  public load_smem_visc_cent${VISC_ORDER}$_x, load_smem_visc_cent${VISC_ORDER}$_y, &
         load_smem_visc_cent${VISC_ORDER}$_z, load_smem_visc_cent${VISC_ORDER}$_z_koff
contains
  ! ... 4 subroutines ...
end module load_smem_visc_cent${VISC_ORDER}$
```

Each subroutine has `integer, parameter :: io_v = ${io_v_val}$`.

Gradient loading conditions use `${D_MIN}$ <= j <= ny-${D_SFX}$` etc.

**BUGS FIXED from load_smem_visc_cent6.f90** — all were wrong `inv_dy(j)`:
- `_x`: `uz` and `wz` multiply by `inv_dz(k)` not `inv_dy(j)`
- `_y`: `ux`, `vx` multiply by `inv_dx(i)` not `inv_dy(j)`; `vz`, `wz` multiply by `inv_dz(k)` not `inv_dy(j)`
- `_z`: `ux`, `wx` multiply by `inv_dx(i)` not `inv_dy(j)`
- `_z_koff`: same fixes as `_z`

Gradient formula template (for `uy`/`vy` in `_x`, y-stencil):
```fortran
      if (1 <= i .and. i <= nx .and. ${D_MIN}$ <= j .and. j <= ny-${D_SFX}$ .and. k <= nz) then
#:if VISC_ORDER == 4
        uy(idx) = (two_third * (-Q(i,2,j-1,k) + Q(i,2,j+1,k)) - one_twelfth * (-Q(i,2,j-2,k) + Q(i,2,j+2,k))) * inv_dy(j)
        vy(idx) = (two_third * (-Q(i,3,j-1,k) + Q(i,3,j+1,k)) - one_twelfth * (-Q(i,3,j-2,k) + Q(i,3,j+2,k))) * inv_dy(j)
#:else
        uy(idx) = (one_sixty * (-Q(i,2,j-3,k) + Q(i,2,j+3,k)) + 0.15d0 * (Q(i,2,j-2,k) - Q(i,2,j+2,k)) + 0.75d0 * (-Q(i,2,j-1,k) + Q(i,2,j+1,k))) * inv_dy(j)
        vy(idx) = (one_sixty * (-Q(i,3,j-3,k) + Q(i,3,j+3,k)) + 0.15d0 * (Q(i,3,j-2,k) - Q(i,3,j+2,k)) + 0.75d0 * (-Q(i,3,j-1,k) + Q(i,3,j+1,k))) * inv_dy(j)
#:endif
      endif
      if (1 <= i .and. i <= nx .and. j <= ny .and. ${D_MIN}$ <= k .and. k <= nz-${D_SFX}$) then
#:if VISC_ORDER == 4
        uz(idx) = (two_third * (-Q(i,2,j,k-1) + Q(i,2,j,k+1)) - one_twelfth * (-Q(i,2,j,k-2) + Q(i,2,j,k+2))) * inv_dz(k)
        wz(idx) = (two_third * (-Q(i,4,j,k-1) + Q(i,4,j,k+1)) - one_twelfth * (-Q(i,4,j,k-2) + Q(i,4,j,k+2))) * inv_dz(k)
#:else
        uz(idx) = (one_sixty * (-Q(i,2,j,k-3) + Q(i,2,j,k+3)) + 0.15d0 * (Q(i,2,j,k-2) - Q(i,2,j,k+2)) + 0.75d0 * (-Q(i,2,j,k-1) + Q(i,2,j,k+1))) * inv_dz(k)
        wz(idx) = (one_sixty * (-Q(i,4,j,k-3) + Q(i,4,j,k+3)) + 0.15d0 * (Q(i,4,j,k-2) - Q(i,4,j,k+2)) + 0.75d0 * (-Q(i,4,j,k-1) + Q(i,4,j,k+1))) * inv_dz(k)
#:endif
      endif
```

Apply same pattern for `_y` (ux/vx: i-stencil × inv_dx; vz/wz: k-stencil × inv_dz) and `_z` (ux/wx: i-stencil × inv_dx; vy/wy: j-stencil × inv_dy).

`_z_koff` differs from `_z` only by `k_base = (blockIdx%z-1)*blockDim%z + k_lo - 1` and extra `k_lo` argument.

Use `load_smem_visc_cent4.f90` as the bug-free reference for the overall subroutine structure (velocity load loop, pipelineCommit, gradient compute loop, pipelineWaitPrior, syncthreads).

---

## File 3 to CREATE: `3D_solver/src/calc_visc_high.f90.fypp`

**Replaces `calc_visc4.f90.fypp` in `_SHARED_FYPP`. Based on that file with the following changes:**

**Header (before `module`):**
```fortran
#:include 'config.fypp'
#:if not defined('VISC_STENCIL')
  #:set VISC_STENCIL = 'ME4Base'
#:endif
#:set is_les = (VISC == 'LES')
#:set base_suf = ""
```

**Module header:**
```fortran
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
  implicit none
  private
  public calc_Ev${VISC_ORDER}$, calc_Fv${VISC_ORDER}$, calc_Gv${VISC_ORDER}$, &
         calc_Ev${VISC_ORDER}$_koff, calc_Fv${VISC_ORDER}$_koff, calc_Gv${VISC_ORDER}$_koff
  real(8), parameter :: one_24 = 1.d0 / 24.d0
contains
  #:if VISC_ORDER == 4
  #:include 'calc_visc_cent.f90.fypp'
  #:elif VISC_STENCIL == 'ME4Base'
  include 'calc_visc_me4_base.f90'
  #:else
  #:include 'calc_visc_cent.f90.fypp'
  #:endif
```

**In each kernel subroutine** — replace `integer, parameter :: io_v = 2` with:
```fortran
    #:if VISC_ORDER == 4
    integer, parameter :: io_v = 1
    #:else
    integer, parameter :: io_v = 2
    #:endif
```

**Subroutine names** — replace `calc_Ev4${suf}$` with `calc_Ev${VISC_ORDER}$${suf}$` (and similarly for F, G).

**load_smem call (E-kernel example):**
```fortran
    #:if VISC_ORDER == 4
    call load_smem_visc_cent4_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    #:elif VISC_STENCIL == 'ME4Base'
    call load_smem_visc_me4_base_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    #:else
    call load_smem_visc_cent6_x(it, jt, kt, j, k, nx, ny, nz, dy, dz, Q, u, v, w, uy, vy, uz, wz)
    #:endif
```

**Interior block (E-kernel, x-direction)** — replace the current `if (3 <= i ...)` block:

```fortran
    #:if VISC_ORDER == 4
    if (2 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) :: mu4(4)
        #:if is_les
        real(8) :: mut4(4)
        mut4(:) = mut(i-1:i+2, j, k)
        #:endif
        mu4(:) = mu(i-1:i+2, j, k)
        kTx = Cp_over_Pr * interp4(mu4) * diff4(dx(i), T(i-1:i+2, j, k))
        #:if is_les
        call calc_tau_straight_LES(mu4, mut4, u(idx-1:idx+2), vy(idx-1:idx+2), wz(idx-1:idx+2), dx(i), txx, utxx)
        call calc_tau_cross_LES(mu4, mut4, v(idx-1:idx+2), uy(idx-1:idx+2), dx(i), txy, vtxy)
        call calc_tau_cross_LES(mu4, mut4, w(idx-1:idx+2), uz(idx-1:idx+2), dx(i), txz, wtxz)
        block
          real(8) :: H(4)
          H(:) = Cp * T(i-1:i+2,j,k) + 0.5d0 * (u(idx-1:idx+2)**2 + v(idx-1:idx+2)**2 + w(idx-1:idx+2)**2) + qc2(i-1:i+2,j,k)
          Hsgs = -interp4(mut4) * diff4(dx(i), H) / Prt
        end block
        #:else
        call calc_tau_straight(mu4, u(idx-1:idx+2), vy(idx-1:idx+2), wz(idx-1:idx+2), dx(i), txx, utxx)
        call calc_tau_cross(mu4, v(idx-1:idx+2), uy(idx-1:idx+2), dx(i), txy, vtxy)
        call calc_tau_cross(mu4, w(idx-1:idx+2), uz(idx-1:idx+2), dx(i), txz, wtxz)
        #:endif
      end block
    #:elif VISC_STENCIL == 'ME4Base'
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        ! UNCHANGED — copy verbatim from calc_visc4.f90.fypp lines 86-113
      end block
    #:else  ! Cent6
    if (3 <= i .and. i <= nx-3 .and. 4 <= j .and. j <= ny-3 .and. 4 <= k .and. k <= nz-3) then
      block
        real(8) :: mu6(6)
        #:if is_les
        real(8) :: mut6(6)
        mut6(:) = mut(i-2:i+3, j, k)
        #:endif
        mu6(:) = mu(i-2:i+3, j, k)
        kTx = Cp_over_Pr * interp6(mu6) * diff6(dx(i), T(i-2:i+3, j, k))
        #:if is_les
        call calc_tau_straight_LES(mu6, mut6, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
        call calc_tau_cross_LES(mu6, mut6, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
        call calc_tau_cross_LES(mu6, mut6, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
        block
          real(8) :: H(6)
          H(:) = Cp * T(i-2:i+3,j,k) + 0.5d0 * (u(idx-2:idx+3)**2 + v(idx-2:idx+3)**2 + w(idx-2:idx+3)**2) + qc2(i-2:i+3,j,k)
          Hsgs = -interp6(mut6) * diff6(dx(i), H) / Prt
        end block
        #:else
        call calc_tau_straight(mu6, u(idx-2:idx+3), vy(idx-2:idx+3), wz(idx-2:idx+3), dx(i), txx, utxx)
        call calc_tau_cross(mu6, v(idx-2:idx+3), uy(idx-2:idx+3), dx(i), txy, vtxy)
        call calc_tau_cross(mu6, w(idx-2:idx+3), uz(idx-2:idx+3), dx(i), txz, wtxz)
        #:endif
      end block
    #:endif
```

The `else` (boundary) block and `E(2,i,j-1,k-1) = ...` output assignments are **UNCHANGED** and shared by all variants.

**Apply same dispatch pattern to F-kernel and G-kernel**, permuting x/y/z roles:
- F-kernel (y): Cent4 `mu4(:) = mu(i,j-1:j+2,k)`, Cent6 `mu6(:) = mu(i,j-2:j+3,k)`, stencils `v(idx±)`, cross `wz/ux(idx±)`
- G-kernel (z): Cent4 `mu4(:) = mu(i,j,k-1:k+2)`, Cent6 `mu6(:) = mu(i,j,k-2:k+3)`, stencils `w(idx±)`, cross `ux/vy(idx±)`

**Interior conditions for all kernels:**

| Kernel | Cent4 (ORDER=4) | ME4-Base (ORDER=6, default) | Cent6 (ORDER=6, VISC_STENCIL='Central') |
|--------|-----------------|-----------------------------|-----------------------------------------|
| E (x)  | `2≤i≤nx-2, 3≤j≤ny-2, 3≤k≤nz-2` | `3≤i≤nx-3, 3≤j≤ny-2, 3≤k≤nz-2` | `3≤i≤nx-3, 4≤j≤ny-3, 4≤k≤nz-3` |
| F (y)  | `3≤i≤nx-2, 2≤j≤ny-2, 3≤k≤nz-2` | `3≤i≤nx-2, 3≤j≤ny-3, 3≤k≤nz-2` | `4≤i≤nx-3, 3≤j≤ny-3, 4≤k≤nz-3` |
| G (z)  | `3≤i≤nx-2, 3≤j≤ny-2, 2≤k≤nz-2` | `3≤i≤nx-2, 3≤j≤ny-2, 3≤k≤nz-3` | `4≤i≤nx-3, 4≤j≤ny-3, 3≤k≤nz-3` |

End with `end module calc_visc_high`.

---

## File 4 to CREATE: `3D_solver/src/calc_visc_high_internal.f90.fypp`

**Same as calc_visc_high.f90.fypp but:**
- Module name: `calc_visc_high_internal`
- Public names append `_in`: `calc_Ev${VISC_ORDER}$_in` etc.
- Subroutine names append `_in`: `calc_Ev${VISC_ORDER}$${suf}$_in`
- Output assignments `E(2,...) = ...` are **inside** the `if` block (not after)
- **No `else` boundary block** — thread exits `if` without writing

See `calc_visc4_internal.f90.fypp` for exact placement of assignments.
End with `end module calc_visc_high_internal`.

---

## Files to MODIFY

### `3D_solver/src/calc_flux_base.f90.fypp`

**Lines 2-4** (VISC_ORDER default):
```fortran
#:if not defined('VISC_ORDER')
  #:set VISC_ORDER = ORDER    ! was: 4 if ORDER >= 4 else 2
#:endif
```

**Lines 19-23** (use statements):
```fortran
  #:if VISC in ('NS', 'LES')
  use calc_visc2
  use calc_visc_high          ! was: calc_visc4
  use calc_visc_high_internal ! was: calc_visc4_internal
  #:endif
```

### `3D_solver/CMakeLists.txt`

**Line 29** (FYPP_FLAGS — add 3D_solver/src/ to search path for `#:include 'calc_visc_cent.f90.fypp'`):
```cmake
  set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${CMAKE_CURRENT_SOURCE_DIR}/../src")
```

**Lines 34-35** in `_SHARED_FYPP`:
```cmake
      "calc_visc_high"              ! was: "calc_visc4"
      "calc_visc_high_internal"     ! was: "calc_visc4_internal"
      "load_smem_visc_cent"         ! NEW (compiled module)
```

**Line 76** in `_BASE` (fix pre-existing bug — file was already renamed in git):
```cmake
      "${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc_me4_base.f90"  ! was: load_smem_visc4.f90
```

---

## Files to DELETE (after all new files are written and verified)

```
3D_solver/src/calc_visc4.f90.fypp
3D_solver/src/calc_visc4_internal.f90.fypp
3D_solver/src/calc_visc_cent4.f90
3D_solver/src/calc_visc_cent6.f90
3D_solver/src/load_smem_visc_cent4.f90
3D_solver/src/load_smem_visc_cent6.f90
```

---

## Reference files (unchanged, use for exact code structure)

- `3D_solver/src/calc_visc_cent4.f90` — bug-free cent4 tau subroutines (copy and parametrize for calc_visc_cent.f90.fypp)
- `3D_solver/src/load_smem_visc_cent4.f90` — bug-free loader reference (base for load_smem_visc_cent.f90.fypp)
- `3D_solver/src/calc_visc4.f90.fypp` — existing kernel for ME4-Base path (copy verbatim for that branch)
- `3D_solver/src/calc_visc4_internal.f90.fypp` — existing internal kernel (same)
- `3D_solver/src/calc_visc_me4_base.f90` — ME4-Base include, keep as-is
- `3D_solver/src/load_smem_visc_me4_base.f90` — ME4-Base loader, keep as-is
