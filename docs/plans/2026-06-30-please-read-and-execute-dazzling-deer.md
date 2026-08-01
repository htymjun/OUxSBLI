# Plan: Unify cent files + rename calc_visc4 → calc_visc_high

## Context

The `calc_visc_cent4.f90` and `calc_visc_cent6.f90` tau-function include files, and their loader counterparts, are near-identical files that only differ in stencil size. This refactor merges them into single fypp-parameterized templates, and renames `calc_visc4` → `calc_visc_high` to support a three-way dispatch at ORDER=4/6: Cent4, ME4-Base (default), or Cent6 (opt-in via `VISC_STENCIL = 'Central'`).

Source plan doc: `docs/plans/2026-06-26-now-this-code-has-stateless-cascade.md`

---

## Step 1 — CREATE `3D_solver/src/calc_visc_cent.f90.fypp`

**fypp include file only** — NOT standalone compiled, NOT in CMakeLists.

Base on `calc_visc_cent4.f90`. Parametrize with `#:set N = VISC_ORDER` (set by caller before `#:include`). Replace:
- Array sizes `4` → `${N}$` in all 4 subroutines and 2 functions
- Function names `interp4`/`diff4` → `interp${N}$`/`diff${N}$`
- interp formula: `#:if N == 4` → cent4 formula; `#:else` → cent6 formula  
- diff formula: `#:if N == 4` → cent4 formula using `one_24`; `#:else` → cent6 formula (hardcoded coefficients from `calc_visc_cent6.f90`)

`one_24` and `two_third` are in scope from the surrounding `calc_visc_high` module. No `#:include 'config.fypp'` — caller has already included it.

---

## Step 2 — CREATE `3D_solver/src/load_smem_visc_cent.f90.fypp`

**Standalone compiled module.** Add `"load_smem_visc_cent"` to `_SHARED_FYPP` in CMakeLists.

Header preamble:
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
```

Public: 4 subroutines with `_cent${VISC_ORDER}$` suffix.

Each subroutine has `integer, parameter :: io_v = ${io_v_val}$`.

**Base on `load_smem_visc_cent4.f90`** (bug-free reference) — copy overall structure verbatim, then:
- Replace hardcoded `3 <= j <= ny-2` with `${D_MIN}$ <= j <= ny-${D_SFX}$` etc.
- Replace 4th-order gradient formulas with `#:if VISC_ORDER == 4` / `#:else` blocks
- Fix the **bugs** from `load_smem_visc_cent6.f90` by using correct `inv_d` axis:

| Subroutine | Gradient var | Correct axis |
|------------|-------------|--------------|
| `_x` | `uy`, `vy` | `inv_dy(j)` ✓ |
| `_x` | `uz`, `wz` | `inv_dz(k)` (was `inv_dy(j)`) |
| `_y` | `ux`, `vx` | `inv_dx(i)` (was `inv_dy(j)`) |
| `_y` | `vz`, `wz` | `inv_dz(k)` (was `inv_dy(j)`) |
| `_z` | `ux`, `wx` | `inv_dx(i)` (was `inv_dy(j)`) |
| `_z` | `vy`, `wy` | `inv_dy(j)` ✓ |
| `_z_koff` | same as `_z` | same fixes |

`_z_koff` differs from `_z` only by `k_base = (blockIdx%z-1)*blockDim%z + k_lo - 1`.

---

## Step 3 — CREATE `3D_solver/src/calc_visc_high.f90.fypp`

**Replaces `calc_visc4.f90.fypp`.** Based on that file with these changes:

**Header (before module):**
```fortran
#:include 'config.fypp'
#:if not defined('VISC_STENCIL')
  #:set VISC_STENCIL = 'ME4Base'
#:endif
#:set is_les = (VISC == 'LES')
#:set base_suf = ""
```

**Module header:** rename to `calc_visc_high`, use `load_smem_visc_cent${VISC_ORDER}$` / `load_smem_visc_me4_base` / `load_smem_visc_cent6` based on dispatch, public names use `${VISC_ORDER}$` instead of `4`, add `one_24` parameter. Use `#:include 'calc_visc_cent.f90.fypp'` or `include 'calc_visc_me4_base.f90'` in `contains` section per dispatch.

**In each kernel (E/F/G, plain and _koff):**

1. Replace `integer, parameter :: io_v = 2` with `#:if VISC_ORDER == 4` → `io_v = 1` / `#:else` → `io_v = 2`

2. Replace load_smem call with three-way `#:if VISC_ORDER == 4` / `#:elif VISC_STENCIL == 'ME4Base'` / `#:else` dispatch.

3. Replace the interior `if` block with three-way dispatch per the interior conditions table:

| Kernel | Cent4 | ME4-Base | Cent6 |
|--------|-------|----------|-------|
| E (x)  | `2≤i≤nx-2, 3≤j≤ny-2, 3≤k≤nz-2` | `3≤i≤nx-3, 3≤j≤ny-2, 3≤k≤nz-2` | `3≤i≤nx-3, 4≤j≤ny-3, 4≤k≤nz-3` |
| F (y)  | `3≤i≤nx-2, 2≤j≤ny-2, 3≤k≤nz-2` | `3≤i≤nx-2, 3≤j≤ny-3, 3≤k≤nz-2` | `4≤i≤nx-3, 3≤j≤ny-3, 4≤k≤nz-3` |
| G (z)  | `3≤i≤nx-2, 3≤j≤ny-2, 2≤k≤nz-2` | `3≤i≤nx-2, 3≤j≤ny-2, 3≤k≤nz-3` | `4≤i≤nx-3, 4≤j≤ny-3, 3≤k≤nz-3` |

- Cent4 interior block: `mu4(4)`, `interp4`/`diff4`, arrays `(idx-1:idx+2)` (for shared mem), `mu(i-1:i+2,j,k)` (for E-kernel)
- ME4-Base interior block: copy verbatim from `calc_visc4.f90.fypp` (uses `mu3(3)`, `flux4()`)
- Cent6 interior block: `mu6(6)`, `interp6`/`diff6`, arrays `(idx-2:idx+3)`, `mu(i-2:i+3,j,k)`

The `else` (boundary) block and output assignments `E(2,i,j-1,k-1) = ...` are **unchanged** from `calc_visc4.f90.fypp` and shared by all three paths.

End with `end module calc_visc_high`.

---

## Step 4 — CREATE `3D_solver/src/calc_visc_high_internal.f90.fypp`

Same as Step 3 but:
- Module name: `calc_visc_high_internal`
- Public names and subroutine names append `_in` (e.g. `calc_Ev${VISC_ORDER}$_in`)
- Output assignments are **inside** the interior `if` block (not after)
- **No `else` boundary block**

Based on `calc_visc4_internal.f90.fypp` structure.

---

## Step 5 — MODIFY `3D_solver/src/calc_flux_base.f90.fypp`

Two changes:

**Lines 2-4** — change VISC_ORDER default:
```fortran
#:if not defined('VISC_ORDER')
  #:set VISC_ORDER = ORDER    ! was: 4 if ORDER >= 4 else 2
#:endif
```

**Lines 19-23** — rename use statements:
```fortran
  use calc_visc_high          ! was: calc_visc4
  use calc_visc_high_internal ! was: calc_visc4_internal
```

---

## Step 6 — MODIFY `3D_solver/CMakeLists.txt`

**Line 29** — add `3D_solver/src/` to FYPP search path:
```cmake
set(FYPP_FLAGS "-I${_CASE_DIR}" "-I${CMAKE_CURRENT_SOURCE_DIR}/../src")
```
(needed so `#:include 'calc_visc_cent.f90.fypp'` inside `calc_visc_high.f90.fypp` resolves)

**Lines 34-35** in `_SHARED_FYPP` — rename visc4 entries and add cent:
```cmake
      "calc_visc_high"
      "calc_visc_high_internal"
      "load_smem_visc_cent"
```

**Line 76** in `_BASE` — fix pre-existing path bug (file was already renamed in git):
```cmake
"${CMAKE_CURRENT_SOURCE_DIR}/../src/load_smem_visc_me4_base.f90"
```
Remove the old `load_smem_visc4.f90` line (which no longer exists).

---

## Step 7 — DELETE obsolete files

After all new files are written and the build succeeds:
```
3D_solver/src/calc_visc4.f90.fypp
3D_solver/src/calc_visc4_internal.f90.fypp
3D_solver/src/calc_visc_cent4.f90
3D_solver/src/calc_visc_cent6.f90
3D_solver/src/load_smem_visc_cent4.f90
3D_solver/src/load_smem_visc_cent6.f90
```

---

## Verification

Build the NSTGV case (ORDER=4, VISC=NS, no VISC_STENCIL set → Cent4 path):
```bash
cd 3D_solver/NSTGV && cmake -B build && cmake --build build -j
```

Build the TBL case (ORDER=6, VISC=LES, Hybrid scheme → ME4-Base default path):
```bash
cd 3D_solver/TBL && cmake -B build && cmake --build build -j
```

To test Cent6 path, temporarily add `#:set VISC_STENCIL = 'Central'` to a case's config.fypp (e.g. TBL) and rebuild. All three paths must compile clean before deleting old files.

The generated `.f90` files in `build/` should be inspected for any obviously wrong fypp expansions (especially the interior condition bounds and the `#:include` resolution for `calc_visc_cent.f90.fypp`).
