# Add a fused conv+visc kernel variant to the 1D solver

## Context

`1D_solver/ST` currently computes the flux in two separate GPU kernel launches per
RK stage: `calc_keep_x` (convective, `1D_solver/src/calc_keep_kernel.f90.fypp`) writes
`E` fresh, then `calc_Ev2` (viscous, `1D_solver/src/calc_visc2.f90.fypp`) re-loads `u`
into its own shared-memory tile and subtracts the viscous stress/heat-flux terms into
the same `E`. Since `1D_solver` exists specifically as a small testbed for studying
GPU occupancy/kernel-launch overhead for the mixed-precision solver work, the next
experiment is a **fused** variant that computes both convective and viscous flux in a
single kernel launch (one shared-memory load pass, one write of `E`), to compare
against the existing split-launch design. The split variant must keep working exactly
as it does today — this is an additive variant, selected at compile time, not a
replacement.

Investigated the rest of the repo for precedent first: `3D_solver_mixed`'s
`OVERLAP_MODE` toggle (`config.fypp`/`calc_flux_base_mixed.f90.fypp`) is *not* kernel
fusion — every mode still launches convective and viscous as separate kernels; it only
changes stream/Green-Context scheduling between them. There is no "fused kernel"
precedent anywhere in this repo (`3D_solver/src/calc_flux_base.f90.fypp` even has a
comment explicitly choosing to keep viscous kernels split, for register-pressure
reasons). What *is* an established pattern is factoring shared flux math into a small
`attributes(device)` helper `#:include`d into more than one kernel module —
`src/calc_visc_cent.f90.fypp` is shared across 4+ 3D-solver kernel modules this way,
and `1D_solver/src/calc_keep_1d.f90` already does this for `KEEP2` (included into
`calc_keep_kernel.f90.fypp`). This plan follows that pattern: factor the viscous
formula out of `calc_Ev2` into a small shared device function, then `include` it from
both the existing split kernel and the new fused kernel, per the user's confirmed
preference (not a zero-touch duplicate).

## Approach

**New config.fypp toggle** (`1D_solver/ST/config.fypp`), compile-time dispatch,
consistent with every other `config.fypp` string toggle (`VISC`, `SCHEME`, ...) and
with the codebase's no-runtime-branching philosophy:
```python
#:set KERNEL_MODE = 'split'   # 'split' (default, unchanged) or 'fused'
```
Default stays `'split'` — behavior is unchanged unless a case deliberately opts in.

**Factor viscous math out for reuse** — new `1D_solver/src/calc_visc_1d.f90` (plain
`.f90`, `include`d, mirroring `calc_keep_1d.f90`'s style exactly):
```fortran
pure attributes(device) function VISC2(u, T, mu) result(Fv)
  real(8), intent(in), dimension(2) :: u, T, mu
  real(8) Fv(2)  ! Fv(1) = txx (viscous stress), Fv(2) = viscous work + heat flux
  real(8) mudx, mux
  mudx  = 0.5d0 * (mu(1) + mu(2)) / dx
  Fv(2) = Cp_over_Pr * mudx * (-T(1) + T(2))
  mux   = mudx * (-u(1) + u(2))
  Fv(1) = two_third * 2.d0 * mux
  Fv(2) = Fv(2) + 0.5d0 * (u(1) + u(2)) * Fv(1)
end function VISC2
```
Relies on `dx`/`Cp_over_Pr`/`two_third` via the including module's own `use`
statements, exactly like `KEEP2` relies on the including module's `R_over_gamma_1` —
no new `use` inside the include file itself.

**Refactor the existing split kernel** (`1D_solver/src/calc_visc2.f90.fypp`) to
`include 'calc_visc_1d.f90'` and call `VISC2` instead of computing `mudx`/`mux`/`txx`
inline in `calc_Ev2` — same numbers, same kernel name/signature/launch config, just
DRY'd. Verify with a rebuild + rerun that `Q.dat` is byte-for-byte (or numerically
identical to solver tolerance) versus the current output before moving on, since this
is the one edit to already-verified code.

**New fused kernel** — `1D_solver/src/calc_fused_kernel.f90.fypp` (new file, mirrors
`calc_keep_kernel.f90.fypp`'s structure): `include`s both `calc_keep_1d.f90` (for
`KEEP2`) and `calc_visc_1d.f90` (for `VISC2`), and defines one
`attributes(global) subroutine calc_fused_x(nx, Q, T, mu, E)` that:
- loads `rho, u, p, tmp(=T)` into shared memory once (identical tile/halo pattern to
  `calc_keep_x`), and reads `mu(i)`/`mu(i+1)` directly from global — matching
  `calc_Ev2`'s existing choice not to tile `mu` (kept as-is for a minimal, faithful
  diff rather than introducing new tiling decisions in this change);
- computes `Fk = KEEP2(...)` and `Fv = VISC2(...)` once `syncthreads()` has completed;
- writes `E(1,i) = Fk(1)`, `E(2,i) = Fk(2) - Fv(1)`, `E(3,i) = Fk(3) - Fv(2)` as a
  single `intent(out)` write (no read-modify-write into `E`, unlike the split path's
  `calc_Ev2` which subtracts into an already-written buffer).

**Dispatch** — `1D_solver/src/calc_flux_base.f90.fypp`'s `calc_E` branches on
`KERNEL_MODE`:
```fortran
#:if KERNEL_MODE == 'fused'
use calc_fused_kernel
#:endif
...
call calc_quantities_T_1D(nx, Q, ruvwp, T, mu)
#:if KERNEL_MODE == 'fused'
call calc_fused_x<<<blocksE,threadsE>>>(nx, ruvwp, T, mu, E)
#:else
call calc_conv(nx, ruvwp, T, E)
call calc_Ev2<<<blocksEv,threadsEv>>>(nx, ruvwp, T, mu, E)
#:endif
```
`calc_conv`/`calc_keep_x`/`calc_Ev2` and their public names are untouched — only
`calc_E`'s body gains the branch. Both `calc_keep_kernel.f90.fypp` and the new
`calc_fused_kernel.f90.fypp` are unconditionally fypp-preprocessed and compiled into
the executable regardless of `KERNEL_MODE` (same precedent as 2D_solver always
compiling both `calc_keep_kernel` and `calc_keep_kernel_internal` regardless of
`BC_X` — only the *call site* is chosen at compile time via the `#:if`, the unused
side just goes uncalled, which is harmless).

**Build wiring** — `1D_solver/CMakeLists.txt`: add `"calc_fused_kernel"` to the
`_SHARED_FYPP` list (so it gets fypp-preprocessed like `calc_keep_kernel`/`calc_visc2`
and its generated `.f90` lands in `_GEN`, already appended to `add_executable`).
`calc_visc_1d.f90` needs **no** CMakeLists change — like `calc_keep_1d.f90` today, it's
a plain-text Fortran `include`, resolved via the compiler's existing
`target_include_directories(... "${CMAKE_CURRENT_SOURCE_DIR}/../src")` search path,
not a separate compilation unit.

## Files touched

```
1D_solver/ST/config.fypp                (edit — add KERNEL_MODE, default 'split')
1D_solver/src/calc_visc_1d.f90           (new — shared VISC2 device function)
1D_solver/src/calc_visc2.f90.fypp        (edit — calc_Ev2 calls VISC2 instead of inlining)
1D_solver/src/calc_fused_kernel.f90.fypp (new — calc_fused_x, single-launch conv+visc)
1D_solver/src/calc_flux_base.f90.fypp    (edit — calc_E branches on KERNEL_MODE)
1D_solver/CMakeLists.txt                 (edit — add calc_fused_kernel to _SHARED_FYPP)
```

No changes needed to `calc_keep_kernel.f90.fypp`, `calc_keep_1d.f90`, `calc_steps.f90`,
`calc_time_dev.f90.fypp`, `preprocess.f90.fypp`, `main.f90`, `set.f90`,
`mod_globals.f90`, or `print_1d.f90` — the RK loop and everything above the flux
dispatcher is unaffected either way.

## Verification

1. **Split path unchanged**: with `config.fypp` left at `KERNEL_MODE = 'split'`,
   rebuild (`cmake --build build -j`) and rerun (`mpirun -n 1 ./a.out`); diff the new
   `Q.dat` against the previously-verified Sod-tube output (same rarefaction/contact/
   shock structure, same values to floating-point tolerance) to confirm the
   `calc_Ev2`→`VISC2` refactor changed nothing numerically.
2. **Fused path correctness**: set `KERNEL_MODE = 'fused'`, rebuild, rerun, and confirm
   `Q.dat` matches the split path's output to floating-point/round-off tolerance (same
   physics, same formula, just fused) — no NaNs, same wave structure via
   `gnuplot ../plot.gnu`.
3. **Generated-code sanity**: inspect `build/calc_flux_base.f90` for whichever mode is
   active to confirm only one call path was emitted (no leftover dead branch text),
   per CLAUDE.md's standing guidance to check fypp output after touching `.fypp`
   templates.
