# CLAUDE.md — 3D_solver_curv/

This file documents the curvilinear O-grid solver. See the repository root `CLAUDE.md` for project overview, shared configuration reference, and tooling rules.

## Architecture

Implements a 2D O-grid curvilinear mesh (ξ, η plane) with a uniform spanwise z direction. Entry point is `main_curv.f90.fypp` instead of `main.f90`.

Key differences from the Cartesian solver:

- **Grid metrics**: `xi_x, xi_y, eta_x, eta_y` (chain-rule coefficients), `n_xi_x, n_xi_y` (area-scaled face normals), and `Jacobian(i,j) = 1/(J_2D·dz)` are stored on device.
- **Conservative variable storage**: `QJ = Q_physical × J_2D × dz`; reconstructed to primitive `Q` via `calc_quantities_3D` / `calc_quantities_T_3D` before kernel dispatch.
- **Flux scaling**: E and F fluxes are area-scaled (include S_ξ or S_η factor); G flux is not (J_2D is folded into `dt_Szeta = dt·J_2D` in `calc_steps_curv.f90`).
- **Computational spacing**: Δξ = Δη = 1 (unit); physical z spacing is the dimensional `dz` passed as an argument.
- **Scheme dispatch**: Only KEEP (`integer(2)`), SLAU (`real(2)`), and Hybrid (`real(8)`) are supported. **LES (`id_visc = integer(8)`) is not implemented** in `calc_flux_base_curv.f90`; only Euler and NS are dispatched.
- **Viscous kernels**: `calc_visc2_curv.f90` (Gaitonde & Visbal, curvilinear); physical gradients use chain rule (∂f/∂x = ξ_x·∂f/∂ξ + η_x·∂f/∂η).
- **Per-case config**: `3D_solver_curv/<CASE>/` containing `mod_globals.f90`, `set.f90`, `config.fypp`, `CMakeLists.txt`, `calc.sh`. Available cases: NACA (O-grid airfoil), CORN (compression corner, M=2, θ=8°).
- `load_smem_visc2.f90` (in `3D_solver/src/`) provides async pipeline shared-memory load helpers (`pipelineMemcpyAsync` / `pipelineCommit` / `pipelineWaitPrior`); requires the `wmma` module.

### Data Flow (Curvilinear)

```
main_curv.f90
  └─ set_grid_curv()         # O-grid generation, metric coefficients, Jacobians
  └─ set() [set.f90]         # IC, BC (case-specific)
  └─ calc_time_dev_curv()    # time-loop entry point
        └─ calc_EFG_curv()   # calc_flux_base_curv: convective + viscous kernels
        └─ calc_R_curv()     # flux divergence
        └─ calc_steps_curv() # RK stage update
        └─ calc_para()       # MPI ghost-cell exchange
        └─ print()           # VTK output
```

## Known issues

### Run-to-run irreproducibility of the Hybrid/SLAU sensor (fixed 2026-11-16)

Until `10219a9` the same binary gave different VTK and restart files on repeated runs:

| case | runs | distinct results | where they differed | size |
|---|---|---|---|---|
| CORN (Euler, Hybrid, ORDER=2, 192x80x4) | 6 | 3 | i = 2..10 on one row (j = 6 or 7), one of the two interior z planes | up to 4e-4 relative |
| NACA (NS, Hybrid, ORDER=2, 513x161x17) | 4 | 4 | both sides of the O-grid cut (i = 1..5, 509..513), j = 1..35, every k | up to O(1) |

Cause: `calc_Ducros_curv` filled the xi ghost columns of the sensor inside the kernel,
`fd(1,j,k) = fd(nx-1,j,k)` and `fd(nx,j,k) = fd(2,j,k)`. The right-hand sides are cells that
*other* threads of the same launch write, so the copy picked up either the new value or the
one left from the previous RK stage (unset memory on the first call), depending on thread
timing. The xi-face fluxes at i = 1 and i = nx-1 compare that sensor with `threshold` to
choose KEEP or SLAU, so a different ghost value switches the scheme on those faces: a
discrete, far-above-round-off difference that starts at the cut (NACA) or at the inlet/outlet
columns (CORN).

Fix: the kernel only fills the eta and z ghosts (copies of the thread's own value, as in the
Cartesian `calc_Ducros`); `calc_conv_curv` fills the xi ghost columns with a
`!$cuf kernel do` loop after the kernel has finished. CORN (9 runs) and NACA (7 runs) are
now bit-identical from run to run, and `test_corn` passes. The fixed CORN result equals one
of the three old outcomes.

Rule: a kernel must not read an output cell that a different thread of the same launch
writes. Ghost fills that are not copies of the thread's own cell go in a separate pass.

### Open: the sensor's xi ghost wraps periodically even when `BC_X = True`

The wrap above is right for the O-grid (NACA, `BC_X = False`). CORN has `BC_X = True`
(inflow/outflow), yet its sensor ghost at the inlet still takes the value of the outlet
column and vice versa, so the faces i = 1 and i = nx-1 choose KEEP/SLAU from the mean of an
inlet and an outlet sensor value. Zero-gradient ghosts (`sensor(1) = sensor(2)`,
`sensor(nx) = sensor(nx-1)`) under `#:if BC_X` would be the consistent choice, but it changes
CORN's results, so it has not been changed.
