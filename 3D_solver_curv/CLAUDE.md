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

### Results are not reproducible from run to run (open, cause not investigated)

Found on 2026-11-15 while checking refactors bit for bit. Running the **same binary**
(built from `a7e56b0`, before any cleanup) twice from the same initial condition gives
different VTK files and a different `recal/Q00001.dat`. The Cartesian solvers do not show
this: 26 configurations of `3D_solver`/`2D_solver` are bit-identical between runs.

Measured with 40 time steps (`nt=20`, `np=2`, 2 MPI ranks, RTX PRO 6000, HPC SDK 26.5),
comparing `recal/Q00001.dat` of different runs:

| case | runs | distinct results | where they differ | size |
|---|---|---|---|---|
| CORN (Euler, Hybrid, ORDER=2, 192x80x4) | 6 | 3 | 18 points per pair of outcomes: i = 2..10 on one row, j = 6 or j = 7, planes k = 1 and 3 only | up to 4e-4 relative |
| NACA (NS, Hybrid, ORDER=2, 513x161x17) | 4 | 4 | 13,891 of 1.4 M points: i = 1..5 and 509..513 (both sides of the O-grid cut) plus a few other columns, j = 1..35, every k | up to O(1) and larger on near-zero quantities |

What this does and does not say:

* The differences are far above round-off and, in CORN, stay in the same few cells after 40
  steps instead of spreading with the stencil. That points to cells that are written
  inconsistently (a race between kernels or thread blocks, or a read of memory nobody
  initialised) rather than to a non-deterministic sum.
* CORN is Euler, so the viscous loaders (`load_smem_visc2_curv.f90`,
  `pipelineMemcpyAsync`) cannot be the only cause.
* NACA differs around the periodic cut of the O-grid (`BC_X=False`), CORN next to the inlet
  (`BC_X=True`); whether the two share a cause is unknown.
* Nothing here has been fixed or narrowed down further.

Consequence: a change to `3D_solver_curv` cannot be validated by comparing run outputs bit
for bit. Compare the fypp-expanded sources instead, and treat small differences in CORN/NACA
results with suspicion until this is understood. To reproduce: build a case, run it twice
with a short `nt`/`np`, and `cmp` the two `recal/Q00001.dat`.
