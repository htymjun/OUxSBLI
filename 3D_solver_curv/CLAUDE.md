# CLAUDE.md — 3D_solver_curv/

This file documents the curvilinear O-grid solver. See the repository root `CLAUDE.md` for project overview, shared configuration reference, and tooling rules.

## Architecture

Implements a 2D O-grid curvilinear mesh (ξ, η plane) with a uniform spanwise z direction. Entry point is `main_curv.f90` instead of `main.f90`.

Key differences from the Cartesian solver:

- **Grid metrics**: `xi_x, xi_y, eta_x, eta_y` (chain-rule coefficients), `n_xi_x, n_xi_y` (area-scaled face normals), and `Jacobian(i,j) = 1/(J_2D·dz)` are stored on device.
- **Conservative variable storage**: `QJ = Q_physical × J_2D × dz`; reconstructed to primitive `Q` via `calc_quantities_3D` / `calc_quantities_T_3D` before kernel dispatch.
- **Flux scaling**: E and F fluxes are area-scaled (include S_ξ or S_η factor); G flux is not (J_2D is folded into `dt_Szeta = dt·J_2D` in `calc_steps_curv.f90`).
- **Computational spacing**: Δξ = Δη = 1 (unit); physical z spacing is the dimensional `dz` passed as an argument.
- **Scheme dispatch**: Only KEEP (`integer(2)`), SLAU (`real(2)`), and Hybrid (`real(8)`) are supported. **LES (`id_visc = integer(8)`) is not implemented** in `calc_flux_base_curv.f90`; only Euler and NS are dispatched.
- **Viscous kernels**: `calc_visc2_curv.f90` (Gaitonde & Visbal, curvilinear); physical gradients use chain rule (∂f/∂x = ξ_x·∂f/∂ξ + η_x·∂f/∂η).
- **Per-case config**: `3D_solver_curv/<CASE>/` containing `mod_globals.f90`, `set.f90`, `Makefile`, `calc.sh`. Available cases: NACA (O-grid airfoil), CORN (compression corner, M=2, θ=8°).
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
