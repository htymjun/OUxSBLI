# Numerical Methods

This document describes the governing equations, discretisation schemes, and time-integration methods implemented in OUxSBLI.

---

## Governing Equations

OUxSBLI solves the compressible Navier-Stokes equations in conservative form. The state vector and flux arrays are:

$$
\mathbf{Q} = \begin{pmatrix} \rho \\ \rho u \\ \rho v \\ \rho w \\ \rho E \end{pmatrix}, \qquad
\frac{\partial \mathbf{Q}}{\partial t} + \frac{\partial \mathbf{E}}{\partial x} + \frac{\partial \mathbf{F}}{\partial y} + \frac{\partial \mathbf{G}}{\partial z} = \frac{\partial \mathbf{E}_v}{\partial x} + \frac{\partial \mathbf{F}_v}{\partial y} + \frac{\partial \mathbf{G}_v}{\partial z}
$$

where ρ is density, (u, v, w) are velocity components, and E is total energy per unit mass.

### Array layout

Conservative variables are stored as `Q(nx, 5, ny, nz)` = [ρ, ρu, ρv, ρw, ρE]. Flux arrays use trimmed index ranges: E-flux drops boundary points in y/z, F-flux in x/z, G-flux in x/y.

### Physics modes

| `VISC` | Equations solved |
|--------|-----------------|
| `'Euler'` | Inviscid compressible Euler (no viscous terms) |
| `'NS'` | Full compressible Navier-Stokes with Sutherland viscosity |
| `'LES'` | NS with selective mixed-scale subgrid-stress model |

---

## Convective Schemes

All four schemes are dispatched from `calc_flux_base.f90.fypp` at compile time. Each scheme has a boundary kernel (handles ghost cells) and an interior-only kernel (avoids redundant boundary checks for better GPU occupancy).

### KEEP — Kinetic Energy and Entropy Preserving

`SCHEME = 'KEEP'`

The KEEP scheme is a split-form finite-difference method that discretely preserves kinetic energy and entropy on uniform grids. It produces no numerical dissipation, making it ideal for smooth vortical flows (Taylor-Green vortex) where artificial dissipation would contaminate the physics.

- **Files:** `calc_keep_kernel.f90.fypp`, `calc_keep_kernel_internal.f90.fypp`
- **Best for:** ETGV
- **Not suitable for:** flows with shocks (no upwinding)

### SLAU — Simple Low-dissipation AUSM

`SCHEME = 'SLAU'` or `SCHEME = 'SLAU'` with `SLAU_VARIANT = 'HRSLAU2'`

SLAU is a pressure-based all-speed Riemann solver from the AUSM family. HRSLAU2 (High Resolution SLAU2) contains less artificial dissipation. Both variants add upwind dissipation that suppresses numerical instabilities at shocks and contact discontinuities.

- **Files:** `calc_slau_kernel.f90.fypp`, `calc_slau_kernel_internal.f90.fypp`
- **Variants:** `SLAU_VARIANT = 'SLAU'` or `'HRSLAU2'` (default)
- **Best for:** SBLI, TBL, BL, OS — any case with shocks

### Hybrid — KEEP ↔ SLAU via Ducros Sensor

`SCHEME = 'Hybrid'`

The Hybrid scheme blends KEEP (smooth regions) and SLAU (near shocks) using the Ducros dilatation-based sensor:

$$
\Phi_D = \frac{(\nabla \cdot \mathbf{u})^2}{(\nabla \cdot \mathbf{u})^2 + |\nabla \times \mathbf{u}|^2 + \epsilon}
$$

When Φ_D ≈ 0 (vorticity dominated, smooth flow) the scheme reduces to KEEP. When Φ_D ≈ 1 (dilatation dominated, near a shock) it switches to SLAU. This makes the Hybrid scheme well-suited to turbulent flows that also contain shock waves (TBL, SBLI at higher resolutions).

- **Files:** `calc_hybrid_kernel.f90.fypp`, `calc_hybrid_kernel_internal.f90.fypp`
- **Best for:** TBL, mixed smooth/shocked regions

---

## Viscous Schemes

### Gaitonde-Visbal 2nd-order (`ORDER = 2` or `VISC_ORDER = 2`)

Central-difference discretisation of viscous fluxes following Gaitonde and Visbal's curvilinear-consistent formulation. Used in both the Cartesian (`calc_visc2.f90.fypp`) and curvilinear (`calc_visc2_curv.f90`) solvers.

### ME4-Base 4th/6th-order (`ORDER = 4` or `6`, or `VISC_ORDER = 4` or `6`)

Higher-order compact viscous stencil from the ME4-Base scheme. Provides higher accuracy for well-resolved DNS/LES at the cost of a wider stencil and additional MPI halo width.

- **3D Cartesian solver (`3D_solver/`):** `calc_visc_high.f90.fypp` / `calc_visc_high_internal.f90.fypp` — one shared, boundary-aware kernel pair for both 4th and 6th order, dispatched on `VISC_ORDER`.
- **2D solver (`2D_solver/`):** `calc_visc4.f90.fypp` / `calc_visc4_internal.f90.fypp` — 4th-order only; the 2D solver has not been refactored to the shared `_high` kernels.

The viscous stencil order can be set independently of the convective order using `VISC_ORDER` in `config.fypp`.

---

## Time Integration

### TVD Runge-Kutta 3 (`RK = 3`)

Shu-Osher 3-stage, 3rd-order total-variation-diminishing Runge-Kutta. Default for all cases. Stability limit: CFL ≤ 1 (in practice CFL ≈ 0.03 for high-order stencils with explicit viscous terms).

### Classical Runge-Kutta 4 (`RK = 4`)

The standard 4-stage, 4th-order classical RK method. Higher accuracy per step but no TVD property. Use for smooth inviscid problems where stability is not the dominant concern.

---

## MPI Decomposition

### z-direction with overlap (`COMMZ = True`)

Enabled for the STZ validation case. Non-blocking MPI halo exchange is posted, interior flux kernels run while MPI is in flight, and halo fluxes are computed after the exchange completes. The overlap depth equals `ORDER // 2` (1, 2, or 3 ghost cells for 2nd, 4th, or 6th order).

> **Constraint:** `COMMZ = True` and `RESCALE = True` cannot be combined.

---

## Curvilinear Solver

`3D_solver_curv/` implements a 2D O-grid (ξ, η) with a uniform spanwise z direction. Key differences from the Cartesian solver:

- Grid metric coefficients (ξ_x, ξ_y, η_x, η_y) and the Jacobian are stored on device.
- Conservative variables are stored as Q × J_2D × Δz; they are converted to primitive form before kernel dispatch.
- E and F fluxes are area-scaled (include the face-area factor S_ξ or S_η); G flux is not.
- Only KEEP, SLAU, and Hybrid are supported; LES is not implemented in `calc_flux_base_curv.f90`.
