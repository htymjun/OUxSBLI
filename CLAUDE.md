# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

OUxSBLI is a GPU-accelerated CFD solver for compressible flows (Euler/Navier-Stokes), written in CUDA Fortran with MPI parallelization. It targets NVIDIA GPUs via the HPC SDK and solves test cases defined in `3D_solver/<CASE>/`. A standalone 2D solver lives under `2D_solver/<CASE>/`. A curvilinear O-grid variant for wing/airfoil cases lives in `3D_solver_curv/<CASE>/`.

## Build & Run

All work happens inside a specific test-case directory. There is no top-level build.

**3D solver cases use CMake:**

```bash
cd 3D_solver/NSTGV    # or DHIT, ETGV, IVST, KHI, SBLI, STZ, TBL

cmake -B build && cmake --build build -j   # fypp preprocess + compile
cd build && mpirun -n 2 ./a.out            # run simulation
```

Profiling (where available):

```bash
cd build && bash ../profile.sh    # nsys/ncu profiling
```

For curvilinear cases (still Makefile-based):

```bash
cd 3D_solver_curv/NACA
make clean && make
bash calc.sh
```

**2D solver cases use CMake:**

```bash
cd 2D_solver/OS    # or BL, DSL, EVC, SBLI, ST
cmake -B build && cmake --build build -j   # fypp preprocess + compile
cd build && mpirun -n 1 ./a.out            # run simulation
```

**Compiler requirement:** NVIDIA HPC SDK (`mpif90` with `-cuda -acc -fast -gpu=ptxinfo,rdc,lto`). Versions 24.* and 25.* are confirmed working.

**Output format:** XML VTK files, readable with ParaView.

## Fypp Preprocessing System

The 3D solver uses [fypp](https://fypp.readthedocs.io/) to generate case-specialised Fortran from templates. CMake drives this automatically.

### How it works

Each case directory contains a **`config.fypp`** file that declares all compile-time configuration variables:

```python
#:set VISC    = 'NS'      # 'Euler', 'NS', 'LES'
#:set SCHEME  = 'SLAU'    # 'KEEP', 'SLAU', 'Roe', 'Hybrid'
#:set ORDER   = 6         # 2, 4, 6  (convective + viscous stencil order)
#:set TVD     = 'hybrid'  # 'none', 'tvd', 'hybrid'
#:set RESCALE = True       # True → SBLI reference-state rescaling
#:set COMMZ   = False      # True → z-direction MPI halo decomposition
#:set RK      = 3          # 3 (TVD-RK3) or 4 (classical RK4)
#:set RESTART = False      # True → read initial condition from file
#:set BC_X    = True       # True → wall/inflow BCs in x; False → periodic
#:set BC_Y    = True       # same for y
#:set BC_Z    = False      # same for z
#:set SLAU_VARIANT = 'HRSLAU2'  # 'SLAU' or 'HRSLAU2'
#:set VISC_ORDER   = 6     # viscous stencil order (defaults to ORDER if omitted)
#:set ORDER_IO = ORDER // 2 - 1  # ghost-cell count for I/O interpolation
```

Source files with the `.f90.fypp` extension are **templates**; CMake runs `fypp -I<case-dir> <template>.f90.fypp <output>.f90` for each one. The generated `.f90` files land in `3D_solver/<CASE>/build/`.

### fypp-preprocessed source files

| Template (in `src/` or `3D_solver/src/`) | Purpose |
|------------------------------------------|---------|
| `src/mod_constant.f90.fypp` | Kind-dispatch Fortran parameters from config values |
| `3D_solver/src/calc_flux_base.f90.fypp` | Top-level flux dispatcher (convective + viscous) |
| `3D_solver/src/calc_time_dev.f90.fypp` | RK time-stepping orchestration |
| `3D_solver/src/calc_keep_kernel.f90.fypp` | KEEP convective kernel |
| `3D_solver/src/calc_slau_kernel.f90.fypp` | SLAU convective kernel |
| `3D_solver/src/calc_roe_kernel.f90.fypp` | Roe convective kernel |
| `3D_solver/src/calc_hybrid_kernel.f90.fypp` | Hybrid KEEP↔SLAU kernel |
| `3D_solver/src/calc_keep_kernel_internal.f90.fypp` | Interior-only KEEP variant |
| `3D_solver/src/calc_slau_kernel_internal.f90.fypp` | Interior-only SLAU variant |
| `3D_solver/src/calc_roe_kernel_internal.f90.fypp` | Interior-only Roe variant |
| `3D_solver/src/calc_hybrid_kernel_internal.f90.fypp` | Interior-only Hybrid variant |
| `3D_solver/src/calc_visc2.f90.fypp` | 2nd-order viscous kernels |
| `3D_solver/src/calc_visc_high.f90.fypp` | 4th/6th-order viscous kernels (`VISC_ORDER`-controlled), boundary-aware (interior/fallback branch per point) |
| `3D_solver/src/calc_visc_high_internal.f90.fypp` | Interior-only 4th/6th-order viscous variant (used where a direction has no BC, i.e. periodic) |
| `3D_solver/src/calc_visc_cent.f90.fypp` | Shared `interp`/`diff`/`calc_tau_straight`/`calc_tau_cross` helpers, `#:include`d into both `calc_visc_high[_internal].f90.fypp` |
| `3D_solver/src/load_smem_visc_cent.f90.fypp` | Shared-memory loaders for `calc_visc_high[_internal]`'s `calc_Ev/Fv/Gv` |
| `3D_solver/src/preprocess.f90.fypp` | Device memory allocation helpers |

The 2D solver uses CMake (like 3D); `2D_solver/src/calc_flux_base.f90.fypp` is preprocessed by the per-case CMakeLists.txt using the same `fypp -I<case-dir>` pattern.

### Changing the scheme or method

Edit the **`config.fypp`** in the case directory and rebuild:

```bash
# e.g. switch NSTGV from KEEP to Hybrid
vim 3D_solver/NSTGV/config.fypp      # change SCHEME = 'Hybrid'
cmake --build build -j               # CMake detects config.fypp changed, re-runs fypp
```

Do **not** manually edit the generated `.f90` files in `build/`; they are overwritten on every build.

## Architecture

Source is split across directories; CMake/vpath merges them at build time:

- `src/` — shared templates and utilities: `mod_constant.f90.fypp`, `main.f90`, `cpu_gpu_mpi.f90`, `calc_muscl.f90`, `calc_physical_quantities.f90`, `print.f90`, `set_coordinate.f90`, `set_compressible_bl.f90`
- `3D_solver/src/` — solver kernels (mostly `.f90.fypp` templates): convective/viscous scheme modules, time-integration, preprocessing
- `3D_solver/<CASE>/` — per-case configuration: `mod_globals.f90` (grid/physical params, thread blocks), `set.f90` (grid/IC/BC), `config.fypp` (build-time scheme/method flags), `CMakeLists.txt`, `calc.sh`

### Curvilinear Solver (`3D_solver_curv/`)

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

### Data Flow (Cartesian)

```
main.f90
  └─ set_coordinate()        # grid metrics & Jacobians
  └─ set() [set.f90]         # grid, IC, BC (case-specific)
  └─ calc_time_dev()         # time-loop entry point  [fypp-generated]
        └─ calc_flux_base()  # dispatches to convective + viscous kernels  [fypp-generated]
        └─ calc_steps()      # RK stage update (Q += dt * RHS)
        └─ calc_para()       # MPI ghost-cell exchange (x) or z-halo (COMMZ)
        └─ print()           # VTK output
```

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

### 2D Solver (`2D_solver/`)

A standalone 2D solver sharing the same convective/viscous kernels as the 3D Cartesian solver. Source layout mirrors the 3D structure:

- `2D_solver/src/` — shared 2D utilities (main, grid, BCs)
- `2D_solver/<CASE>/` — per-case config: `mod_globals.f90`, `set.f90`, `config.fypp`, `CMakeLists.txt`, `calc.sh`

Available cases:

| Case | Description |
|------|-------------|
| BL   | Laminar flat-plate boundary layer (M=0.1, validated vs. Blasius) |
| DSL  | Double shear layer |
| EVC  | Euler vortex convection (grid-convergence study) |
| OS   | 2D oblique shock (M=2, θ=8°, SLAU, Euler) |
| SBLI | Oblique shock / laminar BL interaction (M=2.15, β=30.8°, Re_Xsh≈1e5) |
| ST   | Sod shock tube |

BL and SBLI share the same flat-plate setup: a uniform freestream over a
symmetry wall upstream of the leading edge and a no-slip wall downstream of it,
with the inlet column left frozen at its `set_init` state (the interior-only RK
kernels never touch i=1, and neither case defines an inlet BC). In SBLI that
frozen column is also the oblique-shock generator — above the incident-shock
trace it holds the analytic post-shock state. Both place the leading edge at
x = 0 via the `i_LE` constant in `mod_globals.f90`, so wall quantities can be
compared against theory and reference data without an origin offset.

Digitized reference data for SBLI (Moro et al., Degrez et al., Vila-Perez et
al.) lives in `2D_solver/SBLI/ref/`, together with the `Cf_all.py` / `Cp_all.py`
overlay scripts. Note that `.gitignore` excludes `2D_solver/*/data*` and `*.dat`
globally; `ref/` is kept tracked by an explicit `!2D_solver/*/ref/**` negation.

### 3D Cartesian Cases

| Case | Description |
|------|-------------|
| DHIT | Decaying homogeneous isotropic turbulence |
| ETGV | Supersonic Taylor-Green vortex (entropy-preserving) |
| IVST | Inviscid vortex / smooth test case |
| KHI  | Kelvin-Helmholtz instability |
| NSTGV| NS Taylor-Green vortex |
| SBLI | Shock-boundary layer interaction (RESCALE=True) |
| STZ  | z-direction halo exchange validation (COMMZ=True, BC_Z=True) |
| TBL  | Turbulent boundary layer (LES, Hybrid scheme) |

### Convective Schemes (dispatched from `calc_flux_base.f90.fypp`)

| Scheme | Files | Best for |
|--------|-------|----------|
| KEEP (Kinetic Energy & Entropy Preserving) | `calc_keep_kernel.f90.fypp`, `calc_keep_kernel_internal.f90.fypp` | Smooth vortical flows (TGV, KHI) |
| SLAU (Simple Low-dissipation AUSM) | `calc_slau_kernel.f90.fypp`, `calc_slau_kernel_internal.f90.fypp` | Compressible turbulence with shocks (SBLI, TBL) |
| Roe | `calc_roe_kernel.f90.fypp`, `calc_roe_kernel_internal.f90.fypp` | Supersonic/hypersonic discontinuities |
| Hybrid (KEEP↔SLAU via Ducros sensor) | `calc_hybrid_kernel.f90.fypp`, `calc_hybrid_kernel_internal.f90.fypp` | Mixed smooth/shocked regions |

Viscous discretization: `calc_visc2.f90.fypp` (Gaitonde & Visbal 2nd-order, `VISC_ORDER=2`) or `calc_visc_high.f90.fypp` + `calc_visc_high_internal.f90.fypp` (4th/6th-order, `VISC_ORDER=4` or `6`).

### Cell-Center Velocity Gradients (`ux`, `vy`, `wz`)

For `VISC_ORDER > 2` (NS/LES), `calc_div.f90` precomputes the three diagonal
velocity gradients `ux=∂u/∂x`, `vy=∂v/∂y`, `wz=∂w/∂z` at cell centers once per
RK stage (`calc_div_${VISC_ORDER}$_in`, called from `calc_flux_base.f90.fypp`
right after `calc_quantities_T_3D`) and stores them to DRAM. `calc_Ev/Fv/Gv`
(in both `calc_visc_high.f90.fypp` and `calc_visc_high_internal.f90.fypp`)
then read the two gradients they don't own directly from these arrays instead
of each direction's kernel recomputing the other two independently in shared
memory — e.g. `calc_Ev` needs `vy,wz` (not `ux`, which it derives locally from
its own `u` shared-memory tile); it reads them straight from DRAM rather than
via `load_smem_visc_cent`, which now only loads the two cross-derivatives
unique to that direction (e.g. `uy,uz` for `calc_Ev`). Non-contiguous
(y/z-sweep) reads use scalar `calc_tau_straight_s`/`interp2_${N}$_s` helpers
in `calc_visc_cent.f90.fypp` (same pattern as the existing `mu`/`T` scalar
reads) rather than gathering into a local array.
Each of `ux,vy,wz` only needs a margin in its own direction (see
`calc_div_4_in`/`calc_div_6_in`'s independent per-component gating), which is
why one un-specialized `calc_div` kernel serves every consumer regardless of
that direction's `BC_X/Y/Z` setting. `calc_div_2` is not currently wired to
any consumer (`calc_visc2.f90.fypp` computes its trace terms directly from
`Q` via face differences, not a cell-centered gradient) — it's kept for a
possible future reuse by `calc_hybrid`'s `calc_Ducros` shock sensor at 2nd
order (see the `calc_Ducros` note in Notice below).
**Known limitation:** only wired into `calc_flux_base.f90.fypp`'s non-`COMMZ`
path; no current case combines `COMMZ=True` with `VISC in ('NS','LES')` at
`VISC_ORDER>2` (the only `COMMZ=True` case, STZ, is Euler-only), so
`calc_EFG_halo` still uses the old per-kernel recompute path.

## Configuration

### config.fypp (primary interface — edit this)

All compile-time scheme/method choices live in `<CASE>/config.fypp`. The fypp preprocessor expands `.f90.fypp` templates with these values, generating specialised Fortran with no runtime branching overhead.

| Variable | Values | Effect |
|----------|--------|--------|
| `VISC` | `'Euler'`, `'NS'`, `'LES'` | Physics model |
| `SCHEME` | `'KEEP'`, `'SLAU'`, `'Roe'`, `'Hybrid'` | Convective flux scheme |
| `ORDER` | `2`, `4`, `6` | Spatial accuracy (convective + viscous) |
| `VISC_ORDER` | `2`, `4`, `6` | Override viscous stencil order (defaults to `ORDER`) |
| `TVD` | `'none'`, `'tvd'`, `'hybrid'` | TVD limiter for reconstruction |
| `SLAU_VARIANT` | `'SLAU'`, `'HRSLAU2'` | SLAU flux variant |
| `RESCALE` | `True`, `False` | SBLI reference-state rescaling |
| `COMMZ` | `True`, `False` | z-direction MPI halo decomposition (overlapped comms) |
| `RK` | `3`, `4` | Runge-Kutta stages (TVD-RK3 or classical RK4) |
| `RESTART` | `True`, `False` | Restart from checkpoint file |
| `BC_X` | `True`, `False` | Wall/inflow BCs in x (False → periodic) |
| `BC_Y` | `True`, `False` | Same for y |
| `BC_Z` | `True`, `False` | Same for z |

### mod_globals.f90 (grid, physical parameters, thread blocks)

This file is **not** preprocessed by fypp. It holds:
- Grid size (`nx`, `ny`, `nz`), domain lengths (`Lx`, `Ly`, `Lz`)
- Physical parameters (`gamma`, `R`, `Pr`, `dt`, flow conditions)
- GPU thread-block sizes as `dim3` constants (`threadsE`, `threadsF`, `threadsG`, `threadsEv`, etc.) — must be tuned to the GPU architecture and grid size

### mod_constant.f90 (generated — do not edit directly)

`src/mod_constant.f90.fypp` is preprocessed per-case into `build/mod_constant.f90`. It generates Fortran `parameter` constants that encode config choices via the **type kind** (not the value):

| Parameter | kind=2 | kind=4 | kind=8 |
|-----------|--------|--------|--------|
| `id_visc` (integer) | Euler | NS | LES |
| `id_accuracy` (integer) | 2nd order | 4th order | 6th order |
| `id_tvd` (integer) | no TVD | tvd | hybrid |
| `id_slau` (integer) | SLAU | HRSLAU2 | — |
| `id_rescale` (integer) | off | on | — |
| `id_scheme` | integer(2)=KEEP | real(2)=SLAU | real(4)=Roe / real(8)=Hybrid |
| `id_bc_x/y/z` (integer) | periodic | wall/inflow | — |

The **value** of these parameters is always 0; only the **type kind** matters for compile-time dispatch. `id_recal` (restart flag) is a Fortran `logical` (`.true.`/`.false.`).

## Conservative Variable Layout

`Q(nx, 5, ny, nz)` = `[ρ, ρu, ρv, ρw, ρE]`. Flux arrays use trimmed index ranges: E-flux drops boundary points in y/z, F-flux in x/z, G-flux in x/y. AoS and SoA hybrid memory layout is used.

## MPI Decomposition

**x-direction (default):** 1D decomposition via `calc_para.f90`. Default is 2 MPI ranks (`mpirun -n 2 a.out`), with `mygpu = myrank / 2` (2 ranks per GPU). GPU-aware MPI is optional via `id_gpumpi`.

**z-direction (COMMZ=True):** Enabled for cases like STZ. Overlapped communication: non-blocking z-halo exchange is posted, interior fluxes are computed while MPI is in flight, then halo fluxes are completed. The `overlap` depth equals `ORDER // 2` (1/2/3 for 2nd/4th/6th order). `COMMZ=True` and `RESCALE=True` cannot be combined.

## Python Test Suite

Integration and convergence tests live in `ouxsbli/tests/`. Run with pytest from the repository root:

```bash
pytest ouxsbli/tests/
```

| Test file | What it checks |
|-----------|----------------|
| `test_etgv.py` | Supersonic Taylor-Green vortex (energy decay) |
| `test_st.py` | Sod shock tube (exact Riemann solution) |
| `test_evc.py` | Euler vortex convergence — KEEP 2nd/4th/6th, SLAU 2nd; expected order ≥1.5/3.5 |
| `test_os.py` | 2D oblique shock — pre/post state vs. Rankine-Hugoniot (tol 2%/5%) |
| `test_corn.py` | 3D_solver_curv/CORN — pressure and density ratios vs. θ-β-M theory (tol 5%) |
| `test_bl.py` | 2D laminar BL — Cf and u-profile vs. Blasius (~3 min) |
| `test_sbli.py` | 2D shock/BL interaction — Cp and separation bubble vs. Moro et al. (~8 min, `slow`) |

`test_sbli.py` carries the `slow` marker; deselect it with `pytest -m "not slow"`.
Both new tests patch the case down from its production settings (BL raises dt
~17x to CFL≈0.18; SBLI coarsens to 276×257 and raises dt 10x) and assert a
quasi-steady guard between the last two snapshots before comparing, so a
failure to converge reports itself instead of surfacing as a tolerance miss.
`test_bl.py` normalises by the *local* boundary-layer edge state rather than the
nominal freestream — displacement growth accelerates the flow ~2% in the 12 mm
tall domain, which is case geometry rather than solver error.

Analytical helpers in `ouxsbli/tests/utils/`:
- `oblique_shock.py` — `beta_from_theta()`, `post_shock_state()` via bisection on the θ-β-M relation
- `sod_exact.py` — exact Riemann solver for the Sod shock tube
- `vtk_reader.py` — VTK output reader

Post-processing shared by the tests and the user-facing plotting scripts lives
in `ouxsbli/analysis/`:
- `wall.py` — Sutherland viscosity, one-sided wall derivative, `compute_cf_cp()` / `wall_coeffs_from_vtr()`, `edge_state()`, `find_zero_crossings()`, `load_reference()`
- `blasius.py` — RK4-integrated Blasius similarity solution (`fprime()`, `eta()`, `cf_blasius()`); accurate to ~1e-5, unlike the coarse `fp_tab` in `src/set_compressible_bl.f90`, which is an IC seed only and deviates by up to 7%

## Tutorials

`tutorials/test.py` generates `tutorials/ouxsbli_bl/` — a patched copy of the
`2D_solver/BL` flat-plate case (dimensional parameters, wall-normal grid
stretching, Riemann-invariant top BC) built through the `ouxsbli.Case` API. The
directory is not checked in; run the script to create it. Note that 2D output
requires **two** MPI ranks (the even rank computes, the odd rank writes VTK), so
`mpirun -n 1` produces no files.

## Adding a New Test Case

For 3D Cartesian:
1. `cp -r 3D_solver/NSTGV 3D_solver/MYCASE`
2. Edit `mod_globals.f90` — grid size, physical parameters, thread block sizes
3. Edit `set.f90` — grid generation, initial conditions, boundary condition calls
4. Edit `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
5. Edit `calc.sh` — MPI rank count and runtime args
6. `cmake -B build && cmake --build build -j`

For 2D cases:
1. `cp -r 2D_solver/OS 2D_solver/MYCASE`
2. Edit `mod_globals.f90` — grid size, physical parameters
3. Edit `set.f90` — grid generation, initial conditions, boundary condition calls
4. Edit `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
5. Edit `calc.sh` — runtime args
6. `cmake -B build && cmake --build build -j`

## Notice
* From an occupancy perspective, the subroutines invoked within `calc_flux_base.f90` should not be executed on separate streams.
* `calc_flux_base.f90` and `calc_steps.f90` are the main bottleneck. You should optimize them.
* `calc_flux_base.f90` calls `calc_*_kernel.f90`, `calc_*_kernel_internal.f90`, `calc_visc*.f90`, and `calc_div.f90` (runs once per RK stage, right after `calc_quantities_T_3D`, for `VISC_ORDER>2`). They are the main bottleneck.
* Roe scheme is not used. KEEP, SLAU, Hybrid schemes should be optimized.
* `id_accuracy` controls ghost-cell count for both convective and viscous stencils.
* Do not add `contiguous` and `shared` attributes when passing shared memory as an argument.
* `cpu_gpu_mpi.f90` will be modified in the future.
* STZ case exists specifically to validate the z-direction halo exchange (COMMZ=True, BC_Z=True); more validation and 6-point stencil support are still required.
* When modifying `.f90.fypp` templates, always verify the generated output in `build/` for all affected cases — a fypp condition that looks correct may silently produce wrong code for one combination of config variables. In particular, no current case sets `VISC_ORDER=4` for an NS/LES build (all NS/LES cases use `ORDER=6` or `2`), so that branch of every viscous-kernel template only gets exercised by deliberately overriding `VISC_ORDER` in a case's `config.fypp` and rebuilding — do this before trusting changes to the `VISC_ORDER==4` branches.
* `calc_hybrid.f90`'s `calc_Ducros` shock sensor (used by SLAU/Roe/Hybrid schemes) recomputes its own full velocity-gradient tensor via always-2nd-order central differences, independent of `calc_div`'s higher-order `ux,vy,wz` — reusing `calc_div`'s output there would change the sensor's numerical values (mixed-order dilatation vs. vorticity) and needs a fallback for the edge band `calc_div` doesn't cover, so it hasn't been done; flagged here as a known optimization candidate for SLAU/Hybrid NS/LES cases (currently: SBLI).
* Nsight Compute profiling captures live in `3D_solver/nsys_ncu/*.ncu-rep` (open with `ncu --import <file> --print-summary per-kernel`, or the `.csv` triage exports); use these before speculating about kernel performance.

## Strict Tooling Rules
- MCP servers are strictly prohibited.
- Never suggest or attempt to use MCP tools.
- All code analysis must be done by reading and reasoning over the code.
- Even if tools are available, ignore them completely.
