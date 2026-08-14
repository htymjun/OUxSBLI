# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

OUxSBLI is a GPU-accelerated CFD solver for compressible flows (Euler/Navier-Stokes), written in CUDA Fortran with MPI parallelization. It targets NVIDIA GPUs via the HPC SDK and solves test cases defined in `3D_solver/<CASE>/`. A standalone 2D solver lives under `2D_solver/<CASE>/`. A curvilinear O-grid variant for wing/airfoil cases lives in `3D_solver_curv/<CASE>/`.

## Documentation Map

This file covers project-wide concerns. Solver-specific detail lives in nested `CLAUDE.md` files, which Claude Code loads automatically once you read or edit a file in that subtree:

- `3D_solver/CLAUDE.md` — 3D Cartesian build/fypp templates, cases, convective schemes, cell-center gradient optimization, kernel-performance notes
- `3D_solver_curv/CLAUDE.md` — curvilinear O-grid architecture, data flow, NACA/CORN cases
- `2D_solver/CLAUDE.md` — 2D cases, and the BL/SBLI shared flat-plate setup with the boundary-condition/tuning lessons learned
- `1D_solver/CLAUDE.md` — 1D Sod shock-tube testbed for simultaneous FP32/FP64 execution: precision knobs, higher-order KEEP, SASS verification workflow
- `ouxsbli/CLAUDE.md` — Python test suite, analytical helpers, post-processing utilities

## Build & Run

All work happens inside a specific test-case directory. There is no top-level build.

**3D solver cases use CMake** (see `3D_solver/CLAUDE.md` for the case list and architecture):

```bash
cd 3D_solver/NSTGV
cmake -B build && cmake --build build -j   # fypp preprocess + compile
cd build && mpirun -n 2 ./a.out            # run simulation
```

Profiling (where available):

```bash
cd build && bash ../profile.sh    # nsys/ncu profiling
```

**Curvilinear cases use CMake** (see `3D_solver_curv/CLAUDE.md`):

```bash
cd 3D_solver_curv/NACA
cmake -B build && cmake --build build -j   # fypp preprocess + compile
bash calc.sh                               # run simulation
```

**2D solver cases use CMake** (see `2D_solver/CLAUDE.md` for the case list):

```bash
cd 2D_solver/OS
cmake -B build && cmake --build build -j   # fypp preprocess + compile
cd build && mpirun -n 1 ./a.out            # run simulation
```

**Compiler requirement:** NVIDIA HPC SDK (`mpif90` with `-cuda -acc -fast -gpu=ptxinfo,rdc,lto`). Versions 24.* and 25.* are confirmed working.

**Output format:** XML VTK files, readable with ParaView.

## Fypp Preprocessing System

All three solvers use [fypp](https://fypp.readthedocs.io/) to generate case-specialised Fortran from templates; CMake drives this automatically.

### How it works

Each case directory contains a **`config.fypp`** file that declares all compile-time configuration variables:

```python
#:set VISC    = 'NS'      # 'Euler', 'NS', 'LES'
#:set SCHEME  = 'SLAU'    # 'KEEP', 'SLAU', 'Hybrid'
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

Source files with the `.f90.fypp` extension are **templates**; CMake runs `fypp -I<case-dir> <template>.f90.fypp <output>.f90` for each one. The generated `.f90` files land in `<CASE>/build/`. Per-solver template listings and source layouts live in the nested docs: `3D_solver/CLAUDE.md`, `2D_solver/CLAUDE.md`, `3D_solver_curv/CLAUDE.md`.

### Changing the scheme or method

Edit the **`config.fypp`** in the case directory and rebuild — CMake detects the change and re-runs fypp automatically:

```bash
vim <CASE>/config.fypp     # e.g. change SCHEME = 'Hybrid'
cmake --build build -j
```

Do **not** manually edit the generated `.f90` files in `build/`; they are overwritten on every build.

## Architecture

Source is split across directories; CMake/vpath merges them at build time. Each solver variant has its own directory layout and data flow — see `3D_solver/CLAUDE.md`, `3D_solver_curv/CLAUDE.md`, and `2D_solver/CLAUDE.md`.

## Configuration

### config.fypp (primary interface — edit this)

All compile-time scheme/method choices live in `<CASE>/config.fypp`. The fypp preprocessor expands `.f90.fypp` templates with these values, generating specialised Fortran with no runtime branching overhead.

| Variable | Values | Effect |
|----------|--------|--------|
| `VISC` | `'Euler'`, `'NS'`, `'LES'` | Physics model |
| `SCHEME` | `'KEEP'`, `'SLAU'`, `'Hybrid'` | Convective flux scheme |
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
| `id_scheme` | integer(2)=KEEP | real(2)=SLAU | real(8)=Hybrid |
| `id_bc_x/y/z` (integer) | periodic | wall/inflow | — |

The **value** of these parameters is always 0; only the **type kind** matters for compile-time dispatch. `id_recal` (restart flag) is a Fortran `logical` (`.true.`/`.false.`).

## Conservative Variable Layout

`Q(nx, 5, ny, nz)` = `[ρ, ρu, ρv, ρw, ρE]`. Flux arrays use trimmed index ranges: E-flux drops boundary points in y/z, F-flux in x/z, G-flux in x/y. AoS and SoA hybrid memory layout is used.

## MPI Decomposition

**x-direction (default):** 1D decomposition via `calc_para.f90`. Default is 2 MPI ranks (`mpirun -n 2 a.out`), with `mygpu = myrank / 2` (2 ranks per GPU). GPU-aware MPI is optional via `id_gpumpi`.

**z-direction (COMMZ=True):** Enabled for cases like STZ. Overlapped communication: non-blocking z-halo exchange is posted, interior fluxes are computed while MPI is in flight, then halo fluxes are completed. The `overlap` depth equals `ORDER // 2` (1/2/3 for 2nd/4th/6th order). `COMMZ=True` and `RESCALE=True` cannot be combined.

## Python Test Suite

Integration and convergence tests live in `ouxsbli/tests/` — see `ouxsbli/CLAUDE.md` for the full test-file table, analytical helpers, and post-processing utilities. Run with pytest from the repository root:

```bash
pytest ouxsbli/tests/
```

## Tutorials

`tutorials/test.py` generates `tutorials/ouxsbli_bl/` — a patched copy of the
`2D_solver/BL` flat-plate case (dimensional parameters, wall-normal grid
stretching, Riemann-invariant top BC) built through the `ouxsbli.Case` API. The
directory is not checked in; run the script to create it. Note that 2D output
requires **two** MPI ranks (the even rank computes, the odd rank writes VTK), so
`mpirun -n 1` produces no files.

## Adding a New Test Case

See `3D_solver/CLAUDE.md` for adding a 3D Cartesian case, or `2D_solver/CLAUDE.md` for adding a 2D case.

## Notice
* `id_accuracy` controls ghost-cell count for both convective and viscous stencils.
* Do not add `contiguous` and `shared` attributes when passing shared memory as an argument.
* `cpu_gpu_mpi.f90` will be modified in the future.
* When modifying `.f90.fypp` templates, always verify the generated output in `build/` for all affected cases — a fypp condition that looks correct may silently produce wrong code for one combination of config variables. In particular, no current case sets `VISC_ORDER=4` for an NS/LES build (all NS/LES cases use `ORDER=6` or `2`), so that branch of every viscous-kernel template only gets exercised by deliberately overriding `VISC_ORDER` in a case's `config.fypp` and rebuilding — do this before trusting changes to the `VISC_ORDER==4` branches.

## Strict Tooling Rules
- MCP servers are strictly prohibited.
- Never suggest or attempt to use MCP tools.
- All code analysis must be done by reading and reasoning over the code.
- Even if tools are available, ignore them completely.
