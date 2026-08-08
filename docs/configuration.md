# Configuration Guide

OUxSBLI uses [fypp](https://fypp.readthedocs.io/) to specialise Fortran source at build time. Compile-time choices — physics model, convective scheme, spatial order, boundary conditions — live in `config.fypp` inside each case directory. A separate `mod_globals.f90` holds grid size, physical parameters, and GPU thread-block sizes.

## How It Works

```
config.fypp          ← edit this to change physics/scheme/BCs
     │
     ▼ cmake --build (fypp runs automatically)
src/calc_flux_base.f90.fypp ──► build/calc_flux_base.f90
src/calc_time_dev.f90.fypp  ──► build/calc_time_dev.f90
src/mod_constant.f90.fypp   ──► build/mod_constant.f90
     ...
```

CMake detects changes to `config.fypp` and re-runs fypp on every affected template. The generated `.f90` files in `build/` are overwritten on each build — **do not edit them directly**.

---

## config.fypp Reference

All compile-time options are set here. Below is a full listing with allowed values.

| Variable | Type | Values | Default (NSTGV) | Effect |
|----------|------|--------|-----------------|--------|
| `VISC` | string | `'Euler'`, `'NS'`, `'LES'` | `'NS'` | Physics model |
| `SCHEME` | string | `'KEEP'`, `'SLAU'`, `'Roe'`, `'Hybrid'` | `'SLAU'` | Convective flux scheme |
| `ORDER` | int | `2`, `4`, `6` | `6` | Spatial accuracy (convective + viscous) |
| `VISC_ORDER` | int | `2`, `4`, `6` | *(defaults to `ORDER`)* | Override viscous stencil order |
| `TVD` | string | `'none'`, `'tvd'`, `'hybrid'` | `'hybrid'` | TVD limiter for reconstruction |
| `SLAU_VARIANT` | string | `'SLAU'`, `'HRSLAU2'` | `'HRSLAU2'` | SLAU flux variant |
| `RESCALE` | bool | `True`, `False` | `False` | SBLI reference-state rescaling |
| `COMMZ` | bool | `True`, `False` | `False` | z-direction MPI halo decomposition |
| `RK` | int | `3`, `4` | `3` | Runge-Kutta stages (TVD-RK3 or RK4) |
| `RESTART` | bool | `True`, `False` | `False` | Read initial condition from checkpoint |
| `BC_X` | bool | `True`, `False` | `False` | Wall/inflow BCs in x (False → periodic) |
| `BC_Y` | bool | `True`, `False` | `False` | Wall/inflow BCs in y |
| `BC_Z` | bool | `True`, `False` | `False` | Wall/inflow BCs in z |
| `ORDER_IO` | int | computed | `ORDER // 2 - 1` | Ghost-cell count for I/O interpolation |
| `OUTPUT_PRECISION` | int | `4`, `8` | `4` | VTK output float precision (kind dispatch: real*4 vs real*8). Single precision can't resolve 6th-order convergence trends or tight Cf tolerances — the EVC/BL/OS validation tests set this to `8`. |

> **Note:** `COMMZ = True` and `RESCALE = True` cannot be combined. The Roe scheme is supported but not actively optimised — prefer KEEP, SLAU, or Hybrid.

### Example: NSTGV config.fypp

```python
#:set VISC         = 'NS'
#:set SCHEME       = 'SLAU'
#:set ORDER        = 6
#:set BC_X         = False
#:set BC_Y         = False
#:set BC_Z         = False
#:set RESCALE      = False
#:set RK           = 3
#:set COMMZ        = False
#:set RESTART      = False
#:set TVD          = 'hybrid'
#:set SLAU_VARIANT = 'HRSLAU2'
#:set ORDER_IO     = ORDER // 2 - 1
#:set OUTPUT_PRECISION = 4
```

---

## mod_globals.f90 Reference

This file is **not** preprocessed by fypp. Edit it directly; it is compiled as plain Fortran.

### Grid

| Parameter | Description |
|-----------|-------------|
| `nx`, `ny`, `nz` | Number of grid points (including ghost cells) |
| `Lx`, `Ly`, `Lz` | Domain lengths |

### Physical parameters

| Parameter | Description |
|-----------|-------------|
| `Re` | Reynolds number |
| `M0` | Mach number |
| `gamma` | Ratio of specific heats |
| `Pr` | Prandtl number |
| `R` | Gas constant (J/kg/K) |

### Time stepping

| Parameter | Description |
|-----------|-------------|
| `CFL` | CFL number (sets `dt` implicitly) |
| `dt` | Time-step size |
| `np` | Output interval (steps between VTK writes) |
| `nt` | Total number of output intervals (total steps = `np * nt`) |

### GPU thread-block sizes

| Parameter | Kernel direction | Typical value |
|-----------|-----------------|---------------|
| `threadsE` | x-direction flux | `dim3(32,1,1)` |
| `threadsF` | y-direction flux | `dim3(32,4,1)` |
| `threadsG` | z-direction flux | `dim3(32,1,4)` |
| `threadsEv` | x-direction viscous | `dim3(32,1,1)` |
| `threadsFv` | y-direction viscous | `dim3(32,4,1)` |
| `threadsGv` | z-direction viscous | `dim3(32,1,4)` |

Thread-block sizes must be tuned to the GPU architecture and grid size for best occupancy.

---

## Internal Dispatch Encoding (mod_constant.f90)

`src/mod_constant.f90.fypp` is generated per-case into `build/mod_constant.f90`. It encodes `config.fypp` choices as Fortran `parameter` constants using the **type kind** (not the value) as a dispatch flag. All values are `0`; only the kind matters at compile time.

| Parameter | kind=2 | kind=4 | kind=8 |
|-----------|--------|--------|--------|
| `id_visc` (integer) | Euler | Navier-Stokes | LES |
| `id_accuracy` (integer) | 2nd order | 4th order | 6th order |
| `id_tvd` (integer) | no TVD | tvd | hybrid |
| `id_slau` (integer) | SLAU | HRSLAU2 | — |
| `id_rescale` (integer) | off | on | — |
| `id_bc_x/y/z` (integer) | periodic | wall/inflow | — |
| `id_scheme` | integer(2) = KEEP | real(2) = SLAU | real(4) = Roe / real(8) = Hybrid |

`id_recal` (restart flag) is a Fortran `logical` (`.true.` / `.false.`).

---

## Rebuild Workflow

After editing `config.fypp`:

```bash
cmake --build build -j
```

CMake detects the change and re-runs fypp on every affected template before recompiling. A full rebuild is not needed — only the changed modules are recompiled.

After editing `mod_globals.f90`:

```bash
cmake --build build -j
```

The module is recompiled and all files that `use mod_globals` are relinked.

---

## Worked Example: Switch NSTGV from SLAU to Hybrid Scheme

The Hybrid scheme blends KEEP and SLAU using the Ducros dilatation sensor — it is well-suited to flows with mixed smooth and shocked regions.

1. Open `3D_solver/NSTGV/config.fypp`.
2. Change the `SCHEME` line:

```diff
-#:set SCHEME = 'SLAU'
+#:set SCHEME = 'Hybrid'
```

3. Rebuild:

```bash
cd 3D_solver/NSTGV
cmake --build build -j
```

4. Run:

```bash
cd build && mpirun -n 2 ./a.out
```

CMake automatically selects `calc_hybrid_kernel.f90.fypp` and `calc_hybrid_kernel_internal.f90.fypp` for compilation. No other files need to be changed.

---

## Adding a New Case

### 3D Cartesian

```bash
cp -r 3D_solver/NSTGV 3D_solver/MYCASE
```

Then edit:
1. `mod_globals.f90` — grid size, domain lengths, physical parameters, thread blocks
2. `set.f90` — grid generation, initial conditions, boundary condition calls
3. `config.fypp` — VISC, SCHEME, ORDER, BC_X/Y/Z, etc.
4. `calc.sh` — MPI rank count

Build and run:

```bash
cd 3D_solver/MYCASE
cmake -B build && cmake --build build -j
cd build && mpirun -n 2 ./a.out
```

### 2D

```bash
cp -r 2D_solver/OS 2D_solver/MYCASE
```

Same edit pattern. Run with `mpirun -n 1 ./a.out`.
