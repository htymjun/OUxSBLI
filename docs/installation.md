# Installation & Setup Guide

## System Requirements

| Component | Requirement |
|-----------|-------------|
| OS | Linux (tested on Ubuntu 22.04) |
| GPU | NVIDIA GPU (Ampere or newer recommended) |
| CUDA | 12.x |
| Compiler | NVIDIA HPC SDK 24.x or 25.x (`mpif90`) |
| MPI | OpenMPI bundled with HPC SDK |
| CMake | 3.18 or newer |
| fypp | any version (`pip install fypp`) |
| Python | 3.10 or newer (optional, for Python API and tests) |

---

## 1. Install NVIDIA HPC SDK

The solver requires `mpif90` from the NVIDIA HPC SDK. Download the installer from:

```
https://developer.nvidia.com/hpc-sdk
```

After installation, verify the tools are on your `PATH`:

```bash
mpif90 --version
mpirun --version
```

---

## 2. Install CMake

Most Linux distributions include a recent CMake. If yours does not:

```bash
# Ubuntu 22.04+
sudo apt-get install cmake

# Or download from https://cmake.org/download/
```

Verify:

```bash
cmake --version   # must be 3.18 or newer
```

---

## 3. Install fypp

CMake calls fypp to expand the `.f90.fypp` Fortran templates at build time:

```bash
pip install fypp
```

---

## 4. Clone the Repository

```bash
git clone https://github.com/htymjun/OUxSBLI.git
cd OUxSBLI
```

---

## 5. Install the Python Package (optional)

Install the `ouxsbli` Python package and test dependencies:

```bash
pip install -e ".[dev]"
```

This installs:
- `ouxsbli` — the Python API for parametric runs (editable install)
- `numpy` — array operations
- `scipy` — analytical reference solutions (Blasius, oblique-shock, Sod)
- `vtk` — VTK file reading
- `pytest` — for running the test suite

> **Tip:** Use a virtual environment to keep dependencies isolated:
> ```bash
> python -m venv .venv
> source .venv/bin/activate
> pip install -e ".[dev]"
> ```

---

## 6. Build a Test Case

Each case directory contains a `CMakeLists.txt`. To build the NS Taylor-Green vortex case:

```bash
cd 3D_solver/NSTGV
cmake -B build && cmake --build build -j
```

CMake runs fypp on all `.f90.fypp` templates, then compiles with `mpif90`. The executable `a.out` lands in `build/`.

---

## 7. Run the Simulation

```bash
cd build
mpirun -n 2 ./a.out
```

VTK output files (`Q00000.vtr`, `Q00001.vtr`, …) appear in `data/`. Open them in ParaView as a time series.

---

## Available Cases

### 2D Solver (`2D_solver/`)

| Directory | Physics | Description |
|-----------|---------|-------------|
| `BL/` | Navier-Stokes | Supersonic laminar boundary layer (M=2) |
| `DSL/` | Euler | Double shear layer |
| `EVC/` | Euler | Euler vortex convection (grid-convergence study) |
| `OS/` | Euler | 2D oblique shock (M=2, θ=8°) |
| `SBLI/` | Navier-Stokes | 2D shock-boundary layer interaction |
| `ST/` | Euler | Sod shock tube |

### 3D Cartesian Solver (`3D_solver/`)

| Directory | Physics | Description |
|-----------|---------|-------------|
| `BL/` | Navier-Stokes | Quasi-2D laminar flat-plate boundary layer (extrudes 2D_solver/BL; validated vs. Blasius) |
| `DHIT/` | Navier-Stokes | Decaying homogeneous isotropic turbulence |
| `ETGV/` | Euler | Entropy-preserving Taylor-Green vortex |
| `EVC/` | Euler | Quasi-2D Euler vortex convection (extrudes 2D_solver/EVC; grid-convergence study) |
| `IVST/` | Euler | Inviscid vortex smooth test case |
| `KHI/` | Euler | Kelvin-Helmholtz instability |
| `NSTGV/` | Navier-Stokes | NS Taylor-Green vortex (Re=1600, M=1.25) |
| `OS/` | Euler | Quasi-2D oblique shock + wall reflection (extrudes 2D_solver/OS; validated vs. Rankine-Hugoniot) |
| `SBLI/` | Navier-Stokes | Shock-boundary layer interaction (RESCALE=True) |
| `STZ/` | Navier-Stokes | z-direction MPI halo exchange validation |
| `TBL/` | LES | Turbulent boundary layer (Hybrid scheme) |

### Curvilinear Solver (`3D_solver_curv/`)

| Directory | Physics | Description |
|-----------|---------|-------------|
| `NACA/` | Navier-Stokes | NACA 0012 airfoil on an O-grid |
| `CORN/` | Navier-Stokes | Compression corner (M=2, θ=8°) |

> **Note:** The curvilinear solver also uses CMake. Build the same way: `cmake -B build && cmake --build build -j`

---

## Configuration Reference

Compile-time scheme and physics choices are set in `config.fypp`, not in `mod_globals.f90`. See [Configuration](configuration.md) for details.

At a glance, the dispatch encoding used internally by `mod_constant.f90`:

| Parameter | kind=2 | kind=4 | kind=8 |
|-----------|--------|--------|--------|
| `id_visc` | Euler | Navier-Stokes | LES |
| `id_accuracy` | 2nd order | 4th order | 6th order |
| `id_tvd` | no limiter | tvd | hybrid |

| `id_scheme` | Type | Scheme |
|-------------|------|--------|
| `integer(2)` | KEEP | Kinetic Energy & Entropy Preserving |
| `real(2)` | SLAU | Simple Low-dissipation AUSM |
| `real(8)` | Hybrid | KEEP ↔ SLAU via Ducros sensor |

These Fortran parameters are **generated automatically** from `config.fypp` — you do not edit them directly.
