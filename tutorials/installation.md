# Installation & Setup Guide

## System Requirements

| Component | Requirement |
|-----------|-------------|
| OS | Linux (tested on Ubuntu 22.04) |
| GPU | NVIDIA GPU (Ampere or newer recommended) |
| CUDA | 12.x |
| Compiler | NVIDIA HPC SDK 24.x or 25.x (`mpif90`) |
| MPI | OpenMPI or MPICH (bundled with HPC SDK) |
| Python | 3.10 or newer |

---

## 1. Install NVIDIA HPC SDK

The solver requires `mpif90` from the NVIDIA HPC SDK.

```bash
# Download from https://developer.nvidia.com/hpc-sdk
```

Verify:

```bash
mpif90 --version
mpirun --version
```

---

## 2. Clone the Repository

```bash
git clone <repository-url>
cd OUxSBLI
```

---

## 3. Install the Python Package

Install the `ouxsbli` Python package and its dependencies in one step:

```bash
pip install -e ".[dev]"
```

This installs:
- `ouxsbli` (the Python API, editable install — changes to source take effect immediately)
- `numpy` — array operations
- `pytest` — for running the test suite (dev extra)

> **Tip:** Use a virtual environment to keep dependencies isolated:
> ```bash
> python -m venv .venv
> source .venv/bin/activate
> pip install -e ".[dev]"
> ```

---

## 4. Verify the Installation

```bash
# Run the unit tests (no GPU or build environment required)
python -m pytest ouxsbli/tests/ -v
```

Expected output: all patcher tests pass; I/O tests skip if no simulation data is present yet.

---

## 5. Build a Test Case (Requires GPU)

Each case directory contains a `Makefile`. To build the 2D laminar boundary layer case:

```bash
cd 2D_solver/BL
make clean && make
```

Successful compilation produces `a.out`. Then run:

```bash
mkdir -p data
mpirun -n 2 ./a.out
```

VTK output files (`Q00000.vtr`, `Q00001.vtr`, …) appear in `data/` every `np` timesteps.

---

## Available Cases

### 2D Solver (`2D_solver/`)

| Directory | Physics | Description |
|-----------|---------|-------------|
| `BL/` | Navier-Stokes | Supersonic laminar boundary layer (M=2) |
| `DSL/` | Euler | Double shear layer |
| `EVC/` | Euler | Euler vortex convection |
| `ST/` | Euler | 2D shock tube |

### 3D Solver (`3D_solver/`)

| Directory | Physics | Description |
|-----------|---------|-------------|
| `ETGV/` | Euler | Euler Taylor-Green vortex |
| `NSTGV/` | Navier-Stokes | NS Taylor-Green vortex (Re=1600) |
| `KHI/` | Euler | Kelvin-Helmholtz instability |
| `SBLI/` | Navier-Stokes | Shock-boundary layer interaction |
| `TBL/` | Navier-Stokes | Turbulent boundary layer |
| `DHIT/` | Navier-Stokes | Decaying homogeneous isotropic turbulence |

### Curvilinear Solver (`3D_solver_curv/`)

| Directory | Physics | Description |
|-----------|---------|-------------|
| `NACA/` | Navier-Stokes | NACA 0012 airfoil (O-grid) |
| `CORN/` | Navier-Stokes | Corner geometry |

---

## Scheme / Physics Selection Reference

Parameters in `mod_globals.f90` use the Fortran type **kind** (not value) as a dispatch flag.

| Parameter | kind=2 | kind=4 | kind=8 |
|-----------|--------|--------|--------|
| `id_visc` | Euler | Navier-Stokes | LES |
| `id_accuracy` | 2nd order | 4th order | 6th order |
| `id_tvd` | no limiter | minmod | no limiter/minmod |
| `id_RungeKutta` | TVD-RK3 | RK4 | — |
| `id_recal` | initialize | restart from file | — |

| `id_scheme` | Type | Scheme |
|-------------|------|--------|
| `integer(2)` | KEEP | Kinetic Energy & Entropy Preserving |
| `real(2)` | SLAU | Simple Low-dissipation AUSM |
| `real(8)` | Hybrid | KEEP ↔ SLAU via Ducros sensor |

The Python API handles this automatically via friendly strings:

```python
Case(..., scheme="slau", visc="ns", accuracy=2)
```
