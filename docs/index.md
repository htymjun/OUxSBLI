[![DOI](https://zenodo.org/badge/761651757.svg)](https://doi.org/10.5281/zenodo.19396475)
![CUDA Fortran](https://img.shields.io/badge/CUDA_Fortran-GPU_Accelerated-76B900)
![Modern Fortran](https://img.shields.io/badge/Modern_Fortran-yes-success)
[![Python](https://img.shields.io/badge/Python-3.10+-3776ab?style=flat&logo=python&logoColor=white)](https://python.org)

# OUxSBLI

**OUxSBLI** is a GPU-accelerated CFD solver for compressible flows written in CUDA Fortran, with a Python API for parametric simulation runs.

It solves the compressible Euler / Navier-Stokes / LES equations on structured grids using explicit high-order finite-difference schemes, and runs on NVIDIA GPUs via the HPC SDK with MPI parallelisation.

---

## Quick Start

```bash
# Install fypp (required at build time)
pip install fypp

# Build the NS Taylor-Green vortex case
cd 3D_solver/NSTGV
cmake -B build && cmake --build build -j

# Run with 2 MPI ranks
cd build && mpirun -n 2 ./a.out
```

Or use the Python API:

```python
from ouxsbli import Case

case = Case(
    source   = "3D_solver/NSTGV",
    workdir  = "/tmp/my_run",
    visc     = "NS",
    scheme   = "SLAU",
    accuracy = 6,
    nx = 128, ny = 128, nz = 128,
)
case.build()
case.run(nranks=2)
```

---

## Key Features

| Feature | Description |
|---------|-------------|
| **KEEP scheme** | Kinetic energy and entropy preserving — ideal for smooth vortical flows |
| **SLAU / HRSLAU2** | Low-dissipation AUSM — robust for compressible turbulence with shocks |
| **Hybrid scheme** | KEEP ↔ SLAU blending via Ducros dilatation sensor |
| **High-order stencils** | 2nd, 4th, and 6th-order spatial discretisation |
| **TVD reconstruction** | tvd and hybrid limiters |
| **TVD-RK3 / RK4** | 3-stage and 4-stage time integration |
| **MPI decomposition** | z-direction with overlapped communication |
| **Curvilinear solver** | O-grid for wing/airfoil cases (NACA 0012, compression corner) |
| **Python API** | `Case` class for parametric runs and convergence studies |

---

## Solvers and Cases

### 3D Cartesian (`3D_solver/`)

| Case | Physics | Description |
|------|---------|-------------|
| BL   | NS | Quasi-2D laminar flat-plate boundary layer (extrudes 2D_solver/BL; validated vs. Blasius) |
| DHIT | NS | Decaying homogeneous isotropic turbulence |
| ETGV | Euler | Entropy-preserving Taylor-Green vortex |
| EVC  | Euler | Quasi-2D Euler vortex convection (extrudes 2D_solver/EVC; grid-convergence study) |
| IVST | Euler | Inviscid vortex smooth test case |
| KHI  | Euler | Kelvin-Helmholtz instability |
| NSTGV | NS | NS Taylor-Green vortex (Re=1600, M=1.25) |
| OS   | Euler | Quasi-2D oblique shock + wall reflection (extrudes 2D_solver/OS; validated vs. Rankine-Hugoniot) |
| SBLI | NS | Shock-boundary layer interaction |
| STZ  | NS | z-direction MPI halo exchange validation |
| TBL  | LES | Turbulent boundary layer |

### 2D Solver (`2D_solver/`)

| Case | Physics | Description |
|------|---------|-------------|
| BL   | NS | Supersonic laminar boundary layer |
| DSL  | Euler | Double shear layer |
| EVC  | Euler | Euler vortex convection |
| OS   | Euler | Oblique shock (M=2, θ=8°) |
| SBLI | NS | 2D shock-boundary layer interaction |
| ST   | Euler | Sod shock tube |

### Curvilinear (`3D_solver_curv/`)

| Case | Physics | Description |
|------|---------|-------------|
| NACA | NS | NACA 0012 airfoil on an O-grid |
| CORN | NS | Compression corner (M=2, θ=8°) |

---

## Validations

### Supersonic Taylor-Green Vortex

Results are consistent with Lusher & Sandham (2021). Setup: Re=1600, M=1.25, 512³ grid.

<div align="center">
  <img src="img/Ek.png" alt="TGV kinetic energy" width="420">
  <img src="img/enstrophy.png" alt="TGV total enstrophy" width="420">
</div>

### Shock-Boundary Layer Interaction (SBLI)

<div align="center">
  <img src="img/sbli_2d.png" alt="SBLI" width="800">
</div>

### NACA 0012 Airfoil

<div align="center">
  <img src="img/naca_p.gif" alt="NACA 0012 pressure" width="800">
</div>

---

## Related Publication

Jun Hatayama, Kento Tanaka, and Toshinori Kouchi. "Nonlinear causal relationship between separation bubbles and reflected shock wave in shock wave/turbulent boundary layer interaction based on information theory." *Computers & Fluids* (2026): 107016.

```bibtex
@article{hatayama2026nonlinear,
  title={Nonlinear causal relationship between separation bubbles and reflected shock wave
         in shock wave/turbulent boundary layer interaction based on information theory},
  author={Hatayama, Jun and Tanaka, Kento and Kouchi, Toshinori},
  journal={Computers \& Fluids},
  pages={107016},
  year={2026},
  publisher={Elsevier}
}
```

---

## License

BSD 3-Clause License.
