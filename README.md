[![DOI](https://zenodo.org/badge/761651757.svg)](https://doi.org/10.5281/zenodo.19396475)
![CUDA Fortran](https://img.shields.io/badge/CUDA_Fortran-GPU_Accelerated-76B900)
![Modern Fortran](https://img.shields.io/badge/Modern_Fortran-yes-success)
[![Python](https://img.shields.io/badge/Python-3.10+-3776ab?style=flat&logo=python&logoColor=white)](https://python.org)
[![Docs](https://img.shields.io/badge/docs-online-blue)](https://htymjun.github.io/OUxSBLI/)

<div align="center">
  <img src="./docs/img/OUxSBLI.png" alt="OUxSBLI">  
</div>

OUxSBLI is a GPU-accelerated CFD code with Python-API written in CUDA Fortran. It employs explicit high-order finite-difference schemes on a rectilinear grid (3D and 2D solvers) and a curvilinear grid.

## Documentation

| Guide | Description |
|-------|-------------|
| [Quick Start](docs/quickstart.md) | Run your first simulation in minutes |
| [Installation](docs/installation.md) | Full setup guide for a fresh machine |
| [Configuration](docs/configuration.md) | config.fypp reference and rebuild workflow |
| [Theory](docs/theory.md) | Numerical methods and governing equations |
| [Python API](docs/api.md) | `ouxsbli` package reference |

## Dependency

### CUDA Fortran
* NVIDIA HPC SDK 24.x or 25.x — provides `mpif90` and bundled MPI
* CMake 3.18 or newer
* [fypp](https://fypp.readthedocs.io/) (`pip install fypp`) — Fortran preprocessor
* ParaView — for visualizing VTK output

### Python API
* Python 3.10 or newer
* numpy
* vtk

## Usage

### CUDA Fortran

1. Go to a case directory (e.g. NS Taylor-Green vortex):

```bash
cd 3D_solver/NSTGV
```

2. Edit **`config.fypp`** to choose the physics model, convective scheme, spatial order, boundary conditions, and other compile-time options:

```python
#:set VISC   = 'NS'      # 'Euler', 'NS', or 'LES'
#:set SCHEME = 'SLAU'    # 'KEEP', 'SLAU', or 'Hybrid'
#:set ORDER  = 6         # spatial order: 2, 4, or 6
#:set BC_X   = False     # False → periodic; True → wall/inflow BCs
#:set BC_Y   = False
#:set BC_Z   = False
```

See [docs/configuration.md](docs/configuration.md) for the full reference.

3. Edit **`mod_globals.f90`** to set grid size, domain lengths, physical parameters, and GPU thread-block sizes:

```fortran
integer, parameter :: nx = 513
integer, parameter :: ny = 513
integer, parameter :: nz = 513
real(8), parameter :: Re = 1600.d0
```

4. Edit **`set.f90`** if you need to change the grid geometry, initial conditions, or boundary condition routines.

5. Build with CMake:

```bash
cmake -B build && cmake --build build -j
```

6. Run the simulation:

```bash
cd build && mpirun -n 2 ./a.out
```

VTK output files (`Q00000.vtr`, `Q00001.vtr`, …) appear in the `data/` directory.

**Profiling:** In some case directories, `profile.sh` runs nsys/ncu profiling:

```bash
cd build && bash ../profile.sh
```

### Python API

1. Install the `ouxsbli` package:

```bash
pip install -e ".[dev]"
```

2. Create a `Case`, build, and run:

```python
import pathlib
from ouxsbli import Case

case = Case(
    source  = "3D_solver/NSTGV",
    workdir = "/tmp/my_run",
    # physics (maps to config.fypp)
    visc    = "NS",
    scheme  = "SLAU",
    accuracy = 6,
    # grid (maps to mod_globals.f90)
    nx = 128,
    ny = 128,
    nz = 128,
)

case.build()
case.run(nranks=2)  # mpirun -n 2 ./a.out
```

See [docs/api.md](docs/api.md) for the full API reference.

## Discretization

### Spatial (Convection terms)
* Kinetic energy and entropy preserving (KEEP) scheme
* Simple low-dissipation AUSM (SLAU) scheme
* KEEP / SLAU hybrid scheme

### Spatial (Viscous terms)
* ME4-Base
* Gaitonde and Visbal's 2nd-order scheme

### Spatial SGS
* Selective mixed scale model

### Temporal
* 3-stage TVD Runge-Kutta
* 4-stage classical Runge-Kutta

## Validations and visualizations

### Supersonic Taylor-Green vortex
The results are consistent with Lusher's results.
* Numerical setup

|$Re$   |$1600$ |
| :---: | :---: |
|$Ma$   |$1.25$ |
|$N_x \times N_y \times N_z$|$512\times512\times512$|

~~~bash
@article{lusher2021assessment,
  title={Assessment of low-dissipative shock-capturing schemes for the compressible Taylor--Green vortex},
  author={Lusher, David J and Sandham, Neil D},
  journal={AIAA Journal},
  volume={59},
  number={2},
  pages={533--545},
  year={2021},
  publisher={American Institute of Aeronautics and Astronautics}
}
~~~

<div align="center">
  <img src="./docs/img/Ek.png" alt="TGV_kinetic_energy" width="450">  
</div>

<div align="center">
  <img src="./docs/img/enstrophy.png" alt="TGV_total_enstrophy" width="450">  
</div>

### Shock Boundary Layer Interaction (SBLI)

<div align="center">
  <img src="./docs/img/sbli_2d.png" alt="SBLI" width="900">  
</div>

### 2D Oblique Shock

M=2 freestream with θ=8° flow deflection. Pre- and post-shock states agree with the Rankine-Hugoniot relations within 2% and 5% respectively, verified by `ouxsbli/tests/test_os.py`.

### NACA0012

<div align="center">
  <img src="./docs/img/naca_p.gif" alt="NACA" width="900">
</div>

## Related Publication
This repository contains the implementation used in the following publication:

Jun Hatayama, Kento Tanaka, and Toshinori Kouchi. "Nonlinear causal relationship between separation bubbles and reflected shock wave in shock wave/turbulent boundary layer interaction based on information theory." Computers & Fluids (2026): 107016.

~~~bash
@article{hatayama2026nonlinear,
  title={Nonlinear causal relationship between separation bubbles and reflected shock wave in shock wave/turbulent boundary layer interaction based on information theory},
  author={Hatayama, Jun and Tanaka, Kento and Kouchi, Toshinori},
  journal={Computers \& Fluids},
  pages={107016},
  year={2026},
  publisher={Elsevier}
}
~~~

The repository was made publicly available after publication to improve reproducibility. However, this version may differ slightly from the version used in the paper.

## AI-Assisted Development
Development during 2024 and 2025 was primarily conducted by the project owner.  
Starting in 2026, the project expanded its contributor base and introduced AI-assisted "vibe coding" workflows using Claude Code.

To maintain transparency, we aim to clearly distinguish which parts of the codebase and development workflow involve AI-generated content or AI-assisted modifications. In addition, as part of our effort to share practical knowledge on AI-assisted development in the HPC community, we provide Claude Code plan files under `./docs/plans`.

## License
This project is under BSD 3-Clause License
