# Quick Start

Get OUxSBLI running in under 10 minutes. This guide assumes you already have the NVIDIA HPC SDK installed. If not, start with [Installation](installation.md).

## Prerequisites

| Requirement | Minimum version |
|-------------|-----------------|
| NVIDIA HPC SDK (`mpif90`) | 25.5 or newer |
| CMake | 3.18 |
| fypp | any (`pip install fypp`) |
| GPU | NVIDIA Ampere or newer |
| MPI ranks available | 2 |

## 1. Clone the repository

```bash
git clone https://github.com/htymjun/OUxSBLI.git
cd OUxSBLI
```

## 2. Install fypp

CMake uses fypp to preprocess the Fortran templates at build time:

```bash
pip install fypp
```

## 3. Build the NS Taylor-Green vortex case

`3D_solver/NSTGV` is the recommended starting point — it is a well-validated 3D case with no special boundary conditions.

```bash
cd 3D_solver/NSTGV
cmake -B build && cmake --build build -j
```

CMake automatically runs fypp on all `.f90.fypp` templates and compiles the result. The executable `a.out` is placed in `build/`.

## 4. Run the simulation

```bash
cd build
mpirun -n 2 ./a.out
```

VTK output files (`Q00000.vtr`, `Q00001.vtr`, …) appear in `data/` inside the build directory. The default configuration runs for several hundred timesteps; you can reduce `nt` in `mod_globals.f90` for a quick smoke-test.

## 5. Visualise in ParaView

1. Open ParaView and choose **File → Open**.
2. Navigate to `3D_solver/NSTGV/build/data/` and select `Q00000.vtr` (ParaView will offer to open the whole time series).
3. Click **Apply**, then choose a field (e.g. density or velocity magnitude) from the variable dropdown.

## 6. Next steps

| Goal | Where to look |
|------|---------------|
| Change the convective scheme or spatial order | [Configuration](configuration.md) |
| Run a different test case | [Installation → Available Cases](installation.md#available-cases) |
| Use the Python API for parametric runs | [Python API](api.md) |
| Understand the numerical methods | [Theory](theory.md) |
| Set up on a fresh machine | [Installation](installation.md) |
