# Python API Reference

The `ouxsbli` Python package provides a high-level interface for parametric simulation runs: copy a case, patch `config.fypp` and `mod_globals.f90`, build with CMake, and launch with MPI — all from Python.

## Installation

```bash
pip install -e ".[dev]"
```

This installs `ouxsbli` in editable mode along with `numpy`, `vtk`, and `pytest`.

---

## Case Class

```python
from ouxsbli import Case
```

`Case` manages the full lifecycle of one simulation run.

### Constructor

```python
Case(source: str, workdir: str, **params) -> Case
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `source` | str | Path to the original case directory (e.g. `"3D_solver/NSTGV"`). Never modified. |
| `workdir` | str | Path for the new working directory to be created. Existing directory is removed and recreated on `build()`. |
| `**params` | keyword args | Parameter overrides — see tables below. |

Paths are resolved relative to the repository root (the directory containing `ouxsbli/`).

### config.fypp parameters

These are patched into the copied `config.fypp`:

| Keyword | fypp variable | Type | Example values |
|---------|---------------|------|----------------|
| `visc` | `VISC` | str | `'Euler'`, `'NS'`, `'LES'` |
| `scheme` | `SCHEME` | str | `'KEEP'`, `'SLAU'`, `'Hybrid'` |
| `accuracy` | `ORDER` | int | `2`, `4`, `6` |
| `visc_order` | `VISC_ORDER` | int | `2`, `4` |
| `tvd` | `TVD` | str | `'none'`, `'tvd'`, `'hybdir'` |
| `slau` | `SLAU_VARIANT` | str | `'SLAU'`, `'HRSLAU2'` |
| `rescale` | `RESCALE` | bool | `True`, `False` |
| `recal` | `RESTART` | bool | `True`, `False` |
| `rk` | `RK` | int | `3`, `4` |
| `bc_x` | `BC_X` | bool | `True`, `False` |
| `bc_y` | `BC_Y` | bool | `True`, `False` |
| `commz` | `COMMZ` | bool | `True`, `False` |

String values for `scheme` and `visc` are case-insensitive (`'slau'` → `'SLAU'`, `'ns'` → `'NS'`).

### mod_globals.f90 parameters

Any keyword not found in `config.fypp` is patched into `mod_globals.f90` by matching a Fortran `parameter` declaration:

| Keyword | Description |
|---------|-------------|
| `nx`, `ny`, `nz` | Grid size |
| `Re` | Reynolds number |
| `M0` | Mach number |
| `CFL` | CFL number |
| `np` | Output interval (steps between VTK writes) |
| `nt` | Number of output intervals |
| Any other `parameter` name in `mod_globals.f90` | Matched by name, case-insensitive |

---

### Methods

#### `case.build()`

Copies `source` to `workdir`, patches `config.fypp` and `mod_globals.f90`, then builds with CMake.

```python
case.build()
```

Equivalent shell steps:
```bash
cp -r source workdir
# patch config.fypp and mod_globals.f90
cmake <source> -DCASE_DIR=<workdir>   # from workdir/build/
make
```

Raises `RuntimeError` if CMake or make fails (stdout/stderr included in the message).

#### `case.run(nranks=2)`

Launches the simulation with `mpirun`.

```python
case.run(nranks=2)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nranks` | int | `2` | Number of MPI ranks |

`build()` must be called first. Raises `RuntimeError` if `build()` has not been called or if `mpirun` fails.

VTK output files appear in `workdir/data/`.

---

## Full Example

```python
import pathlib
from ouxsbli import Case

case = Case(
    source  = "3D_solver/NSTGV",
    workdir = "/tmp/nstgv_slau6",
    # config.fypp settings
    visc     = "NS",
    scheme   = "SLAU",
    accuracy = 6,
    tvd      = "hybrid",
    # mod_globals.f90 settings
    nx = 65,
    ny = 65,
    nz = 65,
    Re = 1600.0,
    np = 50,
)

case.build()
case.run(nranks=2)

print("Output written to:", pathlib.Path("/tmp/nstgv_slau6") / "data")
```

---

## Reading Output

The `ouxsbli.tests.utils` module provides a lightweight VTK reader:

```python
from ouxsbli.tests.utils.vtk_reader import read_vtk

data = read_vtk("/tmp/nstgv_slau6/data/Q00010.vtr")
# data is a dict: {"rho": array, "rhou": array, "rhov": array, "rhow": array, "rhoE": array}
```

For full post-processing, open the `.vtr` files in ParaView as a time series.

---

## Test Suite

Integration and convergence tests live in `ouxsbli/tests/`. Run with pytest from the repository root:

```bash
pytest ouxsbli/tests/ -v
```

| Test file | What it checks |
|-----------|----------------|
| `test_etgv.py` | Supersonic Taylor-Green vortex — kinetic energy decay |
| `test_st.py` | Sod shock tube — exact Riemann solution (density, pressure, velocity) |
| `test_evc.py` | Euler vortex convergence — KEEP 2nd/4th/6th order |
| `test_os.py` | 2D oblique shock — pre/post state vs. Rankine-Hugoniot (tol 2% / 5%) |
| `test_corn.py` | 3D_solver_curv/CORN — pressure and density ratios vs. θ-β-M theory (tol 5%) |
| `test_patcher.py` | Unit tests for `patcher.patch()` (no GPU required) |

Integration tests are marked with `@pytest.mark.integration` and are skipped automatically if no GPU or HPC SDK is present.

### Analytical utilities

| Module | Contents |
|--------|----------|
| `ouxsbli/tests/utils/oblique_shock.py` | `beta_from_theta()`, `post_shock_state()` via bisection on the θ-β-M relation |
| `ouxsbli/tests/utils/sod_exact.py` | Exact Riemann solver for the Sod shock tube |
| `ouxsbli/tests/utils/vtk_reader.py` | VTK output reader |

---
