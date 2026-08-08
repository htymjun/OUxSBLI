# Python API Reference

The `ouxsbli` Python package provides a high-level interface for parametric simulation runs: copy a case, patch `config.fypp` and `mod_globals.f90`, build with CMake, and launch with MPI — all from Python.

## Installation

```bash
pip install -e ".[dev]"
```

This installs `ouxsbli` in editable mode along with `numpy`, `scipy`, `vtk`, and `pytest`.

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
| `scheme` | `SCHEME` | str | `'KEEP'`, `'SLAU'`, `'Roe'`, `'Hybrid'` |
| `accuracy` | `ORDER` | int | `2`, `4`, `6` |
| `visc_order` | `VISC_ORDER` | int | `2`, `4`, `6` |
| `tvd` | `TVD` | str | `'none'`, `'tvd'`, `'hybdir'` |
| `slau` | `SLAU_VARIANT` | str | `'SLAU'`, `'HRSLAU2'` |
| `rescale` | `RESCALE` | bool | `True`, `False` |
| `recal` | `RESTART` | bool | `True`, `False` |
| `rk` | `RK` | int | `3`, `4` |
| `bc_x` | `BC_X` | bool | `True`, `False` |
| `bc_y` | `BC_Y` | bool | `True`, `False` |
| `commz` | `COMMZ` | bool | `True`, `False` |
| `gpumpi` | `GPUMPI` | bool | `True`, `False` |
| `output_precision` | `OUTPUT_PRECISION` | int | `4`, `8` |

String values for `scheme` and `visc` are case-insensitive (`'slau'` → `'SLAU'`, `'ns'` → `'NS'`).

> **Note:** `bc_z` has no entry in the alias table above but still works — any keyword not recognised as an alias is matched case-insensitively against the raw `config.fypp` macro name, so `bc_z=True` reaches `BC_Z` the same way an arbitrary `mod_globals.f90` parameter would.

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

### Initial & Boundary Conditions (`ic=`, `bc=`)

`Case` accepts `ic=` and `bc=` keyword arguments that expand into flat `mod_globals.f90`/`config.fypp` parameters, so `set.f90` picks them up via `select case(ic_type)` at runtime:

```python
from ouxsbli import Case, UniformIC, TaylorGreenIC, RiemannIC, PeriodicBC

case = Case(
    source  = "3D_solver/ETGV",
    workdir = "/tmp/tgv_run",
    ic = TaylorGreenIC(Ma=0.1, rho0=1.0, p0=0.71429),
    bc = PeriodicBC(),
)
```

| Class | Purpose | Key fields |
|-------|---------|-----------|
| `UniformIC` | Spatially uniform state | `rho`, `u`, `v`, `w`, `p` |
| `TaylorGreenIC` | 3D Taylor-Green vortex (incompressible-limit) | `Ma`, `rho0`, `p0`; domain should be `[0, 2π]³` |
| `RiemannIC` | 1D Riemann / Sod shock-tube split at `x = x0 * Lx` | `rho_l`, `u_l`, `p_l`, `rho_r`, `u_r`, `p_r`, `x0` |
| `PeriodicBC` | Fully periodic BCs (`BC_X=BC_Y=BC_Z=False`) | — |

`ic=`/`bc=` also accept the dict shorthand (`ic={"type": "riemann", "rho_l": 1.0, ...}`) or, for `bc=`, the string `"periodic"`.

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
| `test_bl.py` | 2D laminar BL — Cf and u-profile vs. Blasius (~4 min) |
| `test_sbli.py` | 2D shock/BL interaction — Cp and separation bubble vs. Moro et al. reference data (~8 min, `slow`) |
| `test_evc_3d.py` | Quasi-2D 3D_solver/EVC — grid convergence, plus z-uniformity (~16.5 min total, `slow`) |
| `test_bl_3d.py` | Quasi-2D 3D_solver/BL — Cf and u-profile vs. Blasius, plus z-uniformity (~13 min, `slow`) |
| `test_os_3d.py` | Quasi-2D 3D_solver/OS — pre/post state vs. Rankine-Hugoniot, plus z-uniformity (~16 min, `slow`) |
| `test_patcher.py` | Unit tests for `patcher.patch()` (no GPU required) |

Integration tests are marked with `@pytest.mark.integration` and are skipped automatically if no GPU or HPC SDK is present. `test_sbli.py` and the three `_3d.py` tests also carry the `slow` marker; deselect both kinds with `pytest -m "not slow"`.

### Analytical utilities

| Module | Contents |
|--------|----------|
| `ouxsbli/tests/utils/oblique_shock.py` | `beta_from_theta()`, `post_shock_state()` via bisection on the θ-β-M relation |
| `ouxsbli/tests/utils/sod_exact.py` | Exact Riemann solver for the Sod shock tube |
| `ouxsbli/tests/utils/vtk_reader.py` | VTK output reader |
| `ouxsbli/tests/utils/bl_common.py` | Shared Blasius tolerances and `test_*` assertions for `test_bl.py`/`test_bl_3d.py` |
| `ouxsbli/tests/utils/evc_common.py` | Shared grid list, physical parameters, `l2_rho_error()`, `assert_grid_convergence()` for `test_evc.py`/`test_evc_3d.py` |
| `ouxsbli/tests/utils/os_common.py` | Shared physical parameters and `assert_pre_post_shock_matches_analytical()` for `test_os.py`/`test_os_3d.py` |

### Post-processing utilities

Shared by the tests and the user-facing plotting scripts, in `ouxsbli/analysis/`:

| Module | Contents |
|--------|----------|
| `wall.py` | Sutherland viscosity, `compute_cf_cp()` / `wall_coeffs_from_vtr()`, `edge_state()`, `find_zero_crossings()`, `load_reference()` |
| `blasius.py` | RK4-integrated Blasius similarity solution (`fprime()`, `eta()`, `cf_blasius()`) |
| `plot_style.py` | Shared matplotlib style settings used by the reference-overlay plotting scripts |
| `dhit_decay_report.py` | Generates the DHIT kinetic-energy decay report |

---
