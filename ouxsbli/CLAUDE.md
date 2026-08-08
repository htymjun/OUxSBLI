# CLAUDE.md — ouxsbli/

This file documents the Python test suite and post-processing helpers. See the repository root `CLAUDE.md` for project overview, shared configuration reference, and tooling rules.

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
| `test_bl.py` | 2D laminar BL — Cf and u-profile vs. Blasius (~4 min) |
| `test_sbli.py` | 2D shock/BL interaction — Cp and separation bubble vs. Moro et al. (~8 min, `slow`) |
| `test_evc_3d.py` | Quasi-2D 3D_solver/EVC — grid convergence, same as `test_evc.py`, plus z-uniformity (~16 min total, `slow`) |
| `test_bl_3d.py` | Quasi-2D 3D_solver/BL — Cf and u-profile vs. Blasius, same as `test_bl.py`, plus z-uniformity (~13 min, `slow`) |
| `test_os_3d.py` | Quasi-2D 3D_solver/OS — pre/post state vs. Rankine-Hugoniot, same as `test_os.py`, plus z-uniformity (~16 min, `slow`) |

`test_sbli.py` and all three `_3d.py` tests carry the `slow` marker; deselect
with `pytest -m "not slow"`. `test_bl.py`/`test_sbli.py` assert a quasi-steady
guard between the last two snapshots before comparing, so a failure to
converge reports itself instead of surfacing as a tolerance miss.
BL's production defaults (`2D_solver/BL/mod_globals.f90`: nx=257, ny=49, Ly=15mm,
dt=4e-8, endT=1e-2 — see the cost note in `2D_solver/CLAUDE.md`) are themselves
cheap (~2.5 min), so `test_bl.py` runs the case unmodified apart from output
cadence/precision.
SBLI still patches its case down (coarsens to 276×257, raises dt 10x) since its
production settings remain expensive. `Case()` can only patch `config.fypp` and
`mod_globals.f90` — anything living inside `set.f90`, e.g. BL's tanh stretch
parameter `s`, is not tunable from a test.

`test_bl.py`/`test_bl_3d.py` share their comparison logic entirely —
`utils/bl_common.py` holds the tolerances, the `wall_profile()` helper, and
the actual `test_*` functions (grid placement, quasi-steady guard, Blasius
Cf/profile match, freestream-not-accelerated guard); both files just `from
.utils.bl_common import test_flow_is_quasi_steady, ...` and pytest collects
the imported names under each file. This works because pytest resolves a
test's fixtures by looking at the *collecting* module (test_bl.py or
test_bl_3d.py), not the module the function was defined in — so each file
only needs to supply its own `bl_run` fixture (build+run, returning
`(prev_snapshot, last_snapshot, k_ref)`; 2D's `k_ref` is always 0) and the
shared bodies resolve against it correctly. `test_bl_3d.py` adds only
`test_field_is_uniform_in_z`, which has no 2D analogue.
Cf and the profile are normalised by the *local* boundary-layer edge state,
which divides out any residual outer-flow acceleration. That is a
convenience, not a licence: it is exactly what hid the Riemann top-BC error
described in `2D_solver/CLAUDE.md`, so `test_freestream_is_not_accelerated`
compares `max(u_e)` against the nominal u0 directly (tolerance 1%).

The three `_3d.py` tests are quasi-2D ports (see `3D_solver/CLAUDE.md` for the
`nz=9`/active-inlet design notes) — same tolerances as their 2D counterparts,
plus an explicit z-uniformity assertion. All three measured over the `slow`
threshold: `test_bl_3d.py` ~13 min (94s build + 685s run), `test_os_3d.py`
~16 min (100s build + 879s run), `test_evc_3d.py` ~16.5 min total across its
6 (scheme, accuracy, nx) combinations. 3D adds whole new kernel launches
(G-flux, `calc_div_wz`) with no 2D analogue, plus the 9x cell count from
`nz=9`, so the per-step cost does not scale down the way the 2D cases' own
cheapness does.

`test_evc.py`/`test_evc_3d.py` and `test_os.py`/`test_os_3d.py` share the same
way, but via helper functions rather than imported `test_*` names — each has
only one (parametrized, for EVC) test per file, so there's no fixture reused
across multiple tests the way `bl_run` is. `utils/evc_common.py` holds
`GRIDS`, the physical parameters, `l2_rho_error()` (ghost-trims x/y and, when
the grid has real z extent, z too — driven purely by `Nz` read back from the
VTK file, not a passed-in flag), and `assert_grid_convergence()`, which both
files' `test_evc*_grid_convergence` call with their own per-grid build/run
callback. `utils/os_common.py` holds the physical parameters and
`assert_pre_post_shock_matches_analytical(rho, u, v, p, ni, nj, k_ref)`; both
files build+run inline (there's only one case config, no grid sweep) and call
it with `k_ref=0` (2D) or the middle z-plane (3D).

## Analytical Helpers

In `ouxsbli/tests/utils/`:
- `oblique_shock.py` — `beta_from_theta()`, `post_shock_state()` via bisection on the θ-β-M relation
- `sod_exact.py` — exact Riemann solver for the Sod shock tube
- `vtk_reader.py` — VTK output reader
- `bl_common.py` — shared Blasius-comparison tolerances and `test_*` assertions for `test_bl.py`/`test_bl_3d.py` (see the Python Test Suite section above for how the sharing works)
- `evc_common.py` — shared grid list, physical parameters, `l2_rho_error()`, and `assert_grid_convergence()` for `test_evc.py`/`test_evc_3d.py`
- `os_common.py` — shared physical parameters and `assert_pre_post_shock_matches_analytical()` for `test_os.py`/`test_os_3d.py`

## Post-Processing

Shared by the tests and the user-facing plotting scripts, in `ouxsbli/analysis/`:
- `wall.py` — Sutherland viscosity, one-sided wall derivative, `compute_cf_cp()` / `wall_coeffs_from_vtr()`, `edge_state()`, `find_zero_crossings()`, `load_reference()`
- `blasius.py` — RK4-integrated Blasius similarity solution (`fprime()`, `eta()`, `cf_blasius()`); accurate to ~1e-5, unlike the coarse `fp_tab` in `src/set_compressible_bl.f90`, which is an IC seed only and deviates by up to 7%
