"""
Quasi-2D 3D-solver analog of test_bl.py (3D_solver/BL).

Extrudes the same laminar flat-plate boundary layer uniformly through a thin,
periodic spanwise `z` (nz=9, fixed) and runs it through the full 3D Cartesian
solver -- the G-flux (z-direction), `calc_div_wz`, 3D grid metrics, and the 3D
VTK writer are code paths 2D_solver/BL structurally cannot exercise. The
Blasius comparison assertions are shared with test_bl.py (see
utils/bl_common.py); this file adds only the fixture (3D build/run) and the
z-uniformity check that has no 2D analogue.

Unlike 2D_solver/BL's "leave i=1 frozen, no inlet BC" idiom, 3D_solver/BL uses
an active Dirichlet inlet (see 3D_solver/BL/set.f90): direct testing showed
the untouched-i=1 idiom corrupts the domain within a single RK step in the 3D
boundary-aware x kernel, reproduced independent of scheme/viscosity/nz, and
never previously exercised since SBLI/TBL (the only other BC_X=True 3D cases)
both actively rewrite their leftmost columns via rescaling.

Production defaults (3D_solver/BL/mod_globals.f90: nx=257, ny=49, nz=9,
Ly=15mm, dt=4e-8, endT=1e-2 -- same as 2D_solver/BL plus the z extrusion) are
themselves what this test runs unmodified apart from output cadence/precision,
mirroring test_bl.py's own convention. Measured build+run ~13 min (build 94s
+ run 685s), over the ~8 min slow threshold -- unlike 2D (whose per-step cost
is launch-overhead dominated), 3D adds whole new kernel launches (G-flux,
calc_div_wz) with no 2D analogue, so the per-step cost does not scale down the
way 2D's cost-reduction pass did; hence `slow`.
"""
import pathlib

import numpy as np
import pytest

from ouxsbli import Case
from .utils.vtk_reader import extract_number, getGrid, getQ
from .utils.bl_common import (  # noqa: F401 -- imported for pytest collection
    test_grid_places_leading_edge_at_origin,
    test_flow_is_quasi_steady,
    test_freestream_is_not_accelerated,
    test_cf_matches_blasius,
    test_velocity_profile_matches_blasius,
)

pytestmark = [pytest.mark.integration, pytest.mark.slow]

# z-uniformity: measured exactly bit-identical (0.0) across every k-plane,
# including in u (the viscous, spanwise-gradient-bearing field) -- ux=du/dx
# and vy=dv/dy are by definition x/y-direction derivatives, so whichever
# formula (calc_div's own "true interior" band vs calc_visc_high's fallback)
# computes them at a given k, it is applied to the SAME (k-invariant) x/y data
# at every k, and produces the SAME result. The "different formula at
# different k" nuance flagged for calc_div (CLAUDE.md) never actually
# produces a numerical difference here, so a single tight tolerance suffices.
Z_UNIFORMITY_RTOL = 1e-10


@pytest.fixture(scope="module")
def bl_run():
    """Build and run the 3D BL case once, returning the last two snapshots and k_ref."""
    workdir = "./tmp/bl3d"

    case = Case(
        source="3D_solver/BL",
        workdir=workdir,
        np=10,
        output_precision=8,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    snaps = sorted(data_dir.glob("Q*.vtr"), key=lambda p: extract_number(p.name))
    assert len(snaps) >= 2, f"expected at least 2 snapshots in {data_dir}, got {len(snaps)}"
    _, _, nk, _, _, _ = getGrid(str(snaps[-1]))
    return snaps[-2], snaps[-1], nk // 2


def test_field_is_uniform_in_z(bl_run):
    """The raw field must be genuinely uniform in z at every k-plane, not just at k_ref."""
    _, last, k_ref = bl_run
    ni, nj, nk, x, y, z = getGrid(str(last))
    rho, u, v, w, p = getQ(str(last), ni, nj, nk)

    assert np.abs(w).max() < 1e-10, f"spanwise velocity should be exactly 0, got max|w|={np.abs(w).max():.3e}"

    for field, name in ((rho, "rho"), (u, "u"), (p, "p")):
        ref = field[k_ref]
        scale = np.abs(ref).mean()
        for k in range(nk):
            rel_err = np.abs(field[k] - ref).max() / scale
            assert rel_err < Z_UNIFORMITY_RTOL, (
                f"{name} is not uniform in z: k={k} vs k_ref={k_ref} relative "
                f"diff {rel_err:.3e} exceeds {Z_UNIFORMITY_RTOL:.0e}")
