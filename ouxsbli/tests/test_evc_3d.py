"""
Quasi-2D 3D-solver analog of test_evc.py (3D_solver/EVC).

Extrudes the same Euler vortex uniformly through a thin, periodic spanwise
`z` (nz=9, fixed) and runs it through the full 3D Cartesian solver -- the
G-flux (z-direction), 3D grid metrics, and the 3D VTK writer are code paths
`2D_solver/EVC` structurally cannot exercise. The convergence machinery and
physical parameters are shared with test_evc.py (see utils/evc_common.py);
this file adds only nz/Lz and the z-uniformity check, which
`assert_grid_convergence` folds in automatically once `l2_rho_error` reports
a real (non-None) z_err.

nz=9 is unchanged across every accuracy tested (4 and 6): verified from
set_grid_cyclic{4,6}_3D's ghost widths, nz=9 gives 3 real z-planes at
accuracy=6 (k=4,5,6) and 5 at accuracy=4 (k=3..7).

nz=9 rather than the theoretical 6th-order-stencil minimum of 7: direct
testing showed nz=7 (1 real interior plane) trips a GPU-kernel edge case that
corrupts the x/y boundary region within a single RK step, even though the
field is perfectly z-uniform (verified: nz=9 runs clean, nz=13 runs clean,
nz=7 does not). See 3D_solver/EVC/mod_globals.f90.

Measured build+run across all 6 (scheme, accuracy, nx) combinations: 992.68s
(~16.5 min total pytest wall-clock), over the ~8 min slow threshold -- same
per-case cause as 3D_solver/BL and 3D_solver/OS (G-flux/z kernels, 9x cell
count from nz=9); hence `slow`.
"""
import pathlib
import os
import pytest
from ouxsbli import Case
from .utils.vtk_reader import initial_vtr, latest_vtr
from .utils.evc_common import (
    SCHEME_CONFIGS,
    Lx, Ly, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0,
    l2_rho_error,
    assert_grid_convergence,
)

pytestmark = pytest.mark.slow

NZ = 9  # fixed; never patched -- see 3D_solver/EVC/mod_globals.f90 notes
Lz = 0.01e0  # inert: VISC='Euler', no viscous kernels, no w-velocity

params = (Lx, Ly, Lz, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0)
SCHEMES = [(scheme, accuracy, min_order, params) for scheme, accuracy, min_order in SCHEME_CONFIGS]


def _run_evc_3d(scheme, accuracy, nx, params):
    """Build and run one quasi-2D 3D EVC configuration; return (L2(rho) error, z_err)."""
    workdir = os.path.join("./tmp", f"evc3d_{scheme}{accuracy}_{nx}")

    Lx, Ly, Lz, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0 = params

    n = nx + accuracy

    case = Case(
        source   = "3D_solver/EVC",
        workdir  = workdir,
        scheme   = scheme,
        accuracy = accuracy,
        visc     = "euler",
        tvd      = "none",
        nx       = n,
        ny       = n,
        nz       = NZ,
        Lx       = Lx,
        Ly       = Ly,
        Lz       = Lz,
        gamma    = gamma,
        Pr       = Pr,
        R        = R,
        M0       = M0,
        beta     = beta,
        theta    = theta,
        Rc       = Rc,
        p0       = p0,
        T0       = T0,
        u0       = u0,
        output_precision=8,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    return l2_rho_error(initial_vtr(data_dir), latest_vtr(data_dir), accuracy)


@pytest.mark.integration
@pytest.mark.parametrize("scheme,accuracy,min_order,params", SCHEMES)
def test_evc_3d_grid_convergence(scheme, accuracy, min_order, params):
    assert_grid_convergence(_run_evc_3d, scheme, accuracy, min_order, params)
