"""
EVC grid-convergence test (2D_solver/EVC).

The Euler vortex convects exactly one full period (t = T = 2*Lx/u0, two
convective periods) and should return to its initial state. The convergence
machinery, physical parameters, and L2-error computation are shared with
test_evc_3d.py -- see utils/evc_common.py.

Schemes tested
--------------
KEEP4 (KEEP, 4th-order accuracy) -> expected order >= 3.6
KEEP6 (KEEP, 6th-order accuracy) -> expected order >= 5.7

Grid levels: nx = ny in [64, 128, 256]
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

params = (Lx, Ly, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0)
SCHEMES = [(scheme, accuracy, min_order, params) for scheme, accuracy, min_order in SCHEME_CONFIGS]


def _run_evc(scheme, accuracy, nx, params):
    """Build and run one EVC configuration; return (L2(rho) error, None)."""
    workdir = os.path.join("./tmp", f"evc_{scheme}{accuracy}_{nx}")

    Lx, Ly, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0 = params

    n = nx + accuracy

    case = Case(
        source   = "2D_solver/EVC",
        workdir  = workdir,
        scheme   = scheme,
        accuracy = accuracy,
        visc     = "euler",
        tvd      = "none",
        nx       = n,
        ny       = n,
        Lx       = Lx,
        Ly       = Ly,
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
        output_precision = 8,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    return l2_rho_error(initial_vtr(data_dir), latest_vtr(data_dir), accuracy)


@pytest.mark.integration
@pytest.mark.parametrize("scheme,accuracy,min_order,params", SCHEMES)
def test_evc_grid_convergence(scheme, accuracy, min_order, params):
    assert_grid_convergence(_run_evc, scheme, accuracy, min_order, params)
