"""
EVC grid-convergence test (2D_solver/EVC).

The Euler vortex convects exactly one full period (t = T = Lx/u0) and
should return to its initial state.  We measure the L2 error of density
at the end time for three grid sizes and check that the convergence rate
matches the expected order of the scheme.

Schemes tested
--------------
KEEP4 (KEEP, 4th-order accuracy) → expected order ≥ 3.6
KEEP6 (KEEP, 6th-order accuracy) → expected order ≥ 5.7

Grid levels: nx = ny ∈ [64, 128, 256]

Note:
Requests output_precision=8 (double-precision VTK output) via Case(...) --
single precision (the default) doesn't have enough dynamic range to resolve
a 4th/6th-order convergence trend.
"""
import math
import pathlib
import numpy as np
import pytest
import os
from ouxsbli import Case
from .utils.vtk_reader import initial_vtr, latest_vtr, getGrid, getScalar, getVector


GRIDS = [64, 128, 256]   # three refinement levels

# mesh
Lx    = 0.1e0
Ly    = 0.1e0
# physical properties
gamma = 1.4e0
Pr    = 0.71e0
R     = 287.15e0
# initial condition
M0    = 0.05e0
beta  = 1.e0 / 50.e0
theta = 0.e0
Rc    = 0.005e0
p0    = 1.e5
T0    = 300.e0
u0    = M0 * np.sqrt(gamma * R * T0)
params = (Lx, Ly, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0)


# (scheme_key, accuracy_order, min_convergence_order)
SCHEMES = [
    ("keep", 4, 3.6, params),
    ("keep", 6, 5.7, params),
]


def _run_evc(scheme, accuracy, nx, params):
    """Build and run one EVC configuration; return L2(rho) error."""
    workdir = os.path.join("./tmp", f"evc_{scheme}{accuracy}_{nx}")

    Lx, Ly, gamma, Pr, R, M0, beta, theta, Rc, p0, T0, u0 = params
    
    n  = nx + accuracy

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

    data_dir  = pathlib.Path(workdir) / "data"
    vtr_path0 = initial_vtr(data_dir)
    vtr_path1 = latest_vtr(data_dir)
    Nx, Ny, Nz, x, y, _ = getGrid(vtr_path0)
    Rho0 = getScalar(vtr_path0, Nx, Ny, Nz, "rho")
    Rho1 = getScalar(vtr_path1, Nx, Ny, Nz, "rho")
    # remove ghost cells
    rho0 = Rho0[0,accuracy//2:-accuracy//2,accuracy//2:-accuracy//2]
    rho1 = Rho1[0,accuracy//2:-accuracy//2,accuracy//2:-accuracy//2]
    l2 = math.sqrt(np.mean((rho1 - rho0)**2))
    print("scheme", scheme, " accuracy", accuracy, " nx", nx, " L2", l2)
    return l2


def _convergence_order(errors, grids):
    """Least-squares estimate of convergence order from error vs grid size."""
    h = np.array([1.0 / n for n in grids], dtype=float)
    e = np.array(errors, dtype=float)
    # log(e) = p*log(h) + C
    log_h = np.log(h)
    log_e = np.log(np.maximum(e, 1e-300))
    p = np.polyfit(log_h, log_e, 1)[0]
    return p


@pytest.mark.integration
@pytest.mark.parametrize("scheme,accuracy,min_order,params", SCHEMES)
def test_evc_grid_convergence(scheme, accuracy, min_order, params):
    errors = []
    for nx in GRIDS:
        l2 = _run_evc(scheme, accuracy, nx, params)
        errors.append(l2)

    order = _convergence_order(errors, GRIDS)

    # Print for visibility in -v output
    for nx, e in zip(GRIDS, errors):
        print(f"  nx={nx:3d}  L2(rho)={e:.3e}")
    print(f"  estimated order = {order:.2f}  (required >= {min_order})")

    assert order >= min_order, (
        f"Convergence order {order:.2f} < required {min_order} "
        f"for scheme={scheme}, accuracy={accuracy}"
    )
