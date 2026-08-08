"""
Shared grid-convergence machinery and physical parameters for test_evc.py
(2D_solver/EVC) and test_evc_3d.py (3D_solver/EVC) -- same vortex, same
tolerances; the 3D file only adds nz/Lz and a z-uniformity check.

The Euler vortex convects and should return to its initial state; each test
measures the L2 error of density at the end time for three grid sizes and
checks that the convergence rate matches the expected order of the scheme.
Requests output_precision=8 in both files -- single precision doesn't have
enough dynamic range to resolve a 4th/6th-order convergence trend.
"""
import math

import numpy as np

from .vtk_reader import getGrid, getScalar

GRIDS = [64, 128, 256]  # three refinement levels

# mesh (x/y only -- z, where applicable, is file-specific)
Lx = 0.1e0
Ly = 0.1e0
# physical properties
gamma = 1.4e0
Pr = 0.71e0
R = 287.15e0
# initial condition
M0 = 0.05e0
beta = 1.e0 / 50.e0
theta = 0.e0
Rc = 0.005e0
p0 = 1.e5
T0 = 300.e0
u0 = M0 * np.sqrt(gamma * R * T0)

# (scheme, accuracy, min_convergence_order)
SCHEME_CONFIGS = [
    ("keep", 4, 3.6),
    ("keep", 6, 5.7),
]

# z-uniformity (3D only): Euler, no viscous kernels -- every k-plane gets
# bit-identical IEEE-754 arithmetic (see 3D_solver/EVC/config.fypp / CLAUDE.md
# notes).
Z_UNIFORMITY_RTOL = 1e-10


def l2_rho_error(vtr_path0, vtr_path1, accuracy):
    """(l2, z_err) between two snapshots, ghost-trimmed in x/y (and z, taking
    the canonical middle plane, when the grid has real z extent).

    z_err is None for a genuinely 2D grid (Nz == 1); otherwise it is the max
    relative difference of rho across all k-planes, which the caller should
    assert stays under Z_UNIFORMITY_RTOL.
    """
    Nx, Ny, Nz, x, y, z = getGrid(vtr_path0)
    Rho0 = getScalar(vtr_path0, Nx, Ny, Nz, "rho")
    Rho1 = getScalar(vtr_path1, Nx, Ny, Nz, "rho")
    g = accuracy // 2

    if Nz == 1:
        rho0 = Rho0[0, g:-g, g:-g]
        rho1 = Rho1[0, g:-g, g:-g]
        z_err = None
    else:
        k_ref = Nz // 2
        z_err = max(
            np.abs(Rho1[k] - Rho1[k_ref]).max() / np.abs(Rho1[k_ref]).mean()
            for k in range(Nz)
        )
        # only one (accuracy=6) or three (accuracy=4) z-planes remain and are
        # bit-identical by construction; take the first rather than averaging.
        rho0 = Rho0[g:Nz - g, g:-g, g:-g][0]
        rho1 = Rho1[g:Nz - g, g:-g, g:-g][0]

    l2 = math.sqrt(np.mean((rho1 - rho0) ** 2))
    return l2, z_err


def convergence_order(errors, grids=GRIDS):
    """Least-squares estimate of convergence order from error vs grid size."""
    h = np.array([1.0 / n for n in grids], dtype=float)
    e = np.array(errors, dtype=float)
    # log(e) = p*log(h) + C
    log_h = np.log(h)
    log_e = np.log(np.maximum(e, 1e-300))
    p = np.polyfit(log_h, log_e, 1)[0]
    return p


def assert_grid_convergence(run_one_grid, scheme, accuracy, min_order, params):
    """Build+run every grid in GRIDS via run_one_grid(scheme, accuracy, nx,
    params) -- which must return (l2, z_err) per `l2_rho_error` -- then assert
    the fitted convergence order meets min_order (and, when z_err is present,
    that the field stayed z-uniform)."""
    errors, z_errs = [], []
    for nx in GRIDS:
        l2, z_err = run_one_grid(scheme, accuracy, nx, params)
        errors.append(l2)
        if z_err is not None:
            z_errs.append(z_err)

    order = convergence_order(errors)

    # Print for visibility in -v output
    if len(z_errs) == len(errors):
        for nx, e, ze in zip(GRIDS, errors, z_errs):
            print(f"  nx={nx:3d}  L2(rho)={e:.3e}  z_err={ze:.3e}")
    else:
        for nx, e in zip(GRIDS, errors):
            print(f"  nx={nx:3d}  L2(rho)={e:.3e}")
    print(f"  estimated order = {order:.2f}  (required >= {min_order})")

    if z_errs:
        assert max(z_errs) < Z_UNIFORMITY_RTOL, (
            f"field is not uniform in z: max relative diff {max(z_errs):.3e} across "
            f"k-planes exceeds {Z_UNIFORMITY_RTOL:.0e} (scheme={scheme}, accuracy={accuracy}); "
            "VISC='Euler' means every k-plane should receive bit-identical arithmetic")

    assert order >= min_order, (
        f"Convergence order {order:.2f} < required {min_order} "
        f"for scheme={scheme}, accuracy={accuracy}"
    )
