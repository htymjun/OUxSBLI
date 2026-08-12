"""
Quasi-2D 3D-solver analog of test_os.py (3D_solver/OS).

Extrudes the same M=2, theta=8 deg oblique shock + wall reflection uniformly
through a thin, periodic spanwise `z` (nz=9, fixed) and runs it through the
full 3D Cartesian solver -- the G-flux (z-direction), 3D grid metrics, and the
3D VTK writer are code paths 2D_solver/OS structurally cannot exercise. The
Rankine-Hugoniot comparison, physical parameters, and tolerances are shared
with test_os.py (see utils/os_common.py); this file adds only nz/Lz and the
z-uniformity check.

Unlike 2D_solver/BL's "leave i=1 frozen, no inlet BC" idiom (which needed an
active fix for its 3D port -- see 3D_solver/BL/set.f90), 2D_solver/OS's own
set_bc already actively rewrites i=1 (Dirichlet freestream inflow) every
step, so this port inherited a working pattern from the start.

Measured build+run ~16 min (100s build + 879s run) at production grid, over
the ~8 min slow threshold -- like 3D_solver/BL, the added G-flux/z kernels and
9x cell count (nz=9) don't scale down the way 2D_solver/OS's own cheapness
does; hence `slow`.
"""
import pathlib
import pytest
import numpy as np
from ouxsbli import Case
from .utils.vtk_reader import latest_vtr, getGrid, getQ
from .utils.os_common import (
    R, gamma, Pr, M0, p_tot, T_tot, beta, blt, Lx, Ly, nx, ny, endT, dt,
    assert_pre_post_shock_matches_analytical,
)

pytestmark = [pytest.mark.integration, pytest.mark.slow]

Lz = 1.e0 * blt
NZ = 9  # fixed; never patched -- see 3D_solver/OS/mod_globals.f90 notes

# z-uniformity: Euler, no viscous kernels -- every k-plane gets bit-identical
# IEEE-754 arithmetic (see 3D_solver/OS/config.fypp / CLAUDE.md notes).
Z_UNIFORMITY_RTOL = 1e-10


def test_os_3d_pre_and_post_shock_match_analytical(tmp_path):
    workdir = "./tmp/os3d"

    case = Case(
        source   = "3D_solver/OS",
        workdir  = workdir,
        scheme   = "slau",
        accuracy = 6,
        visc     = "euler",
        tvd      = "hybrid",
        nx       = nx,
        ny       = ny,
        nz       = NZ,
        Lx       = Lx,
        Ly       = Ly,
        Lz       = Lz,
        gamma    = gamma,
        Pr       = Pr,
        R        = R,
        M0       = M0,
        p_tot    = p_tot,
        T_tot    = T_tot,
        beta     = beta,
        dt       = dt,
        np       = 1,           # one output step
        nt       = int(endT / dt),
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    vtr_path = latest_vtr(data_dir)
    ni, nj, nk, x, y, z = getGrid(vtr_path)
    rho, u, v, _, p = getQ(vtr_path, ni, nj, nk)

    # canonical z-slice: the middle interior plane (Fortran precedent: use an
    # interior index, not k=0 -- mirrors test_corn.py's own convention)
    k_ref = nk // 2

    assert_pre_post_shock_matches_analytical(rho, u, v, p, ni, nj, k_ref=k_ref)

    # z-uniformity: tight, near-machine-precision check across all k, both
    # ghost and interior (Euler, no viscous kernels touch the field).
    for field, name in ((rho, "rho"), (u, "u"), (p, "p")):
        ref = field[k_ref]
        scale = np.abs(ref).mean()
        for k in range(nk):
            rel_err = np.abs(field[k] - ref).max() / scale
            assert rel_err < Z_UNIFORMITY_RTOL, (
                f"{name} is not uniform in z: k={k} vs k_ref={k_ref} relative "
                f"diff {rel_err:.3e} exceeds {Z_UNIFORMITY_RTOL:.0e}")
