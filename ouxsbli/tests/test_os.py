"""
OS: 2D oblique shock (2D_solver/OS).

The case initialises a diagonal oblique shock (M=2, theta=8 deg) across a
symmetric bottom wall; average pre/post-shock states match the analytical
Rankine-Hugoniot values. The comparison, physical parameters, and tolerances
are shared with test_os_3d.py -- see utils/os_common.py.
"""
import pathlib
import pytest
from ouxsbli import Case
from .utils.vtk_reader import latest_vtr, getGrid, getQ
from .utils.os_common import (
    R, gamma, Pr, M0, p_tot, T_tot, beta, Lx, Ly, nx, ny, endT, dt,
    assert_pre_post_shock_matches_analytical,
)


@pytest.mark.integration
def test_os_pre_and_post_shock_match_analytical(tmp_path):
    workdir = "./tmp/os"

    case = Case(
        source   = "2D_solver/OS",
        workdir  = workdir,
        scheme   = "slau",
        accuracy = 6,
        visc     = "euler",
        tvd      = "hybrid",
        nx       = nx,
        ny       = ny,
        Lx       = Lx,
        Ly       = Ly,
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

    assert_pre_post_shock_matches_analytical(rho, u, v, p, ni, nj, k_ref=0)
