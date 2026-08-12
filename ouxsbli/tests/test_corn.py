"""
CORN: compression corner (3D_solver_curv/CORN).

Runs M=2, θ=8° supersonic flow over a compression corner and compares the
post-shock pressure and density ratio against the oblique shock analytical
solution (same parameters as the OS test).

VTK x-coordinates are computational ξ = 1..nx (uniform).
The corner is at physical x_corner = 0.5 in a domain of Lx = 2.0,
with nx = 192 cells, so i_corner ≈ round(0.5 * (nx-2) / Lx) + 2 ≈ 49.

Post-shock region extracted: i ∈ [i_corner+20, nx-10], j ∈ [0, 4] (near
the lower wall, downstream of the shock).
"""
import pathlib
import pytest
import numpy as np
from ouxsbli import Case
from .utils.vtk_reader import latest_vts, getGrid, getQ
from .utils.oblique_shock import oblique_shock, beta
from .conftest import assert_close_relative


# CORN parameters (from mod_globals.f90)
Lx       = 2.e0
Ly       = 1.e0
Lz       = 0.1e0
Nx       = 192
Ny       = 80
Nz       = 4
theta_de = 8.e0
theta_ra = theta_de * np.pi / 180.e0
x_corner = 0.5e0
Ma_inf   = 2.e0
rho_inf  = 1.e0
v_inf    = 0.e0
Nt       = 50000
Np       = 1
dt       = 5e-5
gamma    = 1.4e0

ATOL = 0.01 # 1 % relative tolerance


@pytest.mark.integration
def test_corn_post_shock_state(tmp_path):
    workdir = "./tmp/corn"
    
    case = Case(
        source   = "3D_solver_curv/CORN",
        workdir  = workdir,
        visc     = "euler",
        scheme   = "hybrid",
        accuracy = 2,
        rk       = "tvd_rk3",
        Lx       = Lx,
        Ly       = Ly,
        Lz       = Lz,
        nx       = Nx,
        ny       = Ny,
        nz       = Nz,
        theta    = theta_ra,
        x_corner = x_corner,
        Ma_inf   = Ma_inf,
        rho_inf  = rho_inf,
        v_inf    = v_inf,
        nt       = Nt,
        np       = Np,
        dt       = dt,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    vts_path = latest_vts(data_dir)
    ni, nj, nk, x, y, z = getGrid(vts_path)
    rho, u, v, w, p = getQ(vts_path, ni, nj, nk)
    # Analytical reference
    p_inf = 1.e0 / gamma
    T_inf = p_inf / rho_inf
    rho_ref, p_ref, _ = oblique_shock(Ma_inf, p_inf, T_inf, np.radians(beta(Ma_inf, theta_de)), gamma, 1.e0)
    # Locate the post-shock region in computational index space.
    # i_corner: 1-based interior index of the compression corner.
    # Post-shock region: i ∈ [i_corner+20, ni-10], j ∈ [0, nj//10].
    # VTK data: rho/p shape (nk, nj, ni)
    i_corner = int(x_corner * (ni - 2) / Lx) + 1
    i_lo     = i_corner + 20
    i_hi     = ni - 10
    j_hi     = nj // 10

    rho_post = np.mean(rho[1, :j_hi, i_lo:i_hi])
    p_post   = np.mean(p  [1, :j_hi, i_lo:i_hi])

    assert_close_relative(rho_post, rho_ref, ATOL, f"rho_post (vts={vts_path})")
    assert_close_relative(p_post,   p_ref,   ATOL, "p_post")
