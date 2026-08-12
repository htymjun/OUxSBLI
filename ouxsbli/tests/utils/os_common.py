"""
Shared oblique-shock (Rankine-Hugoniot) comparison and physical parameters
for test_os.py (2D_solver/OS) and test_os_3d.py (3D_solver/OS) -- same case,
same parameters, same tolerances; the 3D file only adds nz/Lz and a
z-uniformity check.

The case initialises a diagonal oblique shock (M=2, theta=8 deg) that
reflects off the bottom wall; both tests average the pre-shock and
post-reflected-shock regions and compare against the analytical
Rankine-Hugoniot state.
"""
import numpy as np

from .oblique_shock import reslected_shock, free_stream
from ..conftest import assert_close_relative

R = 287.03e0
gamma = 1.4e0
Pr = 0.72e0
# free stream
M0 = 2.e0
p_tot = 100.e3
T_tot = 295.e0
# oblique shock
beta = np.pi * 37.2e0 / 180.e0
beta_r = np.pi * 44.1e0 / 180.e0
# space and time
blt = 1.e-3
Lx = 5.e0 * blt
Ly = 2.e0 * blt
nx = 257
ny = 129
endT = 0.1e-3
dt = 3.e-9

PRE_RTOL = 0.03   # 3 % for undisturbed pre-shock
POST_RTOL = 0.03  # 3 % for post-shock (SLAU has some numerical diffusion)


def assert_pre_post_shock_matches_analytical(rho, u, v, p, ni, nj, k_ref=0):
    """Compare the pre-/post-(reflected-)shock mean state at z-plane k_ref
    against Rankine-Hugoniot. rho/u/v/p are the (nk,nj,ni) fields from getQ()."""
    rho3, ux3, uy3, p3 = reslected_shock(M0, gamma, R, p_tot, T_tot, beta, beta_r)

    # free-stream (pre-shock) analytical state
    u0_fs, p0_fs, T0_fs = free_stream(M0, gamma, R, p_tot, T_tot)
    rho0_fs = p0_fs / (R * T0_fs)

    # Pre-shock region: left 10 % of domain (i < ni//10), top half (j >= nj//2).
    # Purely upstream of the incident shock for beta~37 deg in a 5x2 mm domain.
    rho_pre = rho[k_ref, nj // 2 :, : ni // 10].astype(float)
    u_pre   = u[k_ref,   nj // 2 :, : ni // 10].astype(float)
    p_pre   = p[k_ref,   nj // 2 :, : ni // 10].astype(float)

    # Post-reflected-shock region: right 20 % (i >= 4*ni//5), bottom quarter (j < nj//4).
    # At x=4Lx/5=4 mm the reflected shock (from x_hit~2.6 mm) is at y~1.1 mm,
    # so j < nj//4 is well below the reflected shock.
    rho_post = rho[k_ref, : nj // 3, 9 * ni // 10 :].astype(float)
    u_post   = u[k_ref,   : nj // 3, 9 * ni // 10 :].astype(float)
    v_post   = v[k_ref,   : nj // 3, 9 * ni // 10 :].astype(float)
    p_post   = p[k_ref,   : nj // 3, 9 * ni // 10 :].astype(float)

    # pre-shock assertions (undisturbed free stream)
    assert_close_relative(rho_pre.mean(), rho0_fs, PRE_RTOL, "Pre-shock rho")
    assert_close_relative(u_pre.mean(),   u0_fs,   PRE_RTOL, "Pre-shock u")
    assert_close_relative(p_pre.mean(),   p0_fs,   PRE_RTOL, "Pre-shock p")

    # post-shock assertions (after reflected shock)
    assert_close_relative(rho_post.mean(), rho3, POST_RTOL, "Post-shock rho")
    assert_close_relative(u_post.mean(),   ux3,  POST_RTOL, "Post-shock u")
    # After a perfect wall reflection the flow is horizontal; check v~=0 via u scale.
    assert abs(v_post.mean()) / abs(ux3) < POST_RTOL, (
        f"Post-shock v: {v_post.mean():.4f} (expected ~ 0)")
    assert_close_relative(p_post.mean(),   p3,   POST_RTOL, "Post-shock p")
