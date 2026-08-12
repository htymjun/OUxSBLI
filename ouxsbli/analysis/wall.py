"""Wall-quantity post-processing shared by the pytest suite and the
user-facing plotting scripts (2D_solver/BL/Cf.py, 2D_solver/SBLI/ref/*.py).

Everything here is numpy-only; the VTK reader is imported lazily inside
``wall_coeffs_from_vtr`` so the pure-math helpers stay importable without the
``vtk`` package.
"""

import numpy as np

# Sutherland's law constants — must match src/mod_constant.f90.fypp
# (mu0_T0_S_over_T0_2_3 = 1.716d-5 * 384.2d0 * 273.2d0**-1.5).
SUTHERLAND_MU_REF = 1.716e-5  # [Pa s] reference viscosity at T_ref
SUTHERLAND_T_REF = 273.2      # [K]
SUTHERLAND_S = 111.0          # [K]
SUTHERLAND_C = SUTHERLAND_MU_REF * (SUTHERLAND_T_REF + SUTHERLAND_S) / SUTHERLAND_T_REF**1.5


def sutherland_mu(T):
    """Dynamic viscosity mu(T) [Pa s] from Sutherland's law."""
    return SUTHERLAND_C * T**1.5 / (T + SUTHERLAND_S)


def one_sided_deriv_3pt(y0, y1, y2, f0, f1, f2):
    """Derivative at y0 of the quadratic through (y0,f0), (y1,f1), (y2,f2).

    Valid for non-uniform grid spacing (one-sided 3-point finite difference).
    """
    denom0 = (y0 - y1) * (y0 - y2)
    denom1 = (y1 - y0) * (y1 - y2)
    denom2 = (y2 - y0) * (y2 - y1)
    return (
        f0 * ((y0 - y1) + (y0 - y2)) / denom0
        + f1 * (y0 - y2) / denom1
        + f2 * (y0 - y1) / denom2
    )


def wall_dudy(y, u2d):
    """du/dy at the wall (y[0]) for every x column of a (ny, nx) u field.

    Assumes no-slip: the wall value is taken as u = 0, combined with the first
    two off-wall points in a non-uniform one-sided 3-point stencil.
    """
    u_wall = np.zeros_like(u2d[0, :])
    return one_sided_deriv_3pt(y[0], y[1], y[2], u_wall, u2d[1, :], u2d[2, :])


def compute_cf_cp(y, rho2d, u2d, p2d, R, p_inf, q_inf):
    """Bottom-wall skin-friction and pressure coefficients along x.

    Parameters are the (ny, nx) field planes plus the freestream reference
    state. Wall temperature is derived from the wall pressure and density
    (T_w = p_w / (rho_w R)), viscosity from Sutherland's law.

    Returns (cf, cp), each a 1-D array over x.
    """
    rho_w = rho2d[0, :]
    p_w = p2d[0, :]
    T_w = p_w / (rho_w * R)
    mu_w = sutherland_mu(T_w)
    tau_w = mu_w * wall_dudy(y, u2d)
    cf = tau_w / q_inf
    cp = (p_w - p_inf) / q_inf
    return cf, cp


def edge_state(rho2d, u2d, p2d, R):
    """Boundary-layer edge state for every x column: (u_e, rho_e, mu_e).

    The edge is taken at the column maximum of u, which is unambiguous only
    for an attached boundary layer with a monotonic profile — do not use this
    on a separated or shock-laden field.

    Using the local edge state rather than the nominal freestream absorbs the
    mild displacement-driven acceleration that a finite-height domain imposes
    on a growing boundary layer.
    """
    nx = u2d.shape[1]
    j_edge = np.argmax(u2d, axis=0)
    cols = np.arange(nx)
    u_e = u2d[j_edge, cols]
    rho_e = rho2d[j_edge, cols]
    p_e = p2d[j_edge, cols]
    mu_e = sutherland_mu(p_e / (rho_e * R))
    return u_e, rho_e, mu_e


def wall_coeffs_from_vtr(path, R, p_inf, q_inf):
    """Read a VTR snapshot and return (x, cf, cp) along the bottom wall."""
    from ouxsbli.tests.utils.vtk_reader import getGrid, getQ

    ni, nj, nk, x, y, z = getGrid(str(path))
    rho, u, v, w, p = getQ(str(path), ni, nj, nk)
    cf, cp = compute_cf_cp(y, rho[0, :, :], u[0, :, :], p[0, :, :], R, p_inf, q_inf)
    return x, cf, cp


def find_zero_crossings(x, f):
    """Locate f=0 crossings (e.g. Cf separation/reattachment) by linear interpolation."""
    sign_change = np.where(np.diff(np.sign(f)) != 0)[0]
    crossings = []
    for i in sign_change:
        x1, x2 = x[i], x[i + 1]
        f1, f2 = f[i], f[i + 1]
        crossings.append(x1 - f1 * (x2 - x1) / (f2 - f1))
    return crossings


def load_reference(path, x_scale=1.0):
    """Load a 2-column digitized reference dataset (whitespace-separated).

    Digitizer points may arrive out of order, so rows are sorted by the first
    column; ``x_scale`` converts the abscissa (e.g. x/L -> X/Xsh is 1.25 for
    the SBLI datasets digitized from Vila-Perez Fig. 28).

    Returns an (N, 2) array.
    """
    ref = np.loadtxt(path)
    ref = ref[np.argsort(ref[:, 0])]
    ref[:, 0] = ref[:, 0] * x_scale
    return ref
