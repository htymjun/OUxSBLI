"""
Shared Blasius-comparison assertions for test_bl.py (2D_solver/BL) and
test_bl_3d.py (3D_solver/BL) -- same case, same physics, same tolerances; the
3D file only adds a z-uniformity check and a slow build+run. Comparisons use
the *local* boundary-layer edge state rather than the nominal freestream,
which divides out any residual acceleration of the outer flow and is the
standard way to compare a boundary layer in a finite domain against flat-plate
theory.

That normalisation is a convenience, not a licence. An earlier version of
this module's docstring called the ~2% acceleration the 2D case then showed
"a property of the case geometry, not a solver error"; it was neither. The
top boundary was a Riemann-invariant far field, which holds the incoming
invariant fixed and so ties the boundary pressure to the wall-normal velocity
through the acoustic impedance, (p_b - p0)/q_inf = 2*v_b/(u0*M0). At M0 = 0.1
that 1/M0 turned the displacement-induced v_b into a 0.07*q_inf favorable
pressure gradient along the plate. Because edge normalisation hides exactly
that failure, test_freestream_is_not_accelerated below checks the freestream
directly.

Every test here takes a `bl_run` fixture (defined per-file, since the case
source and build cost differ) returning `(prev_snapshot, last_snapshot,
k_ref)`: the last two VTK snapshots and the z-plane index to slice at (always
0 for the 2D case; the middle interior plane for the 3D one).
"""
import numpy as np
import pytest

from ouxsbli.analysis.blasius import fprime
from ouxsbli.analysis.wall import edge_state, one_sided_deriv_3pt, sutherland_mu
from .vtk_reader import getGrid, getQ

# gas constant and freestream; must match {2D,3D}_solver/BL/mod_globals.f90
# (M=0.1 at T0=288.15 K).
R = 287.15
U0 = 0.1 * np.sqrt(1.4 * R * 288.15)

BLASIUS_CF_CONST = 0.664  # Cf = 0.664 / sqrt(Re_x)

# Cf comparison window [m]: starts far enough past the leading edge that the
# under-resolved LE singularity is out, and well short of the 0th-order outlet.
CF_X_MIN = 0.032
CF_X_MAX = 0.072
# Measured (2D, bit-identical across two independent build+run cycles): mean
# 1.212%, max 1.812%. 3D measured within noise of that (mean 1.226%, max
# 1.827%), since its canonical z-plane receives identical treatment to the
# validated 2D case. Tolerances are ~1.4-1.5x the 2D measurement -- tight
# enough to catch a real regression, with headroom for legitimate
# cross-environment drift (different GPU/HPC-SDK version) rather than
# same-run noise, which measured exactly zero here. The old 5%/8% dated from
# when the case carried the Riemann top BC and were loose enough that a -2%
# leading-edge error cancelling a +2.6% pressure-gradient error read as a pass.
CF_MEAN_TOL = 0.018
CF_MAX_TOL = 0.025

# The far field must not do work on the outer flow, checked over the same
# window as Cf. The window matters: the leading-edge singularity produces a
# local ~1.6% overspeed within a few cells of x=0 that is a grid effect,
# whereas a far-field pressure error shows up as a monotone rise peaking at
# the *downstream* end. Measured 0.408% (2D) / 0.407% (3D); the tolerance
# below still rejects the original bug (2.9%) by a wide margin while being
# ~1.5x the current measurement, not ~2.5x.
EDGE_ACCEL_TOL = 0.006

# u-profile station [m] (Re_x ~ 3.2e4) and profile tolerances.
# Measured: 22 points in the window, max|diff| 0.443%/0.447%, rms
# 0.278%/0.276%, |v|/u_e 0.479% (2D/3D).
PROFILE_X = 0.056
PROFILE_ETA_MAX = 6.0
PROFILE_MIN_POINTS = 8
PROFILE_MAX_TOL = 0.0065
PROFILE_RMS_TOL = 0.004
PROFILE_V_TOL = 0.007

# quasi-steady guard between the last two snapshots. Measured 0.0955% (both).
STEADY_TOL = 0.002


def wall_profile(path, k):
    """(x, cf_local, Re_x_local, u_e) at z-plane k using the local edge state."""
    ni, nj, nk, x, y, z = getGrid(str(path))
    rho, u, v, w, p = getQ(str(path), ni, nj, nk)
    rho2d, u2d, p2d = rho[k], u[k], p[k]

    u_e, rho_e, mu_e = edge_state(rho2d, u2d, p2d, R)

    T_w = p2d[0, :] / (rho2d[0, :] * R)
    dudy = one_sided_deriv_3pt(y[0], y[1], y[2], np.zeros(ni), u2d[1, :], u2d[2, :])
    tau_w = sutherland_mu(T_w) * dudy

    cf = tau_w / (0.5 * rho_e * u_e**2)
    Re_x = rho_e * u_e * x / mu_e
    return x, cf, Re_x, u_e


def test_grid_places_leading_edge_at_origin(bl_run):
    """set_grid must put the first no-slip point exactly at x=0 and the wall at y=0."""
    _, last, _ = bl_run
    ni, nj, nk, x, y, z = getGrid(str(last))

    assert y[0] == pytest.approx(0.0, abs=1e-12), f"wall not at y=0: y[0]={y[0]}"
    assert np.abs(x).min() < 1e-9, (
        f"leading edge is off-grid: closest x to 0 is {np.abs(x).min():.3e} m")
    assert x[-1] >= 0.079, f"domain too short for the Cf window: x_max={x[-1]:.4f} m"


def test_flow_is_quasi_steady(bl_run):
    """The window-mean Cf must stop changing between the last two outputs."""
    prev, last, k_ref = bl_run
    x_p, cf_p, _, _ = wall_profile(prev, k_ref)
    x_l, cf_l, _, _ = wall_profile(last, k_ref)

    mask = (x_l >= CF_X_MIN) & (x_l <= CF_X_MAX)
    mean_prev = cf_p[mask].mean()
    mean_last = cf_l[mask].mean()
    change = abs(mean_last - mean_prev) / abs(mean_last)

    assert change < STEADY_TOL, (
        f"boundary layer is still developing: window-mean Cf changed {change:.2%} "
        f"between the last two outputs ({mean_prev:.4e} -> {mean_last:.4e}); "
        "increase endT")


def test_freestream_is_not_accelerated(bl_run):
    """The far field must not do work on the outer flow.

    Every other assertion here normalises by the local edge state and is
    therefore blind to a far-field pressure error -- which is precisely how
    the Riemann-invariant top boundary the 2D case used to carry went
    unnoticed (see module docstring). Guard the freestream directly.
    """
    _, last, k_ref = bl_run
    _, _, _, u_e = wall_profile(last, k_ref)
    ni, nj, nk, x, y, z = getGrid(str(last))

    window = (x >= CF_X_MIN) & (x <= CF_X_MAX)
    accel = u_e[window].max() / U0 - 1.0

    assert accel < EDGE_ACCEL_TOL, (
        f"edge velocity peaks {accel:.2%} above the nominal freestream over "
        f"x = {CF_X_MIN}-{CF_X_MAX} m, more than {EDGE_ACCEL_TOL:.1%}; the top "
        "boundary condition is doing work on the outer flow and the boundary "
        "layer is no longer a zero-pressure-gradient one")


def test_cf_matches_blasius(bl_run):
    """Cf(x) must follow 0.664/sqrt(Re_x) over the comparison window."""
    _, last, k_ref = bl_run
    x, cf, Re_x, _ = wall_profile(last, k_ref)

    mask = (x >= CF_X_MIN) & (x <= CF_X_MAX)
    x_w, cf_w = x[mask], cf[mask]
    assert x_w.size >= 20, f"only {x_w.size} points in the Cf window"

    cf_theory = BLASIUS_CF_CONST / np.sqrt(Re_x[mask])
    err = np.abs(cf_w / cf_theory - 1.0)

    assert (cf_w > 0).all(), "Cf must stay positive on an attached flat plate"
    assert (np.diff(cf_w) < 0).all(), "Cf must decrease monotonically downstream"
    assert err.mean() < CF_MEAN_TOL, (
        f"mean Cf error {err.mean():.3f} exceeds {CF_MEAN_TOL} "
        f"(x = {CF_X_MIN}-{CF_X_MAX} m)")
    assert err.max() < CF_MAX_TOL, (
        f"max Cf error {err.max():.3f} exceeds {CF_MAX_TOL} "
        f"at x = {x_w[np.argmax(err)]:.4f} m")


def test_velocity_profile_matches_blasius(bl_run):
    """u/u_e(eta) at a mid-plate station must match the Blasius f'(eta)."""
    _, last, k_ref = bl_run
    ni, nj, nk, x, y, z = getGrid(str(last))
    rho, u, v, w, p = getQ(str(last), ni, nj, nk)
    rho2d, u2d, p2d = rho[k_ref], u[k_ref], p[k_ref]

    i = int(np.argmin(np.abs(x - PROFILE_X)))
    x_st = x[i]

    u_e, rho_e, mu_e = edge_state(rho2d, u2d, p2d, R)
    nu_e = mu_e[i] / rho_e[i]
    eta = y * np.sqrt(u_e[i] / (nu_e * x_st))

    sel = (eta > 0) & (eta <= PROFILE_ETA_MAX)
    assert sel.sum() >= PROFILE_MIN_POINTS, (
        f"only {sel.sum()} grid points inside eta<={PROFILE_ETA_MAX} at "
        f"x={x_st:.4f} m; the wall-normal grid is too coarse to validate")

    u_sim = u2d[sel, i] / u_e[i]
    u_ref = fprime(eta[sel])
    diff = np.abs(u_sim - u_ref)
    rms = float(np.sqrt(np.mean((u_sim - u_ref) ** 2)))

    assert diff.max() < PROFILE_MAX_TOL, (
        f"max |u/u_e - f'| = {diff.max():.3f} at eta={eta[sel][np.argmax(diff)]:.2f} "
        f"(x={x_st:.4f} m) exceeds {PROFILE_MAX_TOL}")
    assert rms < PROFILE_RMS_TOL, (
        f"profile RMS error {rms:.3f} exceeds {PROFILE_RMS_TOL} at x={x_st:.4f} m")

    v_max = np.abs(v[k_ref, sel, i]).max() / u_e[i]
    assert v_max < PROFILE_V_TOL, (
        f"wall-normal velocity |v|/u_e = {v_max:.3f} is too large at x={x_st:.4f} m")
