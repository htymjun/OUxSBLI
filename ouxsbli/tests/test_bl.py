"""
BL: laminar flat-plate boundary layer (2D_solver/BL).

A uniform M=0.1 freestream develops a boundary layer over a flat plate whose
leading edge sits at x=0 (symmetry BC upstream of it, no-slip downstream).
At this Mach number compressibility is negligible, so both the skin-friction
distribution and the streamwise velocity profile must collapse onto the
incompressible Blasius solution.

Comparisons use the *local* boundary-layer edge state rather than the nominal
freestream: the growing displacement thickness accelerates the flow by ~2% in
the 12 mm tall domain, which is a property of the case geometry, not a solver
error. Normalising by u_e(x) is the standard way to compare a boundary layer
in a finite domain against flat-plate theory.

The committed dt=3e-9 sits ~50x below the stability limit; the test raises it
to 5e-8 (CFL_y ~ 0.18) so ~2.8 flow-through times fit in 200k steps.
"""
import pathlib

import numpy as np
import pytest

from ouxsbli import Case
from ouxsbli.analysis.blasius import fprime
from ouxsbli.analysis.wall import edge_state, one_sided_deriv_3pt, sutherland_mu
from .utils.vtk_reader import extract_number, getGrid, getQ

pytestmark = pytest.mark.integration

# gas constant; must match 2D_solver/BL/mod_globals.f90 (freestream there is
# M=0.1 at T0=288.15 K, but every comparison below uses the local edge state)
R = 287.15

BLASIUS_CF_CONST = 0.664  # Cf = 0.664 / sqrt(Re_x)

# Cf comparison window [m]: starts ~34 cells past the leading edge so the
# startup transient stays out, ends 8 cells short of the outlet.
CF_X_MIN = 0.032
CF_X_MAX = 0.072
CF_MEAN_TOL = 0.05
CF_MAX_TOL = 0.08

# u-profile station [m] (Re_x ~ 3.2e4) and profile tolerances
PROFILE_X = 0.056
PROFILE_ETA_MAX = 6.0
PROFILE_MIN_POINTS = 8
PROFILE_MAX_TOL = 0.02
PROFILE_RMS_TOL = 0.01
PROFILE_V_TOL = 0.02

# quasi-steady guard between the last two snapshots
STEADY_TOL = 0.01


def _wall_profile(path):
    """(x, cf_local, Re_x_local) using the local boundary-layer edge state."""
    ni, nj, nk, x, y, z = getGrid(str(path))
    rho, u, v, w, p = getQ(str(path), ni, nj, nk)
    rho2d, u2d, p2d = rho[0], u[0], p[0]

    u_e, rho_e, mu_e = edge_state(rho2d, u2d, p2d, R)

    T_w = p2d[0, :] / (rho2d[0, :] * R)
    dudy = one_sided_deriv_3pt(y[0], y[1], y[2], np.zeros(ni), u2d[1, :], u2d[2, :])
    tau_w = sutherland_mu(T_w) * dudy

    cf = tau_w / (0.5 * rho_e * u_e**2)
    Re_x = rho_e * u_e * x / mu_e
    return x, cf, Re_x


@pytest.fixture(scope="module")
def bl_run():
    """Build and run the BL case once, returning the last two snapshots."""
    workdir = "./tmp/bl"

    case = Case(
        source="2D_solver/BL",
        workdir=workdir,
        dt=5e-8,      # CFL_y ~ 0.18
        endT=1.0e-2,  # ~2.8 flow-through times (Lx/u0 = 3.5 ms)
        np=10,        # nt derives to 20000; do not pass nt explicitly
        output_precision=8,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    snaps = sorted(data_dir.glob("Q*.vtr"), key=lambda p: extract_number(p.name))
    assert len(snaps) >= 2, f"expected at least 2 snapshots in {data_dir}, got {len(snaps)}"
    return snaps[-2], snaps[-1]


def test_grid_places_leading_edge_at_origin(bl_run):
    """set_grid must put the first no-slip point exactly at x=0 and the wall at y=0."""
    _, last = bl_run
    ni, nj, nk, x, y, z = getGrid(str(last))

    assert y[0] == pytest.approx(0.0, abs=1e-12), f"wall not at y=0: y[0]={y[0]}"
    assert np.abs(x).min() < 1e-9, (
        f"leading edge is off-grid: closest x to 0 is {np.abs(x).min():.3e} m")
    assert x[-1] >= 0.079, f"domain too short for the Cf window: x_max={x[-1]:.4f} m"


def test_flow_is_quasi_steady(bl_run):
    """The window-mean Cf must stop changing between the last two outputs."""
    prev, last = bl_run
    x_p, cf_p, _ = _wall_profile(prev)
    x_l, cf_l, _ = _wall_profile(last)

    mask = (x_l >= CF_X_MIN) & (x_l <= CF_X_MAX)
    mean_prev = cf_p[mask].mean()
    mean_last = cf_l[mask].mean()
    change = abs(mean_last - mean_prev) / abs(mean_last)

    assert change < STEADY_TOL, (
        f"boundary layer is still developing: window-mean Cf changed {change:.2%} "
        f"between the last two outputs ({mean_prev:.4e} -> {mean_last:.4e}); "
        "increase endT")


def test_cf_matches_blasius(bl_run):
    """Cf(x) must follow 0.664/sqrt(Re_x) over the comparison window."""
    _, last = bl_run
    x, cf, Re_x = _wall_profile(last)

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
    _, last = bl_run
    ni, nj, nk, x, y, z = getGrid(str(last))
    rho, u, v, w, p = getQ(str(last), ni, nj, nk)
    rho2d, u2d, p2d = rho[0], u[0], p[0]

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

    v_max = np.abs(v[0, sel, i]).max() / u_e[i]
    assert v_max < PROFILE_V_TOL, (
        f"wall-normal velocity |v|/u_e = {v_max:.3f} is too large at x={x_st:.4f} m")
