"""
SBLI: oblique shock / laminar boundary-layer interaction (2D_solver/SBLI).

An oblique shock (beta = 30.8 deg) generated at the inflow plane impinges on a
laminar flat-plate boundary layer at Xsh = 80 mm (Re_Xsh ~ 1e5), separating it
into a recirculation bubble. Wall pressure and skin friction are compared
against the digitized reference data of Moro et al. (the newest of the three
datasets shipped in 2D_solver/SBLI/ref/).

The committed 564x513 / 5M-step production setup is far too long for CI, so the
test runs a coarsened 276x257 grid with dt raised 10x (CFL_y ~ 0.17). nx=276 is
chosen so dx = 0.64 mm divides the 16 mm inflow run-in exactly, putting the
plate leading edge precisely at x = 0.
"""
import pathlib

import numpy as np
import pytest

from ouxsbli import Case
from ouxsbli.analysis.wall import find_zero_crossings, load_reference, wall_coeffs_from_vtr
from .utils.vtk_reader import extract_number
from .conftest import REPO_ROOT, assert_close_relative

pytestmark = [pytest.mark.integration, pytest.mark.slow]

# must match 2D_solver/SBLI/mod_globals.f90
R = 287.15
gamma = 1.4
M0 = 2.15
p_tot = 25.0e3
T0 = 288.15
blt = 1.0e-3
XSH_MM = 80.0  # shock impingement location [mm]

p0 = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0 = M0 * np.sqrt(gamma * R * T0)
q_inf = 0.5 * rho0 * u0**2

REF_DIR = REPO_ROOT / "2D_solver" / "SBLI" / "ref"
# the published abscissa is x/L (L = 100 mm); X/Xsh = 1.25 x/L
REF_X_SCALE = 1.25

# Cp comparison: window in X/Xsh, and the plateau-normalised tolerances.
# Cp is near zero upstream of the interaction, so a pointwise *relative* error
# is ill-conditioned there; deviations are scaled by the plateau level instead.
CP_X_MIN, CP_X_MAX = 0.4, 1.75
CP_MEAN_TOL = 0.10
CP_MAX_TOL = 0.20
CP_PLATEAU_MIN, CP_PLATEAU_MAX = 1.4, 1.8
CP_PLATEAU_RTOL = 0.12

# Cf: separation-bubble window and location tolerance (in X/Xsh)
CF_X_MIN, CF_X_MAX = 0.3, 1.5
CF_MIN_DEPTH = -1.0e-4
CF_CROSSING_TOL = 0.15

# quasi-steady guards between the last two snapshots
STEADY_CP_TOL = 0.02
STEADY_CROSSING_TOL = 0.05


def _wall(path):
    """(X/Xsh, cf, cp) along the bottom wall for one snapshot."""
    x, cf, cp = wall_coeffs_from_vtr(path, R, p0, q_inf)
    return (x / blt) / XSH_MM, cf, cp


@pytest.fixture(scope="module")
def sbli_run():
    """Build and run the coarsened SBLI case once; return the last two snapshots."""
    workdir = "./tmp/sbli"

    case = Case(
        source="2D_solver/SBLI",
        workdir=workdir,
        nx=276,       # dx = 0.64 mm exactly -> leading edge lands on x = 0
        ny=257,       # y(2) ~ 82 um, enough to resolve the laminar wall layer
        dt=3e-8,      # CFL_y ~ 0.17
        endT=6.0e-3,  # ~25 flow-through times (Lx/u0 = 0.24 ms)
        np=10,        # nt derives to 20000; do not pass nt explicitly
        output_precision=8,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    snaps = sorted(data_dir.glob("Q*.vtr"), key=lambda p: extract_number(p.name))
    assert len(snaps) >= 2, f"expected at least 2 snapshots in {data_dir}, got {len(snaps)}"
    return snaps[-2], snaps[-1]


@pytest.fixture(scope="module")
def cp_reference():
    return load_reference(REF_DIR / "data_Moro_Cp.dat", REF_X_SCALE)


@pytest.fixture(scope="module")
def cf_reference():
    return load_reference(REF_DIR / "data_Moro_Cf.dat", REF_X_SCALE)


def test_interaction_is_quasi_steady(sbli_run, cp_reference):
    """Wall pressure and the bubble extent must settle between the last two outputs."""
    prev, last = sbli_run
    x_p, cf_p, cp_p = _wall(prev)
    x_l, cf_l, cp_l = _wall(last)
    plateau = cp_reference[:, 1].max()

    mask = (x_l >= CP_X_MIN) & (x_l <= CP_X_MAX)
    drift = np.abs(cp_l[mask] - cp_p[mask]).max() / plateau
    assert drift < STEADY_CP_TOL, (
        f"wall pressure still drifting: max |dCp|/Cp_plateau = {drift:.3f} "
        f"between the last two outputs; increase endT")

    bubble = (x_l >= CF_X_MIN) & (x_l <= CF_X_MAX)
    cross_p = find_zero_crossings(x_p[bubble], cf_p[bubble])
    cross_l = find_zero_crossings(x_l[bubble], cf_l[bubble])
    if len(cross_p) == len(cross_l) and cross_l:
        moved = max(abs(a - b) for a, b in zip(cross_p, cross_l))
        assert moved < STEADY_CROSSING_TOL, (
            f"separation/reattachment still moving by {moved:.3f} X/Xsh "
            f"between the last two outputs; increase endT")


def test_cp_matches_reference(sbli_run, cp_reference):
    """Wall Cp(x) must track the Moro et al. distribution through the interaction."""
    _, last = sbli_run
    x, _, cp = _wall(last)

    ref = cp_reference
    sel = (ref[:, 0] >= CP_X_MIN) & (ref[:, 0] <= CP_X_MAX)
    x_ref, cp_ref = ref[sel, 0], ref[sel, 1]
    assert x_ref.size >= 10, f"only {x_ref.size} reference points in the Cp window"

    cp_sim = np.interp(x_ref, x, cp)
    plateau = ref[:, 1].max()
    err = np.abs(cp_sim - cp_ref) / plateau

    assert err.mean() < CP_MEAN_TOL, (
        f"mean Cp deviation {err.mean():.3f} (normalised by the plateau "
        f"Cp={plateau:.3f}) exceeds {CP_MEAN_TOL}")
    assert err.max() < CP_MAX_TOL, (
        f"max Cp deviation {err.max():.3f} exceeds {CP_MAX_TOL} "
        f"at X/Xsh = {x_ref[np.argmax(err)]:.2f}")


def test_cp_plateau_level(sbli_run, cp_reference):
    """The post-interaction pressure plateau must match the reference level."""
    _, last = sbli_run
    x, _, cp = _wall(last)

    ref = cp_reference
    ref_sel = (ref[:, 0] >= CP_PLATEAU_MIN) & (ref[:, 0] <= CP_PLATEAU_MAX)
    sim_sel = (x >= CP_PLATEAU_MIN) & (x <= CP_PLATEAU_MAX)
    assert ref_sel.sum() >= 3 and sim_sel.sum() >= 5

    assert_close_relative(
        cp[sim_sel].mean(), ref[ref_sel, 1].mean(), CP_PLATEAU_RTOL, "Cp plateau")


def test_separation_bubble_matches_reference(sbli_run, cf_reference):
    """A recirculation bubble must form, with separation and reattachment near
    the reference locations."""
    _, last = sbli_run
    x, cf, _ = _wall(last)

    win = (x >= CF_X_MIN) & (x <= CF_X_MAX)
    x_w, cf_w = x[win], cf[win]

    assert cf_w.min() < CF_MIN_DEPTH, (
        f"no separation bubble: min Cf = {cf_w.min():.2e} over "
        f"X/Xsh = {CF_X_MIN}-{CF_X_MAX} (reference reaches "
        f"{cf_reference[:, 1].min():.2e})")

    ref = cf_reference
    ref_win = (ref[:, 0] >= CF_X_MIN) & (ref[:, 0] <= CF_X_MAX)
    ref_cross = find_zero_crossings(ref[ref_win, 0], ref[ref_win, 1])
    sim_cross = find_zero_crossings(x_w, cf_w)

    assert len(ref_cross) == 2, (
        f"expected 2 reference Cf=0 crossings, found {ref_cross}")
    assert len(sim_cross) == 2, (
        f"expected separation and reattachment, found {len(sim_cross)} Cf=0 "
        f"crossings at {[f'{c:.3f}' for c in sim_cross]} "
        f"(reference: {[f'{c:.3f}' for c in ref_cross]})")

    for label, sim_c, ref_c in zip(("separation", "reattachment"), sim_cross, ref_cross):
        assert abs(sim_c - ref_c) < CF_CROSSING_TOL, (
            f"{label} at X/Xsh = {sim_c:.3f} differs from the reference "
            f"{ref_c:.3f} by more than {CF_CROSSING_TOL}")
