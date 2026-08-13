"""Validation of 3D_solver/TBL against Guarini, Moser, Shariff & Wray,
"Direct numerical simulation of a supersonic turbulent boundary layer at
Mach 2.5", JFM 414 (2000) 1-33.

Unlike the other integration tests here, this one does **not** build and run the
solver: a converged supersonic TBL needs ~24 flow-throughs of sampling on top of
~12 of spin-up, which is hours of GPU time and ~13 GB of snapshots. Instead the
expensive part is done once, offline::

    cd 3D_solver/TBL
    # stage A, RESTART=False, endT=3.6d-4, np=10   -> spin-up
    # stage B, RESTART=True,  endT=7.2d-4, np=100  -> sampling
    bash calc.sh
    python -m ouxsbli.analysis.tbl_stats accumulate data --out ouxsbli/tests/data/tbl_acc.npz

and this file asserts against the resulting moment file, which is small enough
to keep in the repo. The tests skip cleanly when it is absent.

Tolerances are deliberately split in two. Integral quantities (Cf, Re_theta, H)
are compared against the paper's published values. Quantities whose *estimator*
is biased at this Reynolds number -- kappa above all, which comes out 0.347
rather than 0.40 even when the analytic profile that went in was built with
0.40 -- are compared against ``reference_station``, i.e. the same estimator run
on the paper's own composite profile on the same grid. Comparing those against
the paper's headline numbers would be measuring the estimator, not the solver.
"""
import pathlib

import numpy as np
import pytest

from ouxsbli.analysis import tbl_stats as T

ACC_PATH = pathlib.Path(__file__).parent / "data" / "tbl_acc.npz"

pytestmark = pytest.mark.skipif(
    not ACC_PATH.exists(),
    reason=f"{ACC_PATH} not present; see this module's docstring to regenerate",
)


@pytest.fixture(scope="module")
def tbl():
    """(station, streamwise sweep, reference station) at the Re_theta=1577 point."""
    acc = dict(np.load(ACC_PATH))
    ix, sw = T.find_station(acc, T.REF["Re_theta"])
    st = T.station(acc, ix)
    return st, sw, T.reference_station(st["y"])


def test_sampling_station_is_inside_the_useful_domain(tbl):
    """The Re_theta=1577 station must sit clear of the inflow and the recycle plane.

    The Lund-Wu-Squires inlet needs several delta to shed its artifacts, and the
    recycle plane at 0.9*Lx is where the profile is copied from; a station
    outside that window is not a statement about the boundary layer.
    """
    st, sw, _ = tbl
    frac = st["ix"] / len(sw["x"])
    assert 0.40 < frac < 0.90, f"station at {frac:.2f} Lx is outside the usable window"
    assert st["x"] / st["d99"] > 4.0, (
        f"station is only {st['x'] / st['d99']:.1f} delta from the inlet")


def test_station_reaches_the_reference_reynolds_number(tbl):
    """Re_theta at the station must actually be the paper's, not merely nearby."""
    st, _, _ = tbl
    assert abs(st["Re_theta"] - T.REF["Re_theta"]) / T.REF["Re_theta"] < 0.05


def test_wall_resolution_is_dns_grade(tbl):
    """First cell below y+ = 1 and at least ~10 points inside y+ = 10."""
    st, _, _ = tbl
    assert st["yp"][1] < 1.0, f"first cell at y+ = {st['yp'][1]:.2f}"
    assert (st["yp"] < 10.0).sum() >= 10


def test_skin_friction_matches(tbl):
    """Cf against the paper's 0.00282."""
    st, _, _ = tbl
    rel = abs(st["Cf"] - T.REF["Cf"]) / T.REF["Cf"]
    assert rel < 0.10, f"Cf = {st['Cf']:.5f} vs {T.REF['Cf']:.5f} ({rel:.1%})"


def test_shape_factor_matches(tbl):
    """H = delta*/theta against the paper's 3.97."""
    st, _, _ = tbl
    rel = abs(st["H"] - T.REF["H"]) / T.REF["H"]
    assert rel < 0.10, f"H = {st['H']:.3f} vs {T.REF['H']:.3f} ({rel:.1%})"


def test_wall_is_adiabatic(tbl):
    """Zero-gradient wall temperature must land on the recovery temperature."""
    st, _, _ = tbl
    T_aw = T.freestream()["T_inf"] * (1.0 + 0.89 * 0.5 * (T.GAMMA - 1.0) * T.M0**2)
    assert abs(st["T_w"] - T_aw) / T_aw < 0.03


def test_van_driest_profile_follows_the_log_law(tbl):
    """U_c+ must track (1/0.40)ln y+ + 4.7 over the paper's 30 <= y+ <= 70.

    Judged against the same measure applied to the reference composite profile,
    which is itself 0.2 off the bare log law because of the wake.
    """
    st, _, ref = tbl
    dev, dev_ref = T.log_law_deviation(st), T.log_law_deviation(ref)
    assert dev < dev_ref + 0.6, f"rms deviation {dev:.3f} vs reference {dev_ref:.3f}"


def test_log_law_constants_match_the_reference_estimator(tbl):
    """kappa and C, compared against the same estimator on the reference profile."""
    st, _, ref = tbl
    k, C = T.log_law_fit(st)
    k_ref, C_ref = T.log_law_fit(ref)
    assert abs(k - k_ref) / k_ref < 0.15, f"kappa {k:.3f} vs reference {k_ref:.3f}"
    assert abs(C - C_ref) < 1.5, f"C {C:.3f} vs reference {C_ref:.3f}"


def test_morkovin_scaling_collapses_the_streamwise_intensity(tbl):
    """sqrt(rho/rho_w) u'/u_tau must peak near Spalart's incompressible ~2.7-2.8.

    Figure 6(b): density-scaled intensities collapse onto the incompressible
    simulations, whereas the raw u'/u_tau of figure 6(a) sits lower.
    """
    st, _, _ = tbl
    scaled = np.sqrt(st["rho"] / st["rho_w"]) * st["urms"] / st["u_tau"]
    raw = st["urms"] / st["u_tau"]
    inner = st["yp"] < 100
    assert 2.2 < scaled[inner].max() < 3.4, f"peak Morkovin u' = {scaled[inner].max():.2f}"
    assert scaled[inner].max() > raw[inner].max(), (
        "density scaling must raise the peak, not lower it")


def test_total_stress_is_constant_near_the_wall(tbl):
    """mu dU/dy - rho<u''v''> should be ~tau_w out to y+ = 30-40 (figure 9)."""
    from ouxsbli.analysis.wall import sutherland_mu

    st, _, _ = tbl
    total = sutherland_mu(st["T"]) * np.gradient(st["u_fav"], st["y"]) - st["Ruv"]
    band = (st["yp"] > 5) & (st["yp"] < 35)
    ratio = total[band] / st["tau_w"]
    assert np.all(np.abs(ratio - 1.0) < 0.20), (
        f"total stress ranges {ratio.min():.2f}-{ratio.max():.2f} of tau_w")


def test_reynolds_shear_stress_is_negative_and_bounded(tbl):
    """-rho<u''v''>/tau_w must rise from 0 at the wall to O(1) and stay bounded."""
    st, _, _ = tbl
    r = -st["Ruv"] / st["tau_w"]
    inner = st["yp"] < 100
    assert 0.5 < r[inner].max() < 1.15, f"peak -rho<u''v''>/tau_w = {r[inner].max():.2f}"
    assert abs(r[0]) < 0.05, "Reynolds stress must vanish at the wall"


def test_pressure_fluctuations_match(tbl):
    """p'_rms/(rho_w u_tau^2): 2.7 at the wall, 3.0 peak, 0.47 in the freestream."""
    st, _, _ = tbl
    q = st["rho_w"] * st["u_tau"] ** 2
    for got, ref, label in (
        (st["prms"][0] / q, T.REF["p_rms_wall"], "wall"),
        (st["prms"].max() / q, T.REF["p_rms_peak"], "peak"),
        (T.freestream_prms(st) / q, T.REF["p_rms_inf"], "freestream"),
    ):
        assert abs(got - ref) / ref < 0.40, f"p'_rms {label} = {got:.2f} vs {ref:.2f}"


def test_freestream_is_not_contaminated(tbl):
    """Turbulence must not leak to the top boundary; the layer has to be resolved.

    A too-shallow domain or a reflecting top condition shows up here first.
    """
    st, _, _ = tbl
    assert st["d99"] < 0.45 * st["y"][-1], (
        f"delta99 = {st['d99'] * 1e3:.2f} mm fills too much of the {st['y'][-1] * 1e3:.2f} mm domain")
    assert st["urms"][-1] / st["u_e"] < 0.02, "freestream streamwise fluctuation too large"


def test_statistics_are_converged(tbl):
    """Split-half check: the two halves of the sample must agree.

    Requires the half-sample files written alongside tbl_acc.npz.
    """
    a, b = ACC_PATH.with_name("tbl_acc_h1.npz"), ACC_PATH.with_name("tbl_acc_h2.npz")
    if not (a.exists() and b.exists()):
        pytest.skip("half-sample accumulations not present")
    st, _, _ = tbl
    sa, sb = T.station(dict(np.load(a)), st["ix"]), T.station(dict(np.load(b)), st["ix"])
    for k in ("Cf", "Re_theta", "u_tau", "d99"):
        rel = abs(sb[k] - sa[k]) / abs(sa[k])
        assert rel < 0.05, f"{k} differs by {rel:.1%} between sample halves"
