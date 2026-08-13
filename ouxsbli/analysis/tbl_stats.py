"""Turbulence statistics for the 3D_solver/TBL supersonic boundary layer.

Averaging is over the homogeneous spanwise direction and over time; the
streamwise direction is inhomogeneous, so every moment is accumulated per
x-column and a sampling station is chosen afterwards.

The reference is Guarini, Moser, Shariff & Wray, "Direct numerical simulation
of a supersonic turbulent boundary layer at Mach 2.5", JFM 414 (2000) 1-33.

Typical use::

    python -m ouxsbli.analysis.tbl_stats accumulate 3D_solver/TBL/data --skip 20
    python -m ouxsbli.analysis.tbl_stats report tbl_acc.npz

``accumulate`` makes a single pass over the snapshots and writes the raw
moments to a .npz; ``report`` turns those into profiles, plots and a summary
table without re-reading the (large) VTK files.
"""

import glob
import os
import re
import sys

import numpy as np

from .wall import one_sided_deriv_3pt, sutherland_mu

# Freestream state of the case, mirroring 3D_solver/TBL/mod_globals.f90.
GAMMA = 1.4
R_GAS = 287.03
PR = 0.72
M0 = 2.5
P_TOT = 100.0e3
T_TOT = 295.0
CP = GAMMA * R_GAS / (GAMMA - 1.0)

# Reference values quoted by Guarini et al. (2000).
REF = {
    "Re_theta": 1577.0,
    "Re_theta0": 849.0,
    "Re_dstar": 6258.0,
    "H": 6258.0 / 1577.0,
    "Cf": 0.00282,
    "kappa": 0.40,
    "C_log": 4.7,
    "p_rms_wall": 2.7,
    "p_rms_peak": 3.0,
    "p_rms_inf": 0.47,
}

# The z-direction carries three periodic ghost planes at each end
# (set_bc_cyclic_z); only k = 4..nz-3 (1-based) are unique.
NGHOST_Z = 3

# Moments accumulated in the single pass, each a (ny, nx) spanwise-and-time
# mean. Reynolds statistics come from the unweighted moments, Favre statistics
# from the rho-weighted ones.
_MOMENTS = (
    "rho", "u", "v", "w", "p", "T", "Tt",
    "rhou", "rhov", "rhow", "rhoT", "rhoTt",
    "rhouu", "rhovv", "rhoww", "rhouv", "rhouT", "rhovT", "rhoTT", "rhoTtTt",
    "uu", "vv", "ww", "uv", "pp", "TT", "TtTt", "uT",
    "ox", "oy", "oz", "oxox", "oyoy", "ozoz",
)


def freestream():
    """Static freestream state derived from the case's total conditions."""
    fac = 1.0 + 0.5 * (GAMMA - 1.0) * M0**2
    T_inf = T_TOT / fac
    p_inf = P_TOT / fac ** (GAMMA / (GAMMA - 1.0))
    rho_inf = p_inf / (R_GAS * T_inf)
    u_inf = M0 * np.sqrt(GAMMA * R_GAS * T_inf)
    mu_inf = sutherland_mu(T_inf)
    return {
        "T_inf": T_inf, "p_inf": p_inf, "rho_inf": rho_inf,
        "u_inf": u_inf, "mu_inf": mu_inf,
        "Re_unit": rho_inf * u_inf / mu_inf,
    }


def _ddx_periodic6(f, d, axis):
    """6th-order central difference along a uniform, periodic axis."""
    r = [np.roll(f, -s, axis=axis) - np.roll(f, s, axis=axis) for s in (1, 2, 3)]
    return (0.75 * r[0] - 0.15 * r[1] + (1.0 / 60.0) * r[2]) / d


def _ddx_uniform6(f, d, axis):
    """6th-order central in the interior, 2nd order in the three edge columns.

    The streamwise direction is not periodic, so the roll-based stencil wraps
    incorrectly near the inlet and outlet; those columns are overwritten and are
    never used as a sampling station anyway.
    """
    out = _ddx_periodic6(f, d, axis)
    edge = np.gradient(f, d, axis=axis)
    sl = [slice(None)] * f.ndim
    for idx in (0, 1, 2, -3, -2, -1):
        sl[axis] = idx
        out[tuple(sl)] = edge[tuple(sl)]
    return out


def _vorticity(u, v, w, x, y, z):
    """(omega_x, omega_y, omega_z) for fields shaped (nz, ny, nx).

    x and z are uniform; y is stretched, so np.gradient's non-uniform form is
    used there.
    """
    dx = x[1] - x[0]
    dz = z[1] - z[0]
    du_dy = np.gradient(u, y, axis=1)
    du_dz = _ddx_periodic6(u, dz, axis=0)
    dv_dx = _ddx_uniform6(v, dx, axis=2)
    dv_dz = _ddx_periodic6(v, dz, axis=0)
    dw_dx = _ddx_uniform6(w, dx, axis=2)
    dw_dy = np.gradient(w, y, axis=1)
    return dw_dy - dv_dz, du_dz - dw_dx, dv_dx - du_dy


def accumulate(data_dir, skip=0, stop=None, stride=1, out=None, verbose=True):
    """Single pass over the snapshots; returns (and optionally saves) moments.

    ``skip``/``stop``/``stride`` slice the snapshot list: ``skip`` drops the
    spin-up, and ``stop`` is what lets the split-half convergence check reuse
    this same routine.

    Reads through :mod:`ouxsbli.analysis.vtr_raw`, which needs only numpy, so
    this can run on a compute node without the ``vtk`` package -- reducing the
    snapshots in place beats copying tens of GB back.
    """
    from .vtr_raw import getGrid, getQ

    def extract_number(p):
        m = re.search(r"Q(\d+)\.vtr$", str(p))
        return int(m.group(1)) if m else float("inf")

    files = sorted(glob.glob(os.path.join(str(data_dir), "Q*.vtr")), key=extract_number)
    files = files[skip:stop:stride]
    if not files:
        raise FileNotFoundError(f"no Q*.vtr under {data_dir} for [{skip}:{stop}:{stride}]")

    ni, nj, nk, x, y, z = getGrid(files[0])
    ks = slice(NGHOST_Z, nk - NGHOST_Z)
    nz_real = nk - 2 * NGHOST_Z

    acc = {k: np.zeros((nj, ni)) for k in _MOMENTS}
    nsamp = 0
    for n, f in enumerate(files):
        rho, u, v, w, p = (a.astype(np.float64) for a in getQ(f, ni, nj, nk))
        T = p / (rho * R_GAS)
        Tt = T + 0.5 * (u**2 + v**2 + w**2) / CP
        ox, oy, oz = _vorticity(u, v, w, x, y, z)

        # drop the periodic ghost planes before averaging
        rho, u, v, w, p, T, Tt, ox, oy, oz = (
            a[ks] for a in (rho, u, v, w, p, T, Tt, ox, oy, oz)
        )

        terms = {
            "rho": rho, "u": u, "v": v, "w": w, "p": p, "T": T, "Tt": Tt,
            "rhou": rho * u, "rhov": rho * v, "rhow": rho * w,
            "rhoT": rho * T, "rhoTt": rho * Tt,
            "rhouu": rho * u * u, "rhovv": rho * v * v, "rhoww": rho * w * w,
            "rhouv": rho * u * v, "rhouT": rho * u * T, "rhovT": rho * v * T,
            "rhoTT": rho * T * T, "rhoTtTt": rho * Tt * Tt,
            "uu": u * u, "vv": v * v, "ww": w * w, "uv": u * v,
            "pp": p * p, "TT": T * T, "TtTt": Tt * Tt, "uT": u * T,
            "ox": ox, "oy": oy, "oz": oz,
            "oxox": ox * ox, "oyoy": oy * oy, "ozoz": oz * oz,
        }
        for key, val in terms.items():
            acc[key] += val.sum(axis=0)
        nsamp += nz_real
        if verbose:
            print(f"  [{n + 1}/{len(files)}] {os.path.basename(f)}", flush=True)

    for key in _MOMENTS:
        acc[key] /= nsamp
    acc.update(x=x, y=y, z=z, nsamp=np.array(nsamp), nfiles=np.array(len(files)))
    if out:
        np.savez_compressed(out, **acc)
        if verbose:
            print(f"wrote {out}  ({len(files)} files, {nsamp} spanwise-time samples)")
    return acc


def _edge(y, u_fav, rho):
    """delta99 plus the edge state, found from the first 0.99*u_e crossing.

    u_e is taken as the median over the top quarter of the domain rather than
    the single top-boundary point, so a slightly non-monotonic freestream does
    not shift the edge.
    """
    top = y >= 0.75 * y[-1]
    u_e = float(np.median(u_fav[top]))
    rho_e = float(np.median(rho[top]))
    tgt = 0.99 * u_e
    j = int(np.argmax(u_fav >= tgt))
    if j < 1:
        raise RuntimeError("no 99% crossing found")
    f = (tgt - u_fav[j - 1]) / (u_fav[j] - u_fav[j - 1])
    d99 = float(y[j - 1] + f * (y[j] - y[j - 1]))
    return d99, u_e, rho_e


def station(acc, ix, half_width=0):
    """Full profile set at streamwise index ``ix``.

    Both decompositions are carried: figures 6 and 10 of the paper use Reynolds
    fluctuations, the stress balance of equation (3.6) uses the Favre form.
    """
    fs = freestream()
    ni = acc["rho"].shape[1]
    lo, hi = max(0, ix - half_width), min(ni, ix + half_width + 1)
    m = {k: acc[k][:, lo:hi].mean(axis=1) for k in _MOMENTS}
    y = acc["y"]

    rho, p, T = m["rho"], m["p"], m["T"]

    # Favre means
    u_fav = m["rhou"] / rho
    v_fav = m["rhov"] / rho
    T_fav = m["rhoT"] / rho

    # Favre stresses:  rho*<a''b''> = <rho a b> - <rho a><rho b>/<rho>
    Ruu = m["rhouu"] - m["rhou"] ** 2 / rho
    Rvv = m["rhovv"] - m["rhov"] ** 2 / rho
    Rww = m["rhoww"] - m["rhow"] ** 2 / rho
    Ruv = m["rhouv"] - m["rhou"] * m["rhov"] / rho
    RuT = m["rhouT"] - m["rhou"] * m["rhoT"] / rho
    RvT = m["rhovT"] - m["rhov"] * m["rhoT"] / rho
    RTT = m["rhoTT"] - m["rhoT"] ** 2 / rho

    # Reynolds fluctuations
    def _rms(sq, mean):
        return np.sqrt(np.maximum(sq - mean**2, 0.0))

    urms = _rms(m["uu"], m["u"])
    vrms = _rms(m["vv"], m["v"])
    wrms = _rms(m["ww"], m["w"])
    prms = _rms(m["pp"], m["p"])
    Trms = _rms(m["TT"], m["T"])
    Ttrms = _rms(m["TtTt"], m["Tt"])
    uv = m["uv"] - m["u"] * m["v"]
    oxrms = _rms(m["oxox"], m["ox"])
    oyrms = _rms(m["oyoy"], m["oy"])
    ozrms = _rms(m["ozoz"], m["oz"])

    # wall quantities: j=0 is the wall (no-slip, zero-gradient rho and p)
    rho_w, T_w, p_w = rho[0], T[0], p[0]
    mu_w = sutherland_mu(T_w)
    nu_w = mu_w / rho_w
    dudy_w = one_sided_deriv_3pt(y[0], y[1], y[2], 0.0, u_fav[1], u_fav[2])
    tau_w = mu_w * dudy_w
    u_tau = np.sqrt(abs(tau_w) / rho_w)
    d_nu = nu_w / u_tau
    Cf = 2.0 * tau_w / (fs["rho_inf"] * fs["u_inf"] ** 2)

    # integral thicknesses, integrated to 1.25*delta99 where the integrands
    # have decayed; truncating at delta99 itself loses ~1-2%
    d99, u_e, rho_e = _edge(y, u_fav, rho)
    jmax = int(np.searchsorted(y, 1.25 * d99))
    sl = slice(0, max(jmax, 4))
    ratio = rho * u_fav / (rho_e * u_e)
    dstar = float(np.trapezoid((1.0 - ratio)[sl], y[sl]))
    theta = float(np.trapezoid((ratio * (1.0 - u_fav / u_e))[sl], y[sl]))

    Re_theta = fs["rho_inf"] * fs["u_inf"] * theta / fs["mu_inf"]
    Re_dstar = fs["rho_inf"] * fs["u_inf"] * dstar / fs["mu_inf"]
    Re_theta0 = fs["rho_inf"] * fs["u_inf"] * theta / mu_w
    Re_tau = d99 / d_nu

    # van Driest transform  U_c^+ = int_0^{u+} sqrt(rho_bar/rho_w) du+
    # (with p ~ const across the layer this is the paper's int sqrt(T_w/T) dU).
    yp = y / d_nu
    up = u_fav / u_tau
    g = np.sqrt(np.maximum(rho / rho_w, 0.0))
    ucp = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(up))))

    # strong Reynolds analogy
    with np.errstate(divide="ignore", invalid="ignore"):
        Ma = u_fav / np.sqrt(GAMMA * R_GAS * T_fav)
        sra = (np.sqrt(np.maximum(RTT / rho, 0.0)) / T_fav) / (
            (GAMMA - 1.0) * Ma**2 * np.sqrt(np.maximum(Ruu / rho, 0.0)) / u_fav
        )
        R_uT = RuT / np.sqrt(np.maximum(Ruu * RTT, 1e-300))
        Pr_t = (Ruv * np.gradient(T_fav, y)) / (RvT * np.gradient(u_fav, y))

    return {
        "ix": ix, "x": float(acc["x"][ix]), "span": (lo, hi), "y": y, "yp": yp,
        "rho": rho, "p": p, "T": T, "T_fav": T_fav, "Tt": m["Tt"],
        "u_rey": m["u"], "u_fav": u_fav, "v_fav": v_fav, "up": up, "ucp": ucp,
        "urms": urms, "vrms": vrms, "wrms": wrms, "prms": prms,
        "Trms": Trms, "Ttrms": Ttrms, "uv": uv,
        "Ruu": Ruu, "Rvv": Rvv, "Rww": Rww, "Ruv": Ruv, "RTT": RTT,
        "oxrms": oxrms, "oyrms": oyrms, "ozrms": ozrms,
        "sra": sra, "R_uT": R_uT, "Pr_t": Pr_t,
        "rho_w": rho_w, "T_w": T_w, "p_w": p_w, "mu_w": mu_w, "nu_w": nu_w,
        "tau_w": tau_w, "u_tau": u_tau, "d_nu": d_nu, "Cf": Cf,
        "u_e": u_e, "rho_e": rho_e, "d99": d99, "dstar": dstar, "theta": theta,
        "H": dstar / theta, "Re_theta": Re_theta, "Re_dstar": Re_dstar,
        "Re_theta0": Re_theta0, "Re_tau": Re_tau, "fs": fs,
    }


def sweep(acc, half_width=0):
    """Cf, delta99, Re_theta, u_tau and Re_tau as functions of x."""
    ni = acc["rho"].shape[1]
    keys = ("Cf", "d99", "Re_theta", "Re_dstar", "H", "u_tau", "Re_tau")
    out = {k: np.full(ni, np.nan) for k in keys}
    for ix in range(3, ni - 3):
        try:
            s = station(acc, ix, half_width)
        except Exception:
            continue
        for k in keys:
            out[k][ix] = s[k]
    out["x"] = acc["x"]
    return out


def find_station(acc, target_re_theta=REF["Re_theta"], lo=0.40, hi=0.92, sw=None):
    """Streamwise index where Re_theta is closest to ``target``, within [lo,hi]*Lx."""
    if sw is None:
        sw = sweep(acc)
    ni = len(sw["x"])
    i0, i1 = int(lo * ni), int(hi * ni)
    seg = sw["Re_theta"][i0:i1]
    if not np.isfinite(seg).any():
        raise RuntimeError("no valid Re_theta in the search window")
    return i0 + int(np.nanargmin(np.abs(seg - target_re_theta))), sw


def reference_profile(y, Cf=REF["Cf"], kappa=REF["kappa"], C=REF["C_log"], Pi=0.25):
    """Guarini's composite mean profile (their eq. 3.3) sampled on grid ``y``.

    Reichardt's inner profile plus Finley's wake, inverted through the Van
    Driest transform and closed with the Crocco-Busemann temperature. Pushing
    this through :func:`station` gives the reference every estimator here
    should reproduce, on the same grid and at the same Reynolds number as the
    DNS -- which is the only fair comparison for quantities like kappa that are
    ill-conditioned at Re_tau ~ 300.
    """
    from scipy.optimize import brentq

    fs = freestream()
    rf = 0.89
    C1 = -np.log(kappa) / kappa + C
    T_inf = fs["T_inf"]
    T_w = T_inf * (1.0 + rf * 0.5 * (GAMMA - 1.0) * M0**2)
    A = np.sqrt(rf * 0.5 * (GAMMA - 1.0) * M0**2 / (T_w / T_inf))
    rho_w = fs["p_inf"] / (R_GAS * T_w)

    def ucp(yp, dp):
        e = np.minimum(np.asarray(yp, float) / dp, 1.0)
        return (np.log(1 + kappa * yp) / kappa
                + C1 * (1 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-0.33 * yp))
                + (e**2 - e**3 + 6 * Pi * e**2 - 4 * Pi * e**3) / kappa)

    u_inf_p = 1.0 / np.sqrt(Cf / 2.0 * (fs["rho_inf"] / rho_w))
    dp = brentq(lambda d: ucp(np.array([d]), d)[0] - np.arcsin(A) / A * u_inf_p, 50, 5000)
    u_tau = fs["u_inf"] / u_inf_p
    d_nu = sutherland_mu(T_w) / rho_w / u_tau

    u = np.clip(np.sin(A * ucp(y / d_nu, dp) / u_inf_p) / A * fs["u_inf"], 0.0, fs["u_inf"])
    u[0] = 0.0
    T = T_w - rf * 0.5 * (GAMMA - 1.0) * M0**2 * T_inf * (u / fs["u_inf"]) ** 2
    rho = fs["p_inf"] / (R_GAS * T)
    return {"y": y, "u": u, "T": T, "rho": rho, "p": np.full_like(y, fs["p_inf"]),
            "delta_plus": dp, "u_tau": u_tau, "d_nu": d_nu}


def reference_station(y, **kw):
    """:func:`station` applied to :func:`reference_profile` on the same grid."""
    r = reference_profile(y, **kw)
    rho, u, T, p = r["rho"], r["u"], r["T"], r["p"]
    ny = len(y)
    col = lambda a: np.repeat(a[:, None], 8, axis=1)  # noqa: E731
    acc = {k: col(np.zeros(ny)) for k in _MOMENTS}
    acc.update(rho=col(rho), u=col(u), p=col(p), T=col(T), rhou=col(rho * u),
               rhoT=col(rho * T), rhouu=col(rho * u * u), uu=col(u * u),
               pp=col(p * p), TT=col(T * T))
    acc["x"] = np.arange(8, dtype=float)
    acc["y"] = y
    acc["z"] = np.arange(8, dtype=float)
    return station(acc, 4)


def freestream_prms(st, y_lo=1.4, y_hi=1.8):
    """p'_rms in the freestream, sampled over y/delta in [y_lo, y_hi].

    Not the top grid point: ``set_bc`` pins the boundary pressure to p0, so
    p'_rms there is identically zero by construction rather than physically.
    Guarini's radiated-noise value of 0.47 is a freestream quantity, so it has
    to be read below the boundary.
    """
    e = st["y"] / st["d99"]
    m = (e >= y_lo) & (e <= y_hi)
    if not m.any():
        m = e >= 0.5 * e.max()
    return float(np.median(st["prms"][m]))


def log_law_deviation(st, yp_lo=30.0, yp_hi=70.0):
    """RMS departure of U_c+ from (1/kappa)ln y+ + C over the paper's band.

    This is the paper's actual claim -- "in the region of 30 <= y+ <= 70, the
    simulation data fall on the log-law curve" -- and unlike a fitted kappa it
    is well conditioned at this Reynolds number.
    """
    yp, ucp = st["yp"], st["ucp"]
    m = (yp >= yp_lo) & (yp <= yp_hi)
    if m.sum() < 2:
        return np.nan
    law = np.log(yp[m]) / REF["kappa"] + REF["C_log"]
    return float(np.sqrt(np.mean((ucp[m] - law) ** 2)))


def log_law_fit(st, yp_lo=20.0, yp_hi=150.0):
    """kappa and C of the van Driest profile, by the paper's own method.

    Guarini et al. determine kappa from the minimum of the diagnostic function
    Xi = y+ dU_c+/dy+ and then back out C at that location. A straight
    least-squares fit over a fixed band does not work here: at Re_tau ~ 300 the
    log region is only a few points wide, and fitting the exact analytic
    composite profile over 30 <= y+ <= 70 returns kappa = 0.33 rather than the
    0.40 that went into it.
    """
    yp, ucp = st["yp"], st["ucp"]
    m = (yp >= yp_lo) & (yp <= yp_hi)
    if m.sum() < 3:
        return np.nan, np.nan
    lg = np.log(yp[m])
    xi = np.gradient(ucp[m], lg)  # y+ dU+/dy+ = dU+/d(ln y+)
    j = int(np.argmin(xi))
    kappa = 1.0 / xi[j]
    return kappa, ucp[m][j] - np.log(yp[m][j]) / kappa


def summary(st):
    """Scalar comparison against Guarini et al.

    Third column is the same estimator applied to the analytic composite
    profile of the paper on this grid, so that quantities which the estimator
    itself biases at Re_tau ~ 300 (kappa, C, Re_tau) can be judged fairly.
    """
    q = st["rho_w"] * st["u_tau"] ** 2
    ref = reference_station(st["y"])
    k_dns, C_dns = log_law_fit(st)
    k_ref, C_ref = log_law_fit(ref)

    rows = [
        ("Re_theta", st["Re_theta"], REF["Re_theta"], ref["Re_theta"]),
        ("Re_theta0 (wall mu)", st["Re_theta0"], REF["Re_theta0"], ref["Re_theta0"]),
        ("Re_delta*", st["Re_dstar"], REF["Re_dstar"], ref["Re_dstar"]),
        ("H = delta*/theta", st["H"], REF["H"], ref["H"]),
        ("Cf", st["Cf"], REF["Cf"], ref["Cf"]),
        ("Re_tau = delta99/delta_nu", st["Re_tau"], np.nan, ref["Re_tau"]),
        ("U_c+ rms dev from log law", log_law_deviation(st), 0.0,
         log_law_deviation(ref)),
        ("kappa (diagnostic minimum)", k_dns, REF["kappa"], k_ref),
        ("C (diagnostic minimum)", C_dns, REF["C_log"], C_ref),
        ("p'_rms/(rho_w u_tau^2) wall", st["prms"][0] / q, REF["p_rms_wall"], np.nan),
        ("p'_rms/(rho_w u_tau^2) peak", st["prms"].max() / q, REF["p_rms_peak"], np.nan),
        ("p'_rms/(rho_w u_tau^2) freestream", freestream_prms(st) / q, REF["p_rms_inf"], np.nan),
    ]
    lines = [f"{'quantity':36s} {'DNS':>11s} {'Guarini':>11s} {'diff':>8s}  "
             f"{'same est. on ref profile':>24s}"]
    for name, got, paper, refval in rows:
        d = "" if not np.isfinite(paper) or paper == 0 else f"{100 * (got - paper) / paper:7.1f}%"
        rv = "" if not np.isfinite(refval) else f"{refval:11.5g}"
        lines.append(f"{name:36s} {got:11.5g} {paper:11.5g} {d:>8s}  {rv:>24s}")
    lines += [
        "",
        f"{'delta99 [mm]':36s} {st['d99'] * 1e3:11.5g}",
        f"{'u_tau [m/s]':36s} {st['u_tau']:11.5g}",
        f"{'wall unit [um]':36s} {st['d_nu'] * 1e6:11.5g}",
        f"{'T_w [K] (adiabatic value 276.97)':36s} {st['T_w']:11.5g}",
        f"{'first cell y+':36s} {st['yp'][1]:11.5g}",
    ]
    return "\n".join(lines), {r[0]: (r[1], r[2], r[3]) for r in rows}


def plot(st, sw=None, out="tbl_stats.png"):
    """Multi-panel figure mirroring Guarini et al.'s figures 4-11."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d99, ut, rw = st["d99"], st["u_tau"], st["rho_w"]
    y_d, yp = st["y"] / d99, st["yp"]
    rr = np.sqrt(st["rho"] / rw)
    inb = y_d <= 1.2

    fig, ax = plt.subplots(2, 4, figsize=(23, 10.5))

    a = ax[0, 0]
    if sw is not None:
        ok = np.isfinite(sw["Cf"])
        a.plot(sw["Re_theta"][ok], sw["Cf"][ok], "-", label="DNS, streamwise sweep")
    a.plot(st["Re_theta"], st["Cf"], "r*", ms=16, label="sampling station")
    a.plot(REF["Re_theta"], REF["Cf"], "ko", ms=10, mfc="none", label="Guarini et al.")
    a.set(xlabel=r"$Re_\theta$", ylabel=r"$C_f$", title="fig 4: skin friction")
    a.legend(); a.grid(alpha=.3)

    a = ax[0, 1]
    m = yp > 0
    a.semilogx(yp[m], st["ucp"][m], "-", lw=2, label=r"$U_c^+$ van Driest")
    a.semilogx(yp[m], st["up"][m], "--", lw=1, label=r"$U^+$ untransformed")
    yl = np.logspace(-1, 0.7, 20)
    a.semilogx(yl, yl, "k:", label=r"$U_c^+=y^+$")
    yl = np.logspace(1, np.log10(max(yp.max(), 20)), 20)
    a.semilogx(yl, np.log(yl) / REF["kappa"] + REF["C_log"], "k-.",
               label=r"$\frac{1}{0.40}\ln y^+ + 4.7$")
    a.set(xlabel=r"$y u_\tau/\nu_w$", ylabel=r"$U_c/u_\tau$", ylim=(0, 25),
          title="fig 5: van Driest transformed velocity")
    a.legend(fontsize=9); a.grid(alpha=.3)

    for a, scale, ttl in (
        (ax[0, 2], np.ones_like(rr), r"fig 6a: intensities $\sqrt{u_i'^2}/u_\tau$"),
        (ax[0, 3], rr, r"fig 6b: $\sqrt{\bar\rho/\rho_w}\,\sqrt{u_i'^2}/u_\tau$"),
    ):
        for q, lab in (("urms", "$u'$"), ("vrms", "$v'$"), ("wrms", "$w'$")):
            a.plot(y_d[inb], (st[q] * scale / ut)[inb], label=lab)
        a.set(xlabel=r"$y/\delta$", ylabel="rms", title=ttl)
        a.legend(); a.grid(alpha=.3)

    a = ax[1, 0]
    mean_stress = sutherland_mu(st["T"]) * np.gradient(st["u_fav"], st["y"])
    a.plot(yp[inb], (mean_stress / st["tau_w"])[inb], label="mean shear")
    a.plot(yp[inb], (-st["Ruv"] / st["tau_w"])[inb], label=r"Reynolds $-\rho\overline{u''v''}$")
    a.plot(yp[inb], ((mean_stress - st["Ruv"]) / st["tau_w"])[inb], "k-", lw=2, label="total")
    a.set(xlabel=r"$y u_\tau/\nu_w$", ylabel=r"$\tau/\tau_w$", xlim=(0, 200), ylim=(0, 1.2),
          title="fig 9: stress balance")
    a.legend(); a.grid(alpha=.3)

    a = ax[1, 1]
    a.plot(y_d[inb], (st["uv"] / ut**2)[inb], label=r"$\overline{u'v'}/u_\tau^2$")
    a.plot(y_d[inb], (st["uv"] * st["rho"] / rw / ut**2)[inb],
           label=r"$(\bar\rho/\rho_w)\overline{u'v'}/u_\tau^2$")
    a.set(xlabel=r"$y/\delta$", ylabel=r"$\overline{u'v'}$",
          title="fig 10: Reynolds shear stress")
    a.legend(); a.grid(alpha=.3)

    a = ax[1, 2]
    a.plot(y_d, st["prms"] / (rw * ut**2), label="DNS")
    for v, lab in ((REF["p_rms_wall"], "wall 2.7"), (REF["p_rms_peak"], "peak 3.0"),
                   (REF["p_rms_inf"], r"freestream 0.47")):
        a.axhline(v, ls=":", c="k")
        a.annotate(lab, (1.35, v * 1.02), fontsize=8)
    a.set(xlabel=r"$y/\delta$", ylabel=r"$p'_{rms}/(\rho_w u_\tau^2)$", xlim=(0, 2),
          title="fig 7: pressure fluctuations")
    a.legend(); a.grid(alpha=.3)

    a = ax[1, 3]
    s = st["nu_w"] / st["u_tau"] ** 2
    for q, lab in (("oxrms", r"$\omega_1$"), ("oyrms", r"$\omega_2$"), ("ozrms", r"$\omega_3$")):
        a.plot(yp, st[q] * s, label=lab)
    a.set(xlabel=r"$y u_\tau/\nu_w$", ylabel=r"$\omega'_{rms}\nu_w/u_\tau^2$",
          xlim=(0, 50), ylim=(0, 0.45), title="fig 8: rms vorticity")
    a.legend(); a.grid(alpha=.3)

    fig.suptitle(
        f"3D_solver/TBL vs Guarini et al. (2000)   |   "
        f"x = {st['x'] * 1e3:.2f} mm,  $Re_\\theta$ = {st['Re_theta']:.0f},  "
        f"$Re_\\tau$ = {st['Re_tau']:.0f},  $C_f$ = {st['Cf']:.5f}", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    return fig


def plot_thermal(st, out="tbl_thermal.png"):
    """Temperature fluctuations and the strong Reynolds analogy (fig 11, sec 4)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y_d = st["y"] / st["d99"]
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))

    a = ax[0]
    a.plot(y_d, st["Trms"] / st["T"], label=r"$T'_{rms}/\bar T$")
    a.plot(y_d, st["Ttrms"] / st["Tt"], label=r"$T'_{t,rms}/\bar T_t$")
    a.plot(y_d, st["Ttrms"] / st["T"], label=r"$T'_{t,rms}/\bar T$")
    a.set(xlabel=r"$y/\delta$", ylabel="rms / mean", xlim=(0, 1.6),
          title="fig 11: temperature fluctuations")
    a.legend(); a.grid(alpha=.3)

    a = ax[1]
    a.plot(y_d, st["sra"], label="DNS")
    a.axhline(1.0, ls=":", c="k", label="SRA prediction (4.9a)")
    a.set(xlabel=r"$y/\delta$", xlim=(0, 1.2), ylim=(0, 3),
          ylabel=r"$\frac{\sqrt{T''^2}/\tilde T}{(\gamma-1)M_a^2\sqrt{u''^2}/\tilde u}$",
          title="strong Reynolds analogy")
    a.legend(); a.grid(alpha=.3)

    a = ax[2]
    a.plot(y_d, -st["R_uT"], label=r"$-R_{u''T''}$")
    a.plot(y_d, st["Pr_t"], label=r"$Pr_t$")
    a.axhline(1.0, ls=":", c="k")
    a.axhline(0.9, ls="--", c="gray", label="$Pr_t=0.9$")
    a.set(xlabel=r"$y/\delta$", xlim=(0, 1.2), ylim=(-1.5, 2),
          title="correlation and turbulent Prandtl number")
    a.legend(); a.grid(alpha=.3)

    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    return fig


def convergence_check(acc_a, acc_b, ix, half_width=0):
    """Compare two half-sample accumulations at the same station."""
    sa, sb = station(acc_a, ix, half_width), station(acc_b, ix, half_width)
    lines = [f"{'quantity':22s} {'first half':>12s} {'second half':>12s} {'diff':>8s}"]
    for k in ("Cf", "Re_theta", "u_tau", "d99", "H"):
        d = 100 * (sb[k] - sa[k]) / sa[k]
        lines.append(f"{k:22s} {sa[k]:12.5g} {sb[k]:12.5g} {d:7.2f}%")
    for lab, fn in (
        ("peak u'/u_tau", lambda s: s["urms"].max() / s["u_tau"]),
        ("peak p'/(rho_w u_tau^2)", lambda s: s["prms"].max() / (s["rho_w"] * s["u_tau"] ** 2)),
        ("peak -rho<u''v''>/tau_w", lambda s: (-s["Ruv"] / s["tau_w"]).max()),
    ):
        pa, pb = fn(sa), fn(sb)
        lines.append(f"{lab:22s} {pa:12.5g} {pb:12.5g} {100 * (pb - pa) / pa:7.2f}%")
    return "\n".join(lines)


def _opt(argv, flag, default):
    if flag in argv:
        return argv[argv.index(flag) + 1]
    return default


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 1
    cmd = argv.pop(0)

    if cmd == "accumulate":
        data_dir = argv.pop(0)
        skip = int(_opt(argv, "--skip", 0))
        stride = int(_opt(argv, "--stride", 1))
        out = _opt(argv, "--out", "tbl_acc.npz")
        accumulate(data_dir, skip=skip, stride=stride, out=out)
        if "--halves" in argv:
            # split-half convergence: the honest test of whether the sampling
            # window was long enough.
            n = len(glob.glob(os.path.join(str(data_dir), "Q*.vtr")))
            mid = skip + ((n - skip) // 2)
            base = out[:-4] if out.endswith(".npz") else out
            accumulate(data_dir, skip=skip, stop=mid, stride=stride,
                       out=f"{base}_h1.npz", verbose=False)
            accumulate(data_dir, skip=mid, stride=stride,
                       out=f"{base}_h2.npz", verbose=False)
        return 0

    if cmd == "report":
        path = argv.pop(0)
        acc = dict(np.load(path))
        target = float(_opt(argv, "--re-theta", REF["Re_theta"]))
        hw = int(_opt(argv, "--half-width", 0))
        ix, sw = find_station(acc, target)
        st = station(acc, ix, hw)
        print(f"sampling station ix={ix}  x={st['x'] * 1e3:.3f} mm "
              f"= {st['x'] / st['d99']:.2f} delta from the inlet, "
              f"window +/-{hw} columns")
        print(f"snapshots={int(acc['nfiles'])}  spanwise-time samples={int(acc['nsamp'])}")
        print()
        text, _ = summary(st)
        print(text)
        out = _opt(argv, "--out", "tbl_stats.png")
        plot(st, sw, out=out)
        plot_thermal(st, out=out.replace(".png", "_thermal.png"))

        base = path[:-4] if path.endswith(".npz") else path
        h1, h2 = f"{base}_h1.npz", f"{base}_h2.npz"
        if os.path.exists(h1) and os.path.exists(h2):
            print()
            print("split-half convergence (was the sampling window long enough?)")
            print(convergence_check(dict(np.load(h1)), dict(np.load(h2)), ix, hw))
        return 0

    print(f"unknown command {cmd!r}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
