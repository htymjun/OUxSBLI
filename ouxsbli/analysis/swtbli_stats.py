"""Post-processing for 3D_solver/SWTBLI (incident shock on a turbulent boundary
layer, two coupled domains). numpy + matplotlib only; .vtr files are read with
:mod:`ouxsbli.analysis.vtr_raw`.

A *stage directory* is what ``3D_solver/SWTBLI/run_stage.sh`` leaves behind,
``runs/<stage>/`` with ``mod_globals.f90`` (the parameters it ran with) and
``data/1/Q*.vtr`` (TBL generator) and ``data/3/Q*.vtr`` (SBLI domain).

Commands (run from the repository root)::

    python -m ouxsbli.analysis.swtbli_stats turb   RUNDIR  [--last-du 60] [--out DIR]
    python -m ouxsbli.analysis.swtbli_stats track  RUNDIR [RUNDIR ...] [--ref NOSHOCK] [--prev TRACK_NPZ] [--out DIR]
    python -m ouxsbli.analysis.swtbli_stats report RUNDIR  [--last-du 100] [--ref NOSHOCK_RUNDIR] [--out DIR]
    python -m ouxsbli.analysis.swtbli_stats zcorr  RUNDIR  [--out DIR]
    python -m ouxsbli.analysis.swtbli_stats wall accumulate RUNDIR [--last-du T] [--out DIR]
    python -m ouxsbli.analysis.swtbli_stats wall report NPZ [NPZ ...] [--compare NPZ[,NPZ]] [--ref-csv CSV] [--out DIR]
    python -m ouxsbli.analysis.swtbli_stats cache-ref RUNDIR [--last-du 100]
    python -m ouxsbli.analysis.swtbli_stats gif    RUNDIR  [--frames 150] [--out DIR]

``turb``   turbulence check of the TBL generator (mean/fluctuation profiles in
           wall units, Cf against Van Driest II, spanwise spectrum, histories).
``track``  per-snapshot separation-shock position, wall-pressure rise, separation
           and reattachment, bubble area, plus the statistical-stationarity test.
           Several stage directories are concatenated in time.
``report`` boundary diagnostics: junction continuity, inlet, top reflection,
           post-shock oscillations, outlet.
``zcorr``  correlation between z and z + Lz/2 (decorrelation of the tiled 4 delta field).
``wall``   ``accumulate`` sums the span-mean wall pressure, tau_w and mean fields of one
           stage into wall_acc.npz; ``report`` merges such files and compares the mean
           C_f and p_w with the free-interaction theory, the inviscid p3, Van Driest II,
           other runs (--compare) and digitised literature data (--ref-csv).
``cache-ref`` stores the mean rows of a no-shock reference run (--ref) as an npz.
``gif``    numerical schlieren, u with the u<0 region, u' at y+~15, and wall
           pressure with the tracked quantities, as GIF animations.
All lengths are printed in units of the inlet delta99 (``blt``) and times in
delta/u0.
"""

import glob
import json
import math
import os
import re
import sys

import numpy as np

from .vtr_raw import read_vtr

GAMMA = 1.4
R_GAS = 287.03
RF = 0.89
# Sutherland law exactly as src/mod_constant.f90.fypp: 1.716e-5*383.6*273.2**-1.5
SUTH_C = 1.716e-5 * 383.6 / 273.2**1.5
SUTH_S = 110.4
NGHOST = 3  # periodic ghost planes on each z end


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def mu_suth(T):
    return SUTH_C * np.asarray(T) ** 1.5 / (np.asarray(T) + SUTH_S)


# ----------------------------------------------------------------------------
# case parameters
# ----------------------------------------------------------------------------
def _fnum(s):
    return float(s.lower().replace("d", "e"))


def _beta_weak(M, th):
    lo, hi = math.asin(1.0 / M) + 1e-10, math.radians(64.0)

    def f(b):
        return math.tan(th) - 2.0 / math.tan(b) * (M * M * math.sin(b) ** 2 - 1.0) / (
            M * M * (GAMMA + math.cos(2.0 * b)) + 2.0)

    flo = f(lo)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0.0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def _oblique(M, b, th, rho, p):
    mn2 = (M * math.sin(b)) ** 2
    p2 = p * (1 + 2 * GAMMA / (GAMMA + 1) * (mn2 - 1))
    rho2 = rho * (GAMMA + 1) * mn2 / ((GAMMA - 1) * mn2 + 2)
    mnb = math.sqrt((1 + 0.5 * (GAMMA - 1) * mn2) / (GAMMA * mn2 - 0.5 * (GAMMA - 1)))
    return rho2, p2, mnb / math.sin(b - th)


def case_params(rundir):
    """Parameters of one stage, parsed from the mod_globals.f90 it ran with."""
    src = open(os.path.join(rundir, "mod_globals.f90")).read()

    def get(name):
        m = re.search(r"::\s*" + name + r"\s*=\s*([-+0-9.dDeE]+)", src)
        if not m:
            raise ValueError(f"{name} not found in {rundir}/mod_globals.f90")
        return _fnum(m.group(1))

    def opt(name, default):
        try:
            return get(name)
        except ValueError:
            return default

    P = {k: get(k) for k in ("blt", "M0", "p_tot", "T_tot", "theta_deg", "ximp_in",
                             "endT", "dt", "dnu0")}
    P["np"] = int(get("np"))
    P["step_offset"] = int(get("step_offset"))
    P["nre1"], P["nre2"] = int(get("nre1")), int(get("nre2"))
    for k in ("nx1", "ny1", "nz1", "nx2", "ny2", "nbuf"):
        P[k] = int(get(k))
    # grid generation 1 (A1..B0) had no gy_top/xdown parameters: 5 %/cell, 18 delta
    P["gy_top"] = opt("gy_top", 1.05)
    P["xdown"] = opt("xdown", 18.0)
    m = re.search(r"dyw0\s*=\s*([0-9.dDeE]+)\s*\*\s*blt", src)
    P["dyw0"] = _fnum(m.group(1)) * P["blt"]
    m = re.search(r"Lx1\s*=\s*([0-9.dDeE]+)\s*\*\s*blt", src)
    P["dx0"] = _fnum(m.group(1)) * P["blt"] / (P["nx1"] - 1)
    m = re.search(r"y_str\s*=\s*([0-9.dDeE]+)\s*\*\s*blt", src)
    P["y_str"] = _fnum(m.group(1)) * P["blt"]
    m = re.search(r"nin\s*=\s*nre2\s*-\s*(\d+)", src)
    P["nin"] = P["nre2"] - int(m.group(1))
    m = re.search(r"Lz1\s*=\s*([0-9.dDeE]+)\s*\*\s*blt", src)
    P["Lz"] = _fnum(m.group(1)) * P["blt"]
    M0 = P["M0"]
    fac = 1 + 0.5 * (GAMMA - 1) * M0**2
    P["T0"] = P["T_tot"] / fac
    P["p0"] = P["p_tot"] / fac ** (GAMMA / (GAMMA - 1))
    P["rho0"] = P["p0"] / (R_GAS * P["T0"])
    P["u0"] = M0 * math.sqrt(GAMMA * R_GAS * P["T0"])
    P["mu0"] = float(mu_suth(P["T0"]))
    P["Taw"] = P["T0"] * (1 + RF * 0.5 * (GAMMA - 1) * M0**2)
    P["tU"] = P["blt"] / P["u0"]
    P["nt"] = int(P["endT"] / (P["np"] * P["dt"]) + 1e-9)
    P["dt_out"] = P["nt"] * P["dt"]
    th = math.radians(P["theta_deg"])
    if th > 0:
        b1 = _beta_weak(M0, th)
        rho2, p2, M2 = _oblique(M0, b1, th, P["rho0"], P["p0"])
        br = _beta_weak(M2, th)
        rho3, p3, M3 = _oblique(M2, br, th, rho2, p2)
    else:
        b1, rho2, p2, rho3, p3 = math.asin(1 / M0), P["rho0"], P["p0"], P["rho0"], P["p0"]
        M2 = M3 = M0
        br = b1
    # beta_r: reflected-shock angle relative to the flow behind the incident shock;
    # beta_r - theta is its angle to the wall
    P.update(beta1=b1, rho2=rho2, p2=p2, rho3=rho3, p3=p3, M2=M2, M3=M3, beta_r=br)
    return P


def make_y(P, ny):
    """Wall-normal grid of set.f90 (make_y), shared by both domains."""
    y = np.zeros(ny)
    dy = P["dyw0"]
    for j in range(1, ny):
        y[j] = y[j - 1] + dy
        dy = min(dy * 1.035, P["dx0"]) if y[j] < P["y_str"] else dy * P["gy_top"]
    return y


def make_x(P, domain):
    """x grid of set.f90 (set_grid): domain 1 = TBL generator, 3 = SBLI."""
    dx0 = P["dx0"]
    if domain == 1:
        return np.arange(P["nx1"]) * dx0
    nx = P["nx2"]
    nu = nx - 1 - P["nbuf"]
    x = np.zeros(nx)
    x[0] = (P["nin"] - 1) * dx0
    x[1:nu + 1] = x[0] + np.arange(1, nu + 1) * dx0
    d = dx0
    for i in range(nu + 1, nx):
        d = min(d * 1.03, 4 * dx0)
        x[i] = x[i - 1] + d
    return x


def snapshots(rundir, rank):
    return sorted(glob.glob(os.path.join(rundir, "data", str(rank), "Q*.vtr")))


def snap_index(path):
    return int(re.search(r"Q(\d+)\.vtr$", path).group(1))


def load(path):
    """Snapshot as float64 arrays shaped (nz_interior, ny, nx); ghost planes dropped."""
    x, y, z, f = read_vtr(path)
    s = slice(NGHOST, -NGHOST)
    vel = f["velocity"]
    return {
        "x": x.astype(float), "y": y.astype(float), "z": z[s].astype(float),
        "rho": f["rho"][s].astype(float), "p": f["p"][s].astype(float),
        "u": vel[s, :, :, 0].astype(float), "v": vel[s, :, :, 1].astype(float),
        "w": vel[s, :, :, 2].astype(float),
    }


def geometry(P, sbli_x, sbli_y):
    """x_in, x_imp, x_top, x_dyn of the SBLI domain (as set.f90 computes them)."""
    blt = P["blt"]
    x_in = sbli_x[0]
    x_imp = x_in + P["ximp_in"] * blt
    x_top = x_imp - sbli_y[-1] / math.tan(P["beta1"])
    return {"x_in": x_in, "x_imp": x_imp, "x_top": x_top, "x_dyn": x_top + 3 * blt}


# ----------------------------------------------------------------------------
# wall and integral quantities
# ----------------------------------------------------------------------------
def _d3(y0, y1, y2, f0, f1, f2):
    d0 = (y0 - y1) * (y0 - y2)
    d1 = (y1 - y0) * (y1 - y2)
    d2 = (y2 - y0) * (y2 - y1)
    return f0 * ((y0 - y1) + (y0 - y2)) / d0 + f1 * (y0 - y2) / d1 + f2 * (y0 - y1) / d2


def wall_tau(y, u2, T_w):
    """tau_w(x) from a (ny, nx) mean u and the wall temperature row."""
    dudy = _d3(y[0], y[1], y[2], u2[0], u2[1], u2[2])
    return mu_suth(T_w) * dudy


def cf_vd2(Re_theta, M0, Tw_Te, Te):
    """Van Driest II skin friction with Karman-Schoenherr; Te = edge temperature [K]."""
    r = RF
    A2 = r * 0.5 * (GAMMA - 1) * M0**2 / Tw_Te
    B = (1 + r * 0.5 * (GAMMA - 1) * M0**2) / Tw_Te - 1
    den = math.sqrt(B * B + 4 * A2)
    al, be = (2 * A2 - B) / den, B / den
    Fc = (Tw_Te - 1) / (math.asin(al) + math.asin(be)) ** 2
    F_rt = float(mu_suth(Te) / mu_suth(Tw_Te * Te))  # mu_e / mu_w
    L = math.log10(F_rt * Re_theta)
    cf_inc = 1.0 / (17.08 * L * L + 25.11 * L + 6.012)
    return cf_inc / Fc


def profile_integrals(y, rho, u, P):
    """delta99, delta*, theta, H from mean profiles (edge = local max of u)."""
    ue = u.max()
    j99 = int(np.argmax(u >= 0.99 * ue))
    if j99 == 0:
        return dict(d99=np.nan, dstar=np.nan, theta=np.nan, H=np.nan, ue=ue)
    f = (0.99 * ue - u[j99 - 1]) / (u[j99] - u[j99 - 1])
    d99 = y[j99 - 1] + f * (y[j99] - y[j99 - 1])
    jj = slice(0, j99 + 1)
    rue = rho[j99] * u[j99]
    dstar = np.trapezoid(1 - rho[jj] * u[jj] / rue, y[jj])
    theta = np.trapezoid(rho[jj] * u[jj] / rue * (1 - u[jj] / u[j99]), y[jj])
    return dict(d99=d99, dstar=dstar, theta=theta, H=dstar / theta, ue=ue,
                rho_e=rho[j99], u_e=u[j99])


# ----------------------------------------------------------------------------
# per-snapshot tracking (SBLI domain)
# ----------------------------------------------------------------------------
def snap_track(P, s, G, pw_base=None):
    """Tracked quantities of one SBLI snapshot (spanwise averages).

    pw_base: time-mean wall pressure of a no-shock run on the same grid (the
    baseline already sits ~2 % above p0 from the boundary layer's displacement);
    p0 is used when it is None.
    """
    blt, p0 = P["blt"], P["p0"]
    x, y = s["x"], s["y"]
    um = s["u"].mean(0)
    pm = s["p"].mean(0)
    Tw = (s["p"][:, 0, :] / (R_GAS * s["rho"][:, 0, :])).mean(0)
    cf = 2 * wall_tau(y, um, Tw) / (P["rho0"] * P["u0"] ** 2)
    pw = pm[0]
    out = {}
    lo = x >= G["x_in"] + 2 * blt
    hi_sep = x <= G["x_imp"] + (P["xdown"] - 3) * blt
    neg = np.where(lo & hi_sep & (cf < 0))[0]
    out["x_sep"] = x[neg[0]] if neg.size else np.nan
    out["x_reat"] = x[neg[-1] + 1] if neg.size else np.nan
    out["L_sep"] = out["x_reat"] - out["x_sep"]
    if P["theta_deg"] > 0:
        # beginning of the interaction: x-smoothed (+-0.3 delta) wall pressure rises
        # by 10 % of (p3 - p0) above the baseline and stays there for 1 delta
        base = pw_base[0] if pw_base is not None else np.full_like(pw, p0)
        nsm = max(1, int(round(0.3 * blt / (x[1] - x[0]))))
        sm = np.convolve(pw - base, np.ones(2 * nsm + 1) / (2 * nsm + 1), mode="same")
        above = lo & (sm >= 0.1 * (P["p3"] - p0))
        n1 = int(round(1.0 * blt / (x[1] - x[0])))
        out["x_0"] = np.nan
        for i in np.where(above)[0]:
            if above[i:i + n1].all():
                out["x_0"] = x[i]
                break
        # separation shock at y = 1.5 delta (just above the incoming layer), upstream
        # of the incident shock's crossing: first point where the x-smoothed pressure
        # exceeds the no-shock baseline by 15 % of (p3 - p0) for 0.5 delta
        j = int(np.argmin(abs(y - 1.5 * blt)))
        x_inc = G["x_imp"] - y[j] / math.tan(P["beta1"])
        base_s = pw_base[2] if pw_base is not None else np.full_like(pw, p0)
        sms = np.convolve(pm[j] - base_s, np.ones(2 * nsm + 1) / (2 * nsm + 1), mode="same")
        ab = lo & (x <= x_inc - 0.5 * blt) & (sms >= 0.15 * (P["p3"] - p0))
        nh = max(1, n1 // 2)
        out["x_sh"] = np.nan
        for i in np.where(ab)[0]:
            if (sms[i:i + nh] >= 0.15 * (P["p3"] - p0)).all():
                out["x_sh"] = x[i]
                break
    else:
        out["x_0"] = np.nan
        out["x_sh"] = np.nan
    # spanwise-averaged reversed-flow area (delta^2)
    dxv = np.gradient(x)
    dyv = np.gradient(y)
    area = ((s["u"] < 0).mean(0) * dyv[:, None] * dxv[None, :])
    jtop = y < 4 * blt
    out["A_bub"] = area[jtop].sum() / blt**2
    out["cf"] = cf
    out["pw"] = pw
    return out


def stage_times(rundirs):
    """[(rundir, P, [(path_tbl, path_sbli, t_global/tU)])] with restart duplicates skipped."""
    out, t0 = [], 0.0
    for k, rd in enumerate(rundirs):
        P = case_params(rd)
        fs = snapshots(rd, 3)
        ft = {snap_index(f): f for f in snapshots(rd, 1)}
        items = []
        for f in fs:
            n = snap_index(f) - P["step_offset"]
            if k > 0 and n == 0:
                continue
            items.append((ft.get(snap_index(f)), f, (t0 + n * P["dt_out"]) / P["tU"]))
        out.append((rd, P, items))
        t0 += P["endT"]
    return out


def stationarity(t, xs, W, tol):
    """Statistical-stationarity test on a tracked position series (units: delta, delta/U).

    Moving mean over windows of width W; over the last 3W the linear-regression
    displacement must be < tol and smaller than the std of the signal, and the
    last three window means must not move monotonically upstream.
    """
    t = np.asarray(t)
    xs = np.asarray(xs)
    ok = np.isfinite(xs)
    t, xs = t[ok], xs[ok]
    res = {"W": W, "tol": tol}
    if t.size < 6 or t[-1] - t[0] < 3 * W:
        res.update(stationary=False, reason="record shorter than 3W")
        return res
    m = t >= t[-1] - 3 * W
    a, b = np.polyfit(t[m], xs[m], 1)
    disp = a * 3 * W
    sd = float(xs[m].std())
    wm = [float(xs[(t >= t[-1] - (q + 1) * W) & (t < t[-1] - q * W + 1e-9)].mean()) for q in (2, 1, 0)]
    mono_up = wm[0] > wm[1] > wm[2]
    res.update(displacement_3W=float(disp), std=sd, window_means=wm,
               monotonic_upstream=bool(mono_up),
               stationary=bool(abs(disp) < tol and abs(disp) < sd and not mono_up))
    return res


def cmd_track(rundirs, out, ref=None, prev=None):
    """prev: track.npz of the preceding stage(s). Its time series is continued, so
    the stationarity test spans several stages without their VTK files (which may
    have been deleted). The first snapshot of a restart duplicates the previous
    stage's last one and is skipped."""
    os.makedirs(out, exist_ok=True)
    pw_base = None
    rows = []
    prev_arr = dict(np.load(prev)) if prev else None
    t_base = float(prev_arr["t"][-1]) if prev_arr is not None else 0.0
    for rd, P, items in stage_times(rundirs):
        for ft, fs, t in items:
            if prev_arr is not None and t == 0.0:
                continue
            t += t_base
            s = load(fs)
            if ref and (pw_base is None or len(pw_base[0]) != len(s["x"])):
                pw_base = _sbli_mean_rows(ref, 1e9, s["x"])
            G = geometry(P, s["x"], s["y"])
            q = snap_track(P, s, G, pw_base)
            blt = P["blt"]
            rows.append(dict(t=t, stage=os.path.basename(rd.rstrip("/")),
                             x_sh=(q["x_sh"] - G["x_imp"]) / blt, x_0=(q["x_0"] - G["x_imp"]) / blt,
                             x_sep=(q["x_sep"] - G["x_imp"]) / blt, x_reat=(q["x_reat"] - G["x_imp"]) / blt,
                             L_sep=q["L_sep"] / blt, A_bub=q["A_bub"],
                             inlet_gap=(q["x_0"] - G["x_in"]) / blt))
    keys = ["t", "x_sh", "x_0", "x_sep", "x_reat", "L_sep", "A_bub", "inlet_gap"]
    arr = {k: np.array([r[k] for r in rows], float) for k in keys}
    arr["stage"] = np.array([r["stage"] for r in rows])
    if prev_arr is not None:
        for k in keys + ["stage"]:
            arr[k] = np.concatenate([prev_arr[k], arr[k]]) if k in prev_arr else arr[k]
    np.savez(os.path.join(out, "track.npz"), **arr)
    Lmean = np.nanmean(arr["L_sep"][-max(5, len(rows) // 4):]) if np.isfinite(arr["L_sep"]).any() else np.nan
    W = 100.0 if not np.isfinite(Lmean) else max(50.0, 10.0 * Lmean)
    st_sh = stationarity(arr["t"], arr["x_sh"], W, 0.25)
    st_0 = stationarity(arr["t"], arr["x_0"], W, 0.25)
    summ = {"n": len(rows), "t_first": float(arr["t"][0]), "t_last": float(arr["t"][-1]),
            "L_sep_recent_mean": float(Lmean), "x_sh_stationarity": st_sh, "x_0_stationarity": st_0,
            "min_inlet_gap": float(np.nanmin(arr["inlet_gap"])) if np.isfinite(arr["inlet_gap"]).any() else None}
    with open(os.path.join(out, "track.json"), "w") as f:
        json.dump(summ, f, indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float))
    _plot_track(arr, out)
    return arr, summ


def _plot_track(arr, out):
    plt = _plt()
    fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for k, c in (("x_sh", "C0"), ("x_0", "C1"), ("x_sep", "C2"), ("x_reat", "C3")):
        ax[0].plot(arr["t"], arr[k], ".-", ms=2, lw=0.8, color=c, label=k)
    ax[0].set_ylabel(r"$(x-x_{imp})/\delta$")
    ax[0].legend(ncol=4, fontsize=8)
    ax[1].plot(arr["t"], arr["L_sep"], ".-", ms=2, lw=0.8)
    ax[1].set_ylabel(r"$L_{sep}/\delta$")
    ax[2].plot(arr["t"], arr["A_bub"], ".-", ms=2, lw=0.8)
    ax[2].set_ylabel(r"reversed-flow area $/\delta^2$")
    ax[2].set_xlabel(r"$t u_0/\delta$")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "track.png"), dpi=110)
    plt.close(fig)


# ----------------------------------------------------------------------------
# turbulence check of the TBL generator
# ----------------------------------------------------------------------------
def cmd_turb(rundir, out, last_du=60.0, x_lo=6.0, x_hi=9.0):
    plt = _plt()
    os.makedirs(out, exist_ok=True)
    (rd, P, items), = stage_times([rundir])
    blt, u0, rho0, p0 = P["blt"], P["u0"], P["rho0"], P["p0"]
    tl = np.array([t for _, _, t in items])
    hist = {"t": [], "cf": [], "urms_peak": []}
    acc = None
    spec = None
    for ft, fs, t in items:
        s = load(ft)
        x, y = s["x"], s["y"]
        ii = np.where((x >= x_lo * blt) & (x <= x_hi * blt))[0]
        um = s["u"][:, :, ii].mean((0, 2))
        Tw = (s["p"][:, 0, ii] / (R_GAS * s["rho"][:, 0, ii])).mean()
        tau = float(wall_tau(y, um[:, None], np.array([Tw]))[0])
        hist["t"].append(t)
        hist["cf"].append(2 * tau / (rho0 * u0**2))
        rho_w = p0 / (R_GAS * Tw)
        ut = math.sqrt(abs(tau) / rho_w)
        up = s["u"][:, :, ii] - s["u"][:, :, ii].mean(0, keepdims=True)
        urms = np.sqrt((up**2).mean((0, 2)))
        hist["urms_peak"].append(float((np.sqrt(s["rho"][:, :, ii].mean((0, 2)) / rho_w) * urms / ut).max()))
        if t < tl[-1] - last_du:
            continue
        r, u, v, w, p = (s[k][:, :, ii] for k in ("rho", "u", "v", "w", "p"))
        fields = {"r": r, "u": u, "v": v, "w": w, "p": p, "ru": r * u, "uu": u * u, "vv": v * v,
                  "ww": w * w, "uv": u * v, "T": p / (R_GAS * r)}
        if acc is None:
            acc = {k: np.zeros(len(y)) for k in fields}
            nacc = 0
        for k, f in fields.items():
            acc[k] += f.mean((0, 2))
        nacc += 1
        # spanwise spectrum of u' at y+ ~ 15 (wall unit of the composite IC)
        j15 = int(np.argmin(abs(y - 15 * P["dnu0"])))
        uz = u[:, j15, :] - u[:, j15, :].mean(0, keepdims=True)
        E = (abs(np.fft.rfft(uz, axis=0)) ** 2).mean(1)
        spec = E if spec is None else spec + E
    m = {k: v / nacc for k, v in acc.items()}
    y = s["y"]
    Tw = m["T"][0]
    rho_w = m["r"][0]
    tau = float(wall_tau(y, m["u"][:, None], np.array([Tw]))[0])
    ut = math.sqrt(tau / rho_w)
    dnu = float(mu_suth(Tw)) / rho_w / ut
    yp = y / dnu
    urms = np.sqrt(np.maximum(m["uu"] - m["u"] ** 2, 0))
    vrms = np.sqrt(np.maximum(m["vv"] - m["v"] ** 2, 0))
    wrms = np.sqrt(np.maximum(m["ww"] - m["w"] ** 2, 0))
    uv = m["uv"] - m["u"] * m["v"]
    mork = np.sqrt(m["r"] / rho_w) / ut
    # Van Driest transformed velocity
    uvd = np.concatenate([[0], np.cumsum(np.sqrt(0.5 * (m["r"][1:] + m["r"][:-1]) / rho_w) * np.diff(m["u"]))]) / ut
    I = profile_integrals(y, m["r"], m["u"], P)
    Re_theta = rho0 * u0 * I["theta"] / P["mu0"]
    cf = 2 * tau / (rho0 * u0**2)
    cf_ref = cf_vd2(Re_theta, P["M0"], Tw / P["T0"], P["T0"])
    jpk = int(np.argmax(urms * mork))
    shear = (-m["r"] * uv / tau)
    spec_n = spec / spec[1:].sum()
    kz_ic = int(round(P["Lz"] / (0.5 * 2.0 * P["blt"])))  # IC streak wavelength = 1 delta -> mode Lz/1delta
    frac_ic = float(spec_n[kz_ic]) if kz_ic < len(spec_n) else float("nan")
    srt = np.sort(spec_n[1:])[::-1]
    n90 = int(np.searchsorted(np.cumsum(srt), 0.9) + 1)
    th = np.array(hist["t"])
    cfh = np.array(hist["cf"])
    mw = th >= th[-1] - 30.0
    res = {
        "station_x_delta": [x_lo, x_hi], "window_du": last_du, "snapshots_averaged": nacc,
        "Cf": cf, "Cf_VanDriestII": cf_ref, "Cf_rel_dev": cf / cf_ref - 1,
        "Re_theta": Re_theta, "Re_tau": I["d99"] / dnu, "delta99_over_blt": I["d99"] / blt,
        "H": I["H"], "u_tau": ut, "delta_nu_um": dnu * 1e6,
        "urms_plus_peak": float((urms * mork)[jpk]), "urms_peak_yplus": float(yp[jpk]),
        "uv_shear_peak": float(shear[1:].max()),
        "vrms_plus_peak": float((vrms * mork).max()), "wrms_plus_peak": float((wrms * mork).max()),
        "spectrum_fraction_in_IC_mode": frac_ic, "modes_for_90pct_energy": n90, "n_modes": len(spec_n) - 1,
        "Cf_last30_mean": float(cfh[mw].mean()), "Cf_last30_relstd": float(cfh[mw].std() / cfh[mw].mean()),
        "Cf_last30_drift": float(np.polyfit(th[mw], cfh[mw], 1)[0] * 30.0 / cfh[mw].mean()) if mw.sum() > 2 else None,
    }
    res["checks"] = {
        "Cf within 10% of Van Driest II": bool(abs(res["Cf_rel_dev"]) < 0.10),
        "u_rms+ peak 2.5-3.2 at y+ 10-20": bool(2.5 <= res["urms_plus_peak"] <= 3.2 and 10 <= res["urms_peak_yplus"] <= 20),
        "-rho u'v'/tau_w peak 0.6-1.0": bool(0.6 <= res["uv_shear_peak"] <= 1.0),
        "spectrum broadband (IC mode < 50%, >=4 modes for 90%)": bool(frac_ic < 0.5 and n90 >= 4),
        "Cf stable over last 30 dU (rel. drift < 3%)": bool(res["Cf_last30_drift"] is not None and abs(res["Cf_last30_drift"]) < 0.03),
    }
    res["turbulent"] = all(res["checks"].values())
    with open(os.path.join(out, "turb.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax[0, 0].semilogx(yp[1:], uvd[1:], "k-", label="DNS (Van Driest)")
    ax[0, 0].semilogx(yp[1:], yp[1:], "k:", lw=0.8)
    ax[0, 0].semilogx(yp[1:], np.log(yp[1:]) / 0.41 + 5.2, "k--", lw=0.8)
    ax[0, 0].set_ylim(0, 25); ax[0, 0].set_xlabel("$y^+$"); ax[0, 0].set_ylabel(r"$u_{vd}^+$")
    ax[0, 1].plot(yp, urms * mork, label="$u'$"); ax[0, 1].plot(yp, vrms * mork, label="$v'$")
    ax[0, 1].plot(yp, wrms * mork, label="$w'$"); ax[0, 1].plot(yp, shear, label=r"$-\bar\rho u'v'/\tau_w$")
    ax[0, 1].set_xlim(0, yp.max() if yp.max() < 400 else 400); ax[0, 1].legend(); ax[0, 1].set_xlabel("$y^+$")
    ax[1, 0].plot(th, cfh * 1e3, "k-"); ax[1, 0].axhline(cf_ref * 1e3, color="r", ls="--", label="Van Driest II")
    ax[1, 0].set_xlabel(r"$t u_0/\delta$"); ax[1, 0].set_ylabel(r"$C_f\times10^3$ (x=6-9$\delta$)"); ax[1, 0].legend()
    lz = P["Lz"] / np.arange(1, len(spec_n))
    ax[1, 1].loglog(lz / dnu, spec_n[1:], "o-"); ax[1, 1].set_xlabel(r"$\lambda_z^+$")
    ax[1, 1].set_ylabel(r"$E_{uu}(k_z)$ at $y^+\approx15$ (normalised)")
    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.suptitle(f"TBL generator: Re_theta={Re_theta:.0f}, Re_tau={res['Re_tau']:.0f}, Cf={cf:.5f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "turb.png"), dpi=110)
    plt.close(fig)
    _plane_figure(rundir, items[-1], P, os.path.join(out, "turb_planes.png"))
    return res


def _plane_figure(rundir, item, P, path):
    """Instantaneous u' at y+~15 (x-z) and u at mid-span (x-y), both domains."""
    plt = _plt()
    ft, fs, t = item
    a, b = load(ft), load(fs)
    blt = P["blt"]
    fig, ax = plt.subplots(2, 1, figsize=(16, 7))
    j15 = int(np.argmin(abs(a["y"] - 15 * P["dnu0"])))
    for s in (a, b):
        up = s["u"][:, j15, :] - s["u"][:, j15, :].mean(0, keepdims=True)
        ax[0].pcolormesh(s["x"] / blt, s["z"] / blt, up / P["u0"], cmap="RdBu_r", vmin=-0.25, vmax=0.25, shading="auto")
        km = len(s["z"]) // 2
        ax[1].pcolormesh(s["x"] / blt, s["y"] / blt, s["u"][km] / P["u0"], cmap="viridis", vmin=0, vmax=1.05, shading="auto")
    ax[0].set_aspect("equal"); ax[0].set_title(f"u'/u0 at y+~15, t={t:.1f}"); ax[0].set_ylabel("z/δ")
    ax[1].set_aspect("equal"); ax[1].set_ylim(0, 4); ax[1].set_title("u/u0 mid-span"); ax[1].set_ylabel("y/δ")
    ax[1].set_xlabel("x/δ (TBL inlet = 0)")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


# ----------------------------------------------------------------------------
# boundary diagnostics
# ----------------------------------------------------------------------------
def _sbli_mean_rows_raw(rundir, last_du):
    """time/span-mean rows of the SBLI domain over the last last_du, on its own grid:
    dict(x, y_top, pw (wall), pt (top), ps (y = 1.5 delta))."""
    (rd, P, items), = stage_times([rundir])
    tl = np.array([t for _, _, t in items])
    pw = pt = ps = None
    n = 0
    for _, fs, t in items:
        if t < tl[-1] - last_du:
            continue
        b = load(fs)
        j = int(np.argmin(abs(b["y"] - 1.5 * P["blt"])))
        a1, a2, a3 = b["p"][:, 0, :].mean(0), b["p"][:, -1, :].mean(0), b["p"][:, j, :].mean(0)
        pw = a1 if pw is None else pw + a1
        pt = a2 if pt is None else pt + a2
        ps = a3 if ps is None else ps + a3
        n += 1
    return {"x": b["x"], "y_top": float(b["y"][-1]), "pw": pw / n, "pt": pt / n, "ps": ps / n}


def _sbli_mean_rows(ref, last_du, x=None, y_top=None):
    """(pw, pt, ps) of a reference run, interpolated in physical x onto ``x``.

    ``ref`` is a stage directory, or an npz written by ``cache-ref`` (the reference
    run's VTK files may have been deleted). Both domains share the physical x of
    the SBLI inlet, so a run on another grid can serve as the reference. ``pt`` is
    returned only when the top boundary sits at the same height (else None).
    """
    if ref.endswith(".npz"):
        R = dict(np.load(ref))
    else:
        cache = os.path.join(ref, "analysis", "ref_rows.npz")
        R = dict(np.load(cache)) if os.path.exists(cache) and not snapshots(ref, 3) else _sbli_mean_rows_raw(ref, last_du)
    if x is None:
        return R["pw"], R["pt"], R["ps"]
    pt = None
    if y_top is not None and abs(float(R["y_top"]) - y_top) < 1e-3 * y_top:
        pt = np.interp(x, R["x"], R["pt"])
    return np.interp(x, R["x"], R["pw"]), pt, np.interp(x, R["x"], R["ps"])


def cmd_cache_ref(rundir, last_du=100.0):
    """Store the mean reference rows of a (no-shock) run so its data/ can be deleted."""
    R = _sbli_mean_rows_raw(rundir, last_du)
    path = os.path.join(rundir, "analysis", "ref_rows.npz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, last_du=last_du, **R)
    print("reference rows written to", path)


def _gridscale(a, n):
    """max |a - (a smoothed over n points)| away from the ends."""
    sm = np.convolve(a, np.ones(n) / n, mode="same")
    h = n // 2
    return float(np.abs(a - sm)[h:-h].max())


def _post_reflected(P, G, xb, yb, pm, prms):
    """Oscillations behind the reflected (separation) shock in the mean field.

    Per row 3 delta < y < y_str - 1 delta the reflected shock is the steepest rise of
    the (0.3 delta smoothed) mean pressure downstream of the incident shock; the
    region 0.5..2.5 delta behind it is then measured like post_shock: grid-scale
    residual of mean p and of p'rms (deviation from a running mean over >= 0.3
    delta, the same physical width on every grid generation), and the mean p'rms.
    Normalised by the inviscid jump of the reflected shock, p3 - p2.
    """
    blt = P["blt"]
    dx = xb[1] - xb[0]
    nw = max(5, int(math.ceil(0.3 * blt / dx)) | 1)
    dpr = P["p3"] - P["p2"]
    rows = []
    for j in np.where((yb > 3 * blt) & (yb < P["y_str"] - blt))[0]:
        x_inc = G["x_imp"] - yb[j] / math.tan(P["beta1"])
        m = (xb > x_inc + blt) & (xb < G["x_imp"] + (P["xdown"] - 1) * blt)
        if m.sum() < 4 * nw:
            continue
        sm = np.convolve(pm[j], np.ones(nw) / nw, mode="same")
        g = np.gradient(sm, xb)
        g[~m] = -np.inf
        xr = xb[int(np.argmax(g))]
        b = (xb > xr + 0.5 * blt) & (xb < xr + 2.5 * blt)
        if b.sum() < 2 * nw:
            continue
        rows.append((yb[j] / blt, (xr - G["x_imp"]) / blt, _gridscale(pm[j, b], nw) / dpr,
                     _gridscale(prms[j, b], nw) / dpr, float(prms[j, b].mean() / dpr)))
    if not rows:
        return None
    a = np.array(rows)
    return {"rows": len(rows), "y_range_delta": [float(a[0, 0]), float(a[-1, 0])],
            "x_reflected_at_y_lo_hi": [float(a[0, 1]), float(a[-1, 1])],
            "smoothing_points": nw,
            "max_gridscale_p_over_dp_r": float(a[:, 2].max()), "mean_gridscale_p_over_dp_r": float(a[:, 2].mean()),
            "max_gridscale_prms_over_dp_r": float(a[:, 3].max()), "mean_gridscale_prms_over_dp_r": float(a[:, 3].mean()),
            "mean_prms_over_dp_r": float(a[:, 4].mean())}


def _top_reflection_landing(P, G, xb, yb, pw):
    """Where a wave reflected at the top by the reflected/separation shock would come
    back to the wall: from the start of the interaction x_0 (10 % of the inviscid wall
    pressure rise) up to the top along the inviscid reflected-shock angle, back down
    along the Mach angle of state 3. It should land behind the uniform x grid."""
    blt = P["blt"]
    p1 = pw[(xb > G["x_in"] + blt) & (xb < G["x_in"] + 5 * blt)].mean()
    rise = np.where((xb > G["x_in"] + 2 * blt) & (pw > p1 + 0.1 * (P["p3"] - P["p0"])))[0]
    if not rise.size:
        return None
    x0 = xb[rise[0]]
    H = yb[-1]
    up = H / math.tan(P["beta_r"] - math.radians(P["theta_deg"]))
    down = H / math.tan(math.asin(1.0 / P["M3"]))
    land = x0 + up + down
    end = G["x_imp"] + P["xdown"] * blt
    return {"x0_delta": float((x0 - G["x_imp"]) / blt), "height_delta": float(H / blt),
            "x_top_crossing_delta": float((x0 + up - G["x_imp"]) / blt),
            "x_landing_delta": float((land - G["x_imp"]) / blt),
            "uniform_end_delta": P["xdown"], "lands_behind_uniform_end": bool(land > end)}


def cmd_report(rundir, out, last_du=100.0, ref=None):
    plt = _plt()
    os.makedirs(out, exist_ok=True)
    (rd, P, items), = stage_times([rundir])
    blt, u0, p0, rho0 = P["blt"], P["u0"], P["p0"], P["rho0"]
    tl = np.array([t for _, _, t in items])
    use = [it for it in items if it[2] >= tl[-1] - last_du]
    A = B = None
    n = 0
    out_neg_inst = []
    for ft, fs, t in use:
        a, b = load(ft), load(fs)
        fa, fb = ({"p": d["p"].mean(0), "pp": (d["p"] ** 2).mean(0), "u": d["u"].mean(0), "r": d["rho"].mean(0),
                   "Tw": (d["p"][:, 0] / (R_GAS * d["rho"][:, 0])).mean(0)} for d in (a, b))
        out_neg_inst.append(float((b["u"][:, 1:-1, -2] < 0).mean()))
        if A is None:
            A = {k: np.zeros_like(v) for k, v in fa.items()}
            B = {k: np.zeros_like(v) for k, v in fb.items()}
        for k in fa:
            A[k] += fa[k]
        for k in fb:
            B[k] += fb[k]
        n += 1
        xa, ya, xb, yb = a["x"], a["y"], b["x"], b["y"]
    A = {k: v / n for k, v in A.items()}
    B = {k: v / n for k, v in B.items()}
    G = geometry(P, xb, yb)
    prms_a = np.sqrt(np.maximum(A["pp"] - A["p"] ** 2, 0))
    prms_b = np.sqrt(np.maximum(B["pp"] - B["p"] ** 2, 0))
    cfa = 2 * wall_tau(ya, A["u"], A["Tw"]) / (rho0 * u0**2)
    cfb = 2 * wall_tau(yb, B["u"], B["Tw"]) / (rho0 * u0**2)
    R = {"window_du": last_du, "snapshots": n}

    # junction: overlap x in [x_in, TBL end) at the same x
    ov = np.where(xa >= xb[0] - 1e-12)[0][:-1]
    ib = np.array([int(np.argmin(abs(xb - xa[i]))) for i in ov])
    jj = slice(1, len(ya) - 1)
    du = np.abs(B["u"][jj][:, ib] - A["u"][jj][:, ov]).max() / u0
    dp = np.abs(B["p"][jj][:, ib] - A["p"][jj][:, ov]).max() / p0
    dcf = np.abs(cfb[ib] / cfa[ov] - 1).max()
    dprms = np.abs(prms_b[0, ib] / prms_a[0, ov] - 1).max()
    R["junction"] = {"overlap_delta": float((xa[ov[-1]] - xa[ov[0]]) / blt), "max_du_over_u0": float(du),
                     "max_dp_over_p0": float(dp), "max_rel_dCf": float(dcf), "max_rel_dprms_wall": float(dprms),
                     "ok": bool(du < 0.03 and dp < 0.03 and dcf < 0.03 and dprms < 0.03)}
    # inlet of the SBLI domain: first 5 delta
    m5 = (xb >= G["x_in"] + 0.5 * blt) & (xb <= G["x_in"] + 5 * blt)
    R["inlet"] = {"pw_mean_rel_dev_max_vs_p0": float(np.abs(B["p"][0, m5] / p0 - 1).max())}
    if ref is not None:
        # upstream influence of the interaction: compare with the no-shock run (the
        # ~2 % rise over p0 there is the boundary layer's own displacement effect)
        pw_ref, pt_ref, _ = _sbli_mean_rows(ref, last_du, xb, yb[-1])
        dv = np.abs(B["p"][0, m5] / pw_ref[m5] - 1).max()
        R["inlet"]["reference_run"] = ref
        R["inlet"]["pw_rel_dev_max_vs_reference"] = float(dv)
        R["inlet"]["ok"] = bool(dv < 0.005)
    # top boundary
    jt = len(yb) - 1
    if P["theta_deg"] > 0:
        up_ = xb < G["x_top"] - 2 * blt
        mid = (xb > G["x_top"] + 2 * blt) & (xb < G["x_dyn"] + 2 * blt)
    else:  # no shock: the whole top boundary should hold the freestream
        up_ = xb < G["x_imp"] + P["xdown"] * blt
        mid = np.zeros_like(up_)
    if ref is not None and pt_ref is not None:
        R["top_vs_reference_upstream_maxdev"] = float(np.abs(B["p"][jt, up_] / pt_ref[up_] - 1).max())
    R["top"] = {"p_top_over_p0_upstream_maxdev": float(np.abs(B["p"][jt, up_] / p0 - 1).max()),
                "p_top_over_p2_behind_incident_maxdev": float(np.abs(B["p"][jt, mid] / P["p2"] - 1).max()) if mid.any() else None}
    xs_unif = (xb > G["x_in"] + blt) & (xb < G["x_imp"] + (P["xdown"] - 1) * blt)
    prof = prms_b[:, xs_unif].mean(1)
    j95 = int(np.argmin(abs(yb - P["y_str"])))
    # p'rms next to the top over p'rms where the y-stretch begins (y_str; 9.5 delta
    # in grid generation 1, where this key was called prms_ratio_top_to_9p5delta)
    R["top"]["prms_ratio_top_to_ystr"] = float(prof[-2] / prof[j95])
    R["top"]["ok"] = bool((R.get("top_vs_reference_upstream_maxdev", R["top"]["p_top_over_p0_upstream_maxdev"])) < 0.01
                          and (R["top"]["p_top_over_p2_behind_incident_maxdev"] or 0) < 0.02
                          and R["top"]["prms_ratio_top_to_ystr"] <= 1.1)
    # post-incident-shock inviscid region: 0.5..2.5 delta behind the shock, 5 < y < 9 delta
    if P["theta_deg"] > 0:
        resid, devs, rms = [], [], []
        for j in np.where((yb > 5 * blt) & (yb < 9 * blt))[0]:
            xs = G["x_top"] + (yb[-1] - yb[j]) / math.tan(P["beta1"])
            m = (xb > xs + 0.5 * blt) & (xb < xs + 2.5 * blt)
            if m.sum() < 9:
                continue
            pr = B["p"][j, m]
            sm = np.convolve(pr, np.ones(9) / 9, mode="same")
            resid.append(np.abs(pr - sm)[4:-4].max() / P["p2"])
            devs.append(np.abs(pr / P["p2"] - 1).mean())
            rms.append((prms_b[j, m] / P["p2"]).mean())
        R["post_shock"] = {"mean_p_over_p2_dev": float(np.mean(devs)), "prms_over_p2": float(np.mean(rms)),
                           "max_gridscale_residual_over_p2": float(np.max(resid)),
                           "ok": bool(np.mean(rms) < 0.01 and np.max(resid) < 0.005)}
    if P["theta_deg"] > 0:
        R["post_reflected"] = _post_reflected(P, G, xb, yb, B["p"], prms_b)
        R["top_reflection"] = _top_reflection_landing(P, G, xb, yb, B["p"][0])
    # outlet
    i_out = len(xb) - 2
    neg = cfb < 0
    x_unif_end = G["x_imp"] + P["xdown"] * blt
    win = neg & (xb < x_unif_end - blt)
    xr = xb[np.where(win)[0][-1] + 1] if win.any() else np.nan
    R["outlet"] = {"min_mean_u_outlet_over_u0": float(B["u"][1:-1, i_out].min() / u0),
                   "max_instantaneous_reverse_fraction": float(max(out_neg_inst)),
                   "reattach_to_uniform_end_delta": float((x_unif_end - xr) / blt) if np.isfinite(xr) else None,
                   "cf_last5delta_slope_positive_or_attached": bool(
                       (cfb[(xb > x_unif_end - 5 * blt) & (xb < x_unif_end)] > 0).all())}
    R["outlet"]["ok"] = bool(R["outlet"]["min_mean_u_outlet_over_u0"] > 0
                             and (R["outlet"]["reattach_to_uniform_end_delta"] is None
                                  or R["outlet"]["reattach_to_uniform_end_delta"] >= 10))
    R["mean_cf_sbli"] = {"x_over_delta_rel_imp": ((xb - G["x_imp"]) / blt).tolist()[::8], "cf": cfb.tolist()[::8]}
    with open(os.path.join(out, "report.json"), "w") as f:
        json.dump(R, f, indent=1, default=float)
    print(json.dumps({k: v for k, v in R.items() if k != "mean_cf_sbli"}, indent=1, default=float))

    fig, ax = plt.subplots(3, 1, figsize=(12, 11))
    ax[0].plot(xa / blt, cfa * 1e3, "b-", label="TBL generator")
    ax[0].plot(xb / blt, cfb * 1e3, "k-", label="SBLI")
    ax[0].axhline(0, color="gray", lw=0.5); ax[0].set_ylabel(r"$C_f\times10^3$"); ax[0].legend()
    ax2 = ax[0].twinx(); ax2.plot(xb / blt, B["p"][0] / p0, "r-"); ax2.set_ylabel("$p_w/p_0$", color="r")
    c = ax[1].pcolormesh(xb / blt, yb / blt, B["u"] / u0, cmap="viridis", shading="auto")
    ax[1].contour(xb / blt, yb / blt, B["u"], [0.0], colors="w", linewidths=0.8)
    ax[1].set_ylim(0, 6); ax[1].set_aspect("equal"); fig.colorbar(c, ax=ax[1], label="mean u/u0")
    c = ax[2].pcolormesh(xb / blt, yb / blt, prms_b / p0, cmap="magma", vmax=np.percentile(prms_b / p0, 99.5), shading="auto")
    ax[2].set_aspect("equal"); fig.colorbar(c, ax=ax[2], label="p'rms/p0"); ax[2].set_xlabel("x/δ")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "report.png"), dpi=100)
    plt.close(fig)
    return R


# ----------------------------------------------------------------------------
# decorrelation of the tiled field
# ----------------------------------------------------------------------------
def cmd_zcorr(rundir, out):
    os.makedirs(out, exist_ok=True)
    (rd, P, items), = stage_times([rundir])
    rows = []
    for ft, fs, t in items:
        r = {"t": t}
        for name, f in (("tbl", ft), ("sbli", fs)):
            s = load(f)
            j = int(np.argmin(abs(s["y"] - 15 * P["dnu0"])))
            u = s["u"][:, j, :]
            u = u - u.mean(0, keepdims=True)
            h = u.shape[0] // 2
            a, b = u[:h].ravel(), u[h:2 * h].ravel()
            r[name] = float((a * b).mean() / math.sqrt((a * a).mean() * (b * b).mean()))
        rows.append(r)
        print(f"t={t:8.2f}  corr(z, z+Lz/2): TBL {r['tbl']:+.3f}  SBLI {r['sbli']:+.3f}")
    with open(os.path.join(out, "zcorr.json"), "w") as fo:
        json.dump(rows, fo, indent=1)
    return rows


# ----------------------------------------------------------------------------
# wall pressure and skin friction (time and span means)
# ----------------------------------------------------------------------------
# Free-interaction theory (Chapman et al. 1957; turbulent constants of Erdos &
# Pallone 1962, as quoted by Delery & Dussauge / Babinsky & Harvey 2011):
#   (p - p1)/q1 = F sqrt(2 Cf0) / (M1^2 - 1)^(1/4),  q1 = gamma p1 M1^2 / 2
# with F = 4.22 at separation and 6.0 at the plateau. Zukoski (1967): p_p/p1 = 1 + M1/2.
FIT_F_SEP = 4.22
FIT_F_PLATEAU = 6.0


def cmd_wall_accumulate(rundir, out, last_du=None):
    """Sums over the SBLI snapshots of one stage (all, or the last last_du), written to
    out/wall_acc.npz: span-mean wall pressure, its square, tau_w, and the span-mean
    rho, u, p, p^2 fields. Segments are merged by ``wall report``, so the VTK files
    of a statistics segment can be deleted after this has run."""
    os.makedirs(out, exist_ok=True)
    (rd, P, items), = stage_times([rundir])
    tl = np.array([t for _, _, t in items])
    if last_du is not None:
        items = [it for it in items if it[2] >= tl[-1] - last_du]
    if P["step_offset"] and len(items) > 1 and snap_index(items[0][1]) == P["step_offset"]:
        items = items[1:]  # snapshot 0 of a restart repeats the previous stage's last one
    S = None
    for _, fs, t in items:
        b = load(fs)
        Tw = (b["p"][:, 0, :] / (R_GAS * b["rho"][:, 0, :])).mean(0)
        f = {"pw": b["p"][:, 0, :].mean(0), "pw2": (b["p"][:, 0, :] ** 2).mean(0),
             "tau": wall_tau(b["y"], b["u"].mean(0), Tw), "Tw": Tw,
             "r": b["rho"].mean(0), "u": b["u"].mean(0), "p": b["p"].mean(0), "pp": (b["p"] ** 2).mean(0)}
        if S is None:
            S = {k: np.zeros_like(v) for k, v in f.items()}
            x, y = b["x"], b["y"]
        for k in f:
            S[k] += f[k]
    G = geometry(P, x, y)
    meta = {k: P[k] for k in ("blt", "p0", "rho0", "u0", "T0", "mu0", "M0", "theta_deg", "p2", "p3", "xdown", "y_str")}
    path = os.path.join(out, "wall_acc.npz")
    np.savez(path, x=x, y=y, n=len(items), t0=items[0][2], t1=items[-1][2], x_in=G["x_in"], x_imp=G["x_imp"],
             stage=os.path.basename(rundir.rstrip("/")), **meta, **{"sum_" + k: v for k, v in S.items()})
    print(f"wall sums of {len(items)} snapshots (t = {items[0][2]:.1f}..{items[-1][2]:.1f}) -> {path}")
    return path


def wall_stats(npzs):
    """Merge wall_acc.npz files (same grid) into mean distributions and the
    interaction scales. Lengths in delta (= blt), x relative to x_imp."""
    A = [dict(np.load(f)) for f in npzs]
    a0 = A[0]
    for a in A[1:]:
        assert a["x"].shape == a0["x"].shape and np.allclose(a["x"], a0["x"]), "wall_acc files from different grids"
    n = sum(int(a["n"]) for a in A)
    m = {k[4:]: sum(a[k] for a in A) / n for k in a0 if k.startswith("sum_")}
    blt, p0, rho0, u0, M0 = (float(a0[k]) for k in ("blt", "p0", "rho0", "u0", "M0"))
    x, y = a0["x"], a0["y"]
    xi = (x - float(a0["x_imp"])) / blt
    xi_in = (float(a0["x_in"]) - float(a0["x_imp"])) / blt
    xdown = float(a0["xdown"])
    cf = 2 * m["tau"] / (rho0 * u0**2)
    pw = m["pw"]
    prms_w = np.sqrt(np.maximum(m["pw2"] - pw**2, 0))
    up = (xi > xi_in + 1) & (xi < xi_in + 5)
    p1 = float(pw[up].mean())
    p3 = float(a0["p3"]) / p0 * p1          # inviscid p3/p1 applied to the measured p1
    # interaction start: 0.3 delta smoothed p_w exceeds p1 by 10 % of (p3 - p1) for 1 delta
    dx = x[1] - x[0]
    nsm = max(1, int(round(0.3 * blt / dx)))
    sm = np.convolve(pw, np.ones(2 * nsm + 1) / (2 * nsm + 1), mode="same")
    n1 = int(round(blt / dx))
    above = (xi > xi_in + 2) & (sm >= p1 + 0.1 * (p3 - p1))
    i0 = next((i for i in np.where(above)[0] if above[i:i + n1].all()), None)
    xi0 = float(xi[i0]) if i0 is not None else np.nan
    neg = np.where((cf < 0) & (xi > xi_in + 2) & (xi < xdown - 1))[0]
    xi_sep = float(xi[neg[0]]) if neg.size else np.nan
    xi_reat = float(xi[neg[-1] + 1]) if neg.size else np.nan
    # upstream reference station 2 delta before the interaction start
    R = {"stages": [str(a["stage"]) for a in A], "snapshots": n,
         "t_range": [float(min(a["t0"] for a in A)), float(max(a["t1"] for a in A))],
         "theta_deg": float(a0["theta_deg"]), "p1_over_p0": p1 / p0}
    ref = np.isfinite(xi0) and (xi > xi0 - 5) & (xi < xi0 - 2)
    if ref is not False and ref.any():
        Cf0 = float(cf[ref].mean())
        i_r = int(np.argmin(abs(xi - (xi0 - 2))))
        I = profile_integrals(y, m["r"][:, i_r], m["u"][:, i_r], None)
        Re_theta = rho0 * u0 * I["theta"] / float(a0["mu0"])
        cf_ref = cf_vd2(Re_theta, M0, float(m["Tw"][i_r]) / float(a0["T0"]), float(a0["T0"]))
        R["upstream"] = {"xi_station": float(xi0 - 2), "Cf0": Cf0, "Cf_VanDriestII": cf_ref,
                         "Cf0_rel_dev_vs_VD2": Cf0 / cf_ref - 1, "Re_theta": Re_theta,
                         "delta99_over_blt": I["d99"] / blt, "H": I["H"]}
    else:
        Cf0 = float(cf[up].mean())
        R["upstream"] = {"Cf0": Cf0}
    d0 = R["upstream"].get("delta99_over_blt", 1.0)
    q1 = 0.5 * GAMMA * p1 * M0**2
    fit = math.sqrt(2 * Cf0) / (M0**2 - 1) ** 0.25 * q1 / p1
    R["interaction"] = {"xi_0": xi0, "xi_sep": xi_sep, "xi_reat": xi_reat, "L_sep": xi_reat - xi_sep,
                        "L_int": -xi0, "L_int_over_delta0": -xi0 / d0, "L_sep_over_delta0": (xi_reat - xi_sep) / d0,
                        "min_Cf": float(cf[(xi > xi_in + 2) & (xi < xdown)].min())}
    th = {"p_sep_over_p1_FIT": 1 + FIT_F_SEP * fit, "p_plateau_over_p1_FIT": 1 + FIT_F_PLATEAU * fit,
          "p_plateau_over_p1_Zukoski": 1 + 0.5 * M0, "p3_over_p1_inviscid": p3 / p1}
    meas = {}
    if neg.size:
        meas["p_sep_over_p1"] = float(np.interp(xi_sep, xi, pw) / p1)
        # plateau: the flattest point of the smoothed p_w inside the bubble
        b = (xi > xi_sep) & (xi < min(xi_reat, 0.0))
        if b.sum() > 3:
            g = np.gradient(sm, x)
            ii = np.where(b)[0]
            k = int(np.argmin(g[ii]))
            # only an interior minimum of dp_w/dx is a plateau; at the ends of the
            # window the pressure simply keeps rising (no plateau, as for theta = 9.5)
            if 2 <= k < len(ii) - 2:
                meas["xi_plateau"] = float(xi[ii[k]])
                meas["p_plateau_over_p1"] = float(sm[ii[k]] / p1)
                meas["plateau_slope_dp_dx_delta_over_p1"] = float(g[ii[k]] * blt / p1)
            else:
                meas["plateau"] = "none (p_w rises monotonically through the bubble)"
    rec = (xi > xdown - 3) & (xi < xdown - 1)
    meas["p_recovery_over_p1"] = float(pw[rec].mean() / p1)
    meas["Cf_recovery_over_Cf0"] = float(cf[rec].mean() / Cf0)
    meas["max_prms_w_over_p1"] = float(prms_w[(xi > xi_in + 1) & (xi < xdown)].max() / p1)
    meas["prms_w_upstream_over_p1"] = float(prms_w[up].mean() / p1)
    R["theory"] = th
    R["measured"] = meas
    R["ratios"] = {k: meas[k1] / th[k2] for k, k1, k2 in (
        ("p_sep_vs_FIT", "p_sep_over_p1", "p_sep_over_p1_FIT"),
        ("p_plateau_vs_FIT", "p_plateau_over_p1", "p_plateau_over_p1_FIT"),
        ("p_plateau_vs_Zukoski", "p_plateau_over_p1", "p_plateau_over_p1_Zukoski"),
        ("p_recovery_vs_inviscid_p3", "p_recovery_over_p1", "p3_over_p1_inviscid")) if k1 in meas}
    D = {"xi": xi, "cf": cf, "pw_p1": pw / p1, "prms_p1": prms_w / p1, "p3_p1": p3 / p1, "Cf0": Cf0,
         "xstar": xi / -xi0 if np.isfinite(xi0) else xi * np.nan,
         "pstar": (pw - p1) / (p3 - p1), "y": y, "u": m["u"], "blt": blt, "u0": u0, "xdown": xdown, "xi_in": xi_in}
    return R, D


def _read_ref_csv(path):
    """Digitised literature data: lines 'kind,x_star,value[,label]' with kind = p
    ((p_w - p1)/(p3 - p1)) or cf (Cf/Cf0), x_star = (x - x_imp)/L_int; '#' comments."""
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = [c.strip() for c in line.split(",")]
        rows.append((f[0].lower(), float(f[1]), float(f[2]), f[3] if len(f) > 3 else os.path.basename(path)))
    return rows


def cmd_wall_report(npzs, out, compare=(), ref_csv=None):
    """Mean C_f and p_w of one run (several wall_acc.npz segments merged) against the
    free-interaction theory, the inviscid p3, Van Driest II upstream, other runs
    (``compare``: wall_acc.npz paths, comma separated, e.g. theta = 9.5 deg) and
    digitised literature data (``ref_csv``) in the IUSTI normalisation
    x* = (x - x_imp)/L_int, p* = (p_w - p1)/(p3 - p1)."""
    plt = _plt()
    os.makedirs(out, exist_ok=True)
    R, D = wall_stats(npzs)
    others = []
    for c in compare:
        Rc, Dc = wall_stats([c])
        others.append((Rc, Dc))
        R.setdefault("compare", []).append({"file": c, "theta_deg": Rc["theta_deg"], "interaction": Rc["interaction"],
                                            "measured": Rc["measured"]})
    with open(os.path.join(out, "wall.json"), "w") as f:
        json.dump(R, f, indent=1, default=float)
    print(json.dumps(R, indent=1, default=float))

    lab = f"θ={R['theta_deg']:g}° ({'+'.join(R['stages'])}, {R['snapshots']} snapshots)"
    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    m = (D["xi"] > D["xi_in"]) & (D["xi"] < D["xdown"] + 2)
    a = ax[0, 0]
    a.plot(D["xi"][m], D["cf"][m] * 1e3, "k-", label=lab)
    for Rc, Dc in others:
        mc = (Dc["xi"] > Dc["xi_in"]) & (Dc["xi"] < Dc["xdown"] + 2)
        a.plot(Dc["xi"][mc], Dc["cf"][mc] * 1e3, "-", lw=0.9, label=f"θ={Rc['theta_deg']:g}° ({'+'.join(Rc['stages'])})")
    if "Cf_VanDriestII" in R["upstream"]:
        a.axhline(R["upstream"]["Cf_VanDriestII"] * 1e3, color="b", ls="--", lw=0.8,
                  label=f"Van Driest II, Re_θ={R['upstream']['Re_theta']:.0f}")
    a.axhline(0, color="gray", lw=0.5)
    for k, c in (("xi_sep", "g"), ("xi_reat", "m")):
        if np.isfinite(R["interaction"][k]):
            a.axvline(R["interaction"][k], color=c, ls=":", lw=0.8)
    a.set_xlabel("(x - x_imp)/δ"); a.set_ylabel(r"$C_f\times10^3$"); a.legend(fontsize=8); a.grid(alpha=0.3)
    a = ax[0, 1]
    a.plot(D["xi"][m], D["pw_p1"][m], "k-", label=lab)
    a.plot(D["xi"][m], 1 + 10 * D["prms_p1"][m], "k:", lw=0.8, label="1 + 10 p'w,rms/p1")
    for Rc, Dc in others:
        mc = (Dc["xi"] > Dc["xi_in"]) & (Dc["xi"] < Dc["xdown"] + 2)
        a.plot(Dc["xi"][mc], Dc["pw_p1"][mc], "-", lw=0.9, label=f"θ={Rc['theta_deg']:g}°")
    th = R["theory"]
    a.plot([D["xi"][m][0], 0, 0, D["xi"][m][-1]], [1, 1, D["p3_p1"], D["p3_p1"]], "b--", lw=0.8, label="inviscid")
    a.axhline(th["p_sep_over_p1_FIT"], color="g", ls="-.", lw=0.8, label="FIT separation (F=4.22)")
    a.axhline(th["p_plateau_over_p1_FIT"], color="c", ls="-.", lw=0.8, label="FIT plateau (F=6.0)")
    a.axhline(th["p_plateau_over_p1_Zukoski"], color="orange", ls="-.", lw=0.8, label="Zukoski plateau 1+M/2")
    a.set_xlabel("(x - x_imp)/δ"); a.set_ylabel("$p_w/p_1$"); a.legend(fontsize=7); a.grid(alpha=0.3)
    refs = _read_ref_csv(ref_csv) if ref_csv else []
    for a, key, ylab, kind in ((ax[1, 0], "pstar", r"$(p_w-p_1)/(p_3-p_1)$", "p"),
                               (ax[1, 1], "cf", r"$C_f/C_{f0}$", "cf")):
        sc = (lambda d: d["cf"] / d["Cf0"]) if key == "cf" else (lambda d: d["pstar"])
        mm = (D["xstar"] > -2.5) & (D["xstar"] < 3)
        a.plot(D["xstar"][mm], sc(D)[mm], "k-", label=lab)
        for Rc, Dc in others:
            mc = (Dc["xstar"] > -2.5) & (Dc["xstar"] < 3)
            a.plot(Dc["xstar"][mc], sc(Dc)[mc], "-", lw=0.9, label=f"θ={Rc['theta_deg']:g}°")
        for lbl in sorted({r[3] for r in refs if r[0] == kind}):
            pts = np.array([(r[1], r[2]) for r in refs if r[0] == kind and r[3] == lbl])
            a.plot(pts[:, 0], pts[:, 1], "o", ms=3, mfc="none", label=lbl)
        a.axhline(0, color="gray", lw=0.5)
        a.set_xlabel(r"$x^* = (x-x_{imp})/L_{int}$"); a.set_ylabel(ylab); a.legend(fontsize=8); a.grid(alpha=0.3)
    I = R["interaction"]
    fig.suptitle(f"Wall pressure and skin friction: L_int = {I['L_int']:.2f}δ, L_sep = {I['L_sep']:.2f}δ "
                 f"(x_sep = {I['xi_sep']:.2f}, x_reat = {I['xi_reat']:.2f})")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "wall.png"), dpi=110)
    plt.close(fig)
    return R


# ----------------------------------------------------------------------------
# GIF animations
# ----------------------------------------------------------------------------
def cmd_gif(rundir, out, frames=150, track_npz=None):
    plt = _plt()
    from matplotlib.animation import PillowWriter
    os.makedirs(out, exist_ok=True)
    (rd, P, items), = stage_times([rundir])
    stride = max(1, int(math.ceil(len(items) / frames)))
    items = items[::stride]
    blt, u0, p0 = P["blt"], P["u0"], P["p0"]
    tr = dict(np.load(track_npz)) if track_npz and os.path.exists(track_npz) else None
    t_off = 0.0
    if tr is not None and "stage" in tr:
        # the track may span several stages (global time). Within it this stage's
        # first entry is local snapshot 1 (snapshot 0 of a restart duplicates the
        # previous stage's last one), except when it is the first stage tracked.
        mine = tr["stage"] == os.path.basename(rundir.rstrip("/"))
        if mine.any():
            t_first = float(tr["t"][mine].min())
            t_off = t_first if t_first == float(tr["t"].min()) else t_first - P["dt_out"] / P["tU"]

    names = ["schlieren", "u_bubble", "uprime_y15", "wall_pressure"]
    figs, writers = {}, {}
    for nm in names:
        figs[nm] = plt.figure(figsize=(16, {"wall_pressure": 10, "schlieren": 6.5}.get(nm, 5.5)))
        writers[nm] = PillowWriter(fps=8)
        writers[nm].setup(figs[nm], os.path.join(out, f"{nm}.gif"), dpi=72)
    G = None
    run_cf = run_pw = None
    nrun = 0
    for ft, fs, t in items:
        a, b = load(ft), load(fs)
        if G is None:
            G = geometry(P, b["x"], b["y"])
            na = int(np.argmin(abs(a["x"] - b["x"][0])))
            X = np.concatenate([a["x"][:na], b["x"]]) / blt
            Yb = b["y"] / blt
            ny1 = len(a["y"])
        km = len(b["z"]) // 2
        # (1) numerical schlieren, mid-span, x-y
        rb = b["rho"][km]
        ra = np.full((len(b["y"]), na), np.nan)
        ra[:ny1] = a["rho"][km][:, :na]
        R2 = np.concatenate([ra, rb], axis=1)
        gy, gx = np.gradient(R2, b["y"], np.concatenate([a["x"][:na], b["x"]]))
        sch = np.log10(np.sqrt(gx**2 + gy**2) * blt / P["rho0"] + 1e-6)
        f = figs["schlieren"]; f.clf(); ax = f.add_subplot(111)
        ax.pcolormesh(X, Yb, sch, cmap="Greys", vmin=-2.5, vmax=0.5, shading="auto")
        ax.set_aspect("equal"); ax.set_xlabel("x/δ"); ax.set_ylabel("y/δ")
        ax.set_title(f"numerical schlieren log10|∇ρ|δ/ρ0, mid-span, t u0/δ = {t:.1f}")
        f.tight_layout(); writers["schlieren"].grab_frame()
        # (2) u with reversed-flow contour
        ua = np.full((len(b["y"]), na), np.nan)
        ua[:ny1] = a["u"][km][:, :na]
        U2 = np.concatenate([ua, b["u"][km]], axis=1) / u0
        f = figs["u_bubble"]; f.clf(); ax = f.add_subplot(111)
        ax.pcolormesh(X, Yb, U2, cmap="viridis", vmin=-0.2, vmax=1.05, shading="auto")
        ax.contour(X, Yb, np.nan_to_num(U2, nan=1.0), [0.0], colors="r", linewidths=0.8)
        ax.set_ylim(0, 6); ax.set_aspect("equal"); ax.set_xlabel("x/δ"); ax.set_ylabel("y/δ")
        ax.set_title(f"u/u0 mid-span (red: u=0), t u0/δ = {t:.1f}")
        f.tight_layout(); writers["u_bubble"].grab_frame()
        # (3) u' at y+ ~ 15
        j15 = int(np.argmin(abs(b["y"] - 15 * P["dnu0"])))
        upa = a["u"][:, j15, :na] - a["u"][:, j15, :na].mean(0, keepdims=True)
        upb = b["u"][:, j15, :] - b["u"][:, j15, :].mean(0, keepdims=True)
        UP = np.concatenate([upa, upb], axis=1) / u0
        f = figs["uprime_y15"]; f.clf(); ax = f.add_subplot(111)
        ax.pcolormesh(X, b["z"] / blt, UP, cmap="RdBu_r", vmin=-0.25, vmax=0.25, shading="auto")
        ax.set_aspect("equal"); ax.set_xlabel("x/δ"); ax.set_ylabel("z/δ")
        ax.set_title(f"u'/u0 at y+≈15, t u0/δ = {t:.1f}")
        f.tight_layout(); writers["uprime_y15"].grab_frame()
        # (4) wall pressure map, span-mean Cf and p_w (instantaneous and running mean
        # over the frames), tracked quantities
        xi = (b["x"] - G["x_imp"]) / blt
        Tw = (b["p"][:, 0, :] / (R_GAS * b["rho"][:, 0, :])).mean(0)
        cf_i = 2 * wall_tau(b["y"], b["u"].mean(0), Tw) / (P["rho0"] * u0**2)
        pw_i = b["p"][:, 0, :].mean(0) / p0
        run_cf = cf_i if run_cf is None else run_cf + cf_i
        run_pw = pw_i if run_pw is None else run_pw + pw_i
        nrun += 1
        f = figs["wall_pressure"]; f.clf()
        ax = f.add_subplot(311)
        ax.pcolormesh(xi, b["z"] / blt, b["p"][:, 0, :] / p0, cmap="inferno",
                      vmin=0.8, vmax=max(1.2, P["p3"] / p0 * 1.1), shading="auto")
        ax.set_aspect("equal"); ax.set_ylabel("z/δ")
        ax.set_title(f"p_w/p0 (SBLI domain), t u0/δ = {t + t_off:.1f}")
        axc = f.add_subplot(312, sharex=ax)
        axc.plot(xi, cf_i * 1e3, color="0.6", lw=0.7, label="Cf, span mean")
        axc.plot(xi, run_cf / nrun * 1e3, "k-", lw=1.2, label=f"Cf, mean of {nrun} frames")
        axc.axhline(0, color="gray", lw=0.5)
        axc.set_ylim(-4, 7); axc.set_ylabel(r"$C_f\times10^3$"); axc.legend(loc="upper left", fontsize=8)
        axp = axc.twinx()
        axp.plot(xi, pw_i, color="salmon", lw=0.7)
        axp.plot(xi, run_pw / nrun, "r-", lw=1.2)
        axp.axhline(P["p3"] / p0, color="r", ls=":", lw=0.8)
        axp.set_ylim(0.8, P["p3"] / p0 * 1.15); axp.set_ylabel("$p_w/p_0$ (red; dotted: inviscid $p_3$)", color="r")
        axc.set_xlabel("(x - x_imp)/δ")
        ax2 = f.add_subplot(313)
        if tr is not None:
            m = tr["t"] <= t + t_off + 1e-9
            for k, c in (("x_sh", "C0"), ("x_0", "C1"), ("x_sep", "C2"), ("x_reat", "C3")):
                ax2.plot(tr["t"][m], tr[k][m], lw=0.8, color=c, label=k)
            ax2.set_xlim(tr["t"][0], tr["t"][-1])
            ax2.legend(ncol=4, fontsize=8, loc="upper left")
        ax2.set_xlabel("t u0/δ"); ax2.set_ylabel("(x - x_imp)/δ"); ax2.grid(alpha=0.3)
        f.tight_layout(); writers["wall_pressure"].grab_frame()
    for nm in names:
        writers[nm].finish()
        plt.close(figs[nm])
    print("GIFs written to", out, "(", len(items), "frames each )")


# ----------------------------------------------------------------------------
def _opt(argv, flag, default):
    if flag in argv:
        i = argv.index(flag)
        v = argv[i + 1]
        del argv[i:i + 2]
        return v
    return default


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 1
    cmd = argv.pop(0)
    out = _opt(argv, "--out", None)
    if cmd == "turb":
        cmd_turb(argv[0], out or os.path.join(argv[0], "analysis"), float(_opt(argv, "--last-du", 60)))
    elif cmd == "track":
        ref = _opt(argv, "--ref", None)
        prev = _opt(argv, "--prev", None)
        cmd_track(argv, out or os.path.join(argv[-1], "analysis"), ref, prev)
    elif cmd == "report":
        ref = _opt(argv, "--ref", None)
        cmd_report(argv[0], out or os.path.join(argv[0], "analysis"), float(_opt(argv, "--last-du", 100)), ref)
    elif cmd == "wall":
        sub = argv.pop(0)
        if sub == "accumulate":
            last = _opt(argv, "--last-du", None)
            cmd_wall_accumulate(argv[0], out or os.path.join(argv[0], "analysis"), float(last) if last else None)
        elif sub == "report":
            cmp_ = _opt(argv, "--compare", None)
            csv = _opt(argv, "--ref-csv", None)
            cmd_wall_report(argv, out or os.path.dirname(os.path.abspath(argv[-1])),
                            cmp_.split(",") if cmp_ else [], csv)
        else:
            print(f"unknown wall subcommand {sub!r}")
            return 1
    elif cmd == "cache-ref":
        cmd_cache_ref(argv[0], float(_opt(argv, "--last-du", 100)))
    elif cmd == "zcorr":
        cmd_zcorr(argv[0], out or os.path.join(argv[0], "analysis"))
    elif cmd == "gif":
        frames = int(_opt(argv, "--frames", 150))
        trk = _opt(argv, "--track", None)
        cmd_gif(argv[0], out or os.path.join(argv[0], "gif"), frames, trk)
    else:
        print(f"unknown command {cmd!r}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
