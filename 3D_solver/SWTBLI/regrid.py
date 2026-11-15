#!/usr/bin/env python3
"""Interpolate SWTBLI restart files onto another grid generation.

    python3 regrid.py OLD_RUNDIR NEW_MOD_GLOBALS DST_RECAL

OLD_RUNDIR is a stage directory (``runs/B0``) with ``mod_globals.f90`` and
``recal/Q0000{1,2}.dat`` + ``Qm.dat``. Its grid is rebuilt from the parameters
with ``swtbli_stats.make_x``/``make_y``, the new one from NEW_MOD_GLOBALS. The
restart files are raw streams of the *physical* conservative variables
Q(nx, ny, nz, 5) (src/sbli.f90 multiplies by the Jacobian before writing), and
Qm.dat holds the five primitive mean profiles (rho, u, v, w, p) over ny1. The
span (nz, Lz) must be the same; only x and y are interpolated (linear).

* TBL generator: the same physical box, plain (x, y) interpolation.
* SBLI domain, in xi = x - x_imp:
  * xi >= XI_KEEP: the old field at the same xi (interaction, recovery, outlet;
    beyond the old outlet the last old column is repeated).
  * the longer inflow region [xi_in_new, XI_KEEP] is mapped linearly onto the old
    one [xi_in_old, XI_KEEP], i.e. the upstream boundary layer is stretched in x
    (factor 1.46 for 25 -> 31 delta) instead of being repeated: the field stays
    continuous, so no artificial front is started at a seam.
  * above the old top (y > y_top_old): the inviscid three-state field of the old
    run, with states 1/2/3 measured on the old top row (span/x means away from
    the shocks) and the incident / reflected shocks continued as straight lines
    (beta1 and beta_r - theta of the old theta) from where they cross the old top.
    The first 1 delta above the old top blends from the old top row into it.
"""
import math
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from ouxsbli.analysis.swtbli_stats import (GAMMA, R_GAS, _beta_weak, _oblique, case_params,  # noqa: E402
                                           make_x, make_y, wall_tau)

XI_KEEP = -12.0  # [delta] identity mapping downstream of this (interaction starts at ~ -7)


def read_q(path, nx, ny, nz):
    q = np.fromfile(path, dtype="<f8")
    assert q.size == nx * ny * nz * 5, f"{path}: {q.size} values, expected {nx}*{ny}*{nz}*5"
    return q.reshape((nx, ny, nz, 5), order="F")


def lin_weights(xo, xn):
    """indices i and weights w with f(xn) = (1-w) f[i] + w f[i+1]; clamped outside xo."""
    xc = np.clip(xn, xo[0], xo[-1])
    i = np.clip(np.searchsorted(xo, xc, side="right") - 1, 0, len(xo) - 2)
    w = (xc - xo[i]) / (xo[i + 1] - xo[i])
    return i, w


def interp_axis(a, axis, i, w):
    sh = [1] * a.ndim
    sh[axis] = len(w)
    w = w.reshape(sh)
    return np.take(a, i, axis=axis) * (1 - w) + np.take(a, i + 1, axis=axis) * w


def regrid_xy(Q, xo, yo, xn_src, yn):
    """Q(nxo, nyo, nz, 5) at (xo, yo) -> at (xn_src, yn) (old-grid coordinates)."""
    i, w = lin_weights(xo, xn_src)
    Qx = interp_axis(Q, 0, i, w)
    j, v = lin_weights(yo, yn)
    return interp_axis(Qx, 1, j, v)


def prim(Q):
    rho = Q[..., 0]
    u, v, w = Q[..., 1] / rho, Q[..., 2] / rho, Q[..., 3] / rho
    p = (GAMMA - 1) * (Q[..., 4] - 0.5 * rho * (u * u + v * v + w * w))
    return rho, u, v, w, p


def cons(rho, u, v, w, p):
    return np.stack([rho, rho * u, rho * v, rho * w, p / (GAMMA - 1) + 0.5 * rho * (u * u + v * v + w * w)], -1)


def wall_rows(Q, y):
    """span-mean Cf and p_w along x (instantaneous)."""
    rho, u, v, w, p = prim(Q[:, :, 3:-3, :])
    um = u.mean(2).T                     # (ny, nx)
    Tw = (p[:, 0] / (R_GAS * rho[:, 0])).mean(1)
    return wall_tau(y, um, Tw), p[:, 0].mean(1)


def main():
    old_rd, new_mg, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    Po = case_params(old_rd)
    tmp = os.path.dirname(os.path.abspath(new_mg))
    Pn = case_params(tmp) if os.path.basename(new_mg) == "mod_globals.f90" else None
    assert Pn is not None, "NEW_MOD_GLOBALS must be a file named mod_globals.f90"
    assert Po["nz1"] == Pn["nz1"] and abs(Po["Lz"] - Pn["Lz"]) < 1e-12, "span must not change"
    blt = Po["blt"]
    nz = Po["nz1"]
    src = os.path.join(old_rd, "recal")
    os.makedirs(dst, exist_ok=True)

    # ---------------- TBL generator ----------------
    xo1, xn1 = make_x(Po, 1), make_x(Pn, 1)
    yo1, yn1 = make_y(Po, Po["ny1"]), make_y(Pn, Pn["ny1"])
    Q1 = read_q(os.path.join(src, "Q00001.dat"), Po["nx1"], Po["ny1"], nz)
    N1 = regrid_xy(Q1, xo1, yo1, xn1, yn1)
    N1[:, 0, :, 1:4] = 0.0
    N1.reshape(-1, order="F").astype("<f8").tofile(os.path.join(dst, "Q00001.dat"))
    print(f"TBL : {Q1.shape[:2]} -> {N1.shape[:2]}")
    del Q1

    # Qm: five primitive profiles over ny1
    qm = np.fromfile(os.path.join(src, "Qm.dat"), dtype="<f8").reshape(5, Po["ny1"])
    j, v = lin_weights(yo1, yn1)
    qmn = qm[:, j] * (1 - v) + qm[:, j + 1] * v
    qmn.astype("<f8").tofile(os.path.join(dst, "Qm.dat"))

    # ---------------- SBLI domain ----------------
    xo, xn = make_x(Po, 3), make_x(Pn, 3)
    yo, yn = make_y(Po, Po["ny2"]), make_y(Pn, Pn["ny2"])
    ximp_o = xo[0] + Po["ximp_in"] * blt
    ximp_n = xn[0] + Pn["ximp_in"] * blt
    xi_n = (xn - ximp_n) / blt
    xi_in_o, xi_in_n = (xo[0] - ximp_o) / blt, (xn[0] - ximp_n) / blt
    xi_src = np.where(xi_n >= XI_KEEP, xi_n,
                      XI_KEEP + (xi_n - XI_KEEP) * (XI_KEEP - xi_in_o) / (XI_KEEP - xi_in_n))
    x_src = ximp_o + xi_src * blt
    Q2 = read_q(os.path.join(src, "Q00002.dat"), Po["nx2"], Po["ny2"], nz)
    N2 = regrid_xy(Q2, xo, yo, x_src, yn)

    # outer three-state field above the old top
    th = math.radians(Po["theta_deg"])
    b1 = _beta_weak(Po["M0"], th)
    rho2, p2, M2 = _oblique(Po["M0"], b1, th, Po["rho0"], Po["p0"])
    br = _beta_weak(M2, th) - th                       # reflected shock angle to the wall
    ytop = yo[-1]
    rho, u, v, w, p = prim(Q2[:, -2, 3:-3, :])         # row ny-1 (interior of the old top)
    rm, um, vm, pm = (a.mean(1) for a in (rho, u, v, p))
    xi_o = (xo - ximp_o) / blt
    xi_inc_top = -ytop / blt / math.tan(b1)
    # reflected shock at the old top: first x behind the incident shock where the row
    # pressure passes half-way from state 2 to the inviscid state 3 level
    _, p3, _ = _oblique(M2, br + th, th, rho2, p2)
    beh = xi_o > xi_inc_top + 3
    p2_meas = pm[beh & (xi_o < xi_inc_top + 6)].mean()
    cross = np.where(beh & (pm > 0.5 * (p2_meas + p2_meas * p3 / p2)))[0]
    xi_ref_top = xi_o[cross[0]] if cross.size else xi_o[-1]
    regions = {1: xi_o < xi_inc_top - 1, 2: (xi_o > xi_inc_top + 1) & (xi_o < xi_ref_top - 1),
               3: (xi_o > xi_ref_top + 2) & (xi_o < Po["xdown"])}
    S = {}
    for k, m in regions.items():
        if not m.any():
            m = np.zeros_like(m); m[-1] = True
        S[k] = np.array([rm[m].mean(), um[m].mean(), vm[m].mean(), pm[m].mean()])
    print(f"outer states over the old top (y = {ytop / blt:.2f} delta): incident crosses at xi = {xi_inc_top:.2f}, "
          f"reflected at xi = {xi_ref_top:.2f}")
    for k in (1, 2, 3):
        print(f"  state {k}: rho/rho0 = {S[k][0] / Po['rho0']:.4f}, u/u0 = {S[k][1] / Po['u0']:.4f}, "
              f"v/u0 = {S[k][2] / Po['u0']:+.4f}, p/p0 = {S[k][3] / Po['p0']:.4f}")
    jj = np.where(yn > ytop)[0]
    ref_row = N2[:, jj[0] - 1:jj[0], :, :]              # last interpolated row, (nx, 1, nz, 5)
    for j in jj:
        dyt = (yn[j] - ytop) / blt
        inc = xi_n > xi_inc_top - dyt / math.tan(b1)
        ref = xi_n > xi_ref_top + dyt / math.tan(br)
        reg = np.where(ref, 3, np.where(inc, 2, 1))
        st = np.stack([S[1], S[2], S[3]])[reg - 1]       # (nx, 4)
        Qs = cons(st[:, 0], st[:, 1], st[:, 2], np.zeros(len(xn)), st[:, 3])[:, None, :]
        wgt = min(dyt, 1.0)
        N2[:, j, :, :] = (1 - wgt) * ref_row[:, 0, :, :] + wgt * Qs
    N2[:, 0, :, 1:4] = 0.0
    N2.reshape(-1, order="F").astype("<f8").tofile(os.path.join(dst, "Q00002.dat"))
    print(f"SBLI: {Q2.shape[:2]} -> {N2.shape[:2]}, inflow stretch {(XI_KEEP - xi_in_n) / (XI_KEEP - xi_in_o):.3f}")

    # ---------------- checks ----------------
    r, uu, vv, ww, pp = prim(N2)
    print(f"SBLI new field: min rho/rho0 = {r.min() / Po['rho0']:.3f}, min p/p0 = {pp.min() / Po['p0']:.3f}, "
          f"max p/p0 = {pp.max() / Po['p0']:.3f}, finite = {np.isfinite(N2).all()}")
    tau_o, pw_o = wall_rows(Q2, yo)
    tau_n, pw_n = wall_rows(N2, yn)
    q0 = 0.5 * Po["rho0"] * Po["u0"] ** 2
    print("   xi    Cf_old*1e3  Cf_new*1e3   pw_old/p0  pw_new/p0   (span-mean, instantaneous)")
    for xi in (-30, -26, -20, -12, -8, -6, -4, -2, 0, 2, 5, 10, 18, 22):
        io = int(np.argmin(abs(xi_o - xi_src[np.argmin(abs(xi_n - xi))])))
        i_n = int(np.argmin(abs(xi_n - xi)))
        print(f"  {xi:5.0f}  {tau_o[io] / q0 * 1e3:10.3f}  {tau_n[i_n] / q0 * 1e3:10.3f}  "
              f"{pw_o[io] / Po['p0']:10.4f}  {pw_n[i_n] / Po['p0']:10.4f}")
    shutil.copy(os.path.join(old_rd, "mod_globals.f90"), os.path.join(dst, "regrid_source_mod_globals.f90"))
    print("written to", dst)


if __name__ == "__main__":
    main()
