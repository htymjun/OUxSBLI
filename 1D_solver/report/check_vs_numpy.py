#!/usr/bin/env python3
"""End-to-end check of the 1D solver against an independent numpy reference.

Reproduces one full TVD-RK3 timestep (3 stages: primitives -> KEEP convective
flux -> viscous flux -> RK update -> BC) exactly as the CUDA Fortran does, and
compares against the solver's Q.dat.

This is the gate that the hand-typed KEEP4/KEEP6 and VISC4/VISC6 coefficients in
src/calc_keep_1d.f90.fypp / src/calc_visc_1d.f90.fypp are what actually runs on
the GPU: check_keep_order.py and check_visc_order.py prove the coefficients are
Nth-order, this proves the Fortran uses them.

There is no boundary order degradation: set_bc fills ng = max(ORDER,VISC_ORDER)/2
ghost cells, so one uniform stencil covers every face the RK update consumes
(Fortran faces ng..nx-ng, cells ng+1..nx-ng).

Usage:  python3 check_vs_numpy.py <Q.dat from an nt=1 run> <ORDER> [VISC_ORDER]
                                  [--press fp32] [--scheme SLAU] [--recon WENO]
        VISC_ORDER defaults to ORDER.

--scheme SLAU checks a SCHEME='SLAU' build instead of the default KEEP:
plain SLAU1 (SLAU_VARIANT='SLAU' -- 1D_solver doesn't implement HRSLAU2),
MUSCL-reconstructed at ORDER>=4 via calc_muscl.f90.fypp's TVD='tvd' (Minmod)
dispatch -- the only TVD value 1D_solver's SCHEME='SLAU' supports. Not
compatible with --press (a KEEP-only/PRESS_PREC concept).

--recon WENO (only meaningful with --scheme SLAU and ORDER=6) checks
SLAU_RECON='WENO' instead of the default MUSCL: src/calc_weno.f90's WENO5-Z
reconstruction (calc_slau_kernel.f90.fypp's delta6_weno) in place of
calc_muscl.f90.fypp's delta6. There is no WENO3, so this is ORDER=6-only,
matching the compile-time guard in calc_slau_kernel.f90.fypp.

--press fp32 mirrors PRESS_PREC='fp32' (now available on KERNEL_MODE in
'seq'/'fused'/'warp_fused' -- see src/calc_keep_1d.f90.fypp's KEEPNP/KEEPP):
the explicit KEEP pressure terms are evaluated in float32. The reference
subtracts the same terms re-evaluated in float64 so everything else stays the
plain fp64 model; the fp64 twin must stay in sync with keep()'s pressure
association or the cancellation degrades (it would still sit orders of
magnitude under the 1e-7 gate).

Scope note: this reference only models a plain float64 KEEP baseline (plus the
optional --press fp32 bias on top of it). It has no KEEP_PREC='fp32'/'term'
model at all, so it must only be run against Q.dat built with KEEP_PREC='fp64'
(--press fp32 or not). For KEEP_PREC='term' or 'fp32' builds -- including
term x PRESS_PREC='fp32' -- use check_sod.py instead: compare against the
exact solution and, separately, against a same-KEEP_PREC PRESS_PREC='fp64'
reference Q.dat (see 1D_solver/CLAUDE.md's Verification section).
"""
import sys
import numpy as np

# --- mod_globals.f90 / mod_constant.f90.fypp ---------------------------------
NX = 4096
GAMMA, RGAS, PR = 1.4, 287.03, 0.72
MU0, RE, TLR, RHO0 = 1.716e-5, 25000.0, 300.0, 1.293
P0 = RHO0 * RGAS * TLR
RHO1, P1 = 0.125 * RHO0, 0.1 * P0
A = np.sqrt(RGAS * TLR)
LX = RE * MU0 / (RHO0 * A)
DX = LX / (NX - 1)
DT = 0.1 * DX / A

R_OVER_GAMMA_1 = RGAS / (GAMMA - 1.0)
CP_OVER_PR = (GAMMA * RGAS / (GAMMA - 1.0)) / PR
TWO_THIRD = 2.0 / 3.0
MU_SUTH = 1.716e-5 * 383.6 * 273.2 ** (-1.5)


def primitives(Q):
    rho = Q[:, 0]
    u = Q[:, 1] / rho
    p = (GAMMA - 1.0) * (Q[:, 2] - 0.5 * rho * u * u)
    T = p / (RGAS * rho)
    mu = MU_SUTH / (T + 111.0) * (T * np.sqrt(T))
    return rho, u, p, T, mu


# Face-centred viscous stencils, mirroring src/calc_visc_1d.f90.fypp.
# interp: coefficients on (a[m+1-j] + a[m+j]); diff: on (a[m+j] - a[m+1-j]).
VISC_INTERP = {2: [1.0 / 2.0],
               4: [9.0 / 16.0, -1.0 / 16.0],
               6: [150.0 / 256.0, -25.0 / 256.0, 3.0 / 256.0]}
VISC_DIFF = {2: [1.0],
             4: [27.0 / 24.0, -1.0 / 24.0],
             6: [2250.0 / 1920.0, -125.0 / 1920.0, 9.0 / 1920.0]}


def visc(u, T, mu, lo, order):
    """Fv(1), Fv(2) at faces `lo` (0-based id of the face's left cell)."""
    m = order // 2
    js = range(1, m + 1)
    mu_f = sum(c * (mu[lo + 1 - j] + mu[lo + j]) for c, j in zip(VISC_INTERP[order], js))
    mudx = mu_f / DX
    dT = sum(c * (T[lo + j] - T[lo + 1 - j]) for c, j in zip(VISC_DIFF[order], js))
    du = sum(c * (u[lo + j] - u[lo + 1 - j]) for c, j in zip(VISC_DIFF[order], js))
    u_f = sum(c * (u[lo + 1 - j] + u[lo + j]) for c, j in zip(VISC_INTERP[order], js))
    fv2 = CP_OVER_PR * mudx * dT
    fv1 = TWO_THIRD * 2.0 * (mudx * du)
    fv2 = fv2 + u_f * fv1
    return fv1, fv2


def keep(rho, u, p, T, lo, order):
    """Flux at faces lo (0-based cell index of the left neighbour of the face),
    using an `order`-point stencil centred on that face."""
    half = order // 2
    idx = [lo - half + 1 + k for k in range(order)]          # 0-based cell ids
    r = [rho[i] for i in idx]
    v = [u[i] for i in idx]
    q = [p[i] for i in idx]
    t = [T[i] for i in idx]

    if order == 2:
        f1 = 0.25 * (r[0] + r[1]) * (v[0] + v[1])
        f2 = 0.5 * (f1 * (v[0] + v[1]) + (q[0] + q[1]))
        f3 = f1 * 0.5 * (t[0] + t[1]) * R_OVER_GAMMA_1
        f3 = f3 + 0.5 * (v[0] * q[1] + v[1] * q[0])
        f3 = f3 + 0.5 * f1 * (v[0] * v[1])
    elif order == 4:
        RV1 = (r[1] + r[2]) * (v[1] + v[2])
        RV2 = (r[1] + r[3]) * (v[1] + v[3])
        RV3 = (r[0] + r[2]) * (v[0] + v[2])
        f1 = RV1 / 3.0 - (RV2 + RV3) / 24.0
        RV1, RV2, RV3 = RV1 / 6.0, RV2 / 48.0, RV3 / 48.0
        a, b = RV1 - RV2, RV1 - RV3
        pres = -(q[0] + q[3]) / 12.0 + 7.0 / 12.0 * (q[1] + q[2])
        f2 = -RV3 * v[0] + a * v[1] + b * v[2] - RV2 * v[3] + pres
        ene = (-RV3 * t[0] + a * t[1] + b * t[2] - RV2 * t[3]) * R_OVER_GAMMA_1
        ene = ene + RV1 * v[1] * v[2] - RV2 * v[1] * v[3] - RV3 * v[0] * v[2]
        f3 = ene + TWO_THIRD * (v[1] * q[2] + v[2] * q[1]) \
            - (v[1] * q[3] + v[3] * q[1] + v[0] * q[2] + v[2] * q[0]) / 12.0
    elif order == 6:
        RV1 = (r[2] + r[3]) * (v[2] + v[3])
        RV2 = (r[2] + r[4]) * (v[2] + v[4])
        RV3 = (r[1] + r[3]) * (v[1] + v[3])
        RV4 = (r[2] + r[5]) * (v[2] + v[5])
        RV5 = (r[1] + r[4]) * (v[1] + v[4])
        RV6 = (r[0] + r[3]) * (v[0] + v[3])
        f1 = 0.375 * RV1 - 0.075 * (RV2 + RV3) + (RV4 + RV5 + RV6) / 120.0
        RV1, RV2, RV3 = 0.1875 * RV1, 0.0375 * RV2, 0.0375 * RV3
        RV4, RV5, RV6 = RV4 / 240.0, RV5 / 240.0, RV6 / 240.0
        c35, c25 = -RV3 + RV5, -RV2 + RV5
        c124, c136 = RV1 - RV2 + RV4, RV1 - RV3 + RV6
        pres = (q[0] + q[5] - 8.0 * (q[1] + q[4]) + 37.0 * (q[2] + q[3])) / 60.0
        f2 = (RV6 * v[0] + c35 * v[1] + c124 * v[2]
              + c136 * v[3] + c25 * v[4] + RV4 * v[5] + pres)
        ene = (RV6 * t[0] + c35 * t[1] + c124 * t[2]
               + c136 * t[3] + c25 * t[4] + RV4 * t[5]) * R_OVER_GAMMA_1
        ene = (ene + RV6 * v[0] * v[3] + RV5 * v[1] * v[4] + RV4 * v[2] * v[5]
               - RV3 * v[1] * v[3] - RV2 * v[2] * v[4] + RV1 * v[2] * v[3])
        f3 = (ene + 0.75 * (v[2] * q[3] + v[3] * q[2])
              - 0.15 * (v[2] * q[4] + v[4] * q[2] + v[1] * q[3] + v[3] * q[1])
              + (v[2] * q[5] + v[5] * q[2] + v[1] * q[4] + v[4] * q[1]
                 + v[0] * q[3] + v[3] * q[0]) / 60.0)
    else:
        raise ValueError(order)
    return np.stack([f1, f2, f3])


def keep_pressure(u, p, lo, order, ft):
    """Pressure part (f2, f3) of the KEEP face flux at dtype ft, transcribing
    KEEPP{2,4,6} of src/calc_warp_kernel.f90.fypp term for term (0-based:
    q[k] is Fortran p(k+1))."""
    half = order // 2
    idx = [lo - half + 1 + k for k in range(order)]
    v = [u[i].astype(ft) for i in idx]
    q = [p[i].astype(ft) for i in idx]
    c = ft
    if order == 2:
        f2 = c(0.5) * (q[0] + q[1])
        f3 = c(0.5) * (v[0] * q[1] + v[1] * q[0])
    elif order == 4:
        f2 = c(7.0 / 12.0) * (q[1] + q[2]) - c(1.0 / 12.0) * (q[0] + q[3])
        f3 = c(2.0 / 3.0) * (v[1] * q[2] + v[2] * q[1]) \
            - c(1.0 / 12.0) * (v[1] * q[3] + v[3] * q[1] + v[0] * q[2] + v[2] * q[0])
    elif order == 6:
        f2 = c(37.0 / 60.0) * (q[2] + q[3]) \
            + (q[0] + q[5] - c(8.0) * (q[1] + q[4])) * c(1.0 / 60.0)
        f3 = c(0.75) * (v[2] * q[3] + v[3] * q[2]) \
            + c(-3.0 / 20.0) * (v[2] * q[4] + v[4] * q[2] + v[1] * q[3] + v[3] * q[1]) \
            + (v[2] * q[5] + v[5] * q[2] + v[1] * q[4] + v[4] * q[1]
               + v[0] * q[3] + v[3] * q[0]) * c(1.0 / 60.0)
    else:
        raise ValueError(order)
    return f2, f3


# --- SLAU / MUSCL, mirroring src/calc_slau_1d.f90.fypp and the repo-root
# src/calc_muscl.f90.fypp's TVD='tvd' (Minmod) dispatch -- the only TVD value
# 1D_solver's SCHEME='SLAU' supports (see calc_flux_base.f90.fypp's guard).
_ONE_THIRD, _ONE_SIXTH = 1.0 / 3.0, 1.0 / 6.0


def minmod2(x, y):
    sgn = np.copysign(1.0, x)
    return sgn * np.maximum(np.minimum(np.abs(x), sgn * y), 0.0)


def minmod3(x, y, z):
    sgn = np.copysign(1.0, x)
    return sgn * np.maximum(np.minimum(np.minimum(np.abs(x), sgn * y), sgn * z), 0.0)


def d33(d1, d2, d3):
    return (minmod3(d1, 2.0 * d2, 2.0 * d3)
            - 2.0 * minmod3(d2, 2.0 * d1, 2.0 * d3)
            + minmod3(d3, 2.0 * d1, 2.0 * d2))


def muscl3rd_minmod(a2, a3, d1, d2, d3):
    """TVD='tvd' dispatch of MUSCL3rd (calc_muscl.f90.fypp's MUSCL3rdMinmod)."""
    b = (3.0 - _ONE_THIRD) / (1.0 - _ONE_THIRD)
    dt1, dt2 = minmod2(d1, b * d2), minmod2(d2, b * d1)
    al = a2 + (2.0 * dt2 + dt1) * _ONE_SIXTH
    dt3, dt4 = minmod2(d3, b * d2), minmod2(d2, b * d3)
    ar = a3 - (2.0 * dt4 + dt3) * _ONE_SIXTH
    return al, ar


def muscl4th_tvd(a2, a3, d1, d2, d3, d4, d5):
    """TVD='tvd' dispatch of MUSCL4th (calc_muscl.f90.fypp's MUSCL4thTVD)."""
    delta2 = d3 - d33(d2, d3, d4) * _ONE_SIXTH
    delta1 = d2 - d33(d1, d2, d3) * _ONE_SIXTH
    dl, dr = minmod2(delta1, 4.0 * delta2), minmod2(delta2, 4.0 * delta1)
    al = a2 + (2.0 * dr + dl) * _ONE_SIXTH
    delta3 = d4 - d33(d3, d4, d5) * _ONE_SIXTH
    dl2, dr2 = minmod2(delta2, 4.0 * delta3), minmod2(delta3, 4.0 * delta2)
    ar = a3 - (2.0 * dl2 + dr2) * _ONE_SIXTH
    return al, ar


def muscl_delta4(a0, a1, a2, a3):
    """calc_muscl.f90.fypp's delta4: 3rd-order MUSCL reconstruction at the
    face between a1 and a2."""
    d1, d2, d3 = -a0 + a1, -a1 + a2, -a2 + a3
    return muscl3rd_minmod(a1, a2, d1, d2, d3)


def muscl_delta6(a0, a1, a2, a3, a4, a5):
    """calc_muscl.f90.fypp's delta6: 4th-order MUSCL reconstruction at the
    face between a2 and a3."""
    d1, d2, d3, d4, d5 = -a0 + a1, -a1 + a2, -a2 + a3, -a3 + a4, -a4 + a5
    return muscl4th_tvd(a2, a3, d1, d2, d3, d4, d5)


# --- WENO5-Z (Borges et al.), mirroring src/calc_weno.f90's
# weno5z_left/weno5z_right -- 1D_solver's alternative to MUSCL at ORDER=6
# (SLAU_RECON='WENO', see calc_slau_kernel.f90.fypp's delta6_weno).
def weno5z_left(v1, v2, v3, v4, v5):
    eps, d0, d1, d2 = 1.0e-20, 1.0 / 10.0, 6.0 / 10.0, 3.0 / 10.0
    p0 = (2.0 * v1 - 7.0 * v2 + 11.0 * v3) / 6.0
    p1 = (-1.0 * v2 + 5.0 * v3 + 2.0 * v4) / 6.0
    p2 = (2.0 * v3 + 5.0 * v4 - 1.0 * v5) / 6.0
    b0 = (13.0 / 12.0) * (v1 - 2.0 * v2 + v3) ** 2 + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) ** 2
    b1 = (13.0 / 12.0) * (v2 - 2.0 * v3 + v4) ** 2 + 0.25 * (v2 - v4) ** 2
    b2 = (13.0 / 12.0) * (v3 - 2.0 * v4 + v5) ** 2 + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) ** 2
    tau5 = np.abs(b0 - b2)
    a0 = d0 * (1.0 + (tau5 / (b0 + eps)) ** 2)
    a1 = d1 * (1.0 + (tau5 / (b1 + eps)) ** 2)
    a2 = d2 * (1.0 + (tau5 / (b2 + eps)) ** 2)
    s = a0 + a1 + a2
    return (a0 * p0 + a1 * p1 + a2 * p2) / s


def weno5z_right(v1, v2, v3, v4, v5):
    eps, d0, d1, d2 = 1.0e-20, 1.0 / 10.0, 6.0 / 10.0, 3.0 / 10.0
    p0 = (-1.0 * v1 + 5.0 * v2 + 2.0 * v3) / 6.0
    p1 = (2.0 * v2 + 5.0 * v3 - 1.0 * v4) / 6.0
    p2 = (11.0 * v3 - 7.0 * v4 + 2.0 * v5) / 6.0
    b0 = (13.0 / 12.0) * (v1 - 2.0 * v2 + v3) ** 2 + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3) ** 2
    b1 = (13.0 / 12.0) * (v2 - 2.0 * v3 + v4) ** 2 + 0.25 * (v2 - v4) ** 2
    b2 = (13.0 / 12.0) * (v3 - 2.0 * v4 + v5) ** 2 + 0.25 * (3.0 * v3 - 4.0 * v4 + v5) ** 2
    tau5 = np.abs(b0 - b2)
    a0 = d0 * (1.0 + (tau5 / (b0 + eps)) ** 2)
    a1 = d1 * (1.0 + (tau5 / (b1 + eps)) ** 2)
    a2 = d2 * (1.0 + (tau5 / (b2 + eps)) ** 2)
    s = a0 + a1 + a2
    return (a0 * p0 + a1 * p1 + a2 * p2) / s


def weno_delta6(a0, a1, a2, a3, a4, a5):
    """calc_weno.f90's delta6_weno: al's window centred on a2, ar's on a3 --
    the face between a2 and a3, matching muscl_delta6's convention."""
    al = weno5z_left(a0, a1, a2, a3, a4)
    ar = weno5z_right(a1, a2, a3, a4, a5)
    return al, ar


def slau_common(rho1, rho2, u1, u2, p1, p2):
    """Mirrors src/calc_slau_1d.f90.fypp's SLAU_common (1D collapse of
    src/calc_scheme_math.f90.fypp: un1==u1, un2==u2, Normal==1)."""
    over_rho1, over_rho2 = 1.0 / rho1, 1.0 / rho2
    c = 0.5 * (np.sqrt(GAMMA * p1 * over_rho1) + np.sqrt(GAMMA * p2 * over_rho2))
    over_c = 1.0 / c
    Mp, Mm = u1 * over_c, u2 * over_c
    g = -np.maximum(np.minimum(Mp, 0.0), -1.0) * np.minimum(np.maximum(Mm, 0.0), 1.0)
    one_g_Vt = (1.0 - g) * (rho1 * np.abs(u1) + rho2 * np.abs(u2)) / (rho1 + rho2)
    Vtp, Vtm = one_g_Vt + g * np.abs(u1), one_g_Vt + g * np.abs(u2)
    bp = np.where(np.abs(Mp) < 1.0, 0.25 * (2.0 - Mp) * (Mp + 1.0) ** 2,
                  0.5 * (1.0 + np.copysign(1.0, Mp)))
    bm = np.where(np.abs(Mm) < 1.0, 0.25 * (2.0 + Mm) * (Mm - 1.0) ** 2,
                  0.5 * (1.0 + np.copysign(1.0, -Mm)))
    dp = -p1 + p2
    return over_rho1, over_rho2, over_c, bp, bm, dp, Vtp, Vtm


def slau_phi(rho, k, p, over_rho):
    return (p * GAMMA / (GAMMA - 1.0) + rho * k) * over_rho


def slau1(rho1, rho2, u1, u2, p1, p2):
    """Mirrors src/calc_slau_1d.f90.fypp's SLAU1 (SLAU_VARIANT='SLAU' --
    1D_solver's only supported SLAU variant, see calc_flux_base.f90.fypp)."""
    over_rho1, over_rho2, over_c, bp, bm, dp, Vtp, Vtm = slau_common(rho1, rho2, u1, u2, p1, p2)
    k1, k2 = 0.5 * u1 * u1, 0.5 * u2 * u2
    M = np.minimum(1.0, np.sqrt(k1 + k2) * over_c)
    x = (1.0 - M) ** 2
    mass = 0.25 * (rho1 * (u1 + Vtp) + rho2 * (u2 - Vtm) - x * dp * over_c)
    mass1, mass2 = mass + np.abs(mass), mass - np.abs(mass)
    pres = 0.5 * ((p1 + p2) * ((1.0 - x) * (bp + bm - 1.0) + 1.0) + (bm - bp) * dp)
    f1 = mass1 + mass2
    f2 = mass1 * u1 + mass2 * u2 + pres
    f3 = mass1 * slau_phi(rho1, k1, p1, over_rho1) + mass2 * slau_phi(rho2, k2, p2, over_rho2)
    return np.stack([f1, f2, f3])


def flux_slau(Q, order, visc_order, ng, recon="MUSCL"):
    """SLAU analogue of flux() below. ORDER==2 is direct upwind (no
    reconstruction), matching calc_slau_kernel.f90.fypp; ORDER>=4
    reconstructs left/right states first, with MUSCL (muscl_delta4/6) or,
    at ORDER==6 only, WENO5-Z (weno_delta6, SLAU_RECON='WENO')."""
    assert recon == "MUSCL" or order == 6, "SLAU_RECON='WENO' needs ORDER=6 (no WENO3)"
    rho, u, p, T, mu = primitives(Q)
    f = np.arange(ng - 1, NX - ng)
    E = np.full((3, NX - 1), np.nan)

    if order == 2:
        rhol, rhor, ul, ur, pl, pr = rho[f], rho[f + 1], u[f], u[f + 1], p[f], p[f + 1]
    else:
        half = order // 2
        idx = [f - half + 1 + k for k in range(order)]
        if order == 4:
            recon_fn = muscl_delta4
        else:
            recon_fn = weno_delta6 if recon == "WENO" else muscl_delta6
        rhol, rhor = recon_fn(*[rho[i] for i in idx])
        ul, ur = recon_fn(*[u[i] for i in idx])
        pl, pr = recon_fn(*[p[i] for i in idx])

    E[:, f] = slau1(rhol, rhor, ul, ur, pl, pr)
    fv1, fv2 = visc(u, T, mu, f, visc_order)
    E[1, f] -= fv1
    E[2, f] -= fv2
    return E


def flux(Q, order, visc_order, ng, press32=False):
    """E(3, nx-1); face f (0-based) sits between cells f and f+1. Only the faces
    the RK update consumes (Fortran ng..nx-ng, i.e. 0-based ng-1..nx-ng-1) are
    written -- the rest are never read, exactly as on the GPU."""
    rho, u, p, T, mu = primitives(Q)
    f = np.arange(ng - 1, NX - ng)           # 0-based face ids actually computed
    E = np.full((3, NX - 1), np.nan)         # NaN so an accidental read is loud

    E[:, f] = keep(rho, u, p, T, f, order)
    if press32:
        f2_32, f3_32 = keep_pressure(u, p, f, order, np.float32)
        f2_64, f3_64 = keep_pressure(u, p, f, order, np.float64)
        E[1, f] += f2_32.astype(np.float64) - f2_64
        E[2, f] += f3_32.astype(np.float64) - f3_64
    fv1, fv2 = visc(u, T, mu, f, visc_order)
    E[1, f] -= fv1
    E[2, f] -= fv2
    return E


def bc(Q, ng):
    """Zero-gradient, ng cells deep; mirrors ST/set.f90."""
    Q[:ng, :] = Q[ng, :]
    Q[NX - ng:, :] = Q[NX - ng - 1, :]


def step(Q1, order, visc_order, press32=False, scheme="KEEP", recon="MUSCL"):
    """One TVD-RK3 timestep, matching calc_time_dev.f90.fypp."""
    ng = max(order, visc_order) // 2
    j = np.arange(ng, NX - ng)               # physical cells (0-based)

    if scheme == "SLAU":
        assert not press32, "press32 is a KEEP-only (PRESS_PREC) concept"
        def calc_flux(Q):
            return flux_slau(Q, order, visc_order, ng, recon)
    else:
        assert recon == "MUSCL", "--recon WENO only applies to --scheme SLAU"
        def calc_flux(Q):
            return flux(Q, order, visc_order, ng, press32)

    def dE(E):
        return (E[:, j] - E[:, j - 1]).T     # (ncell, 3)

    Q2 = Q1.copy()
    Q2[j, :] = Q1[j, :] - (DT / DX) * dE(calc_flux(Q1))
    bc(Q2, ng)

    Q2[j, :] = 0.75 * Q1[j, :] + 0.25 * Q2[j, :] \
        - 0.25 * (DT / DX) * dE(calc_flux(Q2))
    bc(Q2, ng)

    Q1[j, :] = (2.0 * Q2[j, :] + Q1[j, :]
                - 2.0 * (DT / DX) * dE(calc_flux(Q2))) / 3.0
    bc(Q1, ng)
    return Q1


def initial():
    Q = np.empty((NX, 3))
    h = NX // 2
    Q[:h] = [RHO0, 0.0, P0 / (GAMMA - 1.0)]
    Q[h:] = [RHO1, 0.0, P1 / (GAMMA - 1.0)]
    return Q


def main():
    argv = sys.argv[1:]
    press32 = False
    if "--press" in argv:
        i = argv.index("--press")
        press32 = argv[i + 1] == "fp32"
        del argv[i:i + 2]
    scheme = "KEEP"
    if "--scheme" in argv:
        i = argv.index("--scheme")
        scheme = argv[i + 1]
        del argv[i:i + 2]
    recon = "MUSCL"
    if "--recon" in argv:
        i = argv.index("--recon")
        recon = argv[i + 1]
        del argv[i:i + 2]
    qdat, order = argv[0], int(argv[1])
    visc_order = int(argv[2]) if len(argv) > 2 else order
    Q = step(initial(), order, visc_order, press32, scheme, recon)
    rho, u, p, _, _ = primitives(Q)
    ref = np.c_[rho / RHO0, u / A, p / P0]

    got = np.loadtxt(qdat)[:, 1:4]
    print(f"SCHEME = {scheme},  SLAU_RECON = {recon},  ORDER = {order},  VISC_ORDER = {visc_order},"
          f"  1 timestep,  {qdat}"
          + ("  [PRESS_PREC=fp32 reference]" if press32 else ""))
    if not np.isfinite(got).all():
        bad = np.argwhere(~np.isfinite(got))
        i, k = bad[0]
        names = ("rho", "u", "p")
        print(f"  non-finite output: first {names[k]} at row {i + 1} is {got[i, k]}")
        print("\n  worst relative difference = nan  ->  FAIL"
              "   (tolerance 1e-7, set by Q.dat's 8-digit output)")
        return 1
    if not np.isfinite(ref).all():
        bad = np.argwhere(~np.isfinite(ref))
        i, k = bad[0]
        names = ("rho", "u", "p")
        print(f"  non-finite reference: first {names[k]} at row {i + 1} is {ref[i, k]}")
        print("\n  worst relative difference = nan  ->  FAIL"
              "   (tolerance 1e-7, set by Q.dat's 8-digit output)")
        return 1
    worst = 0.0
    for k, name in enumerate(("rho", "u", "p")):
        d = np.abs(got[:, k] - ref[:, k])
        scale = max(np.abs(ref[:, k]).max(), 1e-300)
        print(f"  {name:3s}: max|diff| = {d.max():.3e}   relative = {d.max()/scale:.3e}")
        worst = max(worst, d.max() / scale)
    # Q.dat carries 8 significant digits, so agreement is limited by print width.
    ok = worst < 1e-7
    print(f"\n  worst relative difference = {worst:.3e}  ->  {'PASS' if ok else 'FAIL'}"
          f"   (tolerance 1e-7, set by Q.dat's 8-digit output)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
