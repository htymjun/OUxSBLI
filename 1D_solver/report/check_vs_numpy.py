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
                                  [--press fp32]
        VISC_ORDER defaults to ORDER.

--press fp32 mirrors KERNEL_MODE='warp_fused' with PRESS_PREC='fp32': the
explicit KEEP pressure terms (KEEPP in src/calc_warp_kernel.f90.fypp) are
evaluated in float32. The reference subtracts the same terms re-evaluated in
float64 so everything else stays the plain fp64 model; the fp64 twin must stay
in sync with keep()'s pressure association or the cancellation degrades (it
would still sit orders of magnitude under the 1e-7 gate).
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


def step(Q1, order, visc_order, press32=False):
    """One TVD-RK3 timestep, matching calc_time_dev.f90.fypp."""
    ng = max(order, visc_order) // 2
    j = np.arange(ng, NX - ng)               # physical cells (0-based)

    def dE(E):
        return (E[:, j] - E[:, j - 1]).T     # (ncell, 3)

    Q2 = Q1.copy()
    Q2[j, :] = Q1[j, :] - (DT / DX) * dE(flux(Q1, order, visc_order, ng, press32))
    bc(Q2, ng)

    Q2[j, :] = 0.75 * Q1[j, :] + 0.25 * Q2[j, :] \
        - 0.25 * (DT / DX) * dE(flux(Q2, order, visc_order, ng, press32))
    bc(Q2, ng)

    Q1[j, :] = (2.0 * Q2[j, :] + Q1[j, :]
                - 2.0 * (DT / DX) * dE(flux(Q2, order, visc_order, ng, press32))) / 3.0
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
    qdat, order = argv[0], int(argv[1])
    visc_order = int(argv[2]) if len(argv) > 2 else order
    Q = step(initial(), order, visc_order, press32)
    rho, u, p, _, _ = primitives(Q)
    ref = np.c_[rho / RHO0, u / A, p / P0]

    got = np.loadtxt(qdat)[:, 1:4]
    print(f"ORDER = {order},  VISC_ORDER = {visc_order},  1 timestep,  {qdat}"
          + ("  [PRESS_PREC=fp32 reference]" if press32 else ""))
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
