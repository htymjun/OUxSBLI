#!/usr/bin/env python3
"""Order-of-accuracy check for the 1D KEEP2/KEEP4/KEEP6 stencils.

Re-implements the flux formulas from src/calc_keep_1d.f90.fypp in numpy and
measures the observed convergence rate of the flux divergence against the
analytic derivative on a smooth periodic field. This validates the hand-typed
coefficients independently of the GPU code -- a shock tube cannot, because the
discontinuities are first-order for every scheme.

Expected rates: 2, 4, 6.
"""
import numpy as np

R_OVER_GAMMA_1 = 287.03 / 0.4


def _s(a, k):
    """a shifted by k with periodic wrap: _s(a,k)[i] = a[i+k]."""
    return np.roll(a, -k)


def keep_flux(rho, u, p, T, order):
    """Flux (mass, momentum, energy) at faces i+1/2, periodic. Mirrors the
    Fortran functions term for term; index j there = shift (j-1-order//2+1) here."""
    if order == 2:
        r1, r2 = rho, _s(rho, 1)
        u1, u2 = u, _s(u, 1)
        p1, p2 = p, _s(p, 1)
        T1, T2 = T, _s(T, 1)
        f1 = 0.25 * (r1 + r2) * (u1 + u2)
        f2 = 0.5 * (f1 * (u1 + u2) + (p1 + p2))
        f3 = f1 * 0.5 * (T1 + T2) * R_OVER_GAMMA_1
        f3 = f3 + 0.5 * (u1 * p2 + u2 * p1)
        f3 = f3 + 0.5 * f1 * (u1 * u2)
        return f1, f2, f3

    if order == 4:
        # local indices 1..4 -> shifts -1,0,1,2 ; face between 2 and 3
        rr = [_s(rho, k) for k in (-1, 0, 1, 2)]
        uu = [_s(u, k) for k in (-1, 0, 1, 2)]
        pp = [_s(p, k) for k in (-1, 0, 1, 2)]
        TT = [_s(T, k) for k in (-1, 0, 1, 2)]
        RV1 = (rr[1] + rr[2]) * (uu[1] + uu[2])
        RV2 = (rr[1] + rr[3]) * (uu[1] + uu[3])
        RV3 = (rr[0] + rr[2]) * (uu[0] + uu[2])
        f1 = RV1 / 3.0 - (RV2 + RV3) / 24.0
        RV1, RV2, RV3 = RV1 / 6.0, RV2 / 48.0, RV3 / 48.0
        a, b = RV1 - RV2, RV1 - RV3
        pres = -(pp[0] + pp[3]) / 12.0 + 7.0 / 12.0 * (pp[1] + pp[2])
        f2 = -RV3 * uu[0] + a * uu[1] + b * uu[2] - RV2 * uu[3] + pres
        ene = (-RV3 * TT[0] + a * TT[1] + b * TT[2] - RV2 * TT[3]) * R_OVER_GAMMA_1
        ene = ene + RV1 * uu[1] * uu[2] - RV2 * uu[1] * uu[3] - RV3 * uu[0] * uu[2]
        f3 = ene + 2.0 / 3.0 * (uu[1] * pp[2] + uu[2] * pp[1]) \
            - (uu[1] * pp[3] + uu[3] * pp[1] + uu[0] * pp[2] + uu[2] * pp[0]) / 12.0
        return f1, f2, f3

    if order == 6:
        # local indices 1..6 -> shifts -2..3 ; face between 3 and 4
        sh = (-2, -1, 0, 1, 2, 3)
        rr = [_s(rho, k) for k in sh]
        uu = [_s(u, k) for k in sh]
        pp = [_s(p, k) for k in sh]
        TT = [_s(T, k) for k in sh]
        RV1 = (rr[2] + rr[3]) * (uu[2] + uu[3])
        RV2 = (rr[2] + rr[4]) * (uu[2] + uu[4])
        RV3 = (rr[1] + rr[3]) * (uu[1] + uu[3])
        RV4 = (rr[2] + rr[5]) * (uu[2] + uu[5])
        RV5 = (rr[1] + rr[4]) * (uu[1] + uu[4])
        RV6 = (rr[0] + rr[3]) * (uu[0] + uu[3])
        f1 = 0.375 * RV1 - 0.075 * (RV2 + RV3) + (RV4 + RV5 + RV6) / 120.0
        RV1, RV2, RV3 = 0.1875 * RV1, 0.0375 * RV2, 0.0375 * RV3
        RV4, RV5, RV6 = RV4 / 240.0, RV5 / 240.0, RV6 / 240.0
        c35, c25 = -RV3 + RV5, -RV2 + RV5
        c124, c136 = RV1 - RV2 + RV4, RV1 - RV3 + RV6
        pres = (pp[0] + pp[5] - 8.0 * (pp[1] + pp[4]) + 37.0 * (pp[2] + pp[3])) / 60.0
        f2 = (RV6 * uu[0] + c35 * uu[1] + c124 * uu[2]
              + c136 * uu[3] + c25 * uu[4] + RV4 * uu[5] + pres)
        ene = (RV6 * TT[0] + c35 * TT[1] + c124 * TT[2]
               + c136 * TT[3] + c25 * TT[4] + RV4 * TT[5]) * R_OVER_GAMMA_1
        ene = (ene + RV6 * uu[0] * uu[3] + RV5 * uu[1] * uu[4] + RV4 * uu[2] * uu[5]
               - RV3 * uu[1] * uu[3] - RV2 * uu[2] * uu[4] + RV1 * uu[2] * uu[3])
        f3 = (ene + 0.75 * (uu[2] * pp[3] + uu[3] * pp[2])
              - 0.15 * (uu[2] * pp[4] + uu[4] * pp[2] + uu[1] * pp[3] + uu[3] * pp[1])
              + (uu[2] * pp[5] + uu[5] * pp[2] + uu[1] * pp[4] + uu[4] * pp[1]
                 + uu[0] * pp[3] + uu[3] * pp[0]) / 60.0)
        return f1, f2, f3

    raise ValueError(order)


def field(n):
    """Smooth periodic state on [0,1)."""
    x = np.arange(n) / n
    rho = 1.0 + 0.2 * np.sin(2 * np.pi * x)
    u = 1.0 + 0.1 * np.cos(2 * np.pi * x)
    p = 1.0 + 0.15 * np.sin(4 * np.pi * x)
    T = p / (287.03 * rho)
    dmom_dx = 2 * np.pi * (0.2 * np.cos(2 * np.pi * x) * u
                           - 0.1 * np.sin(2 * np.pi * x) * rho)
    return x, rho, u, p, T, dmom_dx, 1.0 / n


def main():
    print("observed order of accuracy of d(rho u)/dx from the KEEP mass flux\n")
    for order in (2, 4, 6):
        prev_err = prev_n = None
        print(f"  KEEP{order}")
        for n in (64, 128, 256, 512, 1024):
            _, rho, u, p, T, exact, dx = field(n)
            f1, _, _ = keep_flux(rho, u, p, T, order)
            div = (f1 - _s(f1, -1)) / dx          # F_{i+1/2} - F_{i-1/2}
            err = np.abs(div - exact).max()
            rate = ("" if prev_err is None
                    else f"   rate = {np.log2(prev_err / err):.2f}")
            print(f"    n={n:5d}  Linf = {err:.4e}{rate}")
            prev_err, prev_n = err, n
        print()


if __name__ == "__main__":
    main()
