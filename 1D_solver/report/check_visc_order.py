#!/usr/bin/env python3
"""Order-of-accuracy check for the viscous flux stencils.

Re-implements the face interpolation / face derivative from
src/calc_visc_1d.f90.fypp in numpy and measures convergence on a smooth
periodic field. Pure numpy, no GPU.

This is the check that catches a mistyped coefficient. The Sod shock tube
cannot: every scheme is first-order at a discontinuity.

Two tables are printed, and they answer different questions:

  1. FACE FLUX vs exact (4/3) mu du/dx at the face.
     Expect 2 / 4 / 6. This validates the coefficients.

  2. DIVERGENCE (F_{i+1/2} - F_{i-1/2})/dx vs exact d/dx[(4/3) mu du/dx],
     which is what calc_R actually forms.
     Expect 2 / 2 / 2 -- NOT 2 / 4 / 6.

Table 2 is not a bug. A p-th-order accurate flux differenced across one cell
still carries the (h^2/24) F''' term of the difference operator itself, so the
divergence is 2nd order at every stencil width. 3D_solver's calc_visc_high has
exactly the same property; VISC_ORDER is a stencil WIDTH knob, and its purpose
in this testbed is to grow the FP32 instruction count. See 1D_solver/CLAUDE.md.
"""
import numpy as np

# Face-centred stencils, face between indices p/2 and p/2+1 (1-based), written
# as pair sums/differences exactly as calc_visc_1d.f90.fypp writes them.
# interp: coefficients on (a[m+1-j] + a[m+j]),  j = 1..m
# diff:   coefficients on (a[m+j] - a[m+1-j]),  j = 1..m   (times 1/dx)
INTERP = {2: [1.0 / 2.0],
          4: [9.0 / 16.0, -1.0 / 16.0],
          6: [150.0 / 256.0, -25.0 / 256.0, 3.0 / 256.0]}
DIFF = {2: [1.0],
        4: [27.0 / 24.0, -1.0 / 24.0],
        6: [2250.0 / 1920.0, -125.0 / 1920.0, 9.0 / 1920.0]}

TWO_THIRD = 2.0 / 3.0


def _s(a, k):
    """a[i+k] with periodic wrap."""
    return np.roll(a, -k)


def visc_flux(u, mu, dx, order):
    """tau_xx at face i+1/2, i.e. Fv(1) of VISC{order}. Mirrors the Fortran."""
    m = order // 2
    # a(m+1-j) -> offset (1-j), a(m+j) -> offset j, relative to cell i
    mu_f = sum(c * (_s(mu, 1 - j) + _s(mu, j))
               for c, j in zip(INTERP[order], range(1, m + 1)))
    dudx = sum(c * (_s(u, j) - _s(u, 1 - j))
               for c, j in zip(DIFF[order], range(1, m + 1)))
    return TWO_THIRD * 2.0 * (mu_f / dx * dudx)


def field(n):
    """Smooth periodic state on [0,1), with the exact viscous flux and its
    divergence."""
    x = np.arange(n) / n
    tp, fp = 2 * np.pi, 4 * np.pi
    mu = 1.0 + 0.3 * np.sin(tp * x)
    u = np.sin(tp * x) + 0.4 * np.cos(fp * x)

    dmu = 0.3 * tp * np.cos(tp * x)
    du = tp * np.cos(tp * x) - 0.4 * fp * np.sin(fp * x)
    d2u = -tp ** 2 * np.sin(tp * x) - 0.4 * fp ** 2 * np.cos(fp * x)

    # exact tau_xx at the FACE x_{i+1/2}
    dx = 1.0 / n
    xf = x + 0.5 * dx
    mu_f = 1.0 + 0.3 * np.sin(tp * xf)
    du_f = tp * np.cos(tp * xf) - 0.4 * fp * np.sin(fp * xf)
    flux_exact = TWO_THIRD * 2.0 * mu_f * du_f

    # exact d/dx of it at the CELL centre
    div_exact = TWO_THIRD * 2.0 * (dmu * du + mu * d2u)
    return u, mu, dx, flux_exact, div_exact


def table(title, note, err_of):
    print(f"{title}\n  ({note})\n")
    for order in (2, 4, 6):
        prev = None
        print(f"  VISC_ORDER = {order}")
        for n in (64, 128, 256, 512, 1024):
            err = err_of(n, order)
            rate = "" if prev is None else f"   rate = {np.log2(prev / err):.2f}"
            print(f"    n={n:5d}  Linf = {err:.4e}{rate}")
            prev = err
        print()


def main():
    def flux_err(n, order):
        u, mu, dx, flux_exact, _ = field(n)
        return np.abs(visc_flux(u, mu, dx, order) - flux_exact).max()

    def div_err(n, order):
        u, mu, dx, _, div_exact = field(n)
        F = visc_flux(u, mu, dx, order)
        return np.abs((F - _s(F, -1)) / dx - div_exact).max()

    table("1. FACE FLUX  tau_xx = (4/3) mu du/dx  at x_{i+1/2}",
          "expect 2 / 4 / 6 -- this is the coefficient check", flux_err)
    table("2. DIVERGENCE  (F_{i+1/2} - F_{i-1/2})/dx,  what calc_R forms",
          "expect 2 / 2 / 2 at every width -- see this file's docstring", div_err)


if __name__ == "__main__":
    main()
