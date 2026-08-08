"""Incompressible Blasius flat-plate boundary-layer solution.

The similarity ODE  f''' + 0.5 f f'' = 0  with  f(0) = f'(0) = 0  is
integrated with classical RK4 using the well-established shooting constant
f''(0) = 0.33205733622 (Howarth), so no scipy is needed. The resulting dense
table satisfies f'(inf) -> 1 to ~1e-7.

Note: the tabulated ``fp_tab`` in src/set_compressible_bl.f90 is an
approximate profile used only to seed initial conditions — it deviates from
the true Blasius solution by up to ~7% for eta >= 2.4 and must not be used as
validation truth.

Similarity variable: eta = y * sqrt(u_inf / (nu_inf * x)), with x measured
from the plate leading edge.
"""

import numpy as np

_FPP0 = 0.33205733622  # f''(0), Howarth's value
_ETA_MAX = 10.0
_DETA = 0.01


def _integrate():
    n = int(round(_ETA_MAX / _DETA)) + 1
    eta = np.linspace(0.0, _ETA_MAX, n)
    F = np.empty((n, 3))  # columns: f, f', f''
    F[0] = (0.0, 0.0, _FPP0)

    def rhs(s):
        return np.array([s[1], s[2], -0.5 * s[0] * s[2]])

    h = _DETA
    for i in range(n - 1):
        s = F[i]
        k1 = rhs(s)
        k2 = rhs(s + 0.5 * h * k1)
        k3 = rhs(s + 0.5 * h * k2)
        k4 = rhs(s + h * k3)
        F[i + 1] = s + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return eta, F[:, 0], F[:, 1]


ETA_TAB, F_TAB, FP_TAB = _integrate()


def fprime(eta):
    """u/u_inf = f'(eta); clamps to 1 beyond the table (eta > 10)."""
    return np.interp(eta, ETA_TAB, FP_TAB)


def eta(y, x, u_inf, nu_inf):
    """Blasius similarity variable for wall distance y at station x (from the LE)."""
    return y * np.sqrt(u_inf / (nu_inf * x))


def cf_blasius(x, rho_inf, u_inf, mu_inf):
    """Blasius skin friction Cf(x) = 0.664 / sqrt(Re_x), x measured from the LE."""
    Re_x = rho_inf * u_inf * x / mu_inf
    return 2.0 * _FPP0 / np.sqrt(Re_x)
