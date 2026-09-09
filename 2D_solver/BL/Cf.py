"""Plot the DNS skin-friction distribution against the Blasius flat-plate
solution.

Usage:  python Cf.py [snapshot.vtr]
Without an argument the latest snapshot in ./data is used.
"""

import pathlib
import re
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter

CASE_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(CASE_DIR.parents[1]))

from ouxsbli.analysis import plot_style
from ouxsbli.analysis.blasius import cf_blasius
from ouxsbli.analysis.wall import edge_state, sutherland_mu, wall_coeffs_from_vtr
from ouxsbli.tests.utils.vtk_reader import getGrid, getQ, latest_vtr

plot_style.apply()

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif"]
plt.rcParams["mathtext.fontset"] = "stix"

# must match 2D_solver/BL/mod_globals.f90
R = 287.15
gamma = 1.4
M0 = 0.1
p_tot = 25.0e3
T0 = 288.15  # [K] STATIC freestream temperature
blt = 1.0e-3  # length scale [m], used to convert x to mm

p0 = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0 = M0 * np.sqrt(gamma * R * T0)
q_inf = 0.5 * rho0 * u0**2  # freestream dynamic pressure (Cf normalization)
mu0 = sutherland_mu(T0)

# x = 0 is the plate leading edge (set_grid places it there exactly)
# The upper limit stays clear of the 0th-order-extrapolated outlet at x = 100.3 mm.
X_MIN_MM = 5.0
X_MAX_MM = 95.0
Y_MAJOR_INTERVAL = 2.5e-3

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else latest_vtr(CASE_DIR / "data")
    x, cf_dns, _ = wall_coeffs_from_vtr(path, R, p0, q_inf)
    x_mm = x / blt

    mask = (x_mm >= X_MIN_MM) & (x_mm <= X_MAX_MM)
    x_mm_plot = x_mm[mask]
    cf_dns_plot = cf_dns[mask]
    cf_theory = cf_blasius(x_mm_plot * blt, rho0, u0, mu0)

    err = np.abs(cf_dns_plot / cf_theory - 1.0)
    print(f"snapshot: {path}")
    print(f"Cf vs Blasius over x = {X_MIN_MM}-{X_MAX_MM} mm (freestream-normalised): "
          f"mean rel. error {err.mean():.3f}, max {err.max():.3f}")

    # Normalising by the local edge state instead (as ouxsbli/tests/test_bl.py does)
    # divides out any residual acceleration of the outer flow. The gap between the
    # two curves is therefore a direct readout of how much work the far-field
    # boundary is doing on the freestream: it should be near zero.
    ni, nj, nk, xg, yg, zg = getGrid(str(path))
    rho, uu, vv, ww, pp = getQ(str(path), ni, nj, nk)
    u_e, rho_e, mu_e = edge_state(rho[0], uu[0], pp[0], R)
    Re_x_loc = rho_e * u_e * x / mu_e
    cf_loc = cf_dns * (q_inf / (0.5 * rho_e * u_e**2))
    with np.errstate(divide="ignore", invalid="ignore"):
        err_loc = np.abs(cf_loc[mask] / (0.664 / np.sqrt(Re_x_loc[mask])) - 1.0)
    print(f"{'':>4}                                   (edge-normalised): "
          f"mean rel. error {err_loc.mean():.3f}, max {err_loc.max():.3f}")
    print(f"peak edge acceleration u_e/u0 - 1: {u_e[mask].max() / u0 - 1.0:+.4f}")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_mm_plot, cf_dns_plot, color="tab:blue", lw=2.0, label="Present Study")
    #ax.plot(x_mm_plot, cf_loc[mask], color="tab:green", ls=":", lw=2.0,
            #label="Present Study (edge-norm.)")
    ax.plot(x_mm_plot, cf_theory, color="tab:red", ls="--", lw=2.0, label="Blasius")

    ax.set_xlim(0, X_MAX_MM)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$x$ [mm]", fontsize=28)
    ax.set_ylabel(r"$C_f$", fontsize=28)
    ax.tick_params(axis="both", which="major", labelsize=28)

    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{val:g}"))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # y-axis: hide the 0 label, and strip the trailing ".0" from other labels
    # (1.0 -> 1, while 1.5 is left unchanged)
    class ZeroHiddenScalarFormatter(matplotlib.ticker.ScalarFormatter):
        def __call__(self, x, pos=None):
            s = super().__call__(x, pos)
            if np.isclose(x, 0.0):
                return ""
            return re.sub(r"\.0(?=\D*$)", "", s)

    y_formatter = ZeroHiddenScalarFormatter(useMathText=True)
    y_formatter.set_scientific(True)
    y_formatter.set_powerlimits((-3, -3))
    ax.yaxis.set_major_formatter(y_formatter)
    ax.yaxis.get_offset_text().set_fontsize(28)

    if Y_MAJOR_INTERVAL is not None:
        ax.yaxis.set_major_locator(MultipleLocator(Y_MAJOR_INTERVAL))

    ax.legend(frameon=False, fontsize=24, loc="upper right")
    fig.tight_layout()
    fig.savefig(CASE_DIR / "cf_BL_theory_compare.png", dpi=200)
    print("Saved: cf_BL_theory_compare.png")
