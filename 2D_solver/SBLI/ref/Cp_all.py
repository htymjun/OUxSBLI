"""Overlay the DNS wall-pressure distribution on the digitized reference data.

Usage:  python Cp_all.py [snapshot.vtr]
Without an argument the latest snapshot in <case>/data is used.
"""

import pathlib
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter

REF_DIR = pathlib.Path(__file__).resolve().parent
CASE_DIR = REF_DIR.parent
sys.path.insert(0, str(REF_DIR.parents[2]))

from ouxsbli.analysis import plot_style
from ouxsbli.analysis.wall import load_reference, wall_coeffs_from_vtr
from ouxsbli.tests.utils.vtk_reader import latest_vtr

plot_style.apply()

# must match 2D_solver/SBLI/mod_globals.f90
R = 287.15
gamma = 1.4
M0 = 2.15
p_tot = 25.0e3
T0 = 288.15  # [K] STATIC freestream temperature
blt = 1.0e-3  # length scale [m], used to convert x to mm

Xsh_mm = 80.0  # shock impingement location [mm] (Re_Xsh definition point, BL2 spec)

# y-axis zoom range for Cp — matches the ~0 to ~0.17 range implied by
# Vila-Perez Fig. 28(a) / Degrez's ~1.55 overall pressure ratio.
# A small negative margin allows for mild numerical undershoot.
# Set to None to auto-scale instead.
Y_ZOOM_RANGE = (0.0, 0.20)

# font sizes: adjust these to resize the axis numbers (tick labels) and the
# axis titles (Cp, X/Xsh) independently
TICK_LABEL_FONTSIZE = 25
AXIS_LABEL_FONTSIZE = 25
LEGEND_FONTSIZE = 16.5

p0 = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0 = M0 * np.sqrt(gamma * R * T0)
q_inf = 0.5 * rho0 * u0**2  # freestream dynamic pressure (Cp normalization)

# Digitized reference data (columns: x/L, Cp), whitespace-separated, no header;
# x_scale converts the published abscissa to X/Xsh.
REFERENCE_DATA = [
    {
        "path": REF_DIR / "data_Moro_Cp.dat",
        "label": "Moro et al.",
        "style": dict(marker="x", ls="none", color="tab:red", ms=6, mew=2.0),
        "x_scale": 1.25,  # digitized from Vila-Perez Fig. 28(a): x/L -> X/Xsh
    },
    {
        "path": REF_DIR / "data_Degrez_Cp.dat",
        "label": "Degrez et al. (numerical)",
        "style": dict(marker="o", ls="none", mfc="none", color="k", ms=6, mew=2.0),
        "x_scale": 1.25,
    },
    # {
    #     "path": REF_DIR / "data_Vila-Perez_Cp.dat",
    #     "label": "Vila-Perez et al.",
    #     "style": dict(marker="^", ls=":", color="tab:green", lw=1.0, ms=4),
    #     "x_scale": 1.25,
    # },
]

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else latest_vtr(CASE_DIR / "data")

    # --- add more entries here to overlay multiple snapshots (e.g. grid
    #     sensitivity study, or different time steps) ---
    # optional "x_range": (min, max) restricts the plotted range for that
    # case only — reference datasets are unaffected.
    cases = [
        {"path": path, "label": "Present Study", "x_range": (0.125, 1.875)},
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    for case in cases:
        x, _, cp = wall_coeffs_from_vtr(case["path"], R, p0, q_inf)
        x_over_xsh = (x / blt) / Xsh_mm   # X/Xsh normalization

        x_range = case.get("x_range")
        if x_range is not None:
            mask = (x_over_xsh >= x_range[0]) & (x_over_xsh <= x_range[1])
            x_over_xsh = x_over_xsh[mask]
            cp = cp[mask]

        ax.plot(x_over_xsh, cp, label=case["label"], color="tab:blue", lw=2.0)

        print(f"[{case['label']}] Cp range: "
              f"min={cp.min():.3f} at X/Xsh={x_over_xsh[np.argmin(cp)]:.3f}, "
              f"max={cp.max():.3f} at X/Xsh={x_over_xsh[np.argmax(cp)]:.3f}, "
              f"plateau (last pt)={cp[-1]:.3f}")

    # --- overlay all reference / comparison datasets (skip gracefully if a
    #     digitized data file hasn't been created yet) ---
    for ref_entry in REFERENCE_DATA:
        try:
            ref = load_reference(ref_entry["path"], ref_entry.get("x_scale", 1.0))
        except (FileNotFoundError, OSError):
            print(f"[warn] '{ref_entry['path']}' not found — skipping "
                  f"'{ref_entry['label']}' (digitize it first, e.g. with "
                  "WebPlotDigitizer, then re-run).")
            continue

        ax.plot(ref[:, 0], ref[:, 1], label=ref_entry["label"], **ref_entry["style"])

    ax.axhline(0.0, color="k", ls="--", lw=0.6)   # undisturbed-freestream reference level Cp=0
    ax.set_xlim(0, 2.0)
    ax.set_xlabel(r"$X / X_{sh}$", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(r"$C_p$", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)

    if Y_ZOOM_RANGE is not None:
        ax.set_ylim(Y_ZOOM_RANGE[0], Y_ZOOM_RANGE[1])

    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    # "%g" strips unnecessary trailing zeros (e.g. 0.05 stays "0.05", 0.0 -> "0").
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: "" if abs(y) < 1e-12 else f"{y:g}"))

    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE, loc="upper left",
              bbox_to_anchor=(0.0, 1.0), borderaxespad=0.3)
    fig.tight_layout()
    fig.savefig(REF_DIR / "cp_BL2_all.png", dpi=200)
    print("Saved: cp_BL2_all.png")
