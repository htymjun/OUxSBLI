"""Overlay the DNS skin-friction distribution on the digitized reference data.

Usage:  python Cf_all.py [snapshot.vtr]
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
from ouxsbli.analysis.wall import find_zero_crossings, load_reference, wall_coeffs_from_vtr
from ouxsbli.tests.utils.vtk_reader import latest_vtr

plot_style.apply()

# must match 2D_solver/SBLI/mod_globals.f90
R = 287.15
gamma = 1.4
M0 = 2.15
p_tot = 25.0e3
T0 = 288.15  # [K] STATIC freestream temperature
blt = 1.0e-3  # length scale [m], used to convert x to mm

Xsh_mm = 80.0  # shock impingement location [mm]

# y-axis zoom range, in RAW (unscaled) Cf units — matches roughly the
# range shown in Degrez et al. Fig. 15 (-1e-3 to 5e-3). This crops out
# the leading-edge singularity spike so the separation/reattachment
# region is actually visible. Set to None to auto-scale instead.
Y_ZOOM_RANGE = (-1.5e-3, 6.5e-3)

TICK_LABEL_FONTSIZE = 25   # size of the 0,1,2,...,0.25,0.50,... numbers
AXIS_LABEL_FONTSIZE = 25   # size of the "X/Xsh" and "Cf" axis titles
LEGEND_FONTSIZE = 18       # size of the legend entries (Present Study, Moro et al., ...)

p0 = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0 = M0 * np.sqrt(gamma * R * T0)
q_inf = 0.5 * rho0 * u0**2  # freestream dynamic pressure (Cf normalization)

REFERENCE_DATA = [
    {
        "path": REF_DIR / "data_Moro_Cf.dat",
        "label": "Moro et al.",
        "style": dict(marker="x", ls="none", color="tab:red", ms=6, mew=2.0),
        "x_scale": 1.25,  # digitized from Vila-Perez Fig. 28(a): x/L -> X/Xsh
    },
    {
        "path": REF_DIR / "data_Degrez_Cf.dat",
        "label": "Degrez et al. (numerical)",
        "style": dict(marker="o", ls="none", mfc="none", color="k", ms=6, mew=2.0),
        "x_scale": 1.25,
    },
    # {
    #     "path": REF_DIR / "data_Vila-Perez_Cf.dat",
    #     "label": "Vila-Perez et al.",
    #     "style": dict(marker="^", ls=":", color="tab:green", lw=2.0, ms=5),
    #     "x_scale": 1.25,
    # },
]

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else latest_vtr(CASE_DIR / "data")

    cases = [
        {"path": path, "label": "Present Study", "x_range": (0.125, 1.875)},
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    for case in cases:
        x, cf, _ = wall_coeffs_from_vtr(case["path"], R, p0, q_inf)
        x_over_xsh = (x / blt) / Xsh_mm   # X/Xsh normalization

        x_range = case.get("x_range")
        if x_range is not None:
            mask = (x_over_xsh >= x_range[0]) & (x_over_xsh <= x_range[1])
            x_over_xsh = x_over_xsh[mask]
            cf = cf[mask]

        ax.plot(x_over_xsh, cf, label=case["label"], lw=1.5)

        crossings = find_zero_crossings(x_over_xsh, cf)
        if crossings:
            print(f"[{case['label']}] Cf=0 crossings [X/Xsh]: "
                  + ", ".join(f"{xc:.3f}" for xc in crossings))
        else:
            print(f"[{case['label']}] Cf=0 crossings: none found "
                  "(no separation, or out of data range)")
        imin = np.argmin(cf)
        print(f"[{case['label']}] Cf_min = {cf[imin]:.4e} at X/Xsh = {x_over_xsh[imin]:.3f}")

    # --- overlay all reference / comparison datasets ---
    for ref_entry in REFERENCE_DATA:
        ref = load_reference(ref_entry["path"], ref_entry.get("x_scale", 1.0))
        ax.plot(ref[:, 0], ref[:, 1], label=ref_entry["label"], **ref_entry["style"])

        label = ref_entry["label"]
        crossings_ref = find_zero_crossings(ref[:, 0], ref[:, 1])
        imin_ref = np.argmin(ref[:, 1])
        if crossings_ref:
            print(f"[{label}] Cf=0 crossings [X/Xsh]: "
                  + ", ".join(f"{xc:.3f}" for xc in crossings_ref))
        else:
            print(f"[{label}] Cf=0 crossings: none found "
                  "(no separation, or out of data range)")
        print(f"[{label}] Cf_min = {ref[imin_ref, 1]:.4e} "
              f"at X/Xsh = {ref[imin_ref, 0]:.3f}")

    ax.axhline(0.0, color="k", ls="--", lw=0.8)
    ax.set_xlim(0, 2.0)
    ax.set_xlabel(r"$X / X_{sh}$", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(r"$C_f$", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)

    # zoom the y-axis to the interaction region's scale (matches Fig. 15's
    # roughly -1 to 5 x1e-3 range). This clips the leading-edge singularity
    # spike near X/Xsh=0 off the top of the plot, which is expected — that
    # spike isn't the region of interest for the SWBLI comparison.
    if Y_ZOOM_RANGE is not None:
        ax.set_ylim(Y_ZOOM_RANGE[0], Y_ZOOM_RANGE[1])
    # useMathText=True renders the offset text as "x10^-3" with a proper
    # multiplication sign / minus sign instead of the default "1e-3" style.
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3), useMathText=True)

    # major ticks every 1 (x1e-3) so labels read 0,1,2,3,4,5,6 next to the
    # x10^-3 offset text
    ax.yaxis.set_major_locator(MultipleLocator(1e-3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.get_offset_text().set_fontsize(TICK_LABEL_FONTSIZE)  # the "x10^-3" text

    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))

    ax.legend(frameon=False, fontsize=LEGEND_FONTSIZE, loc="upper right",
              bbox_to_anchor=(1.0, 1.0), borderaxespad=0.3)
    fig.tight_layout()
    fig.savefig(REF_DIR / "cf_BL2_all.png", dpi=200)
    print("Saved: cf_BL2_all.png")
