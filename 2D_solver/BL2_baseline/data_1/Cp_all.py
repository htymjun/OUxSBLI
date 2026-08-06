import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
from myvtk import getGrid, getQ, myParams

myParams()

R      = 287.15
gamma  = 1.4
M0     = 2.15
p_tot  = 25.0e3
T0     = 288.15           # [K]  STATIC freestream temperature
blt    = 1.0e-3           # length scale [m], used to convert x to mm

Xsh_mm = 80.0             # shock impingement location [mm] (Re_Xsh definition point, BL2 spec)

# y-axis zoom range for Cp — matches the ~0 to ~0.17 range implied by
# Vila-Pérez Fig. 28(a) / Degrez's ~1.55 overall pressure ratio.
# A small negative margin allows for mild numerical undershoot.
# Set to None to auto-scale instead.
Y_ZOOM_RANGE = (0.0, 0.20)

# font sizes: adjust these to resize the axis numbers (tick labels) and the
# axis titles (Cp, X/Xsh) independently
TICK_LABEL_FONTSIZE = 25
AXIS_LABEL_FONTSIZE = 25
LEGEND_FONTSIZE = 16.5

p0   = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0   = M0 * np.sqrt(gamma * R * T0)
q_inf = 0.5 * rho0 * u0**2   # freestream dynamic pressure (Cp normalization)

def compute_cp(path):
    """
    Compute the bottom-wall (y=0) pressure coefficient
    Cp = (p - p_inf) / (0.5 * rho_inf * v_inf^2)

    Returns
    -------
    x_mm : ndarray  x-coordinate [mm], x=0 assumed at the plate leading edge
    cp   : ndarray  pressure coefficient Cp
    """
    Nx, Ny, Nz, x, y, z = getGrid(path)
    rho, u, v, w, p = getQ(path, Nx, Ny, Nz)

    p2d = p[0, :, :]      # shape (Ny, Nx), extract the 2D (z=0) plane
    p_w = p2d[0, :]       # wall (j=0) pressure

    cp = (p_w - p0) / q_inf
    x_mm = x / blt
    return x_mm, cp

# ============================================================
# 2. Plot (single snapshot, or overlay multiple cases)
# ============================================================

# Digitized reference data (columns: X/Xsh, p/p0 or x/L, p/p0 depending on
# source — see x_scale note in the module docstring above), whitespace-
# separated, no header.

REFERENCE_DATA = [
    {
        "path": "./data_Degrez_Cp.dat",
        "label": "Degrez et al. (numerical)",
        "style": dict(marker="o", ls="none", mfc="none", color="k", ms=5),
        "x_scale": 1.25,  # digitized from Vila-Pérez Fig. 28(a): x/L -> X/Xsh
    },
    {
        "path": "./data_Moro_Cp.dat",
        "label": "Moro et al.",
        "style": dict(marker="x", ls="none", color="tab:red", ms=6, mew=1.2),
        "x_scale": 1.25,  # digitized from Vila-Pérez Fig. 28(a): x/L -> X/Xsh
    },
    {
        "path": "./data_Vila-Pérez_Cp.dat",
        "label": "Vila-Pérez et al.",
        "style": dict(marker="^", ls=":", color="tab:green", lw=1.0, ms=4),
        "x_scale": 1.25,  # digitized from Vila-Pérez Fig. 28(a): x/L -> X/Xsh
    },
]

if __name__ == "__main__":

    # --- add more entries here to overlay multiple files (e.g. grid
    #     sensitivity study, or different time steps) ---
    # optional "x_range": (min, max) restricts the plotted range for that
    # case only — reference datasets are unaffected.
    cases = [
        {"path": "./Q00800.vtr", "label": "Present Study", "x_range": (0.125, 1.875)},
        # {"path": "./Q00600.vtr", "label": "t = step 600"},
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    for case in cases:
        x_mm, cp = compute_cp(case["path"])
        x_over_xsh = x_mm / Xsh_mm   # X/Xsh normalization

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
            ref = np.loadtxt(ref_entry["path"])   # whitespace-separated, no header
        except (FileNotFoundError, OSError):
            print(f"[warn] '{ref_entry['path']}' not found — skipping "
                  f"'{ref_entry['label']}' (digitize it first, e.g. with "
                  "WebPlotDigitizer, then re-run).")
            continue

        order = np.argsort(ref[:, 0])             # digitizer points may be slightly out of order
        ref = ref[order]
        ref[:, 0] = ref[:, 0] * ref_entry.get("x_scale", 1.0)   # X/Xsh
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
    fig.savefig("cp_BL2_all.png", dpi=200)
    print("Saved: cp_BL2_all.png")