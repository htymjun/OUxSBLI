import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
from myvtk import getGrid, getQ, myParams

myParams()

# ============================================================
# 1. BL2 case parameters (must match mod_globals.f90)
# ============================================================
R      = 287.15          # specific gas constant [J/(kg K)]  <- BL2 value
gamma  = 1.4
M0     = 2.15
p_tot  = 25.0e3           # [Pa]
T0     = 288.15           # [K]  STATIC freestream temperature (given directly in mod_globals.f90)
blt    = 1.0e-3           # length scale [m], used to convert x to mm

Xsh_mm = 80.0             # shock impingement location [mm] (Re_Xsh definition point, BL2 spec)

# y-axis zoom range, in RAW (unscaled) Cf units — matches roughly the
# range shown in Degrez et al. Fig. 15 (-1e-3 to 5e-3). This crops out
# the leading-edge singularity spike so the separation/reattachment
# region is actually visible. Set to None to auto-scale instead.
Y_ZOOM_RANGE = (-1.5e-3, 6.5e-3)

# font sizes: adjust these to resize the axis numbers (tick labels) and the
# axis titles (Cf, X/Xsh) independently
TICK_LABEL_FONTSIZE = 25   # size of the 0,1,2,...,0.25,0.50,... numbers
AXIS_LABEL_FONTSIZE = 25   # size of the "X/Xsh" and "Cf" axis titles
LEGEND_FONTSIZE = 18        # size of the legend entries (Present Study, Degrez et al., ...)

p0   = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0   = M0 * np.sqrt(gamma * R * T0)

q_inf = 0.5 * rho0 * u0**2   # freestream dynamic pressure (Cf normalization)

# ------------------------------------------------------------
# Sutherland's law: mu(T) = SUTHERLAND_C * T^1.5 / (T + S)
#   SUTHERLAND_C = mu_ref * (T_ref + S) / T_ref^1.5
# NOTE: verify against mod_constant.f90's mu0_T0_S_over_T0_2_3 (see docstring).
# ------------------------------------------------------------
SUTHERLAND_MU_REF = 1.716e-5   # [Pa s] reference viscosity at T_ref
SUTHERLAND_T_REF  = 273.2       # [K]
SUTHERLAND_S      = 111.0       # [K]  Sutherland constant
SUTHERLAND_C = SUTHERLAND_MU_REF * (SUTHERLAND_T_REF + SUTHERLAND_S) / SUTHERLAND_T_REF**1.5


def sutherland_mu(T):
    """Dynamic viscosity mu(T) [Pa s] from Sutherland's law."""
    return SUTHERLAND_C * T**1.5 / (T + SUTHERLAND_S)


def one_sided_deriv_3pt(y0, y1, y2, f0, f1, f2):
    """
    Derivative at y0 of the quadratic polynomial passing through
    (y0,f0), (y1,f1), (y2,f2). Valid for non-uniform grid spacing
    (one-sided 3-point finite difference).
    """
    denom0 = (y0 - y1) * (y0 - y2)
    denom1 = (y1 - y0) * (y1 - y2)
    denom2 = (y2 - y0) * (y2 - y1)
    return (
        f0 * ((y0 - y1) + (y0 - y2)) / denom0
        + f1 * (y0 - y2) / denom1
        + f2 * (y0 - y1) / denom2
    )


def compute_cf(path):
    """
    Compute the bottom-wall (y=0) skin-friction coefficient Cf(x)
    from a VTR snapshot file.

    Returns
    -------
    x_mm : ndarray  x-coordinate [mm], x=0 assumed at the plate leading edge
    cf   : ndarray  skin-friction coefficient
    """
    Nx, Ny, Nz, x, y, z = getGrid(path)
    rho, u, v, w, p = getQ(path, Nx, Ny, Nz)

    # extract the 2D (z=0) plane
    rho2d = rho[0, :, :]   # shape (Ny, Nx)
    u2d   = u[0, :, :]
    p2d   = p[0, :, :]

    # --- wall (j=0) temperature and viscosity ---
    rho_w = rho2d[0, :]
    p_w   = p2d[0, :]
    T_w   = p_w / (rho_w * R)
    mu_w  = sutherland_mu(T_w)

    # --- wall-normal velocity gradient du/dy|_wall (no-slip: u(y0)=0) ---
    y0, y1, y2 = y[0], y[1], y[2]
    u0_wall = np.zeros_like(u2d[0, :])   # u=0 at the wall
    u1_pt   = u2d[1, :]
    u2_pt   = u2d[2, :]

    dudy_wall = one_sided_deriv_3pt(y0, y1, y2, u0_wall, u1_pt, u2_pt)

    # --- wall shear stress and skin-friction coefficient ---
    tau_w = mu_w * dudy_wall
    cf = tau_w / q_inf

    x_mm = x / blt
    return x_mm, cf


def find_zero_crossings(x_over_xsh, cf):
    """Locate Cf=0 crossings (separation/reattachment points) by linear interpolation."""
    sign_change = np.where(np.diff(np.sign(cf)) != 0)[0]
    crossings = []
    for i in sign_change:
        x1, x2 = x_over_xsh[i], x_over_xsh[i + 1]
        f1, f2 = cf[i], cf[i + 1]
        xc = x1 - f1 * (x2 - x1) / (f2 - f1)
        crossings.append(xc)
    return crossings


# ============================================================
# 2. Plot (single snapshot, or overlay multiple cases)
# ============================================================

# Digitized Degrez et al. 1987 Fig. 15 data (columns: X/Xsh, Cf), whitespace-
# separated, no header. Set to a path if you digitize the curve (e.g. with
# WebPlotDigitizer) and want to overlay it directly. None disables overlay.

REFERENCE_DATA = [
    {
        "path": "./data_Degrez_Cf.dat",
        "label": "Degrez et al. (numerical)",
        "style": dict(marker="o", ls="--", color="k", lw=1.0, ms=3),
        "x_scale": 1.25,
    },
    {
        "path": "./data_Moro_Cf.dat",
        "label": "Moro et al.",
        "style": dict(marker="s", ls="-.", color="tab:red", lw=1.0, ms=3),
        "x_scale": 1.25,
    },
    {
        "path": "./data_Vila-Pérez_Cf.dat",
        "label": "Vila-Pérez et al.",
        "style": dict(marker="^", ls=":", color="tab:green", lw=1.0, ms=3),
        "x_scale": 1.25,
    },
]

if __name__ == "__main__":

    # --- add more entries here to overlay multiple files (e.g. grid
    #     sensitivity study, or different time steps) ---
    # optional "x_range": (min, max) restricts the plotted/printed X/Xsh
    # range for that case only (e.g. to exclude inflow/outflow buffer
    # regions outside the physical plate) — reference datasets are
    # unaffected.
    cases = [
        {"path": "./Q00800.vtr", "label": "Present Study", "x_range": (0.125, 1.875)},
        # {"path": "./Q00600.vtr", "label": "t = step 600"},
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    for case in cases:
        x_mm, cf = compute_cf(case["path"])
        x_over_xsh = x_mm / Xsh_mm   # X/Xsh normalization (matches Fig. 15)

        x_range = case.get("x_range")
        if x_range is not None:
            mask = (x_over_xsh >= x_range[0]) & (x_over_xsh <= x_range[1])
            x_over_xsh = x_over_xsh[mask]
            cf = cf[mask]

        ax.plot(x_over_xsh, cf, label=case["label"], lw=1.5)

        # zero crossings computed on the raw cf
        crossings = find_zero_crossings(x_over_xsh, cf)
        if crossings:
            print(f"[{case['label']}] Cf=0 crossings [X/Xsh]: "
                  + ", ".join(f"{xc:.3f}" for xc in crossings))
        else:
            print(f"[{case['label']}] Cf=0 crossings: none found "
                  "(no separation, or out of data range)")

    # --- overlay all reference / comparison datasets ---
    for ref_entry in REFERENCE_DATA:
        ref = np.loadtxt(ref_entry["path"])       # whitespace-separated, no header
        order = np.argsort(ref[:, 0])             # digitizer points may be slightly out of order
        ref = ref[order]
        ref[:, 0] = ref[:, 0] * ref_entry.get("x_scale", 1.0)   # X/Xsh補正(全データ1.25倍)
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
    ax.set_xlim(0, 2.0)          # matches Fig. 15's axis range
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
    fig.savefig("cf_BL2_all.png", dpi=200)
    print("Saved: cf_BL2_all.png")