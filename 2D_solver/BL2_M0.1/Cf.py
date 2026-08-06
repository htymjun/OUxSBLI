import re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
from myvtk import getGrid, getQ, myParams

myParams()

R      = 287.15
gamma  = 1.4
M0     = 0.1
p_tot  = 25.0e3
T0     = 288.15           # [K]  STATIC freestream temperature
blt    = 1.0e-3           # length scale [m], used to convert x to mm
Pr     = 0.72

Xsh_mm = 80.0             # reference length [mm] for normalization (shock impingement location)

p0   = p_tot / (1.0 + 0.5 * (gamma - 1.0) * M0**2) ** (gamma / (gamma - 1.0))
rho0 = p0 / (R * T0)
u0   = M0 * np.sqrt(gamma * R * T0)

q_inf = 0.5 * rho0 * u0**2   # freestream dynamic pressure (Cf normalization)

# ------------------------------------------------------------
# Sutherland's law: mu(T) = SUTHERLAND_C * T^1.5 / (T + S)
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
    cf   : ndarray  skin-friction coefficient (DNS)
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

# ============================================================
# 2. Theoretical laminar flat-plate skin friction
#    (adiabatic wall, Sutherland viscosity)
# ============================================================

# recovery factor (laminar, adiabatic wall) -- matches rf in mod_globals.f90
rf  = np.sqrt(Pr)
Taw = T0 * (1.0 + rf * 0.5 * (gamma - 1.0) * M0**2)

# --- (a) incompressible Blasius, properties evaluated at freestream T0 ---
mu0 = sutherland_mu(T0)

def cf_blasius_incompressible(x):
    """Cf(x) = 0.664 / sqrt(Re_x), properties at freestream conditions."""
    Re_x = rho0 * u0 * x / mu0
    return 0.664 / np.sqrt(Re_x)

# --- (b) Eckert reference-temperature method (adiabatic wall) ---
#   T* = Te + (0.5 + 0.22 r)(Taw - Te),   here Te = T0 (edge = freestream)
T_star   = T0 + (0.5 + 0.22 * rf) * (Taw - T0)
mu_star  = sutherland_mu(T_star)
rho_star = p0 / (R * T_star)

def cf_eckert(x):
    """Cf(x) via Eckert reference-temperature method (adiabatic wall)."""
    Re_x_star = rho_star * u0 * x / mu_star
    cf_star = 0.664 / np.sqrt(Re_x_star)
    return cf_star * (rho_star / rho0)

print(f"[theory] T0 = {T0:.3f} K,  Taw = {Taw:.3f} K,  T* = {T_star:.3f} K")
print(f"[theory] mu(T0) = {mu0:.4e} Pa s,  mu(T*) = {mu_star:.4e} Pa s")
print(f"[theory] compressibility correction on Cf (T* vs T0): "
      f"{100.0 * (np.sqrt(mu0 / mu_star) * (rho_star / rho0) - 1.0):+.3f} %")

# ============================================================
# 3. Plot: DNS Cf(x) vs. theoretical Cf(x)  (x-axis normalized by Xsh)
# ============================================================

if __name__ == "__main__":

    case = {"path": "./Q03000.vtr", "label": "Present Study"}
    X_MIN_MM = 5.0
    X_MAX_MM = 190.0

    fig, ax = plt.subplots(figsize=(8, 6))

    x_mm, cf_dns = compute_cf(case["path"])
    mask = (x_mm >= X_MIN_MM) & (x_mm <= X_MAX_MM)
    x_mm_plot = x_mm[mask]
    cf_dns_plot = cf_dns[mask]

    x_m = x_mm_plot * blt

    cf_th_incomp = cf_blasius_incompressible(x_m)
    cf_th_eckert = cf_eckert(x_m)

    # --- normalize x-axis by Xsh ---
    x_over_xsh = x_mm_plot / Xsh_mm

    ax.plot(x_over_xsh, cf_dns_plot, color="tab:blue", lw=2.0,
            label=case["label"])
    ax.plot(x_over_xsh, cf_th_eckert, color="k", ls="--", lw=2.0,
            label="Eckert's Reference\nTemperature Method")
    ax.plot(x_over_xsh, cf_th_incomp, color="tab:red", ls=":", lw=2.0,
            label="Incompressible Blasius")  # "Incompressible Blasius (properties at $T_0$)"

    # --- relative error (DNS vs Eckert theory) ---
    rel_err = (cf_dns_plot - cf_th_eckert) / cf_th_eckert * 100.0
    print(f"[compare] max |relative error| DNS vs Eckert theory: "
          f"{np.max(np.abs(rel_err)):.2f} %  "
          f"(mean {np.mean(np.abs(rel_err)):.2f} %) "
          f"over X/Xsh in [{x_over_xsh[0]:.3f}, {x_over_xsh[-1]:.3f}]")

    ax.set_xlim(0, 2)
    ax.set_ylim(bottom=0)

    ax.set_xlabel(r"$X / X_{sh}$", fontsize=25)
    ax.set_ylabel(r"$C_f$", fontsize=25)
    ax.tick_params(axis="both", which="major", labelsize=25)

    # x-axis ticks every 0.5, integer values shown without a decimal point
    # (0.0 -> 0, 0.5 -> 0.5, 1.0 -> 1, 1.5 -> 1.5, 2.0 -> 2)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
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
    ax.yaxis.get_offset_text().set_fontsize(25)

    ax.legend(frameon=False, fontsize=18, loc="upper right")
    fig.tight_layout()
    fig.savefig("cf_BL2_theory_compare.png", dpi=200)
    print("Saved: cf_BL2_theory_compare.png")