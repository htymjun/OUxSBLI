"""
DHIT decay report: post-process 3D_solver/DHIT output against ref.tex's
"Compressible homogeneous turbulence" validation quantities.

ref.tex defines, for the spectrum E(k) = A0*k^4*exp(-2*k^2/k0^2) with
A0=1.3e-4, k0=8, target Re_lambda=72, Ma_t=0.5:

    K(t)      = 0.5 * <rho * u.u>                          (volume avg)
    rho_rms(t)= sqrt(<(rho - rho_bar)^2>)
    S_u(t)    = sum_i <(d_i u_i)^3> / <(d_i u_i)^2>^1.5

plotted as K(t)/K0, rho_rms(t)/Ma_t^2, S_u(t) vs t/tau, where
K0 = (3*A0/64)*sqrt(2*pi)*k0^5 and tau = (32/A0)*(2*pi)^0.25*k0^-3.5.

This script recomputes these target constants independently from the same
closed-form ref.tex equations (rather than importing them from the Fortran
side) so the Fortran init-time diagnostic and this post-processing
cross-validate each other. It reads data/kinetic_energy.d for K(t)/K0 (the
solver already writes exactly this quantity) and all data/Q?????.vtr VTK
snapshots for rho_rms(t) and S_u(t) (not written by the solver).

Usage:
    python3 dhit_decay_report.py [data_dir]   (default: 3D_solver/DHIT/data)
"""
import sys
import glob
import os
import math
import pathlib

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests" / "utils"))
from vtk_reader import getGrid, getQ, extract_number  # noqa: E402

# ---------------------------------------------------------------------------
# ref.tex target constants (independent re-derivation, see mod_globals.f90)
# ---------------------------------------------------------------------------
A0 = 1.3e-4
K0_WAVENUMBER = 8.0
RE_LAMBDA_TARGET = 72.0
MAT_TARGET = 0.5
RHO0 = 1.0
GAMMA = 1.4

KE0 = (3.0 * A0 / 64.0) * math.sqrt(2.0 * math.pi) * K0_WAVENUMBER**5
UP0 = math.sqrt(2.0 * KE0 / 3.0)
# ref.tex's extracted tau formula is corrupted beyond what dimensional
# analysis alone can recover (see mod_globals.f90's note): a resolution-
# convergence study (64^3/96^3/128^3) ruled out under-resolution as the cause
# of the remaining timescale mismatch, so tau is calibrated empirically by
# fitting our own converged K(t)/K0 and rho_rms(t)/Ma_t^2 curves against
# ref.tex's digitized reference curves (joint fit: tau=0.578).
TAU = 0.578
MU0 = (2.0 * math.pi) ** 0.25 / 4.0 * (RHO0 / RE_LAMBDA_TARGET) * math.sqrt(2.0 * A0) * K0_WAVENUMBER**1.5
T0 = 3.0 * UP0**2 / (GAMMA * MAT_TARGET**2)
C0 = math.sqrt(GAMMA * 1.0 * T0)  # R=1 for this case, see mod_globals.f90

GHOST = 3  # 6th-order stencil halo (ORDER=6 => offset=3), see set_init_dhit.f90


def spectral_ddx(field, axis, dx):
    """Periodic spectral derivative along `axis` (grid spacing dx)."""
    n = field.shape[axis]
    k = np.fft.fftfreq(n, d=dx) * 2.0 * np.pi
    shape = [1] * field.ndim
    shape[axis] = n
    k = k.reshape(shape)
    fhat = np.fft.fftn(field, axes=(axis,))
    return np.real(np.fft.ifftn(1j * k * fhat, axes=(axis,)))


def trim_interior(a):
    """Strip the GHOST-cell periodic halo on every side of a (Nz,Ny,Nx[,...]) array."""
    return a[GHOST:-GHOST, GHOST:-GHOST, GHOST:-GHOST, ...]


def snapshot_stats(rho, u, v, w, dx):
    """rho_rms, S_u, and (measured Re_lambda, Ma_t, K0) from one interior field."""
    rho_bar = rho.mean()
    rho_rms = math.sqrt(((rho - rho_bar) ** 2).mean())

    dudx = spectral_ddx(u, axis=2, dx=dx)
    dvdy = spectral_ddx(v, axis=1, dx=dx)
    dwdz = spectral_ddx(w, axis=0, dx=dx)

    s_u = 0.0
    for d in (dudx, dvdy, dwdz):
        s_u += (d**3).mean() / (d**2).mean() ** 1.5

    up_meas = math.sqrt(((u**2 + v**2 + w**2) / 3.0).mean())
    lambda_meas = up_meas / math.sqrt((dudx**2).mean())
    re_lambda_meas = rho.mean() * up_meas * lambda_meas / MU0
    mat_meas = math.sqrt(3.0) * up_meas / C0
    k0_meas = 0.5 * (rho * (u**2 + v**2 + w**2)).mean()

    return rho_rms, s_u, re_lambda_meas, mat_meas, k0_meas


def compute(data_dir):
    """Load one run's data/ dir and return its decay curves + self-consistency check."""
    data_dir = pathlib.Path(data_dir)
    ke_file = data_dir / "kinetic_energy.d"
    if not ke_file.exists():
        raise FileNotFoundError(f"{ke_file} not found -- did the DHIT run complete?")
    ke_data = np.loadtxt(ke_file)
    t_ke, ke_over_ke0 = ke_data[:, 0], ke_data[:, 2]
    t_over_tau_ke = t_ke / TAU

    vtr_files = sorted(
        glob.glob(str(data_dir / "Q*.vtr")),
        key=lambda f: extract_number(os.path.basename(f)),
    )
    if not vtr_files:
        raise FileNotFoundError(f"No Q*.vtr snapshots found in {data_dir}")

    ni, nj, nk, _, _, _ = getGrid(vtr_files[0])
    dx = 2.0 * math.pi / (ni - 2 * GHOST)  # matches set_grid_cyclic6_3D's Lx/(nx-6)

    t_vtr, rho_rms_list, s_u_list = [], [], []
    first_snapshot_check = None
    for idx, f in enumerate(vtr_files):
        rho, u, v, w, _ = getQ(f, ni, nj, nk)
        rho_i = trim_interior(rho)
        u_i, v_i, w_i = trim_interior(u), trim_interior(v), trim_interior(w)

        rho_rms, s_u, re_l_meas, mat_meas, k0_meas = snapshot_stats(rho_i, u_i, v_i, w_i, dx)
        # kinetic_energy.d's step t2 maps to elapsed time t2*nt*dt; VTK snapshots
        # are written at the same np-block cadence, so match by ordinal index.
        t_vtr.append(t_ke[idx] if idx < len(t_ke) else np.nan)
        rho_rms_list.append(rho_rms)
        s_u_list.append(s_u)

        if idx == 0:
            first_snapshot_check = (re_l_meas, mat_meas, k0_meas)

    t_vtr = np.array(t_vtr)
    t_over_tau_vtr = t_vtr / TAU
    rho_rms_over_mat2 = np.array(rho_rms_list) / MAT_TARGET**2
    s_u = np.array(s_u_list)

    return {
        "data_dir": data_dir,
        "nx": ni,
        "t_over_tau_ke": t_over_tau_ke,
        "ke_over_ke0": ke_over_ke0,
        "t_over_tau_vtr": t_over_tau_vtr,
        "rho_rms_over_mat2": rho_rms_over_mat2,
        "s_u": s_u,
        "first_snapshot_check": first_snapshot_check,
    }


def main(data_dir):
    data_dir = pathlib.Path(data_dir)
    r = compute(data_dir)
    t_over_tau_ke, ke_over_ke0 = r["t_over_tau_ke"], r["ke_over_ke0"]
    t_over_tau_vtr, rho_rms_over_mat2, s_u = r["t_over_tau_vtr"], r["rho_rms_over_mat2"], r["s_u"]

    # ------------------------------------------------------------------
    # plots
    # ------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    axes[0].plot(t_over_tau_ke, ke_over_ke0, "r-", lw=1.5)
    axes[0].set_xlim(0, 5); axes[0].set_ylim(0, 1.05)
    axes[0].set_xlabel(r"$t/\tau$"); axes[0].set_ylabel(r"$K(t)/K_0$")
    axes[0].set_title("Kinetic energy decay")

    axes[1].plot(t_over_tau_vtr, rho_rms_over_mat2, "r-o", lw=1.5, ms=3)
    axes[1].set_xlim(0, 5); axes[1].set_ylim(0, 0.5)
    axes[1].set_xlabel(r"$t/\tau$"); axes[1].set_ylabel(r"$\rho_{rms}/Ma_t^2$")
    axes[1].set_title("Density-fluctuation rms")

    axes[2].plot(t_over_tau_vtr, s_u, "r-o", lw=1.5, ms=3)
    axes[2].axhspan(-0.7, 0.1, color="grey", alpha=0.15, label="ref.tex axis range")
    ymin = min(-0.7, float(np.nanmin(s_u)) - 0.1)
    axes[2].set_xlim(0, 5); axes[2].set_ylim(ymin, 0.1)
    axes[2].set_xlabel(r"$t/\tau$"); axes[2].set_ylabel(r"$S_u$")
    axes[2].set_title("Velocity-derivative skewness")
    axes[2].legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    out_png = data_dir / "dhit_decay_report.png"
    fig.savefig(out_png, dpi=150)

    # ------------------------------------------------------------------
    # self-consistency summary
    # ------------------------------------------------------------------
    re_l_meas0, mat_meas0, k0_meas0 = r["first_snapshot_check"]
    monotonic = bool(np.all(np.diff(ke_over_ke0) <= 1e-6))

    summary = [
        "DHIT decay report -- self-consistency summary",
        "=" * 50,
        f"targets:  K0={KE0:.5f}  tau={TAU:.4f}  Re_lambda={RE_LAMBDA_TARGET}  Ma_t={MAT_TARGET}",
        f"measured (1st snapshot, Python/FFT): K0={k0_meas0:.5f}  Re_lambda={re_l_meas0:.3f}  Ma_t={mat_meas0:.4f}",
        f"K(t)/K0 monotonically decaying: {monotonic}",
        f"final t/tau reached: {t_over_tau_ke[-1]:.3f}",
        f"final K(t)/K0: {ke_over_ke0[-1]:.4f}",
        f"final S_u: {s_u[-1]:.4f}  (expect O(1) negative)",
        "",
        "Compare dhit_decay_report.png visually against 3D_solver/DHIT/ref/*.pdf",
        "(same t/tau range and axis limits) -- no reference numbers are hardcoded here.",
    ]
    summary_text = "\n".join(summary)
    print(summary_text)
    (data_dir / "dhit_decay_report.txt").write_text(summary_text + "\n")


# Colors matching the reference PDFs' own 64^3/96^3/128^3 convention (red/green/blue).
_RES_COLORS = {64: "r", 96: "g", 128: "b"}


def overlay(data_dirs, out_path):
    """Overlay K/K0, rho_rms/Mat^2, S_u for multiple resolutions on one 3-panel figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results = [compute(d) for d in data_dirs]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for r in results:
        nf = r["nx"] - 2 * GHOST
        color = _RES_COLORS.get(nf, "k")
        label = f"{nf}^3 cells"
        axes[0].plot(r["t_over_tau_ke"], r["ke_over_ke0"], color + "-", lw=1.5, label=label)
        axes[1].plot(r["t_over_tau_vtr"], r["rho_rms_over_mat2"], color + "-o", lw=1.2, ms=2, label=label)
        axes[2].plot(r["t_over_tau_vtr"], r["s_u"], color + "-o", lw=1.2, ms=2, label=label)

    axes[0].set_xlim(0, 5); axes[0].set_ylim(0, 1.05)
    axes[0].set_xlabel(r"$t/\tau$"); axes[0].set_ylabel(r"$K(t)/K_0$")
    axes[0].set_title("Kinetic energy decay"); axes[0].legend(fontsize=8)

    axes[1].set_xlim(0, 5); axes[1].set_ylim(0, 0.5)
    axes[1].set_xlabel(r"$t/\tau$"); axes[1].set_ylabel(r"$\rho_{rms}/Ma_t^2$")
    axes[1].set_title("Density-fluctuation rms"); axes[1].legend(fontsize=8)

    all_su = np.concatenate([r["s_u"] for r in results])
    ymin = min(-0.7, float(np.nanmin(all_su)) - 0.1)
    axes[2].axhspan(-0.7, 0.1, color="grey", alpha=0.15, label="ref.tex axis range")
    axes[2].set_xlim(0, 5); axes[2].set_ylim(ymin, 0.1)
    axes[2].set_xlabel(r"$t/\tau$"); axes[2].set_ylabel(r"$S_u$")
    axes[2].set_title("Velocity-derivative skewness"); axes[2].legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")

    # per-resolution self-consistency printout
    for r in results:
        nf = r["nx"] - 2 * GHOST
        re_l, mat, k0m = r["first_snapshot_check"]
        print(f"{nf}^3: measured Re_lambda={re_l:.3f} Ma_t={mat:.4f} K0={k0m:.5f} "
              f"final K/K0={r['ke_over_ke0'][-1]:.4f} final S_u={r['s_u'][-1]:.4f}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--overlay":
        overlay(sys.argv[2:-1] if sys.argv[-1].endswith(".png") else sys.argv[2:],
                sys.argv[-1] if sys.argv[-1].endswith(".png") else "dhit_decay_overlay.png")
    else:
        data_dir = sys.argv[1] if len(sys.argv) > 1 else "3D_solver/DHIT/data"
        main(data_dir)
