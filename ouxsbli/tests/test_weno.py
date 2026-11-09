"""
WENO5-Z face reconstruction (RECON='WENO') in FP64 and W32P64, with the SLAU
(HRSLAU2) or local Lax-Friedrichs (RIEMANN='LF') face flux.

WENO5-Z replaces calc_muscl's delta6 (src/calc_weno.f90.fypp), so it runs
inside the SLAU-family kernels (SCHEME='SLAU', ORDER=6).

What is gated, and why the solver gate is not 5th order:
  * test_weno5z_reconstruction_order -- calc_weno alone, on exact cell
    averages: FP64 must be 5th order for BOTH biases (a non-mirrored right
    bias drops v+ to 3rd order and nothing else shows it).
  * test_evc_weno_fp64_lf_convergence -- FP64 WENO-Z + LF on EVC, order >= 3.0.
    The SLAU/LF kernels reconstruct primitive POINT values and evaluate the
    flux from them, which leaves an O(h^2) error for a nonlinear flux whatever
    the reconstruction: on 2D EVC at n=256, FP64 WENO + SLAU (3.9528e-10) and
    MUSCL + SLAU (3.9526e-10) agree to four digits. Measured orders over
    64/128/256: WENO + LF 3.47, WENO + SLAU 3.18, MUSCL + SLAU 1.79
    (2026-11-09, RTX PRO 6000). 3.0 separates WENO from MUSCL with margin.
  * test_os_weno -- Rankine-Hugoniot within 3% (shared helper), boundary-
    aware kernels: WENO inside, MUSCL delta4 / first order next to the walls.
W32P64 has no grid-convergence gate (by decision); it is exercised by the
reconstruction driver (reported, not gated) and by the OS shock case.
"""
import os
import pathlib
import re
import shutil
import subprocess

import numpy as np
import pytest
from ouxsbli import Case
from .utils.vtk_reader import initial_vtr, latest_vtr, getGrid, getQ
from .utils import evc_common as evc
from .utils import os_common as osc

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
ACCURACY = 6
EVC_MIN_ORDER = 3.0
RECON_MIN_ORDER = 4.8
RECON_FIT_GRIDS = (64, 128, 256, 512)


def _weno_params(prec, riemann):
    return dict(scheme="slau", accuracy=ACCURACY, recon="weno",
                weno_prec=prec, riemann=riemann)


# ---------------------------------------------------------------------------
# calc_weno alone
# ---------------------------------------------------------------------------
def _nvfortran():
    exe = shutil.which("nvfortran")
    if exe is None:
        found = sorted(pathlib.Path("/opt/nvidia/hpc_sdk").glob(
            "Linux_x86_64/*/compilers/bin/nvfortran"), reverse=True)
        exe = str(found[0]) if found else None
    if exe is None:
        pytest.skip("nvfortran not found")
    return exe


@pytest.mark.integration
def test_weno5z_reconstruction_order(tmp_path):
    src = REPO_ROOT / "src" / "calc_weno.f90.fypp"
    drv = REPO_ROOT / "ouxsbli" / "tests" / "utils" / "fortran" / "weno_order.cuf"
    gen = tmp_path / "calc_weno.f90"
    exe = tmp_path / "weno_order"
    subprocess.run(["fypp", str(src), str(gen)], check=True)
    # same -cuda -fast as the solver build (-fast is what makes the FP32 path FTZ)
    subprocess.run([_nvfortran(), "-cuda", "-fast", str(gen), str(drv), "-o", str(exe)],
                   check=True, cwd=tmp_path)
    out = subprocess.run([str(exe)], check=True, capture_output=True, text=True).stdout
    rows = [l.split() for l in out.splitlines() if re.match(r"\s*\d+\s+\d+\s", l)]
    table = {(int(r[0]), int(r[1])): [float(x) for x in r[2:]] for r in rows}
    labels = ["FP64 v-", "FP64 v+", "W32P64 v-", "W32P64 v+"]
    for (case, n), errs in sorted(table.items()):
        print(f"  case {case} N={n:5d} " + "  ".join(f"{l}={e:.3e}" for l, e in zip(labels, errs)))
    assert all(np.isfinite(e) for errs in table.values() for e in errs)

    h = np.log([1.0 / n for n in RECON_FIT_GRIDS])
    for col in (0, 1):   # FP64, both biases, O(1) data
        e = np.log([table[(1, n)][col] for n in RECON_FIT_GRIDS])
        order = np.polyfit(h, e, 1)[0]
        print(f"  {labels[col]}: order {order:.2f} over N={RECON_FIT_GRIDS}")
        assert order >= RECON_MIN_ORDER, f"{labels[col]} order {order:.2f} < {RECON_MIN_ORDER}"


# ---------------------------------------------------------------------------
# EVC: FP64 WENO-Z + LF grid convergence
# ---------------------------------------------------------------------------
def _run_evc(dim, prec, riemann, nx):
    """Build and run one EVC grid; return (L2(rho), z_err)."""
    workdir = os.path.join("./tmp", f"evc{dim}_weno_{prec}_{riemann}_{nx}".lower())
    n = nx + ACCURACY
    kw = dict(
        visc="euler", tvd="none", nx=n, ny=n,
        Lx=evc.Lx, Ly=evc.Ly, gamma=evc.gamma, Pr=evc.Pr, R=evc.R, M0=evc.M0,
        beta=evc.beta, theta=evc.theta, Rc=evc.Rc, p0=evc.p0, T0=evc.T0, u0=evc.u0,
        output_precision=8, **_weno_params(prec, riemann),
    )
    if dim == "3d":
        kw.update(nz=9, Lz=0.01e0)   # see test_evc_3d.py for why nz=9
    case = Case(source=f"{dim.upper()}_solver/EVC", workdir=workdir, **kw)
    case.build()
    case.run(nranks=2)
    data_dir = pathlib.Path(workdir) / "data"
    return evc.l2_rho_error(initial_vtr(data_dir), latest_vtr(data_dir), ACCURACY)


@pytest.mark.integration
@pytest.mark.parametrize("dim", [
    "2d",
    pytest.param("3d", marks=pytest.mark.slow),
])
def test_evc_weno_fp64_lf_convergence(dim):
    runner = lambda scheme, accuracy, nx, params: _run_evc(dim, "FP64", "LF", nx)
    evc.assert_grid_convergence(runner, f"{dim} WENO5-Z FP64 + LF", ACCURACY,
                                EVC_MIN_ORDER, None)


# ---------------------------------------------------------------------------
# OS: oblique shock vs Rankine-Hugoniot
# ---------------------------------------------------------------------------
@pytest.mark.integration
@pytest.mark.parametrize("dim,prec,riemann", [
    ("2d", "FP64", "SLAU"),
    ("2d", "W32P64", "LF"),
    pytest.param("3d", "FP64", "LF", marks=pytest.mark.slow),
])
def test_os_weno(dim, prec, riemann):
    workdir = os.path.join("./tmp", f"os{dim}_weno_{prec}_{riemann}".lower())
    kw = dict(
        visc="euler", tvd="hybrid", nx=osc.nx, ny=osc.ny, Lx=osc.Lx, Ly=osc.Ly,
        gamma=osc.gamma, Pr=osc.Pr, R=osc.R, M0=osc.M0, p_tot=osc.p_tot,
        T_tot=osc.T_tot, beta=osc.beta, dt=osc.dt, np=1, nt=int(osc.endT / osc.dt),
        **_weno_params(prec, riemann),
    )
    if dim == "3d":
        kw.update(nz=9)
    case = Case(source=f"{dim.upper()}_solver/OS", workdir=workdir, **kw)
    case.build()
    case.run(nranks=2)
    vtr_path = latest_vtr(pathlib.Path(workdir) / "data")
    ni, nj, nk, x, y, z = getGrid(vtr_path)
    rho, u, v, _, p = getQ(vtr_path, ni, nj, nk)
    osc.assert_pre_post_shock_matches_analytical(rho, u, v, p, ni, nj, k_ref=nk // 2)
