"""
Navier-Stokes energy-budget self-consistency check on 3D_solver/NSTGV's
compressible Taylor-Green vortex.

Unlike every other integration test in this suite (which forces
visc="euler"), this one actually exercises viscous physics: the full
compressible kinetic-energy budget on a periodic domain is exactly

    -dEk/dtau = (epsE + epsD) - pdiv          (tau = t/tc, tc = L0/V0)

where epsE/epsD are the enstrophy/dilatational viscous-dissipation terms and
pdiv = <p*div(u)> is the pressure-dilatation term -- a *reversible* exchange
with internal energy, not part of viscous dissipation, but not negligible
for a compressible (let alone supersonic) TGV. All three terms are computed
independently from the same VTK snapshots via
ouxsbli/analysis/tgv_ke_eps.py's compute(); dEk/dtau is a finite difference
of Ek itself. This also regression-guards the fused convective+viscous
kernels directly: KEEP+NS's calc_keep_visc_kernel and its generalization to
SLAU+NS's calc_slau_visc_kernel both dispatch on NSTGV's own case (SCHEME in
('KEEP','SLAU'), VISC='NS', ORDER=VISC_ORDER=6, periodic).

(An earlier version of this test checked -dEk/dtau against epsE+epsD alone,
omitting pdiv. That holds only near t=0 for a mildly compressible flow --
measured directly on this data, dropping pdiv grows the relative error from
<1% to >100% within the first handful of intervals for the supersonic SLAU
case. Adding pdiv closes the budget to within ~2.5% across the *entire* run
for both schemes, including through SLAU's case where dEk/dtau changes sign
partway through -- so this is the correct budget, not a curve-fit tolerance.)

KEEP requires shock-free flow (per CLAUDE.md's own scheme guidance), so it
runs at a subsonic M0; SLAU is designed for compressible/shocked flow and
runs at NSTGV's actual production M0=1.25 (a supersonic TGV) -- these are
two different flows, not the same flow through two schemes.
"""
import pathlib
import sys

import numpy as np
import pytest

from ouxsbli import Case

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "analysis"))
from tgv_ke_eps import compute  # noqa: E402

# Physical constants shared with 3D_solver/NSTGV/mod_globals.f90's own
# defaults (Re, L0, T, S) -- only M0 varies per scheme below.
RE   = 1600.e0
L0   = 1.524e-3
T    = 530.e0 * 5.e0 / 9.e0
S    = 111.e0

GRID = 129   # nx=ny=nz=129 (128 cells)
NP   = 20    # number of output blocks -> 21 VTK snapshots (incl. t=0)
NT   = 8     # RK steps between each output block (160 total steps)

# Measured max relative error across the whole run is ~2.5% (SLAU) / ~0.4%
# (KEEP) -- this leaves comfortable margin while still being a tight check.
BUDGET_REL_TOL = 0.05e0


def _run_and_check(scheme, M0, workdir):
    case = Case(
        source   = "3D_solver/NSTGV",
        workdir  = workdir,
        scheme   = scheme,   # explicit -- don't rely on config.fypp's default
        visc     = "ns",
        M0       = M0,
        nx       = GRID,
        ny       = GRID,
        nz       = GRID,
        np       = NP,
        nt       = NT,
    )
    case.build()
    case.run(nranks=2)

    data_dir = pathlib.Path(workdir) / "data"
    r = compute(data_dir, Re=RE, M0=M0, L0=L0, T=T, S=S)
    t, Ek, epsE, epsD, pdiv = r["t"], r["Ek"], r["epsE"], r["epsD"], r["pdiv"]
    assert np.all(np.isfinite(Ek)) and np.all(np.isfinite(epsE)) and np.all(np.isfinite(epsD)), (
        "NaN/Inf in Ek or eps -- run likely diverged"
    )
    assert np.all(epsE >= 0.e0) and np.all(epsD >= 0.e0), (
        "epsE/epsD must be non-negative (they are sums of squares) -- a "
        "negative value means the field read from VTK is inconsistent"
    )

    # energy-budget closure: -dEk/dtau == (epsE+epsD) - pdiv, over the whole run.
    V0  = M0 * np.sqrt(1.4e0 * 287.03e0 * T)
    tau = t * V0 / L0
    dEk_dtau  = -np.diff(Ek) / np.diff(tau)
    eps_mid   = 0.5e0 * (epsE[:-1] + epsE[1:]) + 0.5e0 * (epsD[:-1] + epsD[1:])
    pdiv_mid  = 0.5e0 * (pdiv[:-1] + pdiv[1:])
    budget    = eps_mid - pdiv_mid
    rel_err   = np.abs(dEk_dtau - budget) / np.abs(budget)
    assert np.all(rel_err < BUDGET_REL_TOL), (
        f"-dEk/dtau vs (epsE+epsD)-pdiv relative error {rel_err.max():.3f} "
        f">= {BUDGET_REL_TOL} ({scheme}, M0={M0})"
    )


@pytest.mark.integration
def test_nstgv_keep_ns_subsonic_ke_eps():
    """KEEP+NS fused kernel (calc_keep_visc_kernel), subsonic TGV."""
    _run_and_check("keep", 0.4e0, "./tmp/nstgv_keep_ns")


@pytest.mark.integration
def test_nstgv_slau_ns_supersonic_ke_eps():
    """SLAU+NS fused kernel (calc_slau_visc_kernel), supersonic TGV (production M0)."""
    _run_and_check("slau", 1.25e0, "./tmp/nstgv_slau_ns")
