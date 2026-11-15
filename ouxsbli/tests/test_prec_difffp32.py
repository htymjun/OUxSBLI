"""
diffFP32's pre-registered accuracy rules (~/diffFP32/docs/research-plan.md)
applied to the OUxSBLI precision modes of KEEP (Euler, all-periodic):

  FP64        reference
  DF          PREC='DF'              whole state in FmaX double-float
  MIXED       KEEP_PREC='MIXED'      difference-first FP32, face form
  MIXED_CELL  KEEP_PREC='MIXED_CELL' difference-first FP32, cell form (diffFP32)

Rules (a mode keeps FP64 accuracy only if all hold):
  * convergence (3D EVC, KEEP4/KEEP6, n = 64/128/256): L2(rho) within 1.5x of
    FP64 at every n and the observed order within 0.1 of FP64's;
    The EVC vortex is already a small-amplitude case for rho (T' / T0 ~ 2e-7),
    the regime of diffFP32's A = 1e-6 entropy wave where value-first FP32
    fails by orders of magnitude;
  * long-time TGV (ETGV 66^3, RK4, 20000 steps, KEEP2/KEEP6): the time at
    which ke/ke0 or the entropy leaves test_etgv's bounds (or NaN appears)
    is within 10% of FP64's, or neither leaves them.
diffFP32's shock-tube rule concerns WENO/TENO; WENO_PREC='MIXED' is covered
by test_weno.py's reconstruction driver. Every mode is selected through
OUXSBLI_CONFIG, so the cases are the unchanged ones of test_evc_3d/test_etgv.
"""
import os
import pathlib

import numpy as np
import pytest
from ouxsbli import Case
from .test_evc_3d import _run_evc_3d, params as EVC_PARAMS
from .utils.evc_common import GRIDS, convergence_order

pytestmark = pytest.mark.slow

MODES = {"FP64": "", "DF": "prec=df", "MIXED": "keep_prec=mixed", "MIXED_CELL": "keep_prec=mixed_cell"}
L2_FACTOR, ORDER_TOL, COLLAPSE_TOL = 1.5, 0.1, 0.1


def _evc_errors(monkeypatch, mode, accuracy, params):
    monkeypatch.setenv("OUXSBLI_CONFIG", MODES[mode])
    return [_run_evc_3d("keep", accuracy, nx, params)[0] for nx in GRIDS]


def _check_against_fp64(errors, check_order):
    ref = errors["FP64"]
    for mode, e in errors.items():
        print(f"  {mode:10s} " + "  ".join(f"{x:.3e}" for x in e)
              + f"  order {convergence_order(e):.2f}  ratio {max(np.array(e) / ref):.2f}")
    for mode, e in errors.items():
        assert max(np.array(e) / ref) <= L2_FACTOR, f"{mode}: L2 above {L2_FACTOR}x FP64"
        if check_order:
            assert abs(convergence_order(e) - convergence_order(ref)) <= ORDER_TOL, f"{mode}: order"


@pytest.mark.integration
@pytest.mark.parametrize("accuracy", [4, 6])
def test_evc_convergence_vs_fp64(monkeypatch, accuracy):
    errors = {m: _evc_errors(monkeypatch, m, accuracy, EVC_PARAMS) for m in MODES}
    _check_against_fp64(errors, check_order=True)


def _etgv_departure(monkeypatch, mode, accuracy):
    """First time ke/ke0 or the entropy leaves test_etgv's bounds (None if never)."""
    monkeypatch.setenv("OUXSBLI_CONFIG", MODES[mode])
    workdir = pathlib.Path("./tmp") / f"etgv_prec_{mode.lower()}_{accuracy}"
    case = Case(source="3D_solver/ETGV", workdir=str(workdir), visc="euler", scheme="keep",
                accuracy=accuracy, rk="rk4", nx=66, ny=66, nz=66, nt=500, np=40, dt=0.01e0)
    case.build()
    case.run(nranks=2)
    ke = np.loadtxt(workdir / "data" / "kinetic_energy.d")
    ent = np.loadtxt(workdir / "data" / "entropy.d")
    bad = ~np.isfinite(ke[:, -1]) | (ke[:, -1] > 1.01) | ~np.isfinite(ent[:, -1]) | (np.abs(ent[:, -1]) > 3.e-3)
    return (ke[np.argmax(bad), 0] if bad.any() else None), ke[:, -1]


@pytest.mark.integration
@pytest.mark.parametrize("accuracy", [2, 6])
def test_etgv_long_time_vs_fp64(monkeypatch, accuracy):
    res = {m: _etgv_departure(monkeypatch, m, accuracy) for m in MODES}
    t_ref, ke_ref = res["FP64"]
    for mode, (t, ke) in res.items():
        dev = np.nanmax(np.abs(ke - ke_ref)) if len(ke) == len(ke_ref) else np.inf
        print(f"  {mode:10s} departure t = {t}  max|ke/ke0 - FP64| = {dev:.3e}")
    for mode, (t, _) in res.items():
        if t_ref is None:
            assert t is None, f"{mode} leaves the bounds at t = {t}, FP64 never does"
        else:
            assert t is not None and abs(t - t_ref) <= COLLAPSE_TOL * t_ref, f"{mode}: t = {t}, FP64 {t_ref}"
