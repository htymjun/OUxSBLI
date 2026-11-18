"""
Precision modes of 2D_solver against FP64 (PREC = 'FP64' | 'DF' | 'FP32').

Each configuration runs 40 steps in FP64, DF and FP32 from the same build options;
the metric is max |q - q_FP64| / max |q_FP64| of the last snapshot
(OUTPUT_PRECISION = 8) for rho, p and velocity.

  * DF must agree with FP64 to DF rounding (DF_TOL). HR-SLAU2 configurations are
    left out: its wiggle detector is a 0/1 switch on the sign of products of
    pressure differences, which are exactly 0 in FP64 in uniform regions and
    +-1e-30 under DF (docs/precision_modes.md), so the plain-SLAU variant is run.
  * FP32 must agree with FP64 to FP32 rounding (FP32_TOL), and its executable
    must contain no FP64 instruction in any kernel (SASS census).
"""
import glob
import os
import pathlib
import re
import subprocess

import numpy as np
import pytest
from ouxsbli import Case
from ouxsbli.analysis.vtr_raw import read_vtr

pytestmark = pytest.mark.slow

NT = 40
DF_TOL = 1e-11
FP32_TOL = 1e-4
FP64_OP = re.compile(r"\b(DADD|DMUL|DFMA|DSETP|DMNMX)\b|\b[FI]2[FI]\.[^ ;]*64")

CONFIGS = {
    "ST":         ("ST",   {}),                                         # NS KEEP6, walls, calc_div
    "OS_slau1":   ("OS",   dict(slau="slau")),                          # Euler SLAU, shock, TVD hybrid
    "BL_slau1":   ("BL",   dict(slau="slau")),                          # NS SLAU, Riemann far field
    "EVC_hybrid": ("EVC",  dict(scheme="hybrid", slau="slau")),         # Hybrid KEEP/SLAU, periodic
    "EVC_wenolf": ("EVC",  dict(scheme="slau", recon="weno", riemann="lf")),
}


def _run(name, prec):
    source, kw = CONFIGS[name]
    workdir = pathlib.Path("./tmp/prec2d") / f"{name}_{prec}"
    case = Case(source=f"2D_solver/{source}", workdir=str(workdir), nt=NT, np=1,
                output_precision=8, prec=prec, **kw)
    case.build()
    case.run(nranks=2)
    vtr = sorted(glob.glob(str(workdir / "data" / "*.vtr")))[-1]
    return case, workdir, read_vtr(vtr)[3]


def _rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.max(np.abs(a - b)) / np.max(np.abs(a))


@pytest.mark.integration
@pytest.mark.parametrize("name", list(CONFIGS))
def test_prec_2d_vs_fp64(name):
    _, _, ref = _run(name, "fp64")
    _, _, df = _run(name, "df")
    case, workdir, fp32 = _run(name, "fp32")
    for k in ("rho", "p", "velocity"):
        assert _rel(ref[k], df[k]) < DF_TOL, f"{name} DF {k}: {_rel(ref[k], df[k]):.2e}"
        assert np.all(np.isfinite(fp32[k])), f"{name} FP32 {k}: non-finite"
        assert _rel(ref[k], fp32[k]) < FP32_TOL, f"{name} FP32 {k}: {_rel(ref[k], fp32[k]):.2e}"
    # FP32: no FP64 instruction in any kernel
    env = dict(os.environ, PATH=case._build_env()["PATH"])
    sass = subprocess.run(["cuobjdump", "-sass", str(workdir / "build" / "a.out")],
                          capture_output=True, text=True, env=env, check=True).stdout
    bad, kernel = {}, None
    for line in sass.splitlines():
        m = re.search(r"Function : (\S+)", line)
        if m:
            kernel = m.group(1)
        elif kernel and FP64_OP.search(line):
            bad[kernel] = bad.get(kernel, 0) + 1
    assert kernel is not None, "no kernels found by cuobjdump"
    assert not bad, f"{name} FP32: FP64 instructions in {bad}"
