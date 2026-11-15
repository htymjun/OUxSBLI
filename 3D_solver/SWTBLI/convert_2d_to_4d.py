#!/usr/bin/env python3
"""Tile the 2 delta-span restart files of SWTBLI into a 4 delta-span initial field.

    python3 convert_2d_to_4d.py SRC_RECAL DST_RECAL [--amp 0.01] [--seed 1]

Reads SRC_RECAL/Q00001.dat (TBL generator, nx1 x ny1 x nz x 5) and Q00002.dat
(SBLI, nx2 x ny2 x nz x 5), raw stream files of the conservative variables as
written by src/sbli.f90 (Fortran order, Q(nx,ny,nz,5)). The nz-6 interior planes
are written twice (nz' = 2(nz-6)+6) and the 3+3 periodic ghost planes refilled.

A Gaussian-filtered random velocity perturbation of amplitude ``amp*u0``
(u and w, half that on v) is added inside the boundary layer (y < 1.2 delta,
vanishing at the wall) at constant density and pressure. It is required: the
recycling uses a half-span spanwise shift, which maps an exactly Lz/2-periodic
field onto itself, so without it the two halves would stay identical.
Qm.dat depends only on ny and is copied unchanged.
"""
import os
import shutil
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from ouxsbli.analysis.swtbli_stats import case_params, make_y  # noqa: E402


def smooth(a, axis, n=3, periodic=False):
    w = np.exp(-0.5 * (np.arange(-n, n + 1) / (n / 3.0)) ** 2)
    w /= w.sum()
    out = np.zeros_like(a)
    for s, c in zip(range(-n, n + 1), w):
        if periodic:
            out += c * np.roll(a, s, axis=axis)
        else:
            idx = np.clip(np.arange(a.shape[axis]) + s, 0, a.shape[axis] - 1)
            out += c * np.take(a, idx, axis=axis)
    return out


def tile(path_in, path_out, nx, ny, nz, y, blt, u0, amp, rng):
    Q = np.fromfile(path_in, dtype="<f8")
    assert Q.size == nx * ny * nz * 5, f"{path_in}: size {Q.size} != {nx}*{ny}*{nz}*5"
    Q = Q.reshape((nx, ny, nz, 5), order="F")
    ni = nz - 6
    nz2 = 2 * ni + 6
    inner = Q[:, :, 3:3 + ni, :]
    Q2 = np.empty((nx, ny, nz2, 5))
    Q2[:, :, 3:3 + ni, :] = inner
    Q2[:, :, 3 + ni:3 + 2 * ni, :] = inner
    # perturbation on the interior, periodic in z
    fd = np.where(y < 1.2 * blt, np.sin(np.pi * np.minimum(y / (1.2 * blt), 1.0)), 0.0)
    rho = Q2[:, :, 3:3 + 2 * ni, 0]
    for comp, a in ((1, 1.0), (2, 0.5), (3, 1.0)):
        r = rng.uniform(-1, 1, size=(nx, ny, 2 * ni))
        r = smooth(smooth(r, 0), 2, periodic=True)
        r *= amp * u0 * a / np.sqrt((r**2).mean()) * fd[None, :, None]
        mom = Q2[:, :, 3:3 + 2 * ni, comp]
        e = Q2[:, :, 3:3 + 2 * ni, 4]
        vel = mom / rho
        new = vel + r
        e += 0.5 * rho * (new**2 - vel**2)       # keep p and rho
        Q2[:, :, 3:3 + 2 * ni, comp] = rho * new
    Q2[:, 0, :, 1:4] = 0.0                       # wall: no slip
    Q2[:, :, 0:3, :] = Q2[:, :, 2 * ni:2 * ni + 3, :]
    Q2[:, :, 3 + 2 * ni:, :] = Q2[:, :, 3:6, :]
    Q2.reshape(-1, order="F").astype("<f8").tofile(path_out)
    return nz2


def main():
    src, dst = sys.argv[1], sys.argv[2]
    amp = float(sys.argv[sys.argv.index("--amp") + 1]) if "--amp" in sys.argv else 0.01
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 1
    P = case_params(os.path.dirname(os.path.abspath(src)))
    nx1, ny1, nz, nx2, ny2 = (P[k] for k in ("nx1", "ny1", "nz1", "nx2", "ny2"))
    blt = P["blt"]
    y_all = make_y(P, ny2)
    os.makedirs(dst, exist_ok=True)
    rng = np.random.default_rng(seed)
    n1 = tile(os.path.join(src, "Q00001.dat"), os.path.join(dst, "Q00001.dat"), nx1, ny1, nz, y_all[:ny1], blt, P["u0"], amp, rng)
    n2 = tile(os.path.join(src, "Q00002.dat"), os.path.join(dst, "Q00002.dat"), nx2, ny2, nz, y_all, blt, P["u0"], amp, rng)
    if os.path.exists(os.path.join(src, "Qm.dat")):
        shutil.copy(os.path.join(src, "Qm.dat"), os.path.join(dst, "Qm.dat"))
    else:
        print("WARNING: no Qm.dat in", src, "(a RESCALE restart needs it)")
    print(f"tiled nz {nz} -> {n1} (TBL) / {n2} (SBLI), perturbation {amp} u0, written to {dst}")


if __name__ == "__main__":
    main()
