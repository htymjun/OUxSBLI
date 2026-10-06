#!/usr/bin/env python3
"""Mixing statistics from the passive-scalar output of SL_spatial_inlet.

Reads data/<rank>/xiNNNNN.vtr (or data/xiNNNNN.vtr for a single slab), joins the
x slabs (npx in mod_globals.f90), and averages over z and over the selected snapshots.

    python3 scalar_mixing.py --first 200 --last 1000 --out mixing

writes mixing.npz (x, y and the (ny, nx) fields below) and mixing.csv (one row per x):
    mean        <xi>
    rms         sqrt(<xi^2> - <xi>^2)
    mixedness   <4 xi (1 - xi)>   (1 = fully mixed at xi = 0.5, 0 = unmixed)
    y01, y09    y where <xi> crosses 0.1 and 0.9;  delta_xi = y09 - y01
    delta_mix   integral of mixedness over y  (thickness of mixed fluid)
"""
import argparse
import glob
import os
import re

import numpy as np


def read_xi(path):
    """Return (x, y, z, xi[nz, ny, nx]) from a xiNNNNN.vtr written by set.f90 (ghost planes included)."""
    with open(path, "rb") as f:
        head = f.read(8192)
        cut = head.find(b'<AppendedData encoding="raw">')
        xml = head[:cut].decode("ascii", "replace")
        ext = [int(v) for v in re.search(r'WholeExtent="([^"]+)"', xml).group(1).split()]
        n = (ext[1] + 1, ext[3] + 1, ext[5] + 1)
        offs = [int(v) for v in re.findall(r'offset="\s*(\d+)"', xml)]
        base = head.find(b"_", cut) + 1
        out = []
        for off, count in zip(offs, (n[0], n[1], n[2], n[0] * n[1] * n[2])):
            f.seek(base + off + 4)                    # +4 skips the block's own byte count
            out.append(np.fromfile(f, "<f4", count))
    x, y, z, xi = out
    return x.astype(float), y.astype(float), z.astype(float), xi.reshape((n[2], n[1], n[0]))


def slab_dirs(data):
    """Slab directories ordered in x: data/1, data/3, ... or data itself."""
    sub = sorted((d for d in os.listdir(data) if d.isdigit()), key=int)
    sub = [os.path.join(data, d) for d in sub if glob.glob(os.path.join(data, d, "xi*.vtr"))]
    return sub or [data]


def crossing(y, prof, level):
    """First y where a profile crosses `level` (linear interpolation), nan if never."""
    s = prof - level
    idx = np.where(s[:-1] * s[1:] <= 0)[0]
    idx = idx[s[idx] != s[idx + 1]]
    if idx.size == 0:
        return np.nan
    j = idx[0]
    return y[j] + (y[j + 1] - y[j]) * (level - prof[j]) / (prof[j + 1] - prof[j])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data")
    ap.add_argument("--first", type=int, default=0, help="first snapshot number to average")
    ap.add_argument("--last", type=int, default=10**9, help="last snapshot number to average")
    ap.add_argument("--out", default="mixing")
    a = ap.parse_args()

    dirs = slab_dirs(a.data)
    numbers = sorted(int(re.search(r"xi(\d+)\.vtr$", p).group(1)) for p in glob.glob(os.path.join(dirs[0], "xi*.vtr")))
    numbers = [n for n in numbers if a.first <= n <= a.last]
    if not numbers:
        raise SystemExit("no xi snapshots in the requested range")

    s1 = s2 = sm = None
    count = 0
    for n in numbers:
        parts, xs = [], []
        for r, d in enumerate(dirs):
            xl, y, _, xi = read_xi(os.path.join(d, f"xi{n:05d}.vtr"))
            # x slabs overlap by 3 ghost planes on each joined side; keep each plane once
            lo = 0 if r == 0 else 3
            hi = xi.shape[2] if r == len(dirs) - 1 else xi.shape[2] - 3
            parts.append(xi[3:-3, :, lo:hi])          # also drop the 3+3 periodic z ghost planes
            xs.append(xl[lo:hi])
        xi = np.concatenate(parts, axis=2).astype(np.float64)
        x = np.concatenate(xs)
        if s1 is None:
            s1 = np.zeros(xi.shape[1:]); s2 = np.zeros_like(s1); sm = np.zeros_like(s1)
        s1 += xi.sum(axis=0)
        s2 += (xi * xi).sum(axis=0)
        sm += (4.0 * xi * (1.0 - xi)).sum(axis=0)
        count += xi.shape[0]
    mean = s1 / count
    rms = np.sqrt(np.maximum(s2 / count - mean**2, 0.0))
    mixedness = sm / count

    ny, nx = mean.shape

    y01 = np.array([crossing(y, mean[:, i], 0.1) for i in range(nx)])
    y09 = np.array([crossing(y, mean[:, i], 0.9) for i in range(nx)])
    delta_mix = np.trapz(mixedness, y, axis=0)
    rms_max = rms.max(axis=0)

    np.savez(a.out + ".npz", x=x, y=y, mean=mean, rms=rms, mixedness=mixedness,
             y01=y01, y09=y09, delta_xi=y09 - y01, delta_mix=delta_mix, rms_max=rms_max,
             snapshots=np.array(numbers), slabs=len(dirs))
    np.savetxt(a.out + ".csv", np.column_stack([x, y01, y09, y09 - y01, delta_mix, rms_max]),
               delimiter=",", header="x,y01,y09,delta_xi,delta_mix,rms_max", comments="")
    print(f"{len(numbers)} snapshots, {len(dirs)} x slabs joined ({nx} x points), {count} z planes averaged -> {a.out}.npz, {a.out}.csv")


if __name__ == "__main__":
    main()
