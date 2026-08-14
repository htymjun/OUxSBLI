#!/usr/bin/env python3
"""Compare 1D_solver/ST Q.dat against the exact Sod solution, and against a
reference Q.dat when one is given.

Q.dat columns are already normalised: x/Lx, rho/rho0, u/a, p/p0 with
a = sqrt(R*Tlr) = sqrt(p0/rho0), which is exactly the non-dimensionalisation
sod_exact.solve() assumes (rho_l=1, p_l=1 -> velocity scale sqrt(p_l/rho_l)).

Similarity time: t_norm = nt*dt * a/Lx = nt * CFL / (nx-1).
"""
import sys
import pathlib
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[0]))
REPO = pathlib.Path("/home/jhatayama/ouxsbli/mixed/OUxSBLI")
sys.path.insert(0, str(REPO))
from ouxsbli.tests.utils.sod_exact import solve as sod_solve

NX, NT, CFL = 4096, 2500, 0.1


def load(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def main():
    q = sys.argv[1]
    x, rho, u, p = load(q)
    t = NT * CFL / (NX - 1)
    rho_e, u_e, p_e = sod_solve(x, t)

    print(f"file          : {q}")
    print(f"similarity t  : {t:.6f}   (nt={NT}, CFL={CFL}, nx={NX})")
    print(f"finite        : {np.all(np.isfinite(np.c_[rho, u, p]))}")
    print(f"rho range     : [{rho.min():.6f}, {rho.max():.6f}]")
    for name, num, ex in (("rho", rho, rho_e), ("u", u, u_e), ("p", p, p_e)):
        l1 = np.abs(num - ex).mean() / np.abs(ex).mean()
        linf = np.abs(num - ex).max()
        print(f"  {name:3s} vs exact:  L1(rel) = {l1:.4%}   Linf(abs) = {linf:.4e}")

    if len(sys.argv) > 2:
        xr, rr, ur, pr = load(sys.argv[2])
        print(f"\nvs reference  : {sys.argv[2]}")
        for name, num, ref in (("rho", rho, rr), ("u", u, ur), ("p", p, pr)):
            d = np.abs(num - ref)
            scale = np.abs(ref).max()
            print(f"  {name:3s}: max|diff| = {d.max():.6e}   "
                  f"rel to max|ref| = {d.max()/scale:.3e}   "
                  f"L1 = {d.mean():.6e}")


if __name__ == "__main__":
    main()
