"""
TGV kinetic-energy / dissipation report: post-process a compressible
Taylor-Green-vortex run's VTK snapshots into Ek(t) (volume-averaged kinetic
energy) and the enstrophy (epsE) / dilatational (epsD) viscous dissipation
decomposition. The standard compressible-TGV energy-budget self-consistency
check is that -dEk/dt (numerically differentiated from Ek(t)) should match
epsE(t)+epsD(t), both computed independently from the same VTK snapshots --
this validates that a viscous kernel's stress tensor is physically consistent
with the dissipation it should be producing.

Usage (standalone):
    python3 tgv_ke_eps.py [data_dir]   (default: ./Re1600_128_Hybrid)

Also importable: compute(data_dir, Re, M0, L0, T, S) returns
{t, Ek, epsE, epsD} for use from a test -- physical constants are passed in
rather than hardcoded, since different runs (e.g. a subsonic KEEP+NS check
vs. a supersonic SLAU+NS check) use different M0.
"""
import sys
import os
import pathlib

import numpy as np
from numba import njit

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests" / "utils"))
from vtk_reader import getGrid, getVector, getScalar, extract_number  # noqa: E402

RGAS_DEFAULT = 287.03e0


@njit(cache=True, nogil=True, fastmath=True)
def Sutherland(T, S):
  return 1.716e-5 * (273.2e0 + S) / (T + S) * (T / 273.2e0)**1.5e0


@njit(cache=True, nogil=True, fastmath=True)
def calc_ke(rho0, V0, rho, u, v, w, offset):
  Nz, Ny, Nx = rho.shape
  Ek = 0.e0
  for k in range(offset,Nz-offset):
    for j in range(offset,Ny-offset):
      for i in range(offset,Nx-offset):
        Ek += 0.5e0 * rho[k,j,i] * (u[k,j,i]**2 + v[k,j,i]**2 + w[k,j,i]**2)
  Ek /= (rho0 * V0**2 * (Nx - 2 * offset) * (Ny - 2 * offset) * (Nz - 2 * offset))
  return Ek


@njit(cache=True, nogil=True, fastmath=True)
def central2nd(i, j, k, x, y, z, u, v, w):
  ux = (-u[k,j,i-1] + u[k,j,i+1]) / (-x[i-1] + x[i+1])
  uy = (-u[k,j-1,i] + u[k,j+1,i]) / (-y[j-1] + y[j+1])
  uz = (-u[k-1,j,i] + u[k+1,j,i]) / (-z[k-1] + z[k+1])
  vx = (-v[k,j,i-1] + v[k,j,i+1]) / (-x[i-1] + x[i+1])
  vy = (-v[k,j-1,i] + v[k,j+1,i]) / (-y[j-1] + y[j+1])
  vz = (-v[k-1,j,i] + v[k+1,j,i]) / (-z[k-1] + z[k+1])
  wx = (-w[k,j,i-1] + w[k,j,i+1]) / (-x[i-1] + x[i+1])
  wy = (-w[k,j-1,i] + w[k,j+1,i]) / (-y[j-1] + y[j+1])
  wz = (-w[k-1,j,i] + w[k+1,j,i]) / (-z[k-1] + z[k+1])
  return ux, uy, uz, vx, vy, vz, wx, wy, wz


# 4th is better
@njit(cache=True, nogil=True, fastmath=True)
def central4th(i, j, k, x, y, z, u, v, w):
  ux = (u[k,j,i-2] + 8.e0 * (-u[k,j,i-1] + u[k,j,i+1]) - u[k,j,i+2]) / (6.e0 * (-x[i-1] + x[i+1]))
  uy = (u[k,j-2,i] + 8.e0 * (-u[k,j-1,i] + u[k,j+1,i]) - u[k,j+2,i]) / (6.e0 * (-y[j-1] + y[j+1]))
  uz = (u[k-2,j,i] + 8.e0 * (-u[k-1,j,i] + u[k+1,j,i]) - u[k+2,j,i]) / (6.e0 * (-z[k-1] + z[k+1]))
  vx = (v[k,j,i-2] + 8.e0 * (-v[k,j,i-1] + v[k,j,i+1]) - v[k,j,i+2]) / (6.e0 * (-x[i-1] + x[i+1]))
  vy = (v[k,j-2,i] + 8.e0 * (-v[k,j-1,i] + v[k,j+1,i]) - v[k,j+2,i]) / (6.e0 * (-y[j-1] + y[j+1]))
  vz = (v[k-2,j,i] + 8.e0 * (-v[k-1,j,i] + v[k+1,j,i]) - v[k+2,j,i]) / (6.e0 * (-z[k-1] + z[k+1]))
  wx = (w[k,j,i-2] + 8.e0 * (-w[k,j,i-1] + w[k,j,i+1]) - w[k,j,i+2]) / (6.e0 * (-x[i-1] + x[i+1]))
  wy = (w[k,j-2,i] + 8.e0 * (-w[k,j-1,i] + w[k,j+1,i]) - w[k,j+2,i]) / (6.e0 * (-y[j-1] + y[j+1]))
  wz = (w[k-2,j,i] + 8.e0 * (-w[k-1,j,i] + w[k+1,j,i]) - w[k+2,j,i]) / (6.e0 * (-z[k-1] + z[k+1]))
  return ux, uy, uz, vx, vy, vz, wx, wy, wz


@njit(cache=True, nogil=True, fastmath=True)
def calc_total_dissipation(rho0, V0, L0, Rgas, S, x, y, z, rho, u, v, w, p, offset):
  """epsE (enstrophy), epsD (dilatational) viscous dissipation, and pdiv
  (pressure-dilatation, <p*div(u)>) -- all volume-averaged, same
  (rho0*V0^3/L0) normalisation. For a periodic domain, the compressible
  kinetic-energy budget is exactly -d(Ek)/dt = (epsE+epsD) - pdiv: pdiv is
  a *reversible* exchange with internal energy (not part of viscous
  dissipation), significant whenever the flow has non-negligible dilatation
  (higher Mach number, or later in a run as compressible structure
  develops) -- see ouxsbli/tests/test_nstgv_ke_eps.py for how this is used.
  """
  Nz, Ny, Nx = rho.shape
  epsE = 0.e0
  epsD = 0.e0
  pdiv = 0.e0
  for k in range(offset,Nz-offset):
    for j in range(offset,Ny-offset):
      for i in range(offset,Nx-offset):
        T  = p[k,j,i] / (rho[k,j,i] * Rgas)
        mu = Sutherland(T, S)
        ux, uy, uz, vx, vy, vz, wx, wy, wz = central4th(i, j, k, x, y, z, u, v, w)
        omega = (wy - vz)**2 + (uz - wx)**2 + (vx - uy)**2
        div   = ux + vy + wz
        epsE += mu * omega
        epsD += 4.e0 / 3.e0 * mu * div**2
        pdiv += p[k,j,i] * div
  N = (Nx - 2 * offset) * (Ny - 2 * offset) * (Nz - 2 * offset)
  epsE /= (rho0 * V0**3 / L0 * N)
  epsD /= (rho0 * V0**3 / L0 * N)
  pdiv /= (rho0 * V0**3 / L0 * N)
  return epsE, epsD, pdiv


def compute(data_dir, Re, M0, L0, T, S, Rgas=RGAS_DEFAULT, offset=3):
  """Load one TGV run's data/ dir and return {t, Ek, epsE, epsD, pdiv}.

  Time values come from data/kinetic_energy.d's own time column (real
  physical time, `nt*step*dt`), matched by ordinal index to the sorted
  Q*.vtr snapshots -- same convention ouxsbli/analysis/dhit_decay_report.py
  uses, since both are written at the same np-block cadence. Ek/epsE/epsD
  are computed independently here from each snapshot's own rho/u/v/w/p
  fields (not read from the solver's own ke/ke0 log columns, which use a
  different normalisation), so the two sides of the energy-budget check
  stay self-consistent.
  """
  data_dir = pathlib.Path(data_dir)
  mu0  = Sutherland(T, S)
  V0   = M0 * np.sqrt(1.4e0 * Rgas * T)
  rho0 = mu0 * Re / (V0 * L0)

  ke_file = data_dir / "kinetic_energy.d"
  if not ke_file.exists():
    raise FileNotFoundError(f"{ke_file} not found -- did the run complete?")
  ke_data = np.loadtxt(ke_file)
  if ke_data.ndim == 1:
    ke_data = ke_data.reshape(1, -1)
  t_ke = ke_data[:, 0]

  Q_files = [f for f in os.listdir(data_dir) if f.endswith(".vtr")]
  Q_files.sort(key=extract_number)
  if not Q_files:
    raise FileNotFoundError(f"No Q*.vtr snapshots found in {data_dir}")

  first_path = os.path.join(data_dir, Q_files[0])
  Nx, Ny, Nz, _, _, _ = getGrid(first_path)
  Nt = len(Q_files)

  t    = np.full(Nt, np.nan)
  Ek   = np.zeros(Nt)
  epsE = np.zeros(Nt)
  epsD = np.zeros(Nt)
  pdiv = np.zeros(Nt)
  for itr, fname in enumerate(Q_files):
    file_path = os.path.join(data_dir, fname)
    _, _, _, x, y, z = getGrid(file_path)
    rho     = getScalar(file_path, Nx, Ny, Nz, 'rho')
    u, v, w = getVector(file_path, Nx, Ny, Nz, 'velocity')
    p       = getScalar(file_path, Nx, Ny, Nz, 'p')
    t[itr]        = t_ke[itr] if itr < len(t_ke) else np.nan
    Ek[itr]       = calc_ke(rho0, V0, rho, u, v, w, offset)
    epsE[itr], epsD[itr], pdiv[itr] = calc_total_dissipation(rho0, V0, L0, Rgas, S, x, y, z, rho, u, v, w, p, offset)

  return {"t": t, "Ek": Ek, "epsE": epsE, "epsD": epsD, "pdiv": pdiv}


def main(data_dir):
  from tqdm import tqdm  # noqa: F401  (kept for parity with the original script; unused directly)
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt

  # Defaults match the original Re1600_128_Hybrid capture this script was
  # first written for -- 3D_solver/NSTGV/mod_globals.f90's own production
  # constants.
  Re   = 1600.e0
  M0   = 1.25e0
  L0   = 1.524e-3
  T    = 530.e0 * 5.e0 / 9.e0
  S    = 111.e0

  r = compute(data_dir, Re=Re, M0=M0, L0=L0, T=T, S=S)
  t, Ek, epsE, epsD, pdiv = r["t"], r["Ek"], r["epsE"], r["epsD"], r["pdiv"]

  plt.plot(t, Ek)
  plt.xlabel("t"); plt.ylabel("Ek")
  plt.savefig(os.path.join(data_dir, "Ek.png"))
  plt.close()

  save_path = os.path.join(data_dir, "Ek.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# t       Ek", file=f)
    for i in range(len(t)):
      print(f'{t[i]:.3e}', f'{Ek[i]:.3e}', file=f)

  plt.plot(t, epsE, label="epsE")
  plt.plot(t, epsD, label="epsD")
  plt.plot(t, pdiv, label="pdiv")
  plt.xlabel("t"); plt.ylabel("eps"); plt.legend()
  plt.savefig(os.path.join(data_dir, "eps.png"))
  plt.close()

  save_path = os.path.join(data_dir, "eps.d")
  with open(save_path, "w", encoding="UTF-8") as f:
    print("# t       epsE      epsD      pdiv", file=f)
    for i in range(len(t)):
      print(f'{t[i]:.3e}', f'{epsE[i]:.3e}', f'{epsD[i]:.3e}', f'{pdiv[i]:.3e}', file=f)


if __name__ == "__main__":
  data_dir = sys.argv[1] if len(sys.argv) > 1 else "./Re1600_128_Hybrid"
  main(data_dir)
