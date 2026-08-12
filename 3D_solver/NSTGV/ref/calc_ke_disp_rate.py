import numpy as np
import os
from tqdm import tqdm
from numba import njit
import matplotlib.pyplot as plt
from myvtk import getGrid, getVector, getScalar, extract_number


Q_dir   = "."
Q_files = [f for f in os.listdir(Q_dir) if f.endswith(".vtr")]
Q_files.sort(key=extract_number)


@njit(cache=True, nogil=True, fastmath=True)
def Sutherland(T):
  return 1.716e-5 * (273.2e0 + S) / (T + S) * (T / 273.2e0)**1.5e0


Re   = 1600.e0
M0   = 1.25e0
Rgas = 287.03e0
L0   = 1.524e-3
T    = 530.e0 * 5.e0 / 9.e0 
S    = 111.e0
mu0  = Sutherland(T)
V0   = M0 * np.sqrt(1.4e0 * Rgas * T)
rho0 = mu0 * Re / (V0 * L0)
tc   = L0 / V0
endT = 20.e0 * tc
offset = 3


@njit(cache=True, nogil=True, fastmath=True)
def calc_ke(rho0, rho, u, v, w, offset):
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
def calc_total_dissipation(rho0, Re, x, y, z, rho, u, v, w, p, offset):
  Nz, Ny, Nx = rho.shape
  epsE = 0.e0
  epsD = 0.e0
  for k in range(offset,Nz-offset):
    for j in range(offset,Ny-offset):
      for i in range(offset,Nx-offset):
        T  = p[k,j,i] / (rho[k,j,i] * Rgas)
        mu = Sutherland(T)
        ux, uy, uz, vx, vy, vz, wx, wy, wz = central4th(i, j, k, x, y, z, u, v, w)
        omega = (wy - vz)**2 + (uz - wx)**2 + (vx - uy)**2
        div   = ux + vy + wz
        epsE += mu * omega
        epsD += 4.e0 / 3.e0 * mu * div**2
  N = (Nx - 2 * offset) * (Ny - 2 * offset) * (Nz - 2 * offset)
  epsE /= (rho0 * V0**3 / L0 * N)
  epsD /= (rho0 * V0**3 / L0 * N)
  return epsE, epsD


first_path          = os.path.join(Q_dir, Q_files[0])
Nx, Ny, Nz, x, y, z = getGrid(first_path)
Nt = len(Q_files)
dt = endT / Nt
t  = np.linspace(0.e0, 20.e0, Nt)
Ek   = np.zeros(Nt)
epsE = np.zeros(Nt)
epsD = np.zeros(Nt)
for itr, file, in enumerate(Q_files):
  file_path = os.path.join(Q_dir, file)
  rho       = getScalar(file_path, Nx, Ny, Nz, 'rho')
  u, v, w   = getVector(file_path, Nx, Ny, Nz, 'velocity')
  p         = getScalar(file_path, Nx, Ny, Nz, 'p')
  Ek[itr]   = calc_ke(rho0, rho, u, v, w, offset)
  epsE[itr], epsD[itr] = calc_total_dissipation(rho0, Re, x, y, z, rho, u, v, w, p, offset)


plt.plot(t, Ek)
plt.show()

save_path = os.path.join(Q_dir, "Ek.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       Ek", file=f)
  for i in range(Nt):
    print(f'{t[i]:.3e}', f'{Ek[i]:.3e}', file=f)

plt.plot(t, epsE)
plt.show()

plt.plot(t, epsD)
plt.show()

save_path = os.path.join(Q_dir, "eps.d")
with open(save_path, "w", encoding="UTF-8") as f:
  print("# t       eps", file=f)
  for i in range(Nt):
    print(f'{t[i]:.3e}', f'{epsE[i]:.3e}', f'{epsD[i]:.3e}', file=f)

