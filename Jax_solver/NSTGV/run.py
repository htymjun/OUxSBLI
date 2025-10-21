import numpy as np
import jax
import os
import sys
sys.path.append("../src")
from main import main
from set import set_grid, set_init

# grid info
nx    = 66
ny    = 66
nz    = 66
L0    = 1.524e-3
Lx    = 2.e0 * np.pi * L0
Ly    = 2.e0 * np.pi * L0
Lz    = 2.e0 * np.pi * L0

# physical properties
gamma = 1.4e0
Rgas  = 287.03e0
Re    = 1600.e0
M0    = 0.1e0
T0    = 530.e0 * 5.e0 / 9.e0
S     = 111.e0
mu0   = 1.716e-5 * (273.2e0 + S) / (T0 + S) * (T0 / 273.2e0)**1.5
u0    = M0 * np.sqrt(gamma * Rgas * T0)
rho0  = mu0 * Re / (u0 * L0)
p0    = rho0 * Rgas * T0

# time
CFL   = 0.03e0
dt    = CFL * (Lx / float(nx-1)) / u0
dtn   = u0 * dt / L0
Np    = 100
Nt    = int(20.e0 / (float(Np) * dtn))

# device 'cpu' or 'cuda'
device = jax.devices('cpu')[0]

dir = os.path.join(os.getcwd(), 'data')
os.makedirs(dir, exist_ok = True)

main(nx, ny, nz, Lx, Ly, Lz, gamma, Rgas, dt, Nt, Np, M0, rho0, u0, p0, T0, dir, device, set_grid, set_init)

