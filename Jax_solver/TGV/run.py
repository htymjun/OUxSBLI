import numpy as np
import jax
import os
import sys
sys.path.append("../src")
from main import main
from set import set_grid, set_init

jax.config.update("jax_enable_x64", True)

# grid info
nx    = 65
ny    = 65
nz    = 65
Lx    = 2.e0 * np.pi
Ly    = 2.e0 * np.pi
Lz    = 2.e0 * np.pi

# physical properties
gamma = 1.4e0
Rgas  = 287.03e0
M0    = 1.e0
rho0  = 0.4e0 
u0    = 1.e0
p0    = 1.e0
T0    = 1.e0

# time
nt    = 200
np    = 10
dt    = 0.01e0

# device 'cpu' or 'cuda'
device = jax.devices('gpu')[0]

dir = os.path.join(os.getcwd(), 'data')
os.makedirs(dir, exist_ok = True)

main(nx, ny, nz, Lx, Ly, Lz, gamma, Rgas, dt, nt, np, M0, rho0, u0, p0, T0, dir, device, set_grid, set_init)

