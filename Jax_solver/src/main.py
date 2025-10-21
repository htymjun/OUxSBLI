import numpy as np
import jax
import jax.numpy as jnp
from jax import lax
from jax._src import lib
import time
from tqdm import tqdm
from calc_time_dev import Runge_Kutta
from print import print_vtk
from set_coordinate import set_J


def main(nx, ny, nz, Lx, Ly, Lz, gamma, Rgas, dt, Nt, Np, M0, rho0, u0, p0, T0, dir, device, set_grid, set_init):
  start = time.time()

  x, y, z, dx, dy, dz = set_grid(nx, ny, nz, Lx, Ly, Lz)
  J = set_J(nx, ny, nz, dx, dy, dz)
  Q = set_init(nx, ny, nz, x, y, z, gamma, Rgas, M0, rho0, u0, p0, T0)
  print_vtk(x, y, z, gamma, Q, 0, dir)

  # copy on device
  dxj = jnp.array(dx, dtype=jnp.float32, device=device)
  dyj = jnp.array(dy, dtype=jnp.float32, device=device)
  dzj = jnp.array(dz, dtype=jnp.float32, device=device)
  Jj  = jnp.array(J,  dtype=jnp.float32, device=device)
  Qj  = jnp.array(Q / J[:,:,:,None], dtype=jnp.float32, device=device)

  Cp = gamma * Rgas / (gamma - 1.e0)
  Pr = 0.71e0
  x0 = (gamma, Rgas, Cp, Pr, M0, rho0, u0, p0, T0, dt, dxj, dyj, dzj, Jj, Qj)
  for itr in tqdm(range(Np)):
    xs = lax.fori_loop(0, Nt, Runge_Kutta, x0)
    x0 = xs
    Q  = jax.device_get(xs[14])
    Q  = jax.block_until_ready(Q)

    print_vtk(x, y, z, gamma, Q * J[:,:,:,None], itr+1, dir)

  elapsed_time = time.time() - start
  print("elapsed_time:", elapsed_time)

