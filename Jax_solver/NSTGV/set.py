import numpy as np
import jax
from numba import njit


def set_grid(nx, ny, nz, Lx, Ly, Lz):
  dx = np.zeros(nx-1, dtype=np.float32)
  dy = np.zeros(ny-1, dtype=np.float32)
  dz = np.zeros(nz-1, dtype=np.float32)
  dx[:] = Lx / np.float32(nx-2) 
  dy[:] = Ly / np.float32(ny-2) 
  dz[:] = Lz / np.float32(nz-2) 

  x = np.zeros(nx, dtype=np.float32)
  x[1:-1] = np.linspace(0, Lx, nx-2, dtype=np.float32)
  x[0]  = x[1]  - dx[0]
  x[-1] = x[-2] + dx[0]

  y = x
  z = y
  return x, y, z, dx, dy, dz


@njit(cache=True, nogil=True)
def set_init(nx, ny, nz, x, y, z, gamma, Rgas, M0, rho0, u0, p0, T0):
  Q = np.zeros((nz, ny, nx, 5), dtype=np.float32)
  L0 = 1.524e-3
  for k in range(1, nz-1):
    for j in range(1, ny-1):
      for i in range(1, nx-1):
        p = p0 + rho0 * (u0**2) * (np.cos(2.e0 * x[i] / L0) + np.cos(2.e0 * y[j] / L0)) \
                                * (np.cos(2.e0 * z[k] / L0) + 2.e0) / 16.e0
        rho = p / (Rgas * T0)

        Q[k,j,i,0] =  rho
        Q[k,j,i,1] =  rho * u0 * np.sin(x[i] / L0) * np.cos(y[j] / L0) * np.cos(z[k] / L0)
        Q[k,j,i,2] = -rho * u0 * np.cos(x[i] / L0) * np.sin(y[j] / L0) * np.cos(z[k] / L0)
        Q[k,j,i,3] = 0.e0
        Q[k,j,i,4] = p / (gamma - 1.e0) + 0.5e0 * (Q[k,j,i,1]**2 + Q[k,j,i,2]**2) / Q[k,j,i,0]

  for k in range(1, nz-1):
    for j in range(1, ny-1):
      Q[k,j,0,:]  = Q[k,j,-2,:]
      Q[k,j,-1:]  = Q[k,j,1,:]
    
    for i in range(1, nx-1):
      Q[k,0,i,:]  = Q[k,-2,i,:]
      Q[k,-1,i,:] = Q[k,1,i,:]

    Q[k,0,0,:]    = Q[k,-2,-2,:]
    Q[k,0,-1,:]   = Q[k,-2,1,:]
    Q[k,-1,0,:]   = Q[k,1,-2,:]
    Q[k,-1,-1,:]  = Q[k,1,1,:]

  for j in range(ny):
    for i in range(nx):
      Q[0,j,i,:]  = Q[-2,j,i,:]
      Q[-1,j,i,:] = Q[1,j,i,:]
  return Q


@jax.jit
def set_bc(gamma, Rgas, M0, rho0, u0, p0, T0, J, Q):
  Q = Q.at[1:-1,1:-1,0,:].set(Q[1:-1,1:-1,-2,:])
  Q = Q.at[1:-1,1:-1,-1,:].set(Q[1:-1,1:-1,1,:])
  Q = Q.at[1:-1,0,1:-1,:].set(Q[1:-1,-2,1:-1,:])
  Q = Q.at[1:-1,-1,1:-1,:].set(Q[1:-1,1,1:-1,:])

  Q = Q.at[1:-1,0,0,:].set(Q[1:-1,-2,-2,:])
  Q = Q.at[1:-1,0,-1,:].set(Q[1:-1,-2,1,:])
  Q = Q.at[1:-1,-1,0,:].set(Q[1:-1,1,-2,:])
  Q = Q.at[1:-1,-1,-1,:].set(Q[1:-1,1,1,:])

  Q = Q.at[0,:,:,:].set(Q[-2,:,:,:])
  Q = Q.at[-1,:,:,:].set(Q[1,:,:,:])
  return Q

