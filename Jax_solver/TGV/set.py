import numpy as np
import jax

def set_grid(nx, ny, nz, Lx, Ly, Lz):
  dx = np.zeros(nx-1, dtype=np.float64)
  dy = np.zeros(ny-1, dtype=np.float64)
  dz = np.zeros(nz-1, dtype=np.float64)
  dx[:] = Lx / np.float64(nx-2) 
  dy[:] = Ly / np.float64(ny-2) 
  dz[:] = Lz / np.float64(nz-2) 

  x = np.zeros(nx, dtype=np.float64)
  x[1:-1] = np.linspace(0, Lx, nx-2, dtype=np.float64)
  x[0]  = x[1]  - dx[0]
  x[-1] = x[-2] + dx[0]

  y = x
  z = y
  return x, y, z, dx, dy, dz


def set_init(nx, ny, nz, x, y, z, gamma, Rgas, M0, rho0, u0, p0, T0):
  Q = np.zeros((nz, ny, nx, 5), dtype=np.float64)
  for k in range(nz):
    for j in range(ny):
      for i in range(nx):
        Q[k,j,i,0] =  rho0
        Q[k,j,i,1] =  rho0 * M0 * np.sin(x[i]) * np.cos(y[j]) * np.cos(z[k])
        Q[k,j,i,2] = -rho0 * M0 * np.cos(x[i]) * np.sin(y[j]) * np.cos(z[k])
        Q[k,j,i,3] = 0.e0
        Q[k,j,i,4] = (1.e0 / gamma + 0.0625e0 * rho0 * (M0**2) * \
        (np.cos(2.e0 * x[i]) + np.cos(2.e0 * y[j])) * (np.cos(2.e0 * z[k]) + 2.e0) \
        ) / (gamma - 1.e0) + 0.5e0 * (Q[k,j,i,1]**2 + Q[k,j,i,2]**2) / Q[k,j,i,0]

  Q[1:-1,1:-1,0,:] = Q[1:-1,1:-1,-2,:]
  Q[1:-1,1:-1,-1,:] = Q[1:-1,1:-1,1,:]
  Q[1:-1,0,1:-1,:] = Q[1:-1,-2,1:-1,:]
  Q[1:-1,-1,1:-1,:] = Q[1:-1,1,1:-1,:]

  Q[1:-1,0,0,:] = Q[1:-1,-2,-2,:]
  Q[1:-1,0,-1,:] = Q[1:-1,-2,1,:]
  Q[1:-1,-1,0,:] = Q[1:-1,1,-2,:]
  Q[1:-1,-1,-1,:] = Q[1:-1,1,1,:]

  Q[0,1:-1,1:-1,:] = Q[-2,1:-1,1:-1,:]
  Q[-1,1:-1,1:-1,:] = Q[1,1:-1,1:-1,:]
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

  Q = Q.at[0,1:-1,1:-1,:].set(Q[-2,1:-1,1:-1,:])
  Q = Q.at[-1,1:-1,1:-1,:].set(Q[1,1:-1,1:-1,:])
  return Q

