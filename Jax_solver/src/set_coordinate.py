import numpy as np

def set_J(nx, ny, nz, dx, dy, dz):
  J = np.zeros((nz,ny,nx), dtype=np.float64)
  for k in range(1,nz-1):
    for j in range(1,ny-1):
      for i in range(1,nx-1):
        J[k,j,i] = 8.e0 / \
                   ((dx[i-1] + dx[i]) * (dy[j-1] + dy[j]) * (dz[k-1] + dz[k]))
  J[:,:,0]  = J[:,:,1]
  J[:,:,-1] = J[:,:,-2]
  J[:,0,:]  = J[:,1,:]
  J[:,-1,:] = J[:,-2,:]
  J[0,:,:]  = J[1,:,:]
  J[-1,:,:] = J[-2,:,:]
  return J

