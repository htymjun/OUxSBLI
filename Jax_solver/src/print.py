import numpy as np
import os
import vtk

def print_vtk(x, y, z, gamma, Q, num, dir):
  rho1d = np.float32(Q[:,:,:,0].flatten())
  u1d   = np.float32((Q[:,:,:,1] / Q[:,:,:,0]).flatten())
  v1d   = np.float32((Q[:,:,:,2] / Q[:,:,:,0]).flatten())
  w1d   = np.float32((Q[:,:,:,3] / Q[:,:,:,0]).flatten())
  p1d   = np.float32((gamma - 1.e0) * \
          (Q[:,:,:,4] - 0.5e0 * (Q[:,:,:,1]**2 + Q[:,:,:,2]**2 + Q[:,:,:,3]**2) / Q[:,:,:,0]).flatten())
  
  directory = dir
  os.makedirs(directory, exist_ok=True)
  filename  = "Q" + str(num).zfill(5) + ".vtr" 
  filepath  = os.path.join(directory, filename)

  x_coords = vtk.vtkFloatArray()
  y_coords = vtk.vtkFloatArray()
  z_coords = vtk.vtkFloatArray()
  x_coords.SetName("X-Axis")
  y_coords.SetName("Y-Axis")
  z_coords.SetName("Z-Axis")

  nx = len(x)
  ny = len(y)
  nz = len(z)

  for i in range(nx):
    x_coords.InsertNextValue(x[i])
  for j in range(ny):
    y_coords.InsertNextValue(y[j])
  for k in range(nz):
    z_coords.InsertNextValue(z[k])
  
  grid = vtk.vtkRectilinearGrid()
  grid.SetDimensions(nx, ny, nz)
  grid.SetXCoordinates(x_coords)
  grid.SetYCoordinates(y_coords)
  grid.SetZCoordinates(z_coords)

  rho = vtk.vtkFloatArray()
  rho.SetName("rho")
  for i in range(nx * ny * nz):
    rho.InsertNextValue(rho1d[i])
  grid.GetPointData().AddArray(rho)

  velocity = vtk.vtkFloatArray()
  velocity.SetName("velocity")
  velocity.SetNumberOfComponents(3)
  for i in range(nx * ny * nz):
    velocity.InsertNextTuple3(u1d[i], v1d[i], w1d[i])
  grid.GetPointData().SetVectors(velocity)
  
  p = vtk.vtkFloatArray()
  p.SetName("p")
  for i in range(nx * ny * nz):
    p.InsertNextValue(p1d[i])
  grid.GetPointData().AddArray(p)

  writer = vtk.vtkXMLRectilinearGridWriter()
  writer.SetFileName(filepath)
  writer.SetInputData(grid)
  writer.Write()

