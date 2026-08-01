import numpy as np
import os
import re
import vtk
from vtk.util import numpy_support


def extract_number(filename, ext='vtr'):
  match = re.search(r'Q(\d+)\.' + ext + r'$', filename)
  if match:
    return int(match.group(1))
  return float('inf')


def getGrid(file_path):
  if str(file_path).endswith('.vts'):
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(str(file_path))
    reader.Update()
    grid = reader.GetOutput()
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    ni, nj, nk = dims
    bounds = grid.GetBounds()  # (xmin, xmax, ymin, ymax, zmin, zmax)
    x = np.linspace(bounds[0], bounds[1], ni)
    y = np.linspace(bounds[2], bounds[3], nj)
    z = np.linspace(bounds[4], bounds[5], nk)
    return ni, nj, nk, x, y, z
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(str(file_path))
  reader.Update()
  grid = reader.GetOutput()
  x = numpy_support.vtk_to_numpy(grid.GetXCoordinates())
  y = numpy_support.vtk_to_numpy(grid.GetYCoordinates())
  z = numpy_support.vtk_to_numpy(grid.GetZCoordinates())
  return len(x), len(y), len(z), x, y, z


def getGrid_Str(file_path):
  # make VTK Structured Grid Reader
  reader = vtk.vtkXMLStructuredGridReader()
  reader.SetFileName(file_path)
  reader.Update()
  # get grid
  grid = reader.GetOutput()
  dims = [0, 0, 0]
  grid.GetDimensions(dims)
  Nx, Ny, Nz = dims
  points = numpy_support.vtk_to_numpy(grid.GetPoints().GetData())
  points = points.reshape((Nz, Ny, Nx, 3))
  return Nx, Ny, Nz, points[:,:,:,0], points[:,:,:,1], points[:,:,:,2]


def getVector(file_path, Nx, Ny, Nz, name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(str(file_path))
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  Q = reader.GetOutput()
  V = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  V = V.reshape((Nz,Ny,Nx,3))
  u = V[:,:,:,0]
  v = V[:,:,:,1]
  w = V[:,:,:,2]
  return u, v, w


def getScalar(file_path, Nx, Ny, Nz, name):
  reader = vtk.vtkXMLRectilinearGridReader()
  reader.SetFileName(str(file_path))
  reader.GetPointDataArraySelection().DisableAllArrays()
  reader.GetPointDataArraySelection().EnableArray(name)
  reader.Update()
  Q = reader.GetOutput()
  a = numpy_support.vtk_to_numpy(Q.GetPointData().GetArray(name))
  a = a.reshape((Nz,Ny,Nx))
  return a


def _reshape_q(point_data, Nz, Ny, Nx):
  """Extract and reshape rho, velocity, p from a VTK point-data object."""
  rho = numpy_support.vtk_to_numpy(point_data.GetArray("rho")).reshape((Nz, Ny, Nx))
  V   = numpy_support.vtk_to_numpy(point_data.GetArray("velocity")).reshape((Nz, Ny, Nx, 3))
  p   = numpy_support.vtk_to_numpy(point_data.GetArray("p")).reshape((Nz, Ny, Nx))
  return rho, V[:,:,:,0], V[:,:,:,1], V[:,:,:,2], p


def getQ(file_path, Nx, Ny, Nz, reader=None):
  if str(file_path).endswith('.vts'):
    r = vtk.vtkXMLStructuredGridReader()
    r.SetFileName(str(file_path))
    r.Update()
    return _reshape_q(r.GetOutput().GetPointData(), Nz, Ny, Nx)
  if reader is None:
    ext = get_ext(file_path)
    if ext == 'vtr':
      reader = vtk.vtkXMLRectilinearGridReader()
    elif ext == 'vts':
      reader = vtk.vtkXMLStructuredGridReader()
    else:
      raise ValueError("Invalid file type:", file_path)
    reader.GetPointDataArraySelection().DisableAllArrays()
    reader.GetPointDataArraySelection().EnableArray("rho")
    reader.GetPointDataArraySelection().EnableArray("velocity")
    reader.GetPointDataArraySelection().EnableArray("p")
  reader.SetFileName(str(file_path))
  reader.Update()
  return _reshape_q(reader.GetOutput().GetPointData(), Nz, Ny, Nx)


def initial_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vt[rs]"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr or Q*.vts files found in {data_dir}")
  return min(files, key=lambda f: extract_number(os.path.basename(f)))


def latest_vtr(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vt[rs]"))
  if not files:
    raise FileNotFoundError(f"No Q*.vtr or Q*.vts files found in {data_dir}")
  return max(files, key=lambda f: extract_number(os.path.basename(f)))


def latest_vts(data_dir):
  import glob
  files = glob.glob(os.path.join(str(data_dir), "Q*.vts"))
  if not files:
    raise FileNotFoundError(f"No Q*.vts files found in {data_dir}")
  return max(files, key=lambda f: extract_number(os.path.basename(f), ext='vts'))
