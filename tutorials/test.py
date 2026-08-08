import pathlib
import numpy as np
import ouxsbli
from ouxsbli import Case
from ouxsbli.patcher import patch


WORKDIR = "./ouxsbli_bl"

case = Case(
  source  = "2D_solver/BL",    # template — never modified; resolved relative to the repo root
  workdir = WORKDIR,           # new run directory (created fresh each time)
  # Physics
  visc    = "ns",              # Navier-Stokes
  scheme  = "slau",            # SLAU convective scheme
  accuracy = 2,                # 2nd-order spatial
  # Grid (smaller than the template for speed)
  nx = 129,
  ny = 65,
  # Time
  np = 10,                     # output every 10 steps
  nt = 5,                      # total output intervals
)

# Build: copies source → patches mod_globals.f90 → runs cmake + make
# Takes ~1-3 minutes depending on hardware
print("Building...")
case.build()
print("Build complete. a.out:", (pathlib.Path(WORKDIR) / "a.out").exists())

patched_src = (pathlib.Path(WORKDIR) / "mod_globals.f90").read_text()

print("Running simulation...")
case.run(nranks=2)   # mpirun -n 2 ./a.out

vtk_files = sorted((pathlib.Path(WORKDIR) / "data").glob("Q*.vtr"))
print(f"\nSimulation complete. {len(vtk_files)} VTK files produced:")
for f in vtk_files:
  print(f"  {f.name}")

