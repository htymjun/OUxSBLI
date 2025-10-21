# OUxSBLI
OUxSBLI is a GPU-accelerated CFD code written in CUDA Fortran. It employs explicit high-order finite-difference schemes on a rectilinear grid.

## Dependency
1. HPC SDK (version 24.* and 25.* are better)
2. ParaView (for visualization output files are XML VTK format)
3. CUDA & pybind11... (for Python interface. This is under development.)

## Usage
1. Go to a directory (./3D_solver/NSTGV)
   1. 1D_solver and 2D_solver are under development.
   2. pyETGV and pyNSTGV are under development.
2. Edit mod_globals.f90
   1. Choose equation type "id_visc" (Euler or NS)
   2. Choose scheme "id_scheme", "id_accuracy", "id_tvd", and "id_slau"
   3. Choose grid size "nx", "ny", and "nz"
   4. Optimize block-size
   5. Choose time integration method
3. Edit set.f90
   1. Set grid
   2. Set initial conditions
   3. Set boundary conditions
4. $ make
5. $ bash calc.sh
6. In some directory, you can get nsys and ncu iformation by running profile.sh

## Discretization
### Spatial (Convection terms)
   1. Kinetic energy and entropy preserving (KEEP) scheme
   2. Simple low-dissipation AUSM (SLAU) scheme
   3. Roe scheme
### Spatial (Viscous terms)
   1. Sandham's Laplacian form (Recently added. Validation is not enougth)
   2. ME4-Base (./3D_solver/src/old_visc)
   3. Gaitonde and Visbal's 2nd-order scheme (./3D_solver/src/old_visc)
### Spatial SGS
   1. Selective mixed scale model (Under development)
### Temporal
   1. 3-3 TVD Runge-Kutta
   2. 4-4 Runge-Kutta
   3. Gauss Legendre Runge-Kutta (Under development)

## Validation and Verification

### 1. Euler Vortex Convection

#### Grid Convergence
<div align="center">
  <img src="./img/Grid_convergence.png" alt="L2_norm" width="450">  
</div>

Error between numerical solution and theoretical solution at 50 period; solid line for theoretical rate, black square points for 2nd-order KEEP, green square points for 4th-order KEEP, red square points for 5th-order HR-SLAU2 and blue square points for 6th-order KEEP

### 2. Sod Shock Tube

#### Visualization
<div align="center">
  <img src="./img/ST_rho.png" alt="density" width="450">  
</div>

Instantaneous density field at t=0.2; red solid line for 3rd-order SLAU, red dashed line for 4th-order SLAU, blue solid line for 3rd-order HR-SLAU2 and blue dashed line for 4th-order HR-SLAU2

### 3. 3D Viscous Taylor-Green Vortex

#### Evolution of Kinetic Energy and Enstrophy
<div align="center">
  <img src="./img/TGV_ke.png" alt="ke" width="450">  
</div>

Time evolution of the total kinetic energy as a function of the dimensionless time; red solid line for 5th-order HR-SLAU2 with grid points [256, 256, 256], red dashed line for 5th-order HR-SLAU2 with grid points [128, 128, 128], blue solid line for 6th-order KEEP with grid points [256, 256, 256], blue dashed line for 6th-order KEEP with grid points [128, 128, 128] and black solid line for Reference*

<div align="center">
  <img src="./img/TGV_enstrophy.png" alt="enstorophy" width="450">  
</div>

Time evolution of the enstrophy as a function of the dimensionless time; red solid line for 5th-order HR-SLAU2 with grid points [256, 256, 256], red dashed line for 5th-order HR-SLAU2 with grid points [128, 128, 128], blue solid line for 6th-order KEEP with grid points [256, 256, 256], blue dashed line for 6th-order KEEP with grid points [128, 128, 128] and black solid line for Reference*

### 4. M=1.9 Supersonic Turbulent Boundary Layer

#### Log-Law
<div align="center">
  <img src="./img/log_law.png" alt="log-law" width="450">  
</div>

Distribution of Van Driest transformed mean streamwise velocity normalized by wall unit*

#### Reynolds Stress
<div align="center">
  <img src="./img/Reynolds_stress.png" alt="Reynolds-stress" width="450">  
</div>

Distribution of the normalized Reynolds normal stress; blue lines for present simulation; black lines for Reference; solid lines for; dashed lines for; dotted lines for*

<div align="center">
  <img src="./img/stress_balance.png" alt="stress-balance" width="450">  
</div>

Distribution of the normalized shear stress; blue lines for present simulation; black lines for Reference; solid lines for Reynolds stress; dashed lines for viscous stress; dotted lines for total stress*

#### TKE Budget
<div align="center">
  <img src="./img/TKE_budget.png" alt="TKE_budget" width="450">  
</div>

Turbulent kinetic energy budget; solid lines for present simulation, dashed lines for References, dotted lines for Reference*

---

## Reference

1. [Yuichi Kuya, Kosuke Totani, Soshi Kawai, *Kinetic energy and entropy preserving schemes for compressible flows by split convective forms*, Journal of Computational Physics, 2018](https://www.sciencedirect.com/science/article/abs/pii/S0021999118305916)
2. [Yuichi Kuya, Soshi Kawai, *High-order accurate kinetic-enrgy and entropy preserving (KEEP) schemes on curvilinear grids*, Journal of Computational Physics, 2021](https://www.sciencedirect.com/science/article/abs/pii/S0021999121003776)
3. [Christian T. Jacobs, Satya P. Jammy, Neil D. Sandham, *OpenSBLI: A framework for the automated derivation and parallel execution of finite difference solvers on a range of computer architectures*, 2017](https://www.sciencedirect.com/science/article/pii/S187775031630299X?via%3Dihub)
4. [Yves Allaneau, Antony Jameson, *Direct Numerical Simulations of a Two-Dimensional Viscous Flow in a Shocktube Using Kinetic Energy Preserving Scheme*, 2012](https://arc.aiaa.org/doi/abs/10.2514/6.2009-3797)
5. [Yoshiharu Tamaki, Soshi Kawai, *Wall-modeled LES of transonic buffet over NASA-CRM using Cartesian-grid-based flow solver FFVHC-ACE*, 2023](https://arc.aiaa.org/doi/10.2514/6.2023-0429)
6. [CUDA FORTRAN PROGRAMMING GUIDE AND REFERENCE, 2017](https://docs.nvidia.com/hpc-sdk/pgi-compilers/2017/pgi17cudaforug.pdf)
