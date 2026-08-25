! Pantano & Sarkar (2002), JFM 451:329-371.
! Temporally evolving compressible turbulent shear layer.
!
! Default parameters reproduce case A7:
!   Mc = 0.7, s = rho2/rho1 = 1, Re = rho0*DeltaU*delta_theta0/mu0 = 800,
!   Pr = 0.7, gamma = 1.4,
!   Lx x Ly x Lz = 172 x 129 x 86 delta_theta0,
!   Nx x Ny x Nz = 256 x 192 x 128 physical points.
!
! Change case_Mc and density_ratio to obtain A3/A11/B2/B4/B8.
module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0)
  real(sp), parameter        :: threshold   = 0.4_sp
  real(8), parameter         :: blt         = 0.d0

  ! reference nondimensional parameters in the paper
  real(8), parameter :: pi            = acos(-1.d0)
  real(8), parameter :: gamma         = 1.4d0
  real(8), parameter :: Pr            = 0.7d0
  real(8), parameter :: Prt           = 0.9d0
  real(8), parameter :: Re_theta0     = 800.d0
  real(8), parameter :: case_Mc       = 0.7d0
  real(8), parameter :: density_ratio = 1.d0

  ! Solver variables are dimensional.  Set p0=rho0=1 and choose DeltaU from Mc.
  ! T_ref fixes the Sutherland viscosity used by calc_physical_quantities.f90;
  ! delta_theta0 is then chosen so rho0*DeltaU*delta_theta0/mu_ref = 800.
  real(8), parameter :: rho0 = 1.d0
  real(8), parameter :: p0   = 1.d0
  real(8), parameter :: Tref = 300.d0
  real(8), parameter :: R    = p0 / (rho0 * Tref)
  real(8), parameter :: c0   = sqrt(gamma * p0 / rho0)
  real(8), parameter :: du   = 2.d0 * case_Mc * c0
  real(8), parameter :: mu_ref = 1.716d-5 * 383.6d0 * 273.2d0**(-1.5d0) &
                                  / (Tref + 110.4d0) * (Tref * sqrt(Tref))
  real(8), parameter :: delta_theta0 = Re_theta0 * mu_ref / (rho0 * du)

  ! stream labels follow the paper: stream 1 is upper/high-speed side, stream 2 lower.
  real(8), parameter :: rho1 = 2.d0 * rho0 / (1.d0 + density_ratio)
  real(8), parameter :: rho2 = density_ratio * rho1
  real(8), parameter :: u1   = -0.5d0 * du
  real(8), parameter :: u2   =  0.5d0 * du
  real(8), parameter :: T1   = p0 / (R * rho1)
  real(8), parameter :: T2   = p0 / (R * rho2)

  ! mesh: include ghost cells for sixth-order stencils, so array sizes are +6.
  real(8), parameter :: Lx = 172.d0 * delta_theta0
  real(8), parameter :: Ly = 129.d0 * delta_theta0
  real(8), parameter :: Lz = 86.d0  * delta_theta0
  integer, parameter :: nx = 256 + 6
  integer, parameter :: ny = 192 + 6
  integer, parameter :: nz = 128 + 6

  integer, parameter :: nre1   = 1
  integer, parameter :: nre2   = nx
  integer, parameter :: rerank = 0

  ! GPU
  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! initial disturbance: paper uses 10% broadband fluctuations, localized over one
  ! initial shear-layer thickness.  set.f90 builds a deterministic solenoidal proxy.
  real(8), parameter :: disturbance_intensity = 0.10d0
  integer, parameter :: peak_wavelengths_x    = 24

  ! time: normalized paper time is tau = t*DeltaU/delta_theta0.
  integer, parameter :: step_offset   = 0
  integer, parameter :: start_rescale = 0
  real(8), parameter :: CFL           = 0.15d0
  real(8), parameter :: max_speed     = 0.5d0 * du + c0
  real(8), parameter :: dt            = CFL * min(Lx / dble(nx-6), min(Ly / dble(ny-6), Lz / dble(nz-6))) / max_speed
  real(8), parameter :: end_tau       = 600.d0
  real(8), parameter :: endT          = end_tau * delta_theta0 / du
  integer, parameter :: np            = 100
  integer, parameter :: nt            = max(1, int(endT / (dble(np) * dt)))
end module mod_globals
