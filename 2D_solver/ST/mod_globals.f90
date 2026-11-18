module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 2
#ifdef PREC_FP32
  integer, parameter         :: sp          = kind(1.e0) ! PREC = 'FP32' (config.fypp)
#else
  integer, parameter         :: sp          = kind(1.d0) ! single or double
#endif
  real(sp), parameter        :: threshold   = 0.4_sp
  real(8), parameter         :: blt         = 1.d-3
  
  ! mesh
  integer, parameter :: nx = 4097
  integer, parameter :: ny = 257
 
  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksEv, blocksFv, blocks

  ! time
  integer, parameter         :: step_offset   = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0
 
  real(8), parameter :: T0   = 300.d0
  real(8), parameter :: C    = 1.461d-6
  real(8), parameter :: S    = 110.3d0
  real(8), parameter :: mu0  = C * T0**1.5 / (T0 + S)
  real(8), parameter :: Re   = 25000.d0
  real(8), parameter :: rho0 = 1.293d0
  real(8), parameter :: p0   = rho0 * R * T0
  real(8), parameter :: rho1 = 0.125d0 * rho0
  real(8), parameter :: p1   = 0.1d0 * p0
  real(8), parameter :: Lx   = Re * mu0 / sqrt(rho0 * p0)
  real(8), parameter :: Ly   = 0.1d0 * Lx
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2136d0 * Lx / sqrt(p0 / rho0)
  integer, parameter :: np   = 100
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals
