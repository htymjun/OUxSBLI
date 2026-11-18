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
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  integer, parameter :: nx = 513
  integer, parameter :: ny = 513
 

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

  ! HR-SLAU2 and HR-AUSM+-up towards High Resolution Unsteady Aerodynamic Simulations
  ! Keiichi Kitamura, Atsushi Hashimoto
  ! JAXA-SP-14-010

  ! initial condition
  real(8), parameter :: M0   = 0.1d0
  real(8), parameter :: rho0 = 1.d0
  real(8), parameter :: u0   = 1.d0
  real(8), parameter :: d1   = pi / 15.d0
  real(8), parameter :: d2   = 0.05d0
  real(8), parameter :: dt   = 0.25d0 * 0.25d0 * 1.d-3
  real(8), parameter :: endT = 8.d0
  integer, parameter :: np   = 10
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

