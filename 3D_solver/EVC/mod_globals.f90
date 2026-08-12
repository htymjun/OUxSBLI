module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.4_sp
  real(8), parameter         :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 0.1d0
  real(8), parameter :: Ly = 0.1d0
  real(8), parameter :: Lz = 0.01d0 ! inert: VISC='Euler', no viscous kernels, no w-velocity
  integer, parameter :: nx = 258
  integer, parameter :: ny = 258
  ! nz=7 (the theoretical minimum for a 6th-order periodic stencil: 1 real
  ! interior plane) triggers a GPU-kernel edge case -- confirmed by direct
  ! testing -- that corrupts the x/y boundary region within one RK step, even
  ! though the field is perfectly z-uniform. nz=9 (3 interior planes) runs
  ! clean; kept odd so nz/2 lands exactly on the middle interior plane.
  integer, parameter :: nz = 9

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
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

  ! time
  integer, parameter          :: step_offset   = 0
  integer, parameter          :: start_rescale = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.15d0

  ! initial condition
  real(8), parameter :: M0    = 0.05d0
  real(8), parameter :: beta  = 1.d0 / 50.d0
  real(8), parameter :: theta = 0.d0 / 180.d0
  real(8), parameter :: Rc    = 0.005d0
  real(8), parameter :: p0    = 1.d5
  real(8), parameter :: T0    = 300.d0
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: CFL   = 0.03d0
  real(8), parameter :: dt    = CFL * Lx / (dble(nx-1) * u0)
  real(8), parameter :: T     = 2.d0 * Lx / u0
  integer, parameter :: np    = 1
  integer, parameter :: nt    = int(T / (dble(np) * abs(dt)))
end module mod_globals
