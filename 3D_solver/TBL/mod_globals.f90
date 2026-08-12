module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 3
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  real(8), parameter  :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 10.d0 * blt
  real(8), parameter :: Ly = 4.d0 * blt
  real(8), parameter :: Lz = 1.d0 * blt
  ! DNS
  integer, parameter :: nx = 257
  integer, parameter :: ny = 161
  integer, parameter :: nz = 65

  integer, parameter :: nre1 = int(0.5 * nx)
  integer, parameter :: nre2 = int(0.9 * nx)
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,2)
  type(dim3), parameter :: threadsEv = dim3(128,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,2)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  integer, parameter :: step_offset   = 0
  integer, parameter :: start_rescale = 10
  real(8), parameter :: endT  = 0.1d-3
  integer, parameter :: np    = 10
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.5d0
  real(8), parameter :: p_tot = 100.d3
  real(8), parameter :: T_tot = 295.d0
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = T_tot /  (1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 3.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))
  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = 0.89d0
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
end module mod_globals

