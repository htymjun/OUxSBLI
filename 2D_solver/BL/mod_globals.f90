module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 2
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  real(8), parameter  :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 120.d0 * blt
  real(8), parameter :: Ly = 12.d0 * blt
  ! DNS
  integer, parameter :: nx = 129
  integer, parameter :: ny = 33

  ! flat-plate geometry
  integer, parameter :: i_LE = 2 * nx / 12 + 1 ! first no-slip wall point; set_grid places the leading edge at x = 0

  ! GPU thread blocks
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsEv = dim3(64,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksEv, blocksFv, blocks

  ! time
  integer, parameter :: step_offset = 0
  real(8), parameter :: endT  = 1.d-2
  integer, parameter :: np    = 100
  real(8), parameter :: R     = 287.15d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 0.1
  real(8), parameter :: p_tot = 25.d3
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = 288.15d0
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 3.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = dsqrt(Pr)
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
end module mod_globals
