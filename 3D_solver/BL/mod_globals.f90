module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 3
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  real(8), parameter  :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 120.d0 * blt
  real(8), parameter :: Ly = 15.d0 * blt
  ! Lz is physically inert (flow is z-uniform, w=0) but not shrunk to microns:
  ! the z-direction still carries an acoustic eigenvalue +-c0~340 m/s in the
  ! SLAU flux, so CFL_z = c0*dt/dz_wall must stay far looser than the
  ! y-direction's ~0.58 (see nz note below for the actual dz_wall this gives).
  real(8), parameter :: Lz = 1.d0 * blt
  integer, parameter :: nx = 257
  integer, parameter :: ny = 49
  ! nz=7 (the theoretical minimum for a 6th-order periodic stencil: 1 real
  ! interior plane) triggers a GPU-kernel edge case -- confirmed by direct
  ! testing on 3D_solver/EVC -- that corrupts the domain within a single RK
  ! step even for a perfectly z-uniform field. nz=9 (3 interior planes) runs
  ! clean; kept odd so nz/2 lands exactly on the middle interior plane.
  integer, parameter :: nz = 9

  ! flat-plate geometry
  integer, parameter :: i_LE = 2 * nx / 12 + 1 ! first no-slip wall point; set_grid places the leading edge at x = 0

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0

  ! GPU thread blocks (NS, BC_X=True precedent -- 3D_solver/NSTGV/SBLI)
  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  integer, parameter :: step_offset   = 0
  integer, parameter :: start_rescale = 0
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
  real(8), parameter :: dt    = 4.d-8
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = dsqrt(Pr)
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
end module mod_globals
