module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0)
  real(sp), parameter        :: threshold   = 0.4_sp
  real(8), parameter         :: blt         = 1.d-3

  ! Mesh: real x/y domain, periodic-thin z
  real(8), parameter :: Lx = 5.d0 * blt
  real(8), parameter :: Ly = 2.d0 * blt
  real(8), parameter :: Lz = 1.d0 * blt
  integer, parameter :: nx = 257
  integer, parameter :: ny = 129
  ! nz=7 (the theoretical minimum for a 6th-order periodic stencil) triggers a
  ! GPU-kernel edge case -- confirmed by direct testing on 3D_solver/EVC and
  ! 3D_solver/BL -- that corrupts the domain within a single RK step even for
  ! a perfectly z-uniform field. nz=9 (3 interior planes) runs clean.
  integer, parameter :: nz = 9

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0

  ! GPU thread blocks (real-x/y-BC precedent, same as 3D_solver/BL)
  type(dim3), parameter :: threadsE  = dim3(32, 1, 1)
  type(dim3), parameter :: threadsF  = dim3(32, 4, 1)
  type(dim3), parameter :: threadsG  = dim3(32, 1, 4)
  type(dim3), parameter :: threadsEv = dim3(32, 1, 1)
  type(dim3), parameter :: threadsFv = dim3(32, 4, 1)
  type(dim3), parameter :: threadsGv = dim3(32, 1, 4)
  type(dim3), parameter :: threads   = dim3(32, 4, 1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! Time stepping
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 0
  real(8), parameter :: endT  = 0.1d-3
  integer, parameter :: np    = 10
  real(8), parameter :: dt    = 3.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! physical properties
  real(8), parameter :: M0    = 2.d0
  real(8), parameter :: p_tot = 100.d3
  real(8), parameter :: T_tot = 295.d0
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = T_tot /  (1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: rho0  = p0 / (R * T0)
  ! oblique shock
  real(8), parameter :: beta  = dacos(-1.d0) * 37.2d0 / 180.d0
  real(8), parameter :: Ms    = M0 * dsin(beta)
  real(8), parameter :: Ms2   = Ms**2
  real(8), parameter :: theta = datan(2.d0 * (1.d0 / dtan(beta)) * (Ms2 - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: T2    = T0 * (1.d0 + 2.d0 * (gamma - 1.d0) * (Ms2 - 1.d0) * (1.d0 + gamma * Ms2) / (Ms2 * (gamma + 1.d0)**2))
  real(8), parameter :: p2    = p0 * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
  real(8), parameter :: rho2  = p2 / (R * T2)
  real(8), parameter :: u1    = u0 * dsin(beta)
  real(8), parameter :: v1    = u0 * dcos(beta)
  real(8), parameter :: a1    = u0 / M0
  real(8), parameter :: u2    = u1 - 2.d0 * a1 * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
  real(8), parameter :: v2    = u0 * dcos(beta)
  real(8), parameter :: u_magnitude = sqrt(u2**2 + v2**2)
  real(8), parameter :: ux    = u_magnitude * dcos(theta)
  real(8), parameter :: uy    = - u_magnitude * dsin(theta)
end module mod_globals
