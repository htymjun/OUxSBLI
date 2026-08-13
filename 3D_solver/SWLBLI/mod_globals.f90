module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 3
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  real(8), parameter  :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 176.d0 * blt !20
  real(8), parameter :: Ly = 80.d0 * blt !5
  real(8), parameter :: Lz = 2.d0 * blt
  ! DNS
  integer, parameter :: nx = 257
  integer, parameter :: ny = 257
  integer, parameter :: nz = 9
  
  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0

  ! flat-plate geometry
  real(8), parameter :: x_in = -16.d0 * blt ! x(1); the inlet-shock threshold (Xsh - x_in)*tan(beta) depends on this
  real(8), parameter :: Xsh  = 80.d0 * blt  ! inviscid shock impingement point on the wall
  integer, parameter :: i_LE = nint(-x_in * dble(nx-1) / Lx) + 1 ! first no-slip wall point; leading edge stays at x ~ 0 for any nx

  ! RTX 4090
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,2)
  type(dim3), parameter :: threadsEv = dim3(64,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,2)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  integer, parameter :: step_offset = 0
  real(8), parameter :: endT  = 1.d-3
  integer, parameter :: np    = 100
  real(8), parameter :: R     = 287.15d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.15d0
  real(8), parameter :: p_tot = 25.d3
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = 288.15d0
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 5.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = dsqrt(Pr) !0.89d0
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
  ! oblique shock
  real(8), parameter :: beta  = dacos(-1.d0) * 30.8d0 / 180.d0 !theta4, 30.96 !30.8_baseline
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
