module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 2
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  real(8), parameter  :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 150.d0 * blt !140!20
  real(8), parameter :: Ly = 12.5d0 * blt !10!5
  ! DNS
  integer, parameter :: nx = 1921 !257
  integer, parameter :: ny = 321 !129

  ! RTX 4090
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsEv = dim3(64,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksEv, blocksFv, blocks

  ! time
  integer, parameter :: step_offset = 0
  real(8), parameter :: endT  = 8.d-3 !0.3d-3
  integer, parameter :: np    = 800 !30
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.d0
  real(8), parameter :: p_tot = 100.d3
  real(8), parameter :: T_tot = 518.4d0 !295.d0
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = T_tot /  (1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 3.d-9
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = dsqrt(Pr) !0.89d0
  real(8), parameter :: Taw   = 1.676 * T0 !T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2) !M=2
  ! oblique shock
  real(8), parameter :: beta  = dacos(-1.d0) * 31.65d0 / 180.d0 !theta=2 !M2, 40.03(theta=10.65)
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

module mod_shock
  use mod_globals, only : gamma, R, T0, T2, M0, Ms, Ms2, beta, theta
  implicit none

  real(8) :: p0_init
  real(8) :: rho0_init
  real(8) :: p2_init
  real(8) :: rho2_init
  real(8) :: u0_init
  real(8) :: v0_init
  real(8) :: u1_init
  real(8) :: v1_init
  real(8) :: a1_init
  real(8) :: u2_init
  real(8) :: v2_init
  real(8) :: u_magnitude_init
  real(8) :: ux_init
  real(8) :: uy_init

  contains

  subroutine Qin_init(Qp)

    real(8), intent(in) :: Qp(4)

    p0_init   = Qp(4)
    !rho0_init = p0_init / (R * T0)
    rho0_init = Qp(1)
    p2_init   = p0_init * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
    rho2_init = p2_init / (R * T2)

    u0_init    = Qp(2)
    v0_init    = Qp(3)

    u1_init    = u0_init * dsin(beta)
    v1_init    = u0_init * dcos(beta)
    a1_init    = u0_init / M0
    u2_init    = u1_init - 2.d0 * a1_init * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
    v2_init    = u0_init * dcos(beta)
    u_magnitude_init = sqrt(u2_init**2 + v2_init**2)
    ux_init    = u_magnitude_init * dcos(theta)
    uy_init    = - u_magnitude_init * dsin(theta)    

  end subroutine Qin_init

end module mod_shock