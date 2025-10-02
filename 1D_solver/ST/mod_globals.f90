module mod_globals
  use cudafor
  implicit none
  ! mesh
  integer, parameter :: nx = 4097

  ! GPU
  type(dim3) :: blocks  = dim3((nx-1)/128,1,1)
  type(dim3) :: threads = dim3(128,1,1)

  ! time
  integer, parameter :: np = 1

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: R_gas = 287.03d0

  ! initial condition
  real(8), parameter :: T0   = 300.d0
  real(8), parameter :: C    = 1.461d-6
  real(8), parameter :: S    = 110.3d0
  real(8), parameter :: mu0  = C * T0**1.5 / (T0 + S)
  real(8), parameter :: Re   = 25000.d0
  real(8), parameter :: rho0 = 1.293d0
  real(8), parameter :: p0   = rho0 * R_gas * T0
  real(8), parameter :: rho1 = 0.125d0 * rho0
  real(8), parameter :: p1   = 0.1d0 * p0
  real(8), parameter :: Lx   = Re * mu0 / sqrt(rho0 * p0)
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2136d0 * Lx / sqrt(p0 / rho0)
  real(8), parameter :: nt   = endT / (dble(np) * dt)
  real(8), parameter :: dx   = Lx / dble(nx-1)
  real(8), parameter :: dtdx = dt / dx
end module mod_globals

