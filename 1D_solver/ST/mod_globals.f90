module mod_globals
  use cudafor
  implicit none
  integer, parameter :: dimension = 1
  integer, parameter :: sp        = kind(1.d0)

  ! calc_muscl.f90.fypp (repo-root src/) unconditionally compiles
  ! MUSCL3rdThreshold/MUSCL4thThreshold, which reference this constant,
  ! regardless of which TVD dispatch this case actually selects (1D_solver
  ! only ever uses TVD='tvd', the Minmod-limited path -- see
  ! 1D_solver/src/calc_flux_base.f90.fypp's SCHEME='SLAU' guard).
  real(sp), parameter :: threshold = 0.4d0

  ! mesh
  integer, parameter :: nx = 4096

  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsEv = dim3(128,1,1)
  type(dim3), parameter :: threads   = dim3(128,1,1)
  type(dim3) :: blocksE, blocksEv, blocks

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: R     = 287.03d0

  ! Reference viscosity, used only to size the domain length Lx below --
  ! per-point Sutherland's-law mu(T) instead uses the shared
  ! mu0_T0_S_over_T0_2_3 constant in mod_constant.f90.fypp (same physical
  ! value: mu0=1.716d-5 at T0=273.2, S=111).
  real(8), parameter :: mu0 = 1.716d-5

  real(8), parameter :: Re   = 25000.d0
  real(8), parameter :: Tlr  = 300.d0
  real(8), parameter :: rho0 = 1.293d0
  real(8), parameter :: p0   = rho0 * R * Tlr
  real(8), parameter :: rho1 = 0.125d0 * rho0
  real(8), parameter :: p1   = 0.1d0 * p0
  real(8), parameter :: a    = sqrt(R * Tlr)
  real(8), parameter :: Lx   = Re * mu0 / (rho0 * a)
  real(8), parameter :: dx   = Lx / dble(nx-1)
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * dx / a

  integer, parameter :: nt = 2500
end module mod_globals
