module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0)
  real(sp),   parameter      :: threshold   = 0.4_sp
  real(8),    parameter      :: blt         = 1.d-3

  ! mesh — nz is LOCAL per rank; global nz cells = (nranks/2)*(nz-2*ovlp), ovlp = 2 (ORDER=4)
  ! 2D_solver/ST: nx = 4097 (4096 cells) -> 2 z-ranks x 2048 cells
  integer, parameter :: nx = 6
  integer, parameter :: ny = 6
  integer, parameter :: nz = 2052


  integer, parameter :: nre1   = 1
  integer, parameter :: nre2   = nx
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time — RK z-decomposition
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 0

  ! physical properties (same as 2D_solver/ST, tube axis x -> z)
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  real(8), parameter :: T0   = 300.d0
  real(8), parameter :: C    = 1.461d-6
  real(8), parameter :: S    = 110.3d0
  real(8), parameter :: mu0  = C * T0**1.5 / (T0 + S)
  real(8), parameter :: Re   = 25000.d0
  real(8), parameter :: rho0 = 1.293d0
  real(8), parameter :: p0   = rho0 * R * T0
  real(8), parameter :: rho1 = 0.125d0 * rho0
  real(8), parameter :: p1   = 0.1d0 * p0

  real(8), parameter :: rho_L = rho0          ! left state density
  real(8), parameter :: p_L   = p0            ! left state pressure
  real(8), parameter :: rho_R = rho1          ! right state density
  real(8), parameter :: p_R   = p1            ! right state pressure

  real(8), parameter :: Lz   = Re * mu0 / sqrt(rho0 * p0)
  real(8), parameter :: Lx   = 0.1d0 * Lz
  real(8), parameter :: Ly   = 0.1d0 * Lz
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lz / (dble(2*(nz-4)) * sqrt(p0 / rho0))
  real(8), parameter :: endT = 0.2136d0 * Lz / sqrt(p0 / rho0)
  integer, parameter :: np   = 100
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals
