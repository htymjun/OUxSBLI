module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0)
  real(sp),   parameter      :: threshold   = 0.4_sp
  real(8),    parameter      :: blt         = 0.d0

  ! mesh — nz is LOCAL per rank; global nz = nranks*(nz-2)
  real(8), parameter :: Lx = 0.1d0
  real(8), parameter :: Ly = 0.1d0
  real(8), parameter :: Lz = 1.0d0
  integer, parameter :: nx = 513 !6    ! small: 6 interior cells
  integer, parameter :: ny = 513 !6
  integer, parameter :: nz = 129 !128 interior + 1 ghost each end


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

  ! time — RK4 z-decomposition
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: start_rescale = 0

  ! Physical (dimensionless Sod shock tube)
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0  ! unused for Euler, kept for NS switch

  real(8), parameter :: rho_L = 1.0d0          ! left state density
  real(8), parameter :: p_L   = 1.0d0          ! left state pressure
  real(8), parameter :: rho_R = 0.125d0         ! right state density
  real(8), parameter :: p_R   = 0.1d0           ! right state pressure

  real(8), parameter :: CFL  = 0.01d0
  real(8), parameter :: dt   = CFL * Lz / (dble(Nz-1) * sqrt(p_L / rho_L))
  real(8), parameter :: endT = 0.1d0
  integer, parameter :: np   = 1!50
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals
