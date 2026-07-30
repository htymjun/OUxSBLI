
module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0)
  real(sp), parameter        :: threshold   = 0.4_sp
  real(8), parameter         :: blt         = 0.d0

  ! mesh
  real(8), parameter :: L0 = 1.524d-3
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi * L0
  real(8), parameter :: Ly = 2.d0 * pi * L0
  real(8), parameter :: Lz = 2.d0 * pi * L0
  integer, parameter :: nx = 513
  integer, parameter :: ny = 513
  integer, parameter :: nz = 513

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0

  ! GPU
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

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: Re   = 1600.d0
  real(8), parameter :: M0   = 1.25d0
  real(8), parameter :: T    = 530.d0 * 5.d0 / 9.d0
  real(8), parameter :: S    = 111.d0
  real(8), parameter :: mu0  = 1.716d-5 * (273.2d0 + S) / (T + S) * (T / 273.2d0)**1.5d0
  real(8), parameter :: V0   = M0 * sqrt(gamma * R * T)
  real(8), parameter :: RHO0 = mu0 * Re / (V0 * L0)
  real(8), parameter :: p0   = RHO0 * R * T
  real(8), parameter :: CFL  = 0.03d0
  real(8), parameter :: dt   = CFL * (Lx / dble(nx-1)) / V0
  real(8), parameter :: dtn  = V0 * dt / L0
  integer, parameter :: np   = 1!100
  integer, parameter :: nt   = 1!int(20.d0 / (dble(np) * dtn))
end module mod_globals
