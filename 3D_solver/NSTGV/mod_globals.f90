module mod_globals
  use cudafor
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler       !
  !               ! kind4 NS          !
  !               ! kind8 LES         !
  !               ! 1 2nd             !
  !               ! 2 4th             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2)  KEEP    !
  !             ! real(2)     SLAU    !
  !             ! real(4)   Weighted  !
  !             ! real(8)   Threshold !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_sensor   ! 1 Ducros            !
  !             ! 2 Albada            !
  !             ! 3 Ducros + Albada   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_accuracy ! kind2 2nd           !
  !             ! kind4 4th           !
  !             ! kind8 6th           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd      ! kind2 non TVD       !
  !             ! kind4 minmod        !
  !             ! kind8 MUSCL4th      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU          !
  !             ! kind4 HR-SLAU2      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter         :: dimension = 3
  integer(4), parameter      :: id_visc   = 2
  integer(2), parameter      :: id_LL     = 0
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.2d0
  integer(kind=4), parameter :: id_accuracy = 0
  integer(kind=8), parameter :: id_tvd      = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: id_rescale  = 0

  ! mesh
  real(8), parameter :: L0 = 1.524d-3
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi * L0
  real(8), parameter :: Ly = 2.d0 * pi * L0
  real(8), parameter :: Lz = 2.d0 * pi * L0
  integer, parameter :: nx = 129
  integer, parameter :: ny = 129
  integer, parameter :: nz = 129

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0
  real(8), parameter :: blt = 0.d0

  ! GPU
  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,8,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,8)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,8,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,8)
  type(dim3), parameter :: threads   = dim3(32,2,2)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter  :: id_recal = 0
  integer(kind=4), parameter  :: id_RungeKutta = 0
  integer, parameter          :: step_offset = 0
  integer, parameter          :: start_rescale = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: Re   = 1600.d0
  real(8), parameter :: M0   = 0.1d0
  real(8), parameter :: T    = 530.d0 * 5.d0 / 9.d0 
  real(8), parameter :: S    = 111.d0
  real(8), parameter :: mu0  = 1.716d-5 * (273.2d0 + S) / (T + S) * (T / 273.2d0)**1.5d0
  real(8), parameter :: V0   = M0 * sqrt(gamma * R * T)
  real(8), parameter :: RHO0 = mu0 * Re / (V0 * L0)
  real(8), parameter :: p0   = RHO0 * R * T

  real(8), parameter :: CFL = 0.03d0
  real(8), parameter :: dt  = CFL * (Lx / dble(nx-1)) / V0
  real(8), parameter :: dtn = V0 * dt / L0
  integer, parameter :: np  = 100
  integer, parameter :: nt  = int(20.d0 / (dble(np) * dtn))
end module mod_globals

