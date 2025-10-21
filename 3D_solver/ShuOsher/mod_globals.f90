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
  integer(2), parameter      :: id_visc   = 1
  real(2), parameter         :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=8), parameter :: id_tvd      = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: id_rescale  = 0
  real(8), parameter         :: blt         = 2.d-3

  ! exchange
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_exchange ! kind2 off !
  !             ! kind4 on  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_exchange = 0
  ! mesh
  real(8), parameter :: Lx   = 1.d0
  real(8), parameter :: Ly   = 0.1d0 * Lx
  real(8), parameter :: Lz   = 0.1d0 * Lx
  integer, parameter :: nx = 257
  ! 6th-order accuracy
  integer, parameter :: ny = 7
  integer, parameter :: nz = 7

  integer, parameter :: nre1 = int(1.d0 * dble(nx) / 7.d0)
  integer, parameter :: nre2 = int(2.d0 * dble(nx) / 7.d0)
  integer, parameter :: rerank = 0

  ! GPU
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(128,1,1)
  type(dim3), parameter :: threadsG  = dim3(128,1,1)
  type(dim3), parameter :: threadsEv = dim3(128,1,1)
  type(dim3), parameter :: threadsFv = dim3(128,1,1)
  type(dim3), parameter :: threadsGv = dim3(128,1,1)
  type(dim3), parameter :: threads   = dim3(128,1,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! 10step  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_recal      = 0
  integer(kind=4), parameter :: id_RungeKutta = 0
  integer, parameter         :: step_offset   = 0
  integer, parameter         :: np            = 100

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.75d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: rhol = 3.857143d0
  real(8), parameter :: ul   = 2.629369d0
  real(8), parameter :: pl   = 10.3333d0
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = CFL * Lx / (dble(nx-1) * ul)
  real(8), parameter :: endT = 0.2d0
  integer, parameter :: nt   = int(endT / (dble(np) * dt))
end module mod_globals

