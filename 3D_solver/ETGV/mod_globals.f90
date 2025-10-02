module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension = 3
  integer, parameter    :: accuracy  = 2 
  integer, parameter    :: offset    = accuracy / 2
  integer(2), parameter :: id_visc   = 0
  integer(2), parameter :: id_LL     = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc       ! kind2 Euler       !
  !               ! kind4 NS          !
  !               ! kind8 LES         !
  !               ! 1 2nd             !
  !               ! 2 4th             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_av         ! 0 no              !
  !               ! 1 Neumann         !
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
  ! id_keep     ! kind2 KEEP          !
  !             ! kind4 KEEPPE        !
  !             ! kind8 KEP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_slau     ! kind2 SLAU          !
  !             ! kind4 HR-SLAU2      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! slau_wall   ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_rescale  ! kind2 off           !
  !             ! kind4 on            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(2), parameter      :: id_scheme   = 0
  integer, parameter         :: id_sensor   = 1
  real(8), parameter         :: threshold   = 0.9d0
  integer(kind=2), parameter :: id_accuracy = 0
  integer(kind=2), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: slau_wall   = 0
  integer(kind=2), parameter :: id_rescale  = 0
  real(8), parameter         :: blt         = 0.d0

  ! exchange
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_exchange ! kind2 off !
  !             ! kind4 on  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_exchange = 0
  integer(kind=2), parameter :: id_forcing  = 0
  ! mesh
  real(8), parameter :: pi = acos(-1.d0)
  real(8), parameter :: Lx = 2.d0 * pi
  real(8), parameter :: Ly = 2.d0 * pi
  real(8), parameter :: Lz = 2.d0 * pi
  ! 4th-order accuracy
  integer, parameter :: nx = 66!258!130!66
  integer, parameter :: ny = 66!258!130!66
  integer, parameter :: nz = 66!258!130!66
  ! 6th-order accuracy
  !integer, parameter :: nx = 70
  !integer, parameter :: ny = 70
  !integer, parameter :: nz = 70

  integer, parameter :: nre1 = 1
  integer, parameter :: nre2 = nx
  integer, parameter :: rerank = 0

  ! GPU
  ! 4th-order accuracy
  type(dim3), parameter :: threadsE  = dim3(5,8,8)
  type(dim3), parameter :: threadsF  = dim3(8,5,8)
  type(dim3), parameter :: threadsG  = dim3(8,8,5)
  type(dim3), parameter :: threadsEv = dim3(5,8,8)
  type(dim3), parameter :: threadsFv = dim3(8,5,8)
  type(dim3), parameter :: threadsGv = dim3(8,8,5)
  type(dim3), parameter :: threads   = dim3(8,8,8)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks
  ! 6th-order accuracy
  !type(dim3) :: blocksE = dim3((nx-accuracy+1)/23,(ny-accuracy)/4,(nz-accuracy)/4)
  !type(dim3) :: blocksF = dim3((nx-accuracy)/4,(ny-accuracy+1)/23,(nz-accuracy)/4)
  !type(dim3) :: blocksG = dim3((nx-accuracy)/4,(ny-accuracy)/4,(nz-accuracy+1)/23)
  !type(dim3) :: blocks  = dim3((nx-accuracy)/4,(ny-accuracy)/4,(nz-accuracy)/4)
  !type(dim3) :: threadsE = dim3(23,4,4)
  !type(dim3) :: threadsF = dim3(4,23,4)
  !type(dim3) :: threadsG = dim3(4,4,23)
  !type(dim3) :: threads  = dim3(4,4,4)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind=2 ! 3rd_TVD !
  !               ! kind=4 ! 4th     !
  !               ! kind=8 ! Gauss   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind=2 ! set 0   !
  !               ! kind=4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter  :: id_recal = 0
  integer(kind=4), parameter  :: id_RungeKutta = 0
  integer, parameter          :: step_offset = 0
  integer, parameter          :: start_rescale = 0
  integer, parameter          :: nt = 200!1
  integer, parameter          :: np = 50!200
  real(8), parameter          :: dt = 0.01d0!0.02d0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr = 0.71d0
  real(8), parameter :: Prt = 0.9d0
  real(8), parameter :: R = 287.03d0

  ! initial condition
  real(8), parameter :: M0 = 0.4d0
  real(8), parameter :: RHO0 = 1.d0
end module mod_globals

