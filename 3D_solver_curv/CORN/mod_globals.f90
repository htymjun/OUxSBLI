module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension   = 3
  integer, parameter    :: sp          = 4
  real(sp), parameter   :: threshold   = 0.1_sp
  
  ! mesh
  real(8), parameter :: Lx = 2.0d0
  real(8), parameter :: Ly = 1.0d0
  real(8), parameter :: Lz = 0.1d0
  integer, parameter :: nx = 192
  integer, parameter :: ny = 80
  integer, parameter :: nz = 4
  ! Compression corner geometry
  real(8), parameter :: theta    = 8.d0 * acos(-1.d0) / 180.d0 ! 8° in radians
  real(8), parameter :: x_corner = 0.5d0                       ! corner location


  ! GPU thread block dimensions (tuned for nx≈192, ny≈80)
  type(dim3), parameter :: threads   = dim3(32, 8, 1)
  type(dim3), parameter :: threadsE  = dim3(32, 4, 2)
  type(dim3), parameter :: threadsFv = dim3(8,  4, 4)
  type(dim3), parameter :: threadsF  = dim3(8,  16, 2)
  type(dim3), parameter :: threadsEv = dim3(16, 4, 2)
  type(dim3), parameter :: threadsG  = dim3(8,  8, 4)
  type(dim3), parameter :: threadsGv = dim3(8,  8, 4)
  type(dim3) :: blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv

  ! time
  integer, parameter    :: step_offset   = 0

  ! Free-stream flow conditions (non-dimensional)
  real(8), parameter :: Ma_inf  = 2.d0    ! Mach number (supersonic)
  real(8), parameter :: gamma   = 1.4d0   ! heat capacity ratio
  real(8), parameter :: R       = 1.d0    ! gas constant (non-dimensional)
  real(8), parameter :: Pr      = 0.72d0  ! Prandtl number
  real(8), parameter :: Prt     = 0.9d0    ! turbulent Prandtl number (LES branch of calc_visc2_curv)
    
  ! Non-dimensional free-stream
  real(8), parameter :: rho_inf = 1.d0
  real(8), parameter :: u_inf   = Ma_inf
  real(8), parameter :: v_inf   = 0.d0
  real(8), parameter :: p_inf   = 1.d0 / gamma

  ! Time stepping parameters
  real(8), parameter :: dt = 5.d-5
  ! Derived time parameters
  integer, parameter :: nt = 10000
  integer, parameter :: np = 10
  integer, parameter :: rerank = -1    ! no re-scaling by default
end module mod_globals
