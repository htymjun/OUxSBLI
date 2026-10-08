module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension   = 3
  integer, parameter    :: sp          = 4
  real(sp), parameter   :: threshold   = 0.1_sp
  
  ! NACA 0012 O-grid geometry
  real(8), parameter :: chord = 0.05d0
  real(8), parameter :: aoa   = -acos(-1.d0) * 5.d0 / 180.d0 ! angle of attack [radians]
  real(8), parameter :: far_r = 8.d0 * chord                 ! far-field radius [chords]

  ! mesh
  real(8), parameter :: AR = 0.01d0
  real(8), parameter :: Lx = 1.d0       ! dummy (main_curv passes it to set_grid; set_grid ignores it)
  real(8), parameter :: Ly = 1.d0       ! dummy (= far_r, for consistency)
  real(8), parameter :: Lz = AR * chord ! spanwise
  integer, parameter :: nx = 513        ! interior cells (i=2..nx-1) wrapping airfoil
  integer, parameter :: ny = 161        ! wall-normal cells
  integer, parameter :: nz = 17         ! quasi-2D spanwise

  real(8), parameter :: beta  = dacos(-1.d0) * 37.2d0 / 180.d0

  ! GPU thread block dimensions (tuned for nx≈194, ny≈80)
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
  real(8), parameter :: gamma   = 1.4d0    ! heat capacity ratio
  real(8), parameter :: R       = 287.03d0 ! gas constant (non-dimensional)
  real(8), parameter :: Pr      = 0.72d0   ! Prandtl number
  real(8), parameter :: Prt     = 0.9d0    ! turbulent Prandtl number (LES branch of calc_visc2_curv)
  real(8), parameter :: Ma_inf  = 0.8d0    ! Mach number (subsonic)
  real(8), parameter :: p_inf   = 101.3d3
  real(8), parameter :: T_inf   = 288.15d0
  real(8), parameter :: rho_inf = p_inf / (R * T_inf)

  ! Time stepping parameters
  ! dt=1e-3: CFL estimate with stretch=1.05, ny=80, far_r=10 gives Δη_wall≈0.011,
  ! (|u|+a)=1.5 → dt_max≈7e-3; dt=1e-3 is safely below that.
  real(8), parameter :: dt = 2.d-9
  integer, parameter :: nt = 1000
  integer, parameter :: np = 100
  integer, parameter :: rerank = -1
end module mod_globals
