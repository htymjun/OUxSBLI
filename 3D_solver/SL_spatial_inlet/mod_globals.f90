

module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.9_sp !0.9_sp
  real(8), parameter         :: blt         = 0.d0

  ! mesh  
  real(8), parameter :: delta_bl = 2.0d-3 !1.0d-3 !2.0d-3


  real(8), parameter :: Lx_main = 300d-3!600*delta_bl
  real(8), parameter :: Lx_buf  = 100d-3!50*delta_bl 
  real(8), parameter :: Lx = Lx_main + Lx_buf 
  
  real(8), parameter :: Ly_main = 80d-3!120*delta_bl
  real(8), parameter :: Ly_buf  = 40d-3!30*delta_bl  
  real(8), parameter :: Ly = Ly_main + 2.0d0*Ly_buf

  real(8), parameter :: Lz = 70d-3!35*delta_bl
  integer, parameter :: nx_main = 482!4620
  integer, parameter :: nx_buf =  150
  integer, parameter :: nx = nx_main + nx_buf + 2
  
  integer, parameter :: ny_main = 256!536 
  integer, parameter :: ny_buf =  14!24    
  integer, parameter :: ny = ny_main + 2*ny_buf + 2 

  ! x主計算領域は全域一様: 格子幅 = Lx_main/(nx_main-1)。
  ! yのみ中心一様区間と外側伸長区間に分割する（長さの単位: m）。
  ! ny_uniformは両端を含む点数。ny_main-ny_uniform は偶数。
  ! バッファの長さ・点数は主計算領域の生成とは独立。
  real(8), parameter :: Ly_uniform = 40d-3 !30*delta_bl  ! y中心を挟む一様領域の全幅 [m]
  integer, parameter :: ny_uniform = 130 !248  ! y一様領域の格子生成点数

  integer, parameter :: nz = 128!256 !290 
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(32,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,4,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,4)
  type(dim3), parameter :: threadsEv = dim3(32,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,4,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,4)
  type(dim3), parameter :: threads   = dim3(32,4,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  integer, parameter          :: step_offset =  0
  integer, parameter          :: start_rescale = 0

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.15d0

  ! initial condition
  ! real(8), parameter :: p    = 2.72d3
  ! ! M0.5
  ! real(8), parameter :: M1 = 2.0d0!11.0d0/3.0d0
  ! real(8), parameter :: T1   = 162.7908d0
  ! real(8), parameter :: rho1 = p / (R * T1)
  ! real(8), parameter :: u1   = M1 * sqrt(gamma * R * T1)
  ! ! M0.2
  ! real(8), parameter :: M2 = 3.0d0!22.0d0/15.0d0
  ! real(8), parameter :: T2   = 104.6303d0
  ! real(8), parameter :: rho2 = p / (R * T2)
  ! real(8), parameter :: u2   = M2 * sqrt(gamma * R * T2)

  ! initial condition
  real(8), parameter :: p    = 80d3
  ! M0.5
  real(8), parameter :: T1   = 290d0
  real(8), parameter :: rho1 = p / (R * T1)
  real(8), parameter :: u1   = 0.5d0 * sqrt(gamma * R * T1)
  ! M0.2
  real(8), parameter :: T2   = 290d0
  real(8), parameter :: rho2 = p / (R * T2)
  real(8), parameter :: u2   = 0.2d0 * sqrt(gamma * R * T2)
 !Mc = 1.1

  ! The initial field is unperturbed.  Velocity fluctuations are supplied
  ! continuously at the inlet by set_bc instead.
  real(8), parameter :: amp = 0.d0
  integer, parameter :: perturb_nspots = 64
  integer, parameter :: perturb_seed = 20260826
  real(8), parameter :: perturb_sigma_x = 5.d0 * delta_bl
  real(8), parameter :: perturb_sigma_z = 0.10d0 * Lz

  ! Inlet Gaussian velocity fluctuations, localized around the shear layer:
  ! q' = 0.05*abs(U_mean(y))*E(y)*N(0,1), independently for u, v and w,
  ! where E=sech^2(2*(y-Ly/2)/delta_bl).  The centerline RMS remains 5%.
  real(8), parameter :: inlet_fluctuation_rms = 0.05d0
  integer, parameter :: inlet_random_seed = 20260827
  integer, parameter :: inlet_rk_stages = 4
  real(8), parameter :: CFL  = 0.1d0 !0.05d0
  real(8), parameter :: dt =  9.20460E-009!.3727344786232972E-008!3.9757331116907360E-009 != 0.05d0 * (Lx / dble(nx-1)) / (u1 + 340.d0) !CFL * Lx / (dble(nx-1) * abs(u1))
  real(8), parameter :: endT = 10.d0 * Lx / u2
  integer, parameter :: np = 1000 !500 !1000 !2000 !100
  integer, parameter :: nt = 1000!1000 !int(endT / (dble(np) * dt))

end module mod_globals
