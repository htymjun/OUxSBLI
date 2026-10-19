

module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.9_sp !0.9_sp
  real(8), parameter         :: blt         = 0.d0

  real(8), parameter :: delta_bl = 2.0d-3
  ! mesh
  real(8), parameter :: Lx_main = 400*delta_bl
  real(8), parameter :: Lx_buf  = 50*delta_bl !110.d-3 !100.d-3
  real(8), parameter :: Lx = Lx_main + Lx_buf !200.d-3
  
  real(8), parameter :: Ly_main = 80*delta_bl
  real(8), parameter :: Ly_buf  = 30*delta_bl
  real(8), parameter :: Ly = Ly_main + 2.0d0*Ly_buf!40.d-3

  ! z uses the same center-clustered mesh and stretched buffers as y.
  real(8), parameter :: Lz_main = 80*delta_bl
  real(8), parameter :: Lz_buf = Ly_buf
  real(8), parameter :: Lz = Lz_main + 2.d0*Lz_buf


  ! x主計算領域: 流入端の格子幅 dx_inlet から始めて、下流に向かって格子幅が5次smoothstepで滑らかに
  ! 粗くなる(一様区間なし。流入端の傾きは0)。終端の格子幅 = dx_inlet + 振幅で、
  ! 振幅 = (Lx_main - (nx_main-1)*dx_inlet) / ((nx_main-1)/2)。バッファはその幅から続けて伸長する。
  ! Lx_main を (nx_main-1)*dx_inlet より長くすること(現在の終端幅は約 2.4 mm = 流入端の4倍)。
  real(8), parameter :: dx_inlet = 0.14d0*delta_bl
  integer, parameter :: nx_main = 1720
  integer, parameter :: nx_buf =  60!20  !30
  ! x方向MPI分割: npx = xスラブ数(= 計算ランク数 = MPIランク数/2)。npx=1 は分割なし。
  ! nx はスラブ1枚の点数(接続側にゴースト3面)。(nx_global-6) が npx で割り切れること。
  ! 実行は mpiexec -n 2*npx(set.f90 の set_grid が検査する)。
  integer, parameter :: npx = 1
  integer, parameter :: nx_global = nx_main + nx_buf + 2!438
  integer, parameter :: nx = (nx_global - 6)/npx + 6
  
  integer, parameter :: ny_main = 351 !256 !192
  integer, parameter :: ny_buf =  24     !15
  integer, parameter :: ny = ny_main + 2*ny_buf + 2 !130

 integer, parameter :: nz_main = 351
  integer, parameter :: nz_buf = ny_buf
  integer, parameter :: nz = nz_main + 2*nz_buf + 2
  integer, parameter :: rerank = 0

  ! y/z主計算領域は sinh 写像: y = (Ly_main/2)*(1 + sinh(beta_y*s)/sinh(beta_y)), s=-1..1 を ny_main 点で等分。
  ! 層の中心が最も細かく、外側へ滑らかに粗くなる(set.f90 の sinh_axis)。beta が大きいほど中心に集中する。
  !   中心の格子幅 = Ly_main*beta_y/((ny_main-1)*sinh(beta_y))、端の格子幅 = 中心の格子幅*cosh(beta_y)
  real(8), parameter :: beta_y = 3.0d0
  real(8), parameter :: beta_z = beta_y

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
  real(8), parameter :: Sct   = 0.9d0   ! LES: 乱流シュミット数 (スカラーの乱流拡散 = mut/Sct)
  real(8), parameter :: R     = 287.15d0

  ! passive scalar (mixture fraction xi: 1 in stream 1 [u1], 0 in stream 2 [u2]), solved in set.f90
  logical, parameter :: scalar_on = .false.          ! .false. -> no scalar, no extra cost
  real(8), parameter :: Sc        = 0.7d0           ! Schmidt number: rho*D = mu/Sc
  integer, parameter :: scalar_output_every = 1     ! write xi every this many flow outputs (nt steps each)

  ! Stop the whole run (MPI_ABORT) when NaN/Inf shows up in the conserved variables or in xi.
  ! Checked every nan_check_every steps (and at every stage of the first two steps); 0 disables the check.
  ! Cost: about 1.4% of the run per check per step, i.e. 0.14% at 10.
  integer, parameter :: nan_check_every = 10

  ! initial condition
  ! real(8), parameter :: p    = 80d3
  ! ! M2
  ! real(8), parameter :: M1 = 1.0d0
  ! real(8), parameter :: T1   = 290
  ! real(8), parameter :: rho1 = p / (R * T1)
  ! real(8), parameter :: u1   = M1 * sqrt(gamma * R * T1)
  ! ! M0.2
  ! real(8), parameter :: M2 = 0.4d0
  ! real(8), parameter :: T2   = 290
  ! real(8), parameter :: rho2 = p / (R * T2)
  ! real(8), parameter :: u2   = M2 * sqrt(gamma * R * T2)

  real(8), parameter :: p    = 80.0d3
  real(8), parameter :: M1   = 0.6667d0
  real(8), parameter :: T1   = 290d0
  real(8), parameter :: rho1 = p / (R * T1)
  real(8), parameter :: u1   = M1 * sqrt(gamma * R * T1)
  real(8), parameter :: M2   = 0.2667d0
  real(8), parameter :: T2   = 290d0
  real(8), parameter :: rho2 = p / (R * T2)
  real(8), parameter :: u2   = M2 * sqrt(gamma * R * T2)
 !Mc = 1.1
  !4.0d0*0.988d-3 !1.0d-3 !2.0d-3

  ! The initial field is unperturbed.  Velocity fluctuations are supplied
  ! continuously at the inlet by set_bc instead.
  real(8), parameter :: amp = 0.d0
  integer, parameter :: perturb_nspots = 64
  integer, parameter :: perturb_seed = 20260826
  real(8), parameter :: perturb_sigma_x = 5.d0 * delta_bl
  real(8), parameter :: perturb_sigma_z = 0.10d0 * Lz

  ! Inlet Gaussian velocity fluctuations, localized around the shear layer:
  ! q' = inlet_fluctuation_rms*abs(U_mean(y,z))*E(y,z)*N(0,1).
  ! E = 1 - tanh_y**2*tanh_z**2 localizes forcing to both crossing layers.
  real(8), parameter :: inlet_fluctuation_rms = 0.05d0
  integer, parameter :: inlet_random_seed = 20260827
  integer, parameter :: inlet_rk_stages = 4
  ! Non-uniform z on the shared solver: set_bc rescales every RK4 stage update
  ! of plane k by dz(1)/hz(k) so that it becomes the finite-volume update with
  ! the local cell width (see set.f90).  Requires RK=4 and COMMZ=False.
  logical, parameter :: stretched_z_correction = .true.
  real(8), parameter :: CFL  = 0.1d0 !0.05d0
  real(8), parameter :: dt =  1.0E-008!.3727344786232972E-008!3.9757331116907360E-009 != 0.05d0 * (Lx / dble(nx-1)) / (u1 + 340.d0) !CFL * Lx / (dble(nx-1) * abs(u1))
  real(8), parameter :: endT = 10.d0 * Lx / u2
  integer, parameter :: np = 1000!500 !1000 !2000 !100
  integer, parameter :: nt = 1000!1000 !int(endT / (dble(np) * dt))

end module mod_globals
