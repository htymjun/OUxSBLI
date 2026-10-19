

module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.9_sp !0.9_sp
  real(8), parameter         :: blt         = 0.d0

  ! mesh  
  ! delta_bl = 初期渦度厚さ delta_omega0 [m]。tanh(2(y-Ly/2)/delta_bl) の速度分布の最大勾配 = dU/delta_bl。
  real(8), parameter :: delta_bl = 2.0d-3


  real(8), parameter :: Lx_main = 400*delta_bl
  real(8), parameter :: Lx_buf  = 50*delta_bl
  real(8), parameter :: Lx = Lx_main + Lx_buf 
  
  real(8), parameter :: Ly_main = 80*delta_bl          ! y主計算領域 = 100 delta_omega0
  real(8), parameter :: Ly_buf  = 30*delta_bl
  real(8), parameter :: Ly = Ly_main + 2.0d0*Ly_buf

  real(8), parameter :: Lz = 70*delta_bl!40*delta_bl

  !格子のパラメータ
  integer, parameter :: nx_main = 1720
  integer, parameter :: nx_buf =  60
  ! x方向MPI分割: npx = xスラブ数(= 計算ランク数 = MPIランク数/2)。npx=1 は分割なし。
  ! nx はスラブ1枚の点数(接続側にゴースト3面)。(nx_global-6) が npx で割り切れること。
  ! 実行は mpiexec -n 2*npx(set.f90 の set_grid が検査する)。config.fypp は COMMZ=False。
  integer, parameter :: npx = 1
  integer, parameter :: nx_global = nx_main + nx_buf + 2
  integer, parameter :: nx = (nx_global - 6)/npx + 6
  
  integer, parameter :: ny_main = 351!293
  integer, parameter :: ny_buf =  24!14!24    
  integer, parameter :: ny = ny_main + 2*ny_buf + 2 

  ! x主計算領域: 流入端の格子幅 dx_inlet から始めて、下流に向かって格子幅が5次smoothstepで滑らかに
  ! 粗くなる(一様区間なし。流入端の傾きは0)。終端の格子幅 = dx_inlet + 振幅で、
  ! 振幅 = (Lx_main - (nx_main-1)*dx_inlet) / ((nx_main-1)/2)。バッファはその幅から続けて伸長する。
  ! Lx_main を (nx_main-1)*dx_inlet より長くすること(現在の終端幅は約 1.0 delta_bl = 流入端の4倍)。
  real(8), parameter :: dx_inlet = 0.14d0*delta_bl
  ! y主計算領域は sinh 写像: y = (Ly_main/2)*(1 + sinh(beta_y*s)/sinh(beta_y)), s=-1..1 を ny_main 点で等分。
  ! 層の中心が最も細かく、外側へ滑らかに粗くなる(set.f90 の sinh_axis)。beta_y が大きいほど中心に集中する。
  !   中心の格子幅 = Ly_main*beta_y/((ny_main-1)*sinh(beta_y))、端の格子幅 = 中心の格子幅*cosh(beta_y)
  ! バッファの長さ・点数は主計算領域の生成とは独立(主領域の端の格子幅から続けて伸長する)。
  real(8), parameter :: beta_y = 3.0d0

  ! nz は全体のz点数(周期、ゴースト3+3面を含む)。dz = Lz/(nz-6)。
  ! (z方向2分割のときの nz=32 [内部26面x2スラブ] と同じ格子が nz=58)
  integer, parameter :: nz = 342!174!100 !400
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
  real(8), parameter :: Sct   = 0.9d0   ! LES: 乱流シュミット数 (スカラーの乱流拡散 = mut/Sct)
  real(8), parameter :: R     = 287.15d0

  ! passive scalar (mixture fraction xi: 1 in stream 1, 0 in stream 2), solved in set.f90
  logical, parameter :: scalar_on = .false.!.true.          ! .false. -> no scalar, no extra cost
  real(8), parameter :: Sc        = 0.7d0           ! Schmidt number: rho*D = mu/Sc
  integer, parameter :: scalar_output_every = 1     ! write xi every this many flow outputs (nt steps each)

  ! Stop the whole run (MPI_ABORT) when NaN/Inf shows up in the conserved variables or in xi.
  ! Checked every nan_check_every steps (and at every stage of the first two steps); 0 disables the check.
  ! Cost: about 1.4% of the run per check per step, i.e. 0.14% at 10.
  integer, parameter :: nan_check_every = 10

  real(8), parameter :: p    = 80.0d3
  real(8), parameter :: M1   = 0.6667d0
  real(8), parameter :: T1   = 290d0
  real(8), parameter :: rho1 = p / (R * T1)
  real(8), parameter :: u1   = M1 * sqrt(gamma * R * T1)
  real(8), parameter :: M2   = 0.2667d0
  real(8), parameter :: T2   = 290d0
  real(8), parameter :: rho2 = p / (R * T2)
  real(8), parameter :: u2   = M2 * sqrt(gamma * R * T2)

  ! The initial field is unperturbed.  Velocity fluctuations are supplied
  ! continuously at the inlet by set_bc instead.
  real(8), parameter :: amp = 0.d0

  ! Inlet Gaussian velocity fluctuations, localized around the shear layer:
  ! q' = 0.05*abs(U_mean(y))*E(y)*N(0,1), independently for u, v and w,
  ! where E=sech^2(2*(y-Ly/2)/delta_bl).  The centerline RMS remains 5%.
  real(8), parameter :: inlet_fluctuation_rms = 0.05d0
  integer, parameter :: inlet_random_seed = 20260827
  integer, parameter :: inlet_rk_stages = 4
  real(8), parameter :: CFL  = 0.1d0 !0.05d0
  ! dt は set_grid が表示する Suggested dt (対流CFL=0.1, 3方向合算) の少し下。Mc=0.3: 6.0d-10, Mc=0.8: 3.2d-10
  real(8), parameter :: dt = 1.0E-008!2.5E-08!3.0d-10
  real(8), parameter :: endT = 10.d0 * Lx / u2
  ! 出力回数 np、出力間隔 nt ステップ。np*nt*dt が endT (= 高速側の流入速度で10回通過する時間) になる。
  integer, parameter :: np = 1000!200
  integer, parameter :: nt = 1000!int(endT / (dble(np) * dt))

end module mod_globals
