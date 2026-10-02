

module mod_globals
  use cudafor
  implicit none
  integer, parameter         :: dimension   = 3
  integer, parameter         :: sp          = kind(1.d0) ! single or double
  real(sp), parameter        :: threshold   = 0.9_sp !0.9_sp
  real(8), parameter         :: blt         = 0.d0

  
  ! mesh
  real(8), parameter :: Lx_main = 300.d-3!200.d-3 !150.d-3 !200.d-3
  real(8), parameter :: Lx_buf  = 100.d-3 !110.d-3 !100.d-3
  real(8), parameter :: Lx = Lx_main + Lx_buf !200.d-3
  
  real(8), parameter :: Ly_main = 40.d-3!80.d-3 !80
  real(8), parameter :: Ly_buf  = 20.d-3 !80 
  real(8), parameter :: Ly = Ly_main + 2.0d0*Ly_buf!40.d-3

  ! z uses the same center-clustered mesh and stretched buffers as y.
  real(8), parameter :: Lz_main = 40.d-3
  real(8), parameter :: Lz_buf = Ly_buf
  real(8), parameter :: Lz = Lz_main + 2.d0*Lz_buf


  integer, parameter :: nx_main = 492!428*2 !384
  integer, parameter :: nx_buf =  30!20  !30
  integer, parameter :: nx = nx_main + nx_buf + 2!438
  
  integer, parameter :: ny_main = 128 !256 !192
  integer, parameter :: ny_buf =  14     !15
  integer, parameter :: ny = ny_main + 2*ny_buf + 2 !130

 integer, parameter :: nz_main = 128
  integer, parameter :: nz_buf = ny_buf
  integer, parameter :: nz = nz_main + 2*nz_buf + 2
  integer, parameter :: rerank = 0

  ! y/z main domain: uniform central section + stretched outer section.
! *_uniform は両端を含む点数。ny_main-ny_uniform, nz_main-nz_uniform は偶数にすること。
real(8), parameter :: Ly_uniform = 20.d-3   ! y中心の一様区間の全幅 [m]
integer, parameter :: ny_uniform = 80       ! y一様区間の格子生成点数
real(8), parameter :: Lz_uniform = Ly_uniform
integer, parameter :: nz_uniform = ny_uniform

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
  real(8), parameter :: p    = 80d3
  ! M2
  real(8), parameter :: M1 = 1.0d0
  real(8), parameter :: T1   = 290
  real(8), parameter :: rho1 = p / (R * T1)
  real(8), parameter :: u1   = M1 * sqrt(gamma * R * T1)
  ! M0.2
  real(8), parameter :: M2 = 0.4d0
  real(8), parameter :: T2   = 290
  real(8), parameter :: rho2 = p / (R * T2)
  real(8), parameter :: u2   = M2 * sqrt(gamma * R * T2)

 !Mc = 1.1
  real(8), parameter :: delta_bl = 2.0d-3!4.0d0*0.988d-3 !1.0d-3 !2.0d-3

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
  real(8), parameter :: dt =  8.20460E-009!.3727344786232972E-008!3.9757331116907360E-009 != 0.05d0 * (Lx / dble(nx-1)) / (u1 + 340.d0) !CFL * Lx / (dble(nx-1) * abs(u1))
  real(8), parameter :: endT = 10.d0 * Lx / u2
  integer, parameter :: np = 1000!500 !1000 !2000 !100
  integer, parameter :: nt = 1000!1000 !int(endT / (dble(np) * dt))

end module mod_globals
