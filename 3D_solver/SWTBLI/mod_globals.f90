!> SWTBLI: incident oblique shock on a turbulent boundary layer, M0 = 2.28,
!> theta = 12 deg (strong interaction, large separation), two coupled domains run
!> by src/sbli.f90:
!>   ranks 0/1 : TBL generator (recycling-rescaling), 11 delta x 4 delta
!>   ranks 2/3 : SBLI domain, inlet = TBL planes near the recycling plane,
!>               x_imp - 31 delta ... x_imp + 22 delta (+ outlet buffer), 18 delta high
!> All lengths scale with blt (= inlet delta99, the rescaling target).
!> Grid generation 2 (stages G0, C1, C2, P): dx0 = 0.0403 delta. The stages A1..B0
!> ran on generation 1 (dx0 = 0.0604 delta, 12.2 delta high, ximp_in = 25); their
!> parameters are archived in runs/<stage>/mod_globals.f90.
module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension = 3
  integer, parameter  :: sp        = kind(1.d0)
  real(sp), parameter :: threshold = 0.4_sp

  ! inlet delta99 and the composite-profile skin friction that give Re_theta = 500
  ! for this freestream (Re_unit = 5.64e6 /m); see the header of set.f90.
  ! Re_tau(delta99) = 108, wall unit = 9.95 um
  real(8), parameter :: blt     = 1.0724d-3
  real(8), parameter :: Cf_init = 3.783d-3
  real(8), parameter :: dnu0    = 9.952d-6    ! inlet viscous length (composite profile)

  ! unused since src/sbli.f90 selects the GPU node-locally; kept for its interface
  integer, parameter :: mygpu1 = 0
  integer, parameter :: mygpu2 = 1

  ! ---- spanwise (shared by both domains); nz counts 3+3 periodic ghost planes
  real(8), parameter :: Lz1 = 4.d0 * blt
  integer, parameter :: nz1 = 114
  real(8), parameter :: Lz2 = Lz1
  integer, parameter :: nz2 = nz1

  ! ---- TBL generator: uniform dx = 11 delta / 273 = 0.0403 delta (dx+ = 4.3)
  integer, parameter :: nx1 = 274
  real(8), parameter :: Lx1 = 11.d0 * blt
  real(8), parameter :: dx0 = Lx1 / dble(nx1-1)
  integer, parameter :: ny1 = 128              ! y(ny1) = 4.07 delta
  real(8), parameter :: Ly1 = 4.d0 * blt       ! nominal; the grid sets the exact value

  ! ---- SBLI domain: 1316 uniform cells (53 delta) + 80 buffer cells (+3 %/cell up
  ! to 4 dx0, 9.4 delta). x_top = x_imp - 24.2 delta for theta = 12 deg, so the top
  ! upstream of the incident shock still spans 6.8 delta; the uniform part ends
  ! xdown = 22 delta behind x_imp, at least 10 delta behind the reattachment.
  integer, parameter :: nx2   = 1397
  integer, parameter :: nbuf  = 80
  real(8), parameter :: Lx2   = 62.4d0 * blt   ! nominal
  integer, parameter :: ny2   = 455            ! y(ny2) = 18.08 delta
  real(8), parameter :: Ly2   = 18.d0 * blt    ! nominal
  real(8), parameter :: ximp_in = 31.d0        ! (x_imp - x(SBLI inlet)) / blt
  real(8), parameter :: xdown   = 22.d0        ! (end of the uniform x grid - x_imp) / blt

  ! wall-normal grid: first cell y+ = 0.8, +3.5 %/cell up to dy = dx0, uniform
  ! to y_str, then +gy_top per cell to the top. Generation 1 stretched by 5 %/cell
  ! above 9.5 delta; the reflected shock entered that coarse layer and left a fan
  ! of grid-scale waves behind it, so the uniform part now reaches 16 delta and the
  ! stretch is milder.
  real(8), parameter :: dyw0  = 0.00742d0 * blt
  real(8), parameter :: y_str = 16.d0 * blt
  real(8), parameter :: gy_top = 1.03d0

  ! set_bc_common imports these names; they describe the SBLI domain
  integer, parameter :: nx = nx2, ny = ny2, nz = nz2

  ! recycling plane at 9.95 delta (generation 1: 9.97). The reference profile Qm is the running mean over
  ! x in [9.35, 9.95] delta, i.e. AT the recycling plane: a mean over a wide window
  ! upstream (the former [6, 10] delta) makes the outer scaling use a delta99 smaller
  ! than the recycled plane's, and the inlet delta then grows ~3 % per pass (seen in
  ! the first A1 run: delta99 1.0 -> 1.24 blt at the inlet within 150 delta/U).
  ! The SBLI inlet receives TBL planes nre2-5 .. nre2-3 (calc_para flatten_rescale).
  integer, parameter :: nre1 = 233
  integer, parameter :: nre2 = 248
  integer, parameter :: nin  = nre2 - 5
  integer, parameter :: rerank = 0

  ! RTX PRO 6000 / GH200
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,2)
  type(dim3), parameter :: threadsEv = dim3(128,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,2)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! ---- time (endT, np, step_offset, theta_deg are rewritten by stage.sh)
  ! delta/u0 = 1.968 us = 182 steps
  integer, parameter :: step_offset   = 6400
  integer, parameter :: start_rescale = 10
  real(8), parameter :: endT  = 5.9034d-4
  integer, parameter :: np    = 30
  real(8), parameter :: dt    = 1.08d-8   ! acoustic CFL 0.45 on the first wall cell
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! ---- freestream
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.28d0
  real(8), parameter :: p_tot = 50.d3
  real(8), parameter :: T_tot = 290.d0
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = T_tot /  (1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)

  ! ---- physical properties
  real(8), parameter :: Pr  = 0.72d0
  real(8), parameter :: Prt = 0.9d0
  real(8), parameter :: rf  = 0.89d0
  real(8), parameter :: Taw = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)

  ! ---- incident shock: flow deflection [deg]; 0 switches the shock off.
  ! Generation 1 ramped 4 -> 7 -> 9.5; generation 2 continues 9.5 -> 11 -> 12
  ! (stage.sh). States 2/3 are derived in set.f90.
  real(8), parameter :: theta_deg = 12.d0
end module mod_globals
