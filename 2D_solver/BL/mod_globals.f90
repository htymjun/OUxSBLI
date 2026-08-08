module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 2
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  real(8), parameter  :: blt         = 1.d-3

  ! mesh
  real(8), parameter :: Lx = 120.d0 * blt
  ! Ly must keep the top boundary out of the boundary layer: 15 mm is ~7.2 * delta99
  ! at the trailing edge, where delta* / Ly = 4.7%. Cheaper than 30 mm and still
  ! generous now that the top/outlet impose p = p0 directly instead of coupling
  ! pressure to v -- that coupling (now removed) was the reason a tall domain used
  ! to matter.
  real(8), parameter :: Ly = 15.d0 * blt
  ! ny is halved along with Ly, so with the s = 2.4 tanh stretch in set_grid this
  ! keeps dy_wall ~= 26 um (was 25.3 um at 30mm/97) -- near-wall resolution, where
  ! Cf is measured, is essentially unchanged. Point count halves: 24,929 -> 12,593.
  integer, parameter :: nx = 257
  integer, parameter :: ny = 49

  ! flat-plate geometry
  integer, parameter :: i_LE = 2 * nx / 12 + 1 ! first no-slip wall point; set_grid places the leading edge at x = 0

  ! GPU thread blocks
  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsEv = dim3(64,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksEv, blocksFv, blocks

  ! time
  integer, parameter :: step_offset = 0
  ! 2.84 flow-through times (Lx / u0 = 3.526 ms); the production run's own
  ! snapshot history showed full convergence by 2.27 flow-throughs (t=8ms of a
  ! 20ms/5.67-FT run) with no change out to 5.67, so this keeps margin at half
  ! the step count. This is also what ouxsbli/tests/test_bl.py already uses.
  real(8), parameter :: endT  = 1.d-2
  integer, parameter :: np    = 100
  real(8), parameter :: R     = 287.15d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 0.1
  real(8), parameter :: p_tot = 25.d3
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = 288.15d0
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  ! acoustic CFL_y = (u0+c0)*dt/dy_wall ~= 0.58 at the new dy_wall (~26 um).
  ! Short probes (0.4ms transient, 257x49/Ly=15mm) ran clean with no NaN and no
  ! leading-edge ringing all the way through dt=8e-8 (CFL~1.16); this keeps a
  ! ~2x margin below the highest value actually tested rather than running at
  ! the edge. SLAU's upwind dissipation is what buys this margin -- KEEP does
  ! not have it and diverges on this case regardless of dt (see config.fypp).
  real(8), parameter :: dt    = 4.d-8
  integer, parameter :: nt    = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = dsqrt(Pr)
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
end module mod_globals
