module mod_globals
  use cudafor
  implicit none
  integer, parameter  :: dimension   = 3
  integer, parameter  :: sp          = kind(1.d0) ! single or double
  real(sp), parameter :: threshold   = 0.4_sp
  ! target delta99 at the rescaling inlet. Re_theta = 751 * delta99[mm] for this
  ! freestream, so this sets Re_theta ~ 1430 at the inlet; it grows at 13.9 /mm,
  ! crossing Guarini et al. (2000)'s Re_theta = 1577 at x ~ 0.63 Lx
  real(8), parameter  :: blt         = 1.9d-3

  ! mesh: 2271+ x 858+ x 606+ against the reference's 2269+ x 875+ x 1134+
  real(8), parameter :: Lx = 9.0d0 * blt
  real(8), parameter :: Ly = 3.4d0 * blt
  real(8), parameter :: Lz = 2.5d0 * blt
  ! DNS: dx+ = 8.9, dz+ = 4.9, dy_min+ = 0.60 (reference 8.9, 5.9, 0.48)
  integer, parameter :: nx = 257
  integer, parameter :: ny = 193
  integer, parameter :: nz = 129

  integer, parameter :: nre1 = int(0.5 * nx)
  integer, parameter :: nre2 = int(0.9 * nx)
  integer, parameter :: rerank = 0

  type(dim3), parameter :: threadsE  = dim3(128,1,1)
  type(dim3), parameter :: threadsF  = dim3(32,2,1)
  type(dim3), parameter :: threadsG  = dim3(32,1,2)
  type(dim3), parameter :: threadsEv = dim3(128,1,1)
  type(dim3), parameter :: threadsFv = dim3(32,2,1)
  type(dim3), parameter :: threadsGv = dim3(32,1,2)
  type(dim3), parameter :: threads   = dim3(32,2,1)
  type(dim3) :: blocksE, blocksF, blocksG, blocksEv, blocksFv, blocksGv, blocks

  ! time
  integer, parameter :: step_offset   = 0
  integer, parameter :: start_rescale = 10
  ! one flow-through Lx/u0 = 2.98d-5 s = 4967 steps
  ! stage A (spin-up)  : endT = 3.6d-4, np =  10, RESTART=False -> 12.1 FT
  ! stage B (sampling) : endT = 7.2d-4, np = 100, RESTART=True  -> 24.2 FT
  real(8), parameter :: endT  = 3.6d-4
  integer, parameter :: np    = 10
  real(8), parameter :: R     = 287.03d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.5d0
  real(8), parameter :: p_tot = 100.d3
  real(8), parameter :: T_tot = 295.d0
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = T_tot /  (1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)
  real(8), parameter :: dt    = 6.d-9 ! acoustic CFL 0.44 on the first wall cell
  integer, parameter :: nt    = int(endT / (dble(np) * dt))
  ! physical properties
  real(8), parameter :: Pr    = 0.72d0
  real(8), parameter :: Prt   = 0.9d0
  ! wall temperature
  real(8), parameter :: rf    = 0.89d0
  real(8), parameter :: Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
end module mod_globals
