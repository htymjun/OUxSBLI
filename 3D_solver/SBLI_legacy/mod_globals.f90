module mod_globals
  use cudafor
  implicit none
  integer, parameter    :: dimension = 3
  integer, parameter    :: accuracy  = 2
  integer, parameter    :: offset    = accuracy / 2
  integer(4), parameter :: id_visc   = 2
  integer, parameter    :: id_turbulence = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_visc     ! kind2 Euler         !
  !             ! kind4 NS            !
  !             ! kind8 LES           !
  !             ! 1 2nd               !
  !             ! 2 4th               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_turbulence ! 0 laminar         !
  !               ! 1 SMS             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_scheme   ! integer(2)  KEEP    !
  !             ! real(2)     SLAU    !
  !             ! real(4)   Weighted  !
  !             ! real(8)   Threshold !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_accuracy ! kind2 2nd           !
  !             ! kind4 4th           !
  !             ! kind8 6th           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_tvd      ! kind2 non TVD       !
  !             ! kind4 minmod        !
  !             ! kind8 switch        !
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
  real(2), parameter         :: id_scheme   = 0
  real(8), parameter         :: threshold   = 0.4d0
  integer(kind=8), parameter :: id_accuracy = 0
  integer(kind=8), parameter :: id_tvd      = 0
  integer(kind=2), parameter :: id_keep     = 0
  integer(kind=4), parameter :: id_slau     = 0
  integer(kind=2), parameter :: slau_wall   = 0
  integer(kind=4), parameter :: id_rescale  = 0
  real(8), parameter         :: blt         = 2.d-3

  ! mesh
  real(8), parameter :: Lx = 6.25d-3  ! * 8 25 delta
  real(8), parameter :: Ly = 8d-3     !   4 delta
  real(8), parameter :: Lz = 32d-3    !  16 delta
  integer, parameter :: nx = 161      ! * 8
  integer, parameter :: ny = 321
  integer, parameter :: nz = 2049
  integer, parameter :: rerank = 4
  integer, parameter :: nre1 = int(0.24d0 * nx) ! 4 delta for Lz = delta, 2 delta for Lz = 4 delta, 1 delta for Lz = 16 delta
  integer, parameter :: nre2 = int(0.56d0 * nx)
  integer, parameter :: overlap = 3

  ! RTX 4090
  type(dim3) :: blocksE   = dim3((nx-accuracy+1)/32,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: blocksF   = dim3((nx-accuracy)/1,(ny-accuracy+1)/64,(nz-accuracy)/1)
  type(dim3) :: blocksG   = dim3((nx-accuracy)/1,(ny-accuracy)/1,(nz-accuracy+1)/64)
  type(dim3) :: blocksEv  = dim3((nx-accuracy+1)/32,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: blocksFv  = dim3((nx-accuracy)/1,(ny-accuracy+1)/64,(nz-accuracy)/1)
  type(dim3) :: blocksGv  = dim3((nx-accuracy)/1,(ny-accuracy)/1,(nz-accuracy+1)/64)
  type(dim3) :: blocks    = dim3((nx-accuracy)/159,(ny-accuracy)/1,(nz-accuracy)/1)
  type(dim3) :: threadsE  = dim3(32,1,1)
  type(dim3) :: threadsF  = dim3(1,64,1)
  type(dim3) :: threadsG  = dim3(1,1,64)
  type(dim3) :: threadsEv = dim3(32,1,1)
  type(dim3) :: threadsFv = dim3(1,64,1)
  type(dim3) :: threadsGv = dim3(1,1,64)
  type(dim3) :: threads   = dim3(159,1,1)

  ! time
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_RungeKutta ! kind2 ! 3rd_TVD !
  !               ! kind4 ! 4th     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! id_recal      ! kind2 ! set 0   !
  !               ! kind4 ! recal   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=2), parameter :: id_RungeKutta = 0
  integer(kind=4), parameter :: id_recal      = 0
  integer, parameter         :: step_offset   = 118
  integer, parameter         :: start_rescale = 3 ! this parameter must be greater than 3
  real(8), parameter :: endT = 0.005d-3
  integer, parameter :: np   = 1
  real(8), parameter :: u0   = 506.8d0
  real(8), parameter :: CFL  = 0.1d0
  real(8), parameter :: dt   = 8d-9!7.5d-9!CFL * Lx / (dble(nx-1) * u0) ! 7.7e-9
  integer, parameter :: nt   = int(endT / (dble(np) * dt))

  ! physical properties
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: Pr    = 0.71d0
  real(8), parameter :: Prt   = 0.9d0
  real(8), parameter :: R     = 287.03d0

  ! initial condition
  real(8), parameter :: M0    = 1.9d0
  real(8), parameter :: p0    = 14924.d0
  real(8), parameter :: T0    = 171.31d0
  real(8), parameter :: beta  = dacos(-1.d0) * 39.27d0 / 180.d0
  real(8), parameter :: Ms    = M0 * dsin(beta)
  real(8), parameter :: Ms2   = Ms**2
  real(8), parameter :: theta = datan(2.d0 * (1.d0 / dtan(beta)) * (Ms2 - 1.d0) / (M0**2 * (gamma + dcos(2.d0 * beta)) + 2.d0))
  real(8), parameter :: T2    = T0 * (1.d0 + 2.d0 * (gamma - 1.d0) * (Ms2 - 1.d0) * (1.d0 + gamma * Ms2) / (Ms2 * (gamma + 1.d0)**2))
  real(8), parameter :: p2    = p0 * (1.d0 + 2.d0 * gamma * (Ms2 - 1.d0) / (gamma + 1.d0))
  real(8), parameter :: rho2  = p2 / (R * T2)
  real(8), parameter :: u1    = u0 * dsin(beta)
  real(8), parameter :: v1    = u0 * dcos(beta)
  real(8), parameter :: a1    = u0 / M0
  real(8), parameter :: u2    = u1 - 2.d0 * a1 * (Ms - 1.d0 / Ms) / (gamma + 1.d0)
  real(8), parameter :: v2    = u0 * dcos(beta)
  real(8), parameter :: u_magnitude = sqrt(u2**2 + v2**2)
  real(8), parameter :: ux    = u_magnitude * dcos(theta)
  real(8), parameter :: uy    = - u_magnitude * dsin(theta)
end module mod_globals

