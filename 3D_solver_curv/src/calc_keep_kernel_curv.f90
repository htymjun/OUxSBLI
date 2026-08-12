!> Curvilinear KEEP flux kernels (2nd-order).
!> E and F fluxes are area-scaled: flux = KEEP2(..., unit_normal) * face_area.
!> G flux uses uniform Normal_z (z is Cartesian); caller scales by J_2D via dt_Szeta.
module calc_keep_kernel_curv
  use libm
  use mod_globals, only : gamma, threadsE, threadsF, threadsG
  use mod_constant, only : R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  implicit none
  private
  public calc_keep_xi_curv, calc_keep_eta_curv, calc_keep_z_curv
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0
contains
  ! KEEP2/4/6 AND SLAU_common/SLAU1/HRSLAU2/phi are compiled here too
  ! (calc_scheme_math.f90 covers both schemes) even though this module only
  ! calls KEEP2 -- gamma above is the one SLAU_common's body references.
  include 'calc_scheme_math.f90'

  !> KEEP 2nd-order flux at xi-faces (i+1/2, j, k).
  !> n_xi_x/y(nx-1, ny-2): area-scaled face normals at interior eta cells.
  !> E is multiplied by face area |S_xi| so the step kernel uses scalar dt*dz.
  attributes(global) subroutine calc_keep_xi_curv(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q_1, Q_2, Q_3, Q_4, Q_5, T, E)
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8),    intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: T(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: E(nx-1,ny-2,nz-2,5)
    integer  :: i, j, k
    real(8)  :: nxx, nxy, S, Normal(5)
    real(8)  :: rho(2), u(2), v(2), w(2), uu(2), p(2), Tv(2)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_xi_x(i, j-1); nxy = n_xi_y(i, j-1)
    S   = sqrt(nxx*nxx + nxy*nxy)
    Normal = (/ 0.d0, nxx/S, nxy/S, 0.d0, 0.d0 /)
    rho(1) = Q_1(i,j,k); rho(2) = Q_1(i+1,j,k)
    u(1)   = Q_2(i,j,k); u(2)   = Q_2(i+1,j,k)
    v(1)   = Q_3(i,j,k); v(2)   = Q_3(i+1,j,k)
    w(1)   = Q_4(i,j,k); w(2)   = Q_4(i+1,j,k)
    p(1)   = Q_5(i,j,k); p(2)   = Q_5(i+1,j,k)
    Tv(1)  = T(i,  j,k); Tv(2)  = T(i+1,  j,k)
    uu(1)  = u(1)*Normal(2) + v(1)*Normal(3)
    uu(2)  = u(2)*Normal(2) + v(2)*Normal(3)
    E(i,j-1,k-1,:) = KEEP2(id_accuracy, rho, u, v, w, uu, p, Tv, Normal) * S
  end subroutine calc_keep_xi_curv


  !> KEEP 2nd-order flux at eta-faces (i, j+1/2, k).
  !> n_eta_x/y(nx-2, ny-1): area-scaled face normals at interior xi cells.
  attributes(global) subroutine calc_keep_eta_curv(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q_1, Q_2, Q_3, Q_4, Q_5, T, F)
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8),    intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: T(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: F(nx-2,ny-1,nz-2,5)
    integer  :: i, j, k
    real(8)  :: nex, ney, S, Normal(5)
    real(8)  :: rho(2), u(2), v(2), w(2), uu(2), p(2), Tv(2)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nex = n_eta_x(i-1, j); ney = n_eta_y(i-1, j)
    S   = sqrt(nex*nex + ney*ney)
    Normal = (/ 0.d0, nex/S, ney/S, 0.d0, 0.d0 /)
    rho(1) = Q_1(i,j,k); rho(2) = Q_1(i,j+1,k)
    u(1)   = Q_2(i,j,k); u(2)   = Q_2(i,j+1,k)
    v(1)   = Q_3(i,j,k); v(2)   = Q_3(i,j+1,k)
    w(1)   = Q_4(i,j,k); w(2)   = Q_4(i,j+1,k)
    p(1)   = Q_5(i,j,k); p(2)   = Q_5(i,j+1,k)
    Tv(1)  = T(i,  j,k); Tv(2)  = T(i,  j+1,k)
    uu(1)  = u(1)*Normal(2) + v(1)*Normal(3)
    uu(2)  = u(2)*Normal(2) + v(2)*Normal(3)
    F(i-1,j,k-1,:) = KEEP2(id_accuracy, rho, u, v, w, uu, p, Tv, Normal) * S
  end subroutine calc_keep_eta_curv


  !> KEEP 2nd-order flux at zeta-faces (i, j, k+1/2).
  !> Normal_z = (0,0,0,1,0); G is NOT area-scaled (J_2D accounted in dt_Szeta).
  attributes(global) subroutine calc_keep_z_curv(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5, T, G)
    use mod_constant, only : Normal_z
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: T(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: G(nx-2,ny-2,nz-1,5)
    integer :: i, j, k
    real(8) :: rho(2), u(2), v(2), w(2), uu(2), p(2), Tv(2)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(1) = Q_1(i,j,k); rho(2) = Q_1(i,j,k+1)
    u(1)   = Q_2(i,j,k); u(2)   = Q_2(i,j,k+1)
    v(1)   = Q_3(i,j,k); v(2)   = Q_3(i,j,k+1)
    w(1)   = Q_4(i,j,k); w(2)   = Q_4(i,j,k+1)
    p(1)   = Q_5(i,j,k); p(2)   = Q_5(i,j,k+1)
    Tv(1)  = T(i,  j,k); Tv(2)  = T(i,  j,k+1)
    uu(1)  = w(1);       uu(2)  = w(2)   ! Normal_z(4)=1 → uu = w
    G(i-1,j-1,k,:) = KEEP2(id_accuracy, rho, u, v, w, uu, p, Tv, Normal_z)
  end subroutine calc_keep_z_curv
end module calc_keep_kernel_curv
