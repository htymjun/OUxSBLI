!> Curvilinear SLAU flux kernels (2nd-order, no MUSCL reconstruction).
!> E and F are area-scaled; G uses Normal_z (not area-scaled).
module calc_slau_kernel_curv
  use libm
  use mod_globals, only : threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1, id_slau, R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  use calc_hybrid
  implicit none
  private
  public calc_slau_xi_curv, calc_slau_eta_curv, calc_slau_z_curv
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU
  ! KEEP2/4/6 are compiled here too (calc_scheme_math.f90 include below covers
  ! both schemes) even though this module never calls them -- these constants
  ! are the ones their bodies reference. gamma/sp (needed by SLAU_common) are
  ! already in scope transitively via `use calc_hybrid` above.
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0
contains
  include 'calc_scheme_math.f90'

  !> SLAU 2nd-order flux at xi-faces (i+1/2, j, k). Area-scaled by |S_xi|.
  attributes(global) subroutine calc_slau_xi_curv(id_accuracy, nx, ny, nz, n_xi_x, n_xi_y, Q_1, Q_2, Q_3, Q_4, Q_5, sensor, E)
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8),    intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(sp),   intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: E(nx-1,ny-2,nz-2,5)
    integer :: i, j, k
    real(8) :: nxx, nxy, S, Normal(5)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr, unl, unr
    real(sp) :: fdx
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_xi_x(i,j-1); nxy = n_xi_y(i,j-1)
    S   = sqrt(nxx*nxx + nxy*nxy)
    Normal = (/ 0.d0, nxx/S, nxy/S, 0.d0, 0.d0 /)
    rhol = Q_1(i,j,k); rhor = Q_1(i+1,j,k)
    ul   = Q_2(i,j,k); ur   = Q_2(i+1,j,k)
    vl   = Q_3(i,j,k); vr   = Q_3(i+1,j,k)
    wl   = Q_4(i,j,k); wr   = Q_4(i+1,j,k)
    pl   = Q_5(i,j,k); pr   = Q_5(i+1,j,k)
    unl  = ul*Normal(2) + vl*Normal(3)
    unr  = ur*Normal(2) + vr*Normal(3)
    fdx  = 0.5_sp * (sensor(i,j,k) + sensor(i+1,j,k))
    call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, unl, unr, pl, pr, Normal, fdx, &
              E(i,j-1,k-1,1), E(i,j-1,k-1,2), E(i,j-1,k-1,3), E(i,j-1,k-1,4), E(i,j-1,k-1,5))
    E(i,j-1,k-1,:) = E(i,j-1,k-1,:) * S
  end subroutine calc_slau_xi_curv


  !> SLAU 2nd-order flux at eta-faces (i, j+1/2, k). Area-scaled by |S_eta|.
  attributes(global) subroutine calc_slau_eta_curv(id_accuracy, nx, ny, nz, n_eta_x, n_eta_y, Q_1, Q_2, Q_3, Q_4, Q_5, sensor, F)
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8),    intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(sp),   intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: F(nx-2,ny-1,nz-2,5)
    integer :: i, j, k
    real(8) :: nex, ney, S, Normal(5)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr, unl, unr
    real(sp) :: fdy
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nex = n_eta_x(i-1,j); ney = n_eta_y(i-1,j)
    S   = sqrt(nex*nex + ney*ney)
    Normal = (/ 0.d0, nex/S, ney/S, 0.d0, 0.d0 /)
    rhol = Q_1(i,j,k); rhor = Q_1(i,j+1,k)
    ul   = Q_2(i,j,k); ur   = Q_2(i,j+1,k)
    vl   = Q_3(i,j,k); vr   = Q_3(i,j+1,k)
    wl   = Q_4(i,j,k); wr   = Q_4(i,j+1,k)
    pl   = Q_5(i,j,k); pr   = Q_5(i,j+1,k)
    unl  = ul*Normal(2) + vl*Normal(3)
    unr  = ur*Normal(2) + vr*Normal(3)
    fdy  = 0.5_sp * (sensor(i,j,k) + sensor(i,j+1,k))
    call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, unl, unr, pl, pr, Normal, fdy, &
              F(i-1,j,k-1,1), F(i-1,j,k-1,2), F(i-1,j,k-1,3), F(i-1,j,k-1,4), F(i-1,j,k-1,5))
    F(i-1,j,k-1,:) = F(i-1,j,k-1,:) * S
  end subroutine calc_slau_eta_curv


  !> SLAU 2nd-order flux at zeta-faces. Normal_z; G not area-scaled.
  attributes(global) subroutine calc_slau_z_curv(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5, sensor, G)
    use mod_constant, only : Normal_z
    integer(2), intent(in), value               :: id_accuracy
    integer,    intent(in), value               :: nx, ny, nz
    real(8),    intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),    intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(sp),   intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),    intent(out), device, contiguous :: G(nx-2,ny-2,nz-1,5)
    integer :: i, j, k
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr
    real(sp) :: fdz
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rhol = Q_1(i,j,k); rhor = Q_1(i,j,k+1)
    ul   = Q_2(i,j,k); ur   = Q_2(i,j,k+1)
    vl   = Q_3(i,j,k); vr   = Q_3(i,j,k+1)
    wl   = Q_4(i,j,k); wr   = Q_4(i,j,k+1)
    pl   = Q_5(i,j,k); pr   = Q_5(i,j,k+1)
    fdz  = 0.5_sp * (sensor(i,j,k) + sensor(i,j,k+1))
    call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, wl, wr, pl, pr, Normal_z, fdz, &
              G(i-1,j-1,k,1), G(i-1,j-1,k,2), G(i-1,j-1,k,3), G(i-1,j-1,k,4), G(i-1,j-1,k,5))
  end subroutine calc_slau_z_curv
end module calc_slau_kernel_curv
