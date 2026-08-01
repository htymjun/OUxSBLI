!> Curvilinear Hybrid flux kernels with Ducros-based scheme blending
!> Automatically switches between KEEP (smooth) and SLAU (shock) schemes
module calc_hybrid_kernel_curv
  use libm
  use mod_globals, only : gamma, threshold, threadsE, threadsF, threadsG
  use mod_constant, only : over_gamma_1, R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third, id_accuracy, id_slau
  use calc_muscl
  use calc_hybrid_curv
  implicit none
  private
  public calc_hybrid_xi_curv, calc_hybrid_eta_curv, calc_hybrid_z_curv

  real(8), parameter :: one_24  = 1.d0 / 24.d0
  real(8), parameter :: one_48  = 1.d0 / 48.d0
  real(8), parameter :: one_60  = 1.d0 / 60.d0
  real(8), parameter :: one_120 = 1.d0 / 120.d0
  real(8), parameter :: one_240 = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0

  interface KEEP
    module procedure KEEP2
  end interface KEEP

  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU

contains
  include '../../3D_solver/src/calc_keep_3d.f90'
  include '../../3D_solver/src/calc_slau_3d.f90'


  !> Hybrid flux at xi-faces (i+1/2, j, k). Area-scaled.
  !> Blends KEEP and SLAU based on Ducros sensor threshold
  attributes(global) subroutine calc_hybrid_xi_curv(id_accuracy, nx, ny, nz, &
                                                      n_xi_x, n_xi_y, Q_1, Q_2, Q_3, Q_4, Q_5, T, sensor, E)
    integer(2), intent(in), value             :: id_accuracy
    integer,  intent(in), value               :: nx, ny, nz
    real(8),  intent(in), device, contiguous  :: n_xi_x(nx-1,ny-2), n_xi_y(nx-1,ny-2)
    real(8),  intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: T(nx,ny,nz)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),  intent(out), device, contiguous :: E(nx-1,ny-2,nz-2,5)
    integer :: i, j, k
    real(8) :: nxx, nxy, S, Normal(5)
    real(8) :: rho(2), u(2), v(2), w(2), uu(2), p(2), Tv(2)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr
    real(sp) :: fdx
    real(8) :: EKeep(5), ESLAU(5)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_xi_x(i, j-1);  nxy = n_xi_y(i, j-1)
    S   = sqrt(nxx*nxx + nxy*nxy)
    Normal = (/ 0.d0, nxx/S, nxy/S, 0.d0, 0.d0 /)

    rho(1) = Q_1(i,j,k);  rho(2) = Q_1(i+1,j,k)
    u(1)   = Q_2(i,j,k);  u(2)   = Q_2(i+1,j,k)
    v(1)   = Q_3(i,j,k);  v(2)   = Q_3(i+1,j,k)
    w(1)   = Q_4(i,j,k);  w(2)   = Q_4(i+1,j,k)
    p(1)   = Q_5(i,j,k);  p(2)   = Q_5(i+1,j,k)
    Tv(1)  = T(i,  j,k);  Tv(2)  = T(i+1,  j,k)
    uu(1)  = u(1)*Normal(2) + v(1)*Normal(3)
    uu(2)  = u(2)*Normal(2) + v(2)*Normal(3)

    fdx = 0.5_sp * (sensor(i,j,k) + sensor(i+1,j,k))
    if (fdx <= threshold) then
      ! KEEP: smooth flow region
      EKeep = KEEP2(id_accuracy, rho, u, v, w, uu, p, Tv, Normal) * S
      E(i,j-1,k-1,:) = EKeep
    else
      ! SLAU: shocked region
      rhol = rho(1);  rhor = rho(2)
      ul   = u(1);    ur   = u(2)
      vl   = v(1);    vr   = v(2)
      wl   = w(1);    wr   = w(2)
      pl   = p(1);    pr   = p(2)
      call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, uu(1), uu(2), pl, pr, &
                Normal, fdx, ESLAU(1), ESLAU(2), ESLAU(3), ESLAU(4), ESLAU(5))
      E(i,j-1,k-1,:) = ESLAU * S
    endif
  end subroutine calc_hybrid_xi_curv


  !> Hybrid flux at eta-faces (i, j+1/2, k). Area-scaled.
  attributes(global) subroutine calc_hybrid_eta_curv(id_accuracy, nx, ny, nz, &
                                                       n_eta_x, n_eta_y, Q_1, Q_2, Q_3, Q_4, Q_5, T, sensor, F)
    integer(2), intent(in), value             :: id_accuracy
    integer,  intent(in), value               :: nx, ny, nz
    real(8),  intent(in), device, contiguous  :: n_eta_x(nx-2,ny-1), n_eta_y(nx-2,ny-1)
    real(8),  intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: T(nx,ny,nz)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),  intent(out), device, contiguous :: F(nx-2,ny-1,nz-2,5)
    integer :: i, j, k
    real(8) :: nxx, nxy, S, Normal(5)
    real(8) :: rho(2), u(2), v(2), w(2), uu(2), p(2), Tv(2)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr
    real(sp) :: fdy
    real(8) :: FKeep(5), FSLAU(5)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_eta_x(i-1, j);  nxy = n_eta_y(i-1, j)
    S   = sqrt(nxx*nxx + nxy*nxy)
    Normal = (/ 0.d0, nxx/S, nxy/S, 0.d0, 0.d0 /)

    rho(1) = Q_1(i,j,k);  rho(2) = Q_1(i,j+1,k)
    u(1)   = Q_2(i,j,k);  u(2)   = Q_2(i,j+1,k)
    v(1)   = Q_3(i,j,k);  v(2)   = Q_3(i,j+1,k)
    w(1)   = Q_4(i,j,k);  w(2)   = Q_4(i,j+1,k)
    p(1)   = Q_5(i,j,k);  p(2)   = Q_5(i,j+1,k)
    Tv(1)  = T(i,j,k);    Tv(2)  = T(i,j+1,k)
    uu(1)  = u(1)*Normal(2) + v(1)*Normal(3)
    uu(2)  = u(2)*Normal(2) + v(2)*Normal(3)

    fdy = 0.5_sp * (sensor(i,j,k) + sensor(i,j+1,k))

    if (fdy <= threshold) then
      ! KEEP: smooth flow region
      FKeep = KEEP2(id_accuracy, rho, u, v, w, uu, p, Tv, Normal) * S
      F(i-1,j,k-1,:) = FKeep
    else
      ! SLAU: shocked region
      rhol = rho(1);  rhor = rho(2)
      ul   = u(1);    ur   = u(2)
      vl   = v(1);    vr   = v(2)
      wl   = w(1);    wr   = w(2)
      pl   = p(1);    pr   = p(2)
      call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, uu(1), uu(2), pl, pr, &
                Normal, fdy, FSLAU(1), FSLAU(2), FSLAU(3), FSLAU(4), FSLAU(5))
      F(i-1,j,k-1,:) = FSLAU * S
    endif
  end subroutine calc_hybrid_eta_curv


  !> Hybrid flux at z-faces (i, j, k+1/2). NOT area-scaled (z is uniform Cartesian).
  !> Caller will scale by dt_Szeta = dt * J_2D(i,j)
  attributes(global) subroutine calc_hybrid_z_curv(id_accuracy, nx, ny, nz, &
                                                     Q_1, Q_2, Q_3, Q_4, Q_5, T, sensor, G)
    use mod_constant, only : Normal_z
    integer(2), intent(in), value             :: id_accuracy
    integer,  intent(in), value               :: nx, ny, nz
    real(8),  intent(in), device, contiguous  :: Q_1(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_2(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_3(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_4(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: Q_5(nx,ny,nz)
    real(8),  intent(in), device, contiguous  :: T(nx,ny,nz)
    real(sp), intent(in), device, contiguous  :: sensor(nx,ny,nz)
    real(8),  intent(out), device, contiguous :: G(nx-2,ny-2,nz-1,5)
    integer :: i, j, k
    real(8) :: rho(2), u(2), v(2), w(2), uu(2), p(2), Tv(2)
    real(8) :: rhol, rhor, ul, ur, vl, vr, wl, wr, pl, pr
    real(sp) :: fdz
    real(8) :: GKeep(5), GSLAU(5)
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(1) = Q_1(i,j,k);    rho(2) = Q_1(i,j,k+1)
    u(1)   = Q_2(i,j,k);    u(2)   = Q_2(i,j,k+1)
    v(1)   = Q_3(i,j,k);    v(2)   = Q_3(i,j,k+1)
    w(1)   = Q_4(i,j,k);    w(2)   = Q_4(i,j,k+1)
    p(1)   = Q_5(i,j,k);    p(2)   = Q_5(i,j,k+1)
    Tv(1)  = T(i,j,k);      Tv(2)  = T(i,j,k+1)
    uu(1)  = w(1)           ! Normal_z = (0,0,1)
    uu(2)  = w(2)

    fdz = 0.5_sp * (sensor(i,j,k) + sensor(i,j,k+1))

    if (fdz <= threshold) then
      ! KEEP: smooth flow region
      GKeep = KEEP2(id_accuracy, rho, u, v, w, uu, p, Tv, Normal_z)
      G(i-1,j-1,k,:) = GKeep
    else
      ! SLAU: shocked region
      rhol = rho(1);  rhor = rho(2)
      ul   = u(1);    ur   = u(2)
      vl   = v(1);    vr   = v(2)
      wl   = w(1);    wr   = w(2)
      pl   = p(1);    pr   = p(2)
      call SLAU(id_slau, rhol, rhor, ul, ur, vl, vr, wl, wr, uu(1), uu(2), pl, pr, &
                Normal_z, fdz, GSLAU(1), GSLAU(2), GSLAU(3), GSLAU(4), GSLAU(5))
      G(i-1,j-1,k,:) = GSLAU
    endif
  end subroutine calc_hybrid_z_curv
end module calc_hybrid_kernel_curv
