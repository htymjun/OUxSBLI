module set
  use cudafor
  use mod_globals, only : rho0, p0, u0, rho2, p2, ux, uy
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1
  implicit none
contains
  ! Uniform Cartesian grid.
  subroutine set_grid(myrank, nx, ny, Lx, Ly, xc, yc, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: xc(nx), yc(ny), dx(nx-1), dy(ny-1)
    real(8) :: dx0, dy0
    integer :: i, j
    dx0 = Lx / dble(nx - 1)
    dy0 = Ly / dble(ny - 1)
    dx  = dx0
    dy  = dy0
    xc(1) = 0.0d0
    do i = 1, nx - 1
      xc(i + 1) = xc(i) + dx0
    end do
    yc(1) = 0.0d0
    do j = 1, ny - 1
      yc(j + 1) = yc(j) + dy0
    end do
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, x, y, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(out) :: Q(nx,ny,4)
    integer :: i, j
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1) = rho0
        Q(i,j,2) = rho0 * u0
        Q(i,j,3) = 0.d0
        Q(i,j,4) = p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2
    enddo;enddo
    do i = int(0.1d0 * nx), nx
      Q(i,ny,1) = rho2
      Q(i,ny,2) = rho2 * ux
      Q(i,ny,3) = rho2 * uy
      Q(i,ny,4) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo
  end subroutine set_init

  subroutine set_bc(myrank, nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    real(8) Jacobian_tmp
    integer :: i, j
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      Jacobian_tmp = 1.d0 / Jacobian(1,ny-1)
      if (i < int(0.1d0 * nx)) then
        QJ_1(i,ny) = rho0 * Jacobian_tmp
        QJ_2(i,ny) = rho0 * u0 * Jacobian_tmp
        QJ_3(i,ny) = 0.d0
        QJ_4(i,ny) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      else
        QJ_1(i,ny) = rho2 * Jacobian_tmp
        QJ_2(i,ny) = rho2 * ux * Jacobian_tmp
        QJ_3(i,ny) = rho2 * uy * Jacobian_tmp
        QJ_4(i,ny) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
      endif
      ! Slip
      QJ_1(i,1) =  QJ_1(i,2)
      QJ_2(i,1) =  QJ_2(i,2)
      QJ_3(i,1) = -QJ_3(i,2)
      QJ_4(i,1) =  QJ_4(i,2)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do j = 1, ny
      Jacobian_tmp = 1.d0 / Jacobian(1,j)
      QJ_1(1,j)  = rho0 * Jacobian_tmp
      QJ_2(1,j)  = rho0 * u0 * Jacobian_tmp
      QJ_3(1,j)  = 0.d0
      QJ_4(1,j)  = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      QJ_1(nx,j) = QJ_1(nx-1,j)
      QJ_2(nx,j) = QJ_2(nx-1,j)
      QJ_3(nx,j) = QJ_3(nx-1,j)
      QJ_4(nx,j) = QJ_4(nx-1,j)
    enddo
  end subroutine set_bc
end module set
