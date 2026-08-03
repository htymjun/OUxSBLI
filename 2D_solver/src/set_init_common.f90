module set_init_common
  use mod_globals, only : gamma, R, Taw, rf
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  use set_compressible_bl
  implicit none
contains
  subroutine set_init_tbl(nx, ny, x, y, blt0, blt, u0, p0, T0, M0, Q)
    integer, intent(in)  :: nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(in)  :: blt0, blt, u0, p0, T0, M0
    real(8), intent(out) :: Q(nx,ny,4)
    integer i, j
    real(8) :: p_wall
    real(8), allocatable :: rho(:), u(:), v(:), T(:)
    allocate(rho(ny), u(ny), v(ny), T(ny))
    call calc_HD_Blasius(ny, y, blt0, u0, T0, p0, M0, rho, u, v, T)
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1) = p0 / (R * T(j))
        Q(i,j,2) = Q(i,j,1) * u(j)
        Q(i,j,3) = Q(i,j,1) * v(j)
        Q(i,j,4) = p0 * over_gamma_1 + 0.5d0 * (Q(i,j,2)**2 + Q(i,j,3)**2) / Q(i,j,1)
    enddo;enddo
    deallocate(rho, u, v, T)
    ! bottom
    Q(:,1,1) = Q(:,2,1)
    Q(:,1,2) = 0.d0
    Q(:,1,3) = 0.d0
    p_wall = gamma_1 * (Q(2,2,4) - 0.5d0 * (Q(2,2,2)**2 + Q(2,2,3)**2) / Q(2,2,1))
    Q(:,1,4) = p_wall * over_gamma_1
  end subroutine set_init_tbl
end module set_init_common

