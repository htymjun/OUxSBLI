module set
  use mod_globals, only : nx, ny, gamma, R
  use mod_constant, only : id_accuracy
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, xc, yc, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: xc(nx), yc(ny), dx(nx-1), dy(ny-1)
    call set_grid_cyclic(id_accuracy, nx, ny, Lx, Ly, xc, yc, dx, dy)
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, x, y, Q)
    use mod_globals, only : pi, M0, rho0, u0, d1, d2
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(out) :: Q(nx,ny,4)
    integer i, j
    real(8) :: p = rho0 * u0**2 / (gamma * M0**2)
    do j = 1, ny
      do i = 1, nx
        if (y(j) <= pi) then
          Q(i,j,1) = rho0
          Q(i,j,2) = rho0 * u0 * tanh((y(j) - 0.5d0 * pi) / d1)
          Q(i,j,3) = rho0 * u0 * d2 * sin(x(i))
          Q(i,j,4) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,j,2)**2 + Q(i,j,3)**2) / Q(i,j,1)
        else
          Q(i,j,1) = rho0
          Q(i,j,2) = rho0 * u0 * tanh((1.5d0 * pi - y(j)) / d1)
          Q(i,j,3) = rho0 * u0 * d2 * sin(x(i))
          Q(i,j,4) = p / (gamma - 1.d0) + 0.5d0 * (Q(i,j,2)**2 + Q(i,j,3)**2) / Q(i,j,1)
        endif
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, Q_1, Q_2, Q_3, Q_4)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q_1(nx,ny), Q_2(nx,ny), Q_3(nx,ny), Q_4(nx,ny)
    call set_bc_cyclic(id_accuracy, nx, ny, Q_1, Q_2, Q_3, Q_4)
  end subroutine set_bc
end module set

