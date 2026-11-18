module set
#ifdef PREC_DF
  use mod_prec
#endif
  use cudafor
  use mod_globals, only : nx, ny, gamma, R, rho0, p0, rho1, p1
  use mod_constant, only : id_accuracy, gamma_1, over_gamma_1
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    real(8) dx0, dy0
    integer i, j
    dx0 = Lx / dble(nx-1)
    dy0 = Ly / dble(ny-1)
    dx  = dx0
    dy  = dy0
    x(1) = 0.d0
    do i = 1, nx-1
      x(i+1) = x(i) + dx0
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      y(j+1) = y(j) + dy0
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, x, y, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(8), intent(out) :: Q(nx,ny,4)
    integer i, j
    do j = 1, ny
      do i = 1, nx / 2
        Q(i,j,1) = rho0
        Q(i,j,2) = 0.d0
        Q(i,j,3) = 0.d0
        Q(i,j,4) = p0 / (gamma - 1.d0)
      enddo
      do i = nx / 2 + 1, nx
        Q(i,j,1) = rho1
        Q(i,j,2) = 0.d0
        Q(i,j,3) = 0.d0
        Q(i,j,4) = p1 / (gamma - 1.d0)
      enddo
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: myrank, nx, ny
    REAL_W, intent(in), device     :: Jacobian(nx,ny)
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    call bc_inout<<<(ny+127)/128, 128>>>(nx, ny, QJ_1, QJ_2, QJ_3, QJ_4)
    call bc_wall_sym<<<(nx+127)/128, 128>>>(nx, ny, QJ_1, QJ_2, QJ_3, QJ_4)
  end subroutine set_bc


  !> inlet / outlet: zero gradient in x (j = 2 .. ny-1)
  attributes(global) subroutine bc_inout(nx, ny, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: nx, ny
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    integer j
    j = (blockIdx%x - 1) * blockDim%x + threadIdx%x + 1
    if (j > ny-1) return
    QJ_1(1,j)  = QJ_1(2,j);    QJ_2(1,j)  = QJ_2(2,j);    QJ_3(1,j)  = QJ_3(2,j);    QJ_4(1,j)  = QJ_4(2,j)
    QJ_1(nx,j) = QJ_1(nx-1,j); QJ_2(nx,j) = QJ_2(nx-1,j); QJ_3(nx,j) = QJ_3(nx-1,j); QJ_4(nx,j) = QJ_4(nx-1,j)
  end subroutine bc_inout


  !> no-slip wall at j = 1, symmetric boundary at j = ny (i = 1 .. nx)
  attributes(global) subroutine bc_wall_sym(nx, ny, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: nx, ny
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    REAL_W, parameter :: zero = 0.d0, half = 0.5d0
    REAL_T p_wall
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx) return
    ! no-slip wall
    p_wall = gamma_1 * (QJ_4(i,2) - half * (QJ_2(i,2)**2 + QJ_3(i,2)**2) / QJ_1(i,2))
    QJ_1(i,1) = QJ_1(i,2)
    QJ_2(i,1) = zero
    QJ_3(i,1) = zero
    QJ_4(i,1) = p_wall * over_gamma_1
    ! symetric boundary condition
    QJ_1(i,ny) =  QJ_1(i,ny-1)
    QJ_2(i,ny) =  QJ_2(i,ny-1)
    QJ_3(i,ny) = -QJ_3(i,ny-1)
    QJ_4(i,ny) =  QJ_4(i,ny-1)
  end subroutine bc_wall_sym
end module set
