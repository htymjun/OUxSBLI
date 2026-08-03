module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy, rf, Taw
  use mod_constant, only : id_rescale, Cp, gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  use set_init_common
  implicit none
  real(8) ptbl
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    integer i, j, ny_b
    real(8) dx1, dy1
    dx1 = Lx / dble(nx-1)
    dy1 = dx1

    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
        ny_b  = j
      else
        dy(j) = dy1 * (1.d0 + 0.75d0 * dble(j-ny_b) / dble(ny-ny_b))
      endif
      y(j+1) = y(j) + dy(j)
    enddo
  end subroutine set_grid


  subroutine interpolate_Q(nx, nyi, ny, y, Qi, Qp)
    integer, intent(in)  :: nx, nyi, ny
    real(8), intent(in)  :: y(ny)
    real(8), intent(in)  :: Qi(5,nyi)
    real(8), intent(out) :: Qp(4,ny)
    integer :: i, j
    real(8) :: y1, y2, alpha
    do j = 1, ny
      if (y(j) <= Qi(1,1)) then
        Qp(:,j) = Qi(2:5,1)
      elseif (y(j) >= Qi(1,nyi)) then
        Qp(:,j) = Qi(2:5,nyi)
      else
        do i = 1, nyi-1
          y1 = Qi(1,i)
          y2 = Qi(1,i+1)
          if (y1 <= y(j) .and. y(j) <= y2) then
            alpha = (y(j) - y1) / (y2 - y1)
            Qp(:,j) = (1.0d0 - alpha) * Qi(2:5,i) + alpha * Qi(2:5,i+1)
            exit
          endif
        enddo
      endif
    enddo
  end subroutine interpolate_Q


  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,ny,4)
    character(len=40) filename
    integer j, nyi, filesize, ios
    real(8) rho, u, v, p
    real(8), allocatable :: Qi(:,:), Qp(:,:)
    write(filename, "(a)") "./Qin.dat"
    open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file Qin.dat"
    endif
    inquire(10, size=filesize)
    nyi = filesize / (8 * 5)
    allocate(Qi(5,nyi))
    read(10) Qi
    close(10)
    allocate(Qp(4,ny))
    call interpolate_Q(nx, nyi, ny, ys, Qi, Qp)
    ptbl = 0.d0
    do j = 1, nyi
      ptbl = ptbl + Qp(4,j)
    enddo
    ptbl = ptbl / dble(nyi)
    do j = 1, ny
      rho = Qp(1,j)
      u   = Qp(2,j)
      v   = Qp(3,j)
      Q(:,j,1) = rho
      Q(:,j,2) = rho * u
      Q(:,j,3) = rho * v
      Q(:,j,4) = ptbl * over_gamma_1 + 0.5d0 * rho * (u**2 + v**2)
    enddo
    deallocate(Qi, Qp)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny) ! Q / Jacobian
    integer i, j, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      ! outlet
      QJ_1(nx,j) = QJ_1(nx-1,j); QJ_2(nx,j) = QJ_2(nx-1,j); QJ_3(nx,j) = QJ_3(nx-1,j); QJ_4(nx,j) = QJ_4(nx-1,j)
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! Neumann
      QJ_1(i,ny) = QJ_1(i,ny-1)
      QJ_2(i,ny) = QJ_2(i,ny-1)
      QJ_3(i,ny) = QJ_3(i,ny-1)
      QJ_4(i,ny) = QJ_4(i,ny-1)
      ! NoSlip
      QJ_1(i,1) = QJ_1(i,2)
      QJ_2(i,1) = 0.d0
      QJ_3(i,1) = 0.d0
      p_wall = gamma_1 * (QJ_4(i,2) - 0.5d0 * (QJ_2(i,2)**2 + QJ_3(i,2)**2) / QJ_1(i,2))
      QJ_4(i,1) = p_wall * over_gamma_1
    enddo
  end subroutine set_bc
end module set
