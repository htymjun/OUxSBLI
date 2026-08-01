module set
  use mod_globals, only : nx, ny, nz, gamma, rho_L, p_L, rho_R, p_R, Lx, Ly, Lz
  use set_bc_common
  use set_coordinate
  implicit none
contains

  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    use mpi
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1, xf(nx+1), yf(ny+1), zf(nz+1)
    integer nranks, ierr, nz_int, iz_offset, i, j, k
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! interior per rank = nz - 2*overlap_fb; z-rank index = myrank/2
    nz_int    = nz - 2*(kind(id_accuracy)/3 + 1)
    iz_offset = (myrank/2) * nz_int
    dx1 = Lx / dble(nx - 2)
    dy1 = Ly / dble(ny - 2)
    dz1 = Lz / dble((nranks/2) * nz_int)
    dx(:) = dx1;  dy(:) = dy1;  dz(:) = dz1
    ! x: cell-face positions then cell centres
    do i = 2, nx;  xf(i) = dx1 * dble(i-2);  enddo
    xf(1) = xf(2) - dx1;  xf(nx+1) = xf(nx) + dx1
    do i = 1, nx;  x(i) = 0.5d0*(xf(i)+xf(i+1));  enddo
    ! y
    do j = 2, ny;  yf(j) = dy1 * dble(j-2);  enddo
    yf(1) = yf(2) - dy1;  yf(ny+1) = yf(ny) + dy1
    do j = 1, ny;  y(j) = 0.5d0*(yf(j)+yf(j+1));  enddo
    ! z: rank-local slab, with ghost face positions
    do k = 2, nz;  zf(k) = dz1 * dble(iz_offset + k - 2);  enddo
    zf(1) = zf(2) - dz1;  zf(nz+1) = zf(nz) + dz1
    do k = 1, nz;  z(k) = 0.5d0*(zf(k)+zf(k+1));  enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (myrank < 2) then
            Q(i,j,k,1) = rho_L
            Q(i,j,k,5) = p_L / (gamma - 1.d0)
          else
            Q(i,j,k,1) = rho_R
            Q(i,j,k,5) = p_R / (gamma - 1.d0)
          endif
          Q(i,j,k,2) = 0.d0   ! rho*u
          Q(i,j,k,3) = 0.d0   ! rho*v
          Q(i,j,k,4) = 0.d0   ! rho*w
        enddo
      enddo
    enddo
  end subroutine set_init


  ! Custom BC for z-decomposition:
  !   - periodic in x and y (for all z including ghost layers)
  !   - zero-gradient in z only for outermost ranks
  !   - does NOT call set_bc_cyclic (which would corrupt z ghost cells)
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    use mpi
    use mod_constant, only : id_accuracy
    integer, intent(in), value               :: myrank, nx, ny, nz
    real(8), intent(in), device              :: Jacobian(nx,ny)
    real(8), intent(inout), device           :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(in), device, optional    :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    integer nranks, ierr, i, j, k, m, ovlp
    ovlp = kind(id_accuracy)/3 + 1
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! x: Dirichlet slip wall — reflect ρu (l=2) to enforce u=0 at face; extrapolate rest
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do m = 0, ovlp-1
          QJ_1(ovlp-m,      j,k) =  QJ_1(ovlp+1+m,  j,k)
          QJ_2(ovlp-m,      j,k) = -QJ_2(ovlp+1+m,  j,k)
          QJ_3(ovlp-m,      j,k) =  QJ_3(ovlp+1+m,  j,k)
          QJ_4(ovlp-m,      j,k) =  QJ_4(ovlp+1+m,  j,k)
          QJ_5(ovlp-m,      j,k) =  QJ_5(ovlp+1+m,  j,k)
          QJ_1(nx-ovlp+1+m, j,k) =  QJ_1(nx-ovlp-m, j,k)
          QJ_2(nx-ovlp+1+m, j,k) = -QJ_2(nx-ovlp-m, j,k)
          QJ_3(nx-ovlp+1+m, j,k) =  QJ_3(nx-ovlp-m, j,k)
          QJ_4(nx-ovlp+1+m, j,k) =  QJ_4(nx-ovlp-m, j,k)
          QJ_5(nx-ovlp+1+m, j,k) =  QJ_5(nx-ovlp-m, j,k)
        enddo
      enddo
    enddo
    ! y: Dirichlet slip wall — reflect ρv (l=3) to enforce v=0 at face; extrapolate rest
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        do m = 0, ovlp-1
          QJ_1(i,ovlp-m,     k) =  QJ_1(i,ovlp+1+m,  k)
          QJ_2(i,ovlp-m,     k) =  QJ_2(i,ovlp+1+m,  k)
          QJ_3(i,ovlp-m,     k) = -QJ_3(i,ovlp+1+m,  k)
          QJ_4(i,ovlp-m,     k) =  QJ_4(i,ovlp+1+m,  k)
          QJ_5(i,ovlp-m,     k) =  QJ_5(i,ovlp+1+m,  k)
          QJ_1(i,ny-ovlp+1+m,k) =  QJ_1(i,ny-ovlp-m, k)
          QJ_2(i,ny-ovlp+1+m,k) =  QJ_2(i,ny-ovlp-m, k)
          QJ_3(i,ny-ovlp+1+m,k) = -QJ_3(i,ny-ovlp-m, k)
          QJ_4(i,ny-ovlp+1+m,k) =  QJ_4(i,ny-ovlp-m, k)
          QJ_5(i,ny-ovlp+1+m,k) =  QJ_5(i,ny-ovlp-m, k)
        enddo
      enddo
    enddo
    ! z lo boundary: zero-gradient on rank 0 — set all ovlp ghost cells
    if (myrank == 0) then
      !$cuf kernel do(3)<<<*,*>>>
      do k = 1, ovlp
        do j = 1, ny
          do i = 1, nx
            QJ_1(i,j,k) = QJ_1(i,j,ovlp+1)
            QJ_2(i,j,k) = QJ_2(i,j,ovlp+1)
            QJ_3(i,j,k) = QJ_3(i,j,ovlp+1)
            QJ_4(i,j,k) = QJ_4(i,j,ovlp+1)
            QJ_5(i,j,k) = QJ_5(i,j,ovlp+1)
          enddo
        enddo
      enddo
    endif
    ! z hi boundary: zero-gradient on last compute rank (nranks-2 in even/odd pattern)
    if (myrank == nranks - 2) then
      !$cuf kernel do(3)<<<*,*>>>
      do k = 1, ovlp
        do j = 1, ny
          do i = 1, nx
            QJ_1(i,j,nz-ovlp+k) = QJ_1(i,j,nz-ovlp)
            QJ_2(i,j,nz-ovlp+k) = QJ_2(i,j,nz-ovlp)
            QJ_3(i,j,nz-ovlp+k) = QJ_3(i,j,nz-ovlp)
            QJ_4(i,j,nz-ovlp+k) = QJ_4(i,j,nz-ovlp)
            QJ_5(i,j,nz-ovlp+k) = QJ_5(i,j,nz-ovlp)
          enddo
        enddo
      enddo
    endif
  end subroutine set_bc


  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set
