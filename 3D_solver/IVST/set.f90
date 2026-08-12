module set
  use mod_globals, only : nx, ny, nz, gamma, R, rhol => rho0, rhor => rho1, pl => p0, pr => p1
  use set_bc_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx), dy(ny), dz(nz)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-1)
    dy1 = Ly / dble(ny-1)
    dz1 = Lz / dble(nz-1)
    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      dy(j) = dy1
      y(j+1) = y(j) + dy(j)
    enddo
    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    do i = 1, nx
      if (i < int(0.5*nx)) then
        Q(i,:,:,1) = rhol
        Q(i,:,:,2) = 0.d0
        Q(i,:,:,3) = 0.d0
        Q(i,:,:,4) = 0.d0
        Q(i,:,:,5) = pl / (gamma  - 1.d0)
      else
        Q(i,:,:,1) = rhor
        Q(i,:,:,2) = 0.d0
        Q(i,:,:,3) = 0.d0
        Q(i,:,:,4) = 0.d0
        Q(i,:,:,5) = pr / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q_1, Q_2, Q_3, Q_4, Q_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in), value            :: myrank, nx, ny, nz
    real(8), intent(in), device           :: Jacobian(ny)
    real(8), intent(inout), device        :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz) ! Q / J
    real(8), intent(in), device, optional :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    integer :: i, j, k, jc = 4, kc = 4
    real(8), device :: Qc_1(nx), Qc_2(nx), Qc_3(nx), Qc_4(nx), Qc_5(nx)
    ! inlet and outlet
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, 4
      do j = 4, 4
        Q_1(1,j,k)    = rhol / Jacobian(j)
        Q_2(1,j,k)    = 0.d0
        Q_3(1,j,k)    = 0.d0
        Q_4(1,j,k)    = 0.d0
        Q_5(1,j,k)    = pl / (gamma - 1.d0) / Jacobian(j)
        Q_1(2,j,k)    = rhol / Jacobian(j)
        Q_2(2,j,k)    = 0.d0
        Q_3(2,j,k)    = 0.d0
        Q_4(2,j,k)    = 0.d0
        Q_5(2,j,k)    = pl / (gamma - 1.d0) / Jacobian(j)
        Q_1(3,j,k)    = rhol / Jacobian(j)
        Q_2(3,j,k)    = 0.d0
        Q_3(3,j,k)    = 0.d0
        Q_4(3,j,k)    = 0.d0
        Q_5(3,j,k)    = pl / (gamma - 1.d0) / Jacobian(j)
        Q_1(nx-2,j,k) = rhor / Jacobian(j)
        Q_2(nx-2,j,k) = 0.d0
        Q_3(nx-2,j,k) = 0.d0
        Q_4(nx-2,j,k) = 0.d0
        Q_5(nx-2,j,k) = pr / (gamma - 1.d0) / Jacobian(j)
        Q_1(nx-1,j,k) = rhor / Jacobian(j)
        Q_2(nx-1,j,k) = 0.d0
        Q_3(nx-1,j,k) = 0.d0
        Q_4(nx-1,j,k) = 0.d0
        Q_5(nx-1,j,k) = pr / (gamma - 1.d0) / Jacobian(j)
        Q_1(nx,j,k)   = rhor / Jacobian(j)
        Q_2(nx,j,k)   = 0.d0
        Q_3(nx,j,k)   = 0.d0
        Q_4(nx,j,k)   = 0.d0
        Q_5(nx,j,k)   = pr / (gamma - 1.d0) / Jacobian(j)
    enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      Qc_1(i) = Q_1(i,jc,kc)
      Qc_2(i) = Q_2(i,jc,kc)
      Qc_3(i) = Q_3(i,jc,kc)
      Qc_4(i) = Q_4(i,jc,kc)
      Qc_5(i) = Q_5(i,jc,kc)
    enddo
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q_1(i,j,k) = Qc_1(i)
          Q_2(i,j,k) = Qc_2(i)
          Q_3(i,j,k) = Qc_3(i)
          Q_4(i,j,k) = Qc_4(i)
          Q_5(i,j,k) = Qc_5(i)
    enddo;enddo;enddo
  end subroutine set_bc


  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set

