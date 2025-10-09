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
    real(8), intent(out) :: Q(5,nx,ny,nz)
    integer i, j, k
    do i = 1, nx
      if (i < int(0.5*nx)) then
        Q(1,i,:,:) = rhol
        Q(2,i,:,:) = 0.d0
        Q(3,i,:,:) = 0.d0
        Q(4,i,:,:) = 0.d0
        Q(5,i,:,:) = pl / (gamma  - 1.d0)
      else
        Q(1,i,:,:) = rhor
        Q(2,i,:,:) = 0.d0
        Q(3,i,:,:) = 0.d0
        Q(4,i,:,:) = 0.d0
        Q(5,i,:,:) = pr / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: Q(5,nx,ny,nz) ! Q / J
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    integer :: i, j, k, l, jc = 4, kc = 4
    real(8), device :: Qc(5,nx)
  
    ! inlet and outlet
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, 4
      do j = 4, 4
        Q(1,1,j,k)    = rhol / Jacobian(j)
        Q(2,1,j,k)    = 0.d0
        Q(3,1,j,k)    = 0.d0
        Q(4,1,j,k)    = 0.d0
        Q(5,1,j,k)    = pl / (gamma - 1.d0) / Jacobian(j)
        Q(1,2,j,k)    = rhol / Jacobian(j)
        Q(2,2,j,k)    = 0.d0
        Q(3,2,j,k)    = 0.d0
        Q(4,2,j,k)    = 0.d0
        Q(5,2,j,k)    = pl / (gamma - 1.d0) / Jacobian(j)
        Q(1,3,j,k)    = rhol / Jacobian(j)
        Q(2,3,j,k)    = 0.d0
        Q(3,3,j,k)    = 0.d0
        Q(4,3,j,k)    = 0.d0
        Q(5,3,j,k)    = pl / (gamma - 1.d0) / Jacobian(j)
        Q(1,nx-2,j,k) = rhor / Jacobian(j)
        Q(2,nx-2,j,k) = 0.d0
        Q(3,nx-2,j,k) = 0.d0
        Q(4,nx-2,j,k) = 0.d0
        Q(5,nx-2,j,k) = pr / (gamma - 1.d0) / Jacobian(j)
        Q(1,nx-1,j,k) = rhor / Jacobian(j)
        Q(2,nx-1,j,k) = 0.d0
        Q(3,nx-1,j,k) = 0.d0
        Q(4,nx-1,j,k) = 0.d0
        Q(5,nx-1,j,k) = pr / (gamma - 1.d0) / Jacobian(j)
        Q(1,nx,j,k)   = rhor / Jacobian(j)
        Q(2,nx,j,k)   = 0.d0
        Q(3,nx,j,k)   = 0.d0
        Q(4,nx,j,k)   = 0.d0
        Q(5,nx,j,k)   = pr / (gamma - 1.d0) / Jacobian(j)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do i = 1, nx
      do l = 1, 5
        Qc(l,i) = Q(l,i,jc,kc)
    enddo;enddo

    !$cuf kernel do(4)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            Q(l,i,j,k) = Qc(l,i)
    enddo;enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut

  subroutine calc_forcing(nx, ny, nz, dx, dy, dz, rho, u, v, w, p, fx, fy, fz)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: rho(nx,ny,nz), u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), p(nx,ny,nz)
    real(8), intent(out), device :: fx(nx-2,ny-2,nz-2), fy(nx-2,ny-2,nz-2), fz(nx-2,ny-2,nz-2)
  end subroutine calc_forcing
end module set

