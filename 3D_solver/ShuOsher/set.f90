module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, rhol, ul, pl
  use set_bc_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
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
    real(8) :: pi = acos(-1.d0)
    do i = 1, nx
      if (x(i) < 0.125d0) then
        Q(1,i,:,:) = rhol
        Q(2,i,:,:) = rhol * ul
        Q(3,i,:,:) = 0.d0
        Q(4,i,:,:) = 0.d0
        Q(5,i,:,:) = pl / (gamma - 1.d0) + 0.5d0 * rhol * ul**2
      else
        Q(1,i,:,:) = 1.d0 + 0.2d0 * sin(8.d0 * 2.d0 * pi * x(i))
        Q(2,i,:,:) = 0.d0
        Q(3,i,:,:) = 0.d0
        Q(4,i,:,:) = 0.d0
        Q(5,i,:,:) = 1.d0 / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: QJ(5,nx,ny,nz) ! Q / J
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    integer i, j, k, l, jc, kc, offset
    real(8), device :: Qc(5,nx)
    real(8) pw
    jc     = 4
    kc     = 4
    offset = 3

    ! inlet and outlet
    !$cuf kernel do(2) <<<*,*>>>
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        QJ(1,1,j,k)    = rhol / Jacobian(j)
        QJ(2,1,j,k)    = QJ(1,1,j,k) * ul
        QJ(3,1,j,k)    = 0.d0
        QJ(4,1,j,k)    = 0.d0
        QJ(5,1,j,k)    = (pl / (gamma - 1.d0) + 0.5d0 * rhol * ul**2) / Jacobian(j)

        QJ(1,nx,j,k)   = QJ(1,nx-1,j,k)
        QJ(2,nx,j,k)   = 0.d0
        QJ(3,nx,j,k)   = 0.d0
        QJ(4,nx,j,k)   = 0.d0
        pw             = (gamma - 1.d0) * (QJ(5,nx-1,j,k) - 0.5d0 * (QJ(2,nx-1,j,k)) / QJ(1,nx-1,j,k))
        QJ(5,nx,j,k)   = pw / (gamma - 1.d0)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do i = 1, nx
      do l = 1, 5
        Qc(l,i) = QJ(l,i,jc,kc)
    enddo;enddo

    !$cuf kernel do(4)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            QJ(l,i,j,k) = Qc(l,i)
    enddo;enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set

