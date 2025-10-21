module set
  use mod_globals, only : id_accuracy, nx, ny, nz, Lx, Ly, Lz, gamma, R
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    call set_grid_cyclic(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  end subroutine set_grid
  
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : id_accuracy, u1, rho1, u2, rho2, p, amp
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(5,nx,ny,nz)
    integer i, j, k, offset
    real(8) :: v, w, pi = acos(-1.d0)
    offset = 3
    do k = 1+offset, nz-offset
      w = amp * sin(2.d0 * pi * z(k) / Lz)
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          v = amp * sin(2.d0 * pi * x(i) / Lx)
          if (y(j) > 0.75d0 * Lx .or. y(j) < 0.25d0 * Lx) then
            Q(1,i,j,k) = rho1
            Q(2,i,j,k) = rho1 * u1
            Q(3,i,j,k) = rho1 * v
            Q(4,i,j,k) = rho1 * w
            Q(5,i,j,k) = p / (gamma - 1.d0) + 0.5d0 * rho1 * (u1**2 + v**2 + w**2)
          else
            Q(1,i,j,k) = rho2
            Q(2,i,j,k) = rho2 * u2
            Q(3,i,j,k) = rho2 * v
            Q(4,i,j,k) = rho2 * w
            Q(5,i,j,k) = p / (gamma - 1.d0) + 0.5d0 * rho2 * (u2**2 + v**2 + w**2)
          endif
    enddo;enddo;enddo
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny,nz)
    real(8), intent(inout), device :: Q(5,nx,ny,nz)
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_bc

  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set

