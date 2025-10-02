module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, RHO0, M0, V0, p0, T, dtn
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    call set_grid_cyclic(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  end subroutine set_grid
  
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(5,nx,ny,nz)
    integer i, j, k, offset
    if (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    endif
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          ! rho
          Q(1,i,j,k) =  RHO0
          ! rho u
          Q(2,i,j,k) =  RHO0 * M0 * sin(x(i)) * cos(y(j)) * cos(z(k))
          ! rho v
          Q(3,i,j,k) = -RHO0 * M0 * cos(x(i)) * sin(y(j)) * cos(z(k))
          ! rho w0
          Q(4,i,j,k) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(5,i,j,k) = (1.d0/gamma+0.0625d0*RHO0*(M0**2)*(cos(2.d0*x(i))+cos(2.d0*y(j)))*(cos(2.d0*z(k))+2.d0))&
                       / (gamma - 1.d0) + 0.5d0 * (Q(2,i,j,k)**2 + Q(3,i,j,k)**2 + Q(4,i,j,k)**2) / Q(1,i,j,k)
    enddo;enddo;enddo
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q(5,nx,ny,nz)
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    integer i, j, k, l
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut

  subroutine calc_forcing(nx, ny, nz, dx, dy, dz, Q, fx, fy, fz)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: fx(nx-2,ny-2,nz-2), fy(nx-2,ny-2,nz-2), fz(nx-2,ny-2,nz-2)
  end subroutine calc_forcing
end module set

