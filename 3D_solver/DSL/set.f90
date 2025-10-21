module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R
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
    use mod_globals, only : id_accuracy, rho0, u0, pi, d1, d2, M0
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(5,nx,ny,nz)
    real(8) :: Cp = R * gamma / (gamma - 1.d0), p = rho0 * u0**2 / (gamma * M0**2)
    integer i, j, k, offset
    if (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    endif
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (y(j) <= pi) then
            Q(1,i,j,k) = rho0
            Q(2,i,j,k) = rho0 * u0 * tanh((y(j) - 0.5d0 * pi) / d1)
            Q(3,i,j,k) = rho0 * u0 * d2 * sin(x(i))
            Q(4,i,j,k) = rho0 * u0 * d2 * sin(z(k))
            Q(5,i,j,k) = p / (gamma - 1.d0) + 0.5d0 * (Q(2,i,j,k)**2 + Q(3,i,j,k)**2 + Q(4,i,j,k)**2) / Q(1,i,j,k)
          else
            Q(1,i,j,k) = rho0
            Q(2,i,j,k) = rho0 * u0 * tanh((1.5d0 * pi - y(j)) / d1)
            Q(3,i,j,k) = rho0 * u0 * d2 * sin(x(i))
            Q(4,i,j,k) = rho0 * u0 * d2 * sin(z(k))
            Q(5,i,j,k) = p / (gamma - 1.d0) + 0.5d0 * (Q(2,i,j,k)**2 + Q(3,i,j,k)**2 + Q(4,i,j,k)**2) / Q(1,i,j,k)
          endif
    enddo;enddo;enddo
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
    use mod_globals, only : id_accuracy
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: Qre(ny*(nz-6)*5)
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_bc

  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set

