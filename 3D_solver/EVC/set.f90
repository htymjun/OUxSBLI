module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R, theta
  use set_bc_common
  use set_coordinate
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    call set_grid_cyclic(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : M0, rho0, p0, T0, u0, Rc, beta
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    real(8) xc, yc, ex, T, rho, u, v, du, dv, p
    real(8) :: Cp = R * gamma / (gamma - 1.d0)
    xc = x(int(nx/2))
    yc = y(int(ny/2))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ex  = exp(-0.5d0 * ((x(i) - xc)**2 + (y(j) - yc)**2) / (Rc**2))
          T   = T0 - 0.5d0 * (u0 * beta)**2 / Cp * ex**2
          rho = rho0 * (T / T0)**(1.d0  / (gamma - 1.d0))
          du  = -u0 * beta * (y(j) - yc) / Rc * ex
          dv  =  u0 * beta * (x(i) - xc) / Rc * ex
          u   = u0 * cos(theta) + du * cos(theta) - dv * sin(theta)
          v   = u0 * sin(theta) + du * sin(theta) + dv * cos(theta)
          p   = rho * R * T
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = rho * u
          Q(i,j,k,3) = rho * v
          ! no spanwise velocity: the vortex is uniform in z by construction
          Q(i,j,k,4) = 0.d0
          ! p / (gamma - 1) + 0.5 * (rhou ** 2 + rhov ** 2 ) / rho
          Q(i,j,k,5) = p / (gamma - 1.d0) + 0.5d0 * rho * (u**2 + v**2)
    enddo;enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)
    use mod_constant, only : id_accuracy
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
  end subroutine set_bc


  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set
