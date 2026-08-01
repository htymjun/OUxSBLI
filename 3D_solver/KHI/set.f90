module set
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, gamma, R
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
    use mod_globals, only : u1, rho1, u2, rho2, p, amp
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k, offset
    real(8) :: v, w, pi = acos(-1.d0)
    if (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    endif
    do k = 1+offset, nz-offset
      w = amp * u1 * sin(2.d0 * pi * z(k) / Lz)
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          v = amp * u1 * sin(2.d0 * pi * x(i) / Lx)
          if (y(j) > 0.75d0 * Lx .or. y(j) < 0.25d0 * Lx) then
            Q(i,j,k,1) = rho1
            Q(i,j,k,2) = rho1 * u1
            Q(i,j,k,3) = rho1 * v
            Q(i,j,k,4) = rho1 * w
            Q(i,j,k,5) = p / (gamma - 1.d0) + 0.5d0 * rho1 * (u1**2 + v**2 + w**2)
          else
            Q(i,j,k,1) = rho2
            Q(i,j,k,2) = rho2 * u2
            Q(i,j,k,3) = rho2 * v
            Q(i,j,k,4) = rho2 * w
            Q(i,j,k,5) = p / (gamma - 1.d0) + 0.5d0 * rho2 * (u2**2 + v**2 + w**2)
          endif
    enddo;enddo;enddo
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q_1, Q_2, Q_3, Q_4, Q_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    use mod_constant, only : id_accuracy
    integer, intent(in), value            :: myrank, nx, ny, nz
    real(8), intent(in), device           :: Jacobian(nx,ny,nz)
    real(8), intent(inout), device        :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(in), device, optional :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    call set_bc_cyclic(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
  end subroutine set_bc


  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set

