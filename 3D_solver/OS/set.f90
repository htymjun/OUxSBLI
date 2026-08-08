module set
  use cudafor
  use mod_globals, only : rho0, p0, u0, rho2, p2, ux, uy
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  implicit none
contains
  ! Uniform Cartesian x/y grid; z uniform and periodic.
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) :: dx0, dy0, dz1
    integer :: i, j, k
    dx0 = Lx / dble(nx - 1)
    dy0 = Ly / dble(ny - 1)
    dx  = dx0
    dy  = dy0
    xc(1) = 0.0d0
    do i = 1, nx - 1
      xc(i + 1) = xc(i) + dx0
    end do
    yc(1) = 0.0d0
    do j = 1, ny - 1
      yc(j + 1) = yc(j) + dy0
    end do
    ! z: uniform and periodic, matching set_grid_cyclic6_3D's own convention
    ! (dz1 = Lz/(nz-6), zc(k) = dz1*(k-4)); no absolute origin to respect here
    ! (unlike x), so the same formula is valid for every k, ghost or interior.
    dz1 = Lz / dble(nz-6)
    do k = 1, nz
      zc(k) = dz1 * dble(k-4)
    enddo
    dz(:) = dz1
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer :: i, j, k
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,j,k,1) = rho0
          Q(i,j,k,2) = rho0 * u0
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0 ! no spanwise velocity: uniform in z by construction
          Q(i,j,k,5) = p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2
      enddo;enddo;enddo
    ! incident-shock generator: post-shock state on the top row, downstream
    ! of the shock's origin at 10% of the domain
    do k = 1, nz
      do i = int(0.1d0 * nx), nx
        Q(i,ny,k,1) = rho2
        Q(i,ny,k,2) = rho2 * ux
        Q(i,ny,k,3) = rho2 * uy
        Q(i,ny,k,4) = 0.d0
        Q(i,ny,k,5) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8) Jacobian_tmp
    integer :: i, j, k

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        Jacobian_tmp = 1.d0 / Jacobian(1,ny-1)
        if (i < int(0.1d0 * nx)) then
          QJ_1(i,ny,k) = rho0 * Jacobian_tmp
          QJ_2(i,ny,k) = rho0 * u0 * Jacobian_tmp
          QJ_3(i,ny,k) = 0.d0
          QJ_4(i,ny,k) = 0.d0
          QJ_5(i,ny,k) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
        else
          QJ_1(i,ny,k) = rho2 * Jacobian_tmp
          QJ_2(i,ny,k) = rho2 * ux * Jacobian_tmp
          QJ_3(i,ny,k) = rho2 * uy * Jacobian_tmp
          QJ_4(i,ny,k) = 0.d0
          QJ_5(i,ny,k) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
        endif
        ! Slip (reflecting) wall at the bottom
        QJ_1(i,1,k) =  QJ_1(i,2,k)
        QJ_2(i,1,k) =  QJ_2(i,2,k)
        QJ_3(i,1,k) = -QJ_3(i,2,k)
        QJ_4(i,1,k) =  QJ_4(i,2,k)
        QJ_5(i,1,k) =  QJ_5(i,2,k)
    enddo;enddo

    ! left: active Dirichlet freestream inlet (already the 2D_solver/OS
    ! pattern -- unlike 2D_solver/BL's frozen idiom, this actively rewrites
    ! i=1 every step, which is what the 3D boundary-aware x kernel needs; see
    ! 3D_solver/BL/set.f90 for why). right: zero-gradient outlet.
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        Jacobian_tmp = 1.d0 / Jacobian(1,j)
        QJ_1(1,j,k)  = rho0 * Jacobian_tmp
        QJ_2(1,j,k)  = rho0 * u0 * Jacobian_tmp
        QJ_3(1,j,k)  = 0.d0
        QJ_4(1,j,k)  = 0.d0
        QJ_5(1,j,k)  = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
        QJ_1(nx,j,k) = QJ_1(nx-1,j,k)
        QJ_2(nx,j,k) = QJ_2(nx-1,j,k)
        QJ_3(nx,j,k) = QJ_3(nx-1,j,k)
        QJ_4(nx,j,k) = QJ_4(nx-1,j,k)
        QJ_5(nx,j,k) = QJ_5(nx-1,j,k)
    enddo;enddo

    call set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
  end subroutine set_bc


  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set
