module set_init_common
  use mod_globals, only : gamma, R
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  implicit none
contains
  subroutine calc_diff_eq(nx, ny, nz, x, y, z, L, endT, phi)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(in)    :: x(nx), y(ny), z(nz), L, endT
    real(8), intent(inout) :: phi(nx,ny,nz)
    real(8), allocatable :: phi_new(:,:,:)
    real(8) D, dt, dx, dy, dz, dx2, dy2, dz2, d2x, d2y, d2z
    real(8) :: pi = acos(-1.d0)
    integer itr, endItr, i, j, k
    allocate(phi_new(nx,ny,nz))
    phi_new = phi
    D   = L**2 / (800.d0 * pi)
    dx2 = (-x(1) + x(2))**2
    dy2 = (-y(1) + y(2))**2
    dz2 = (-z(1) + z(2))**2
    dt  = 1.d0 / (4.d0 * D * (1.d0 / dx2 + 1.d0 / dy2 + 1.d0 / dz2))
    endItr = int(endT / dt)
    do itr = 1, endItr
      do k = 2, nz-1
        dz2 = (0.5d0 * (-z(k-1) + z(k+1)))**2
        do j = 2, ny-1
          dy2 = (0.5d0 * (-y(j-1) + y(j+1)))**2
          do i = 2, nx-1
            dx2 = (0.5d0 * (-x(i-1) + x(i+1)))**2
            d2x = (phi(i-1,j,k) - 2.d0 * phi(i,j,k) + phi(i+1,j,k)) / dx2
            d2y = (phi(i,j-1,k) - 2.d0 * phi(i,j,k) + phi(i,j+1,k)) / dy2
            d2z = (phi(i,j,k-1) - 2.d0 * phi(i,j,k) + phi(i,j,k+1)) / dz2
            phi_new(i,j,k) = phi(i,j,k) + dt * D * (d2x + d2y + d2z)
      enddo;enddo;enddo
      ! boundary conditions
      phi_new(1,:,:)  = phi_new(2,:,:)
      phi_new(nx,:,:) = phi_new(nx-1,:,:)
      phi_new(:,1,:)  = 0.d0
      phi_new(:,ny,:) = phi_new(:,ny-1,:)
      phi_new(:,:,1)  = phi_new(:,:,nz-1)
      phi_new(:,:,nz) = phi_new(:,:,2)
      ! update
      phi = phi_new
    enddo
    deallocate(phi_new)
  end subroutine calc_diff_eq


  subroutine calc_Gaussian_filter_x(nx, ny, nz, n, x, phi)
    integer, intent(in)    :: nx, ny, nz, n
    real(8), intent(in)    :: x(nx)
    real(8), intent(inout) :: phi(nx,ny,nz)
    real(8) sigma, d, tmp
    integer :: i, j, k, ix
    real(8), allocatable :: phi_tmp(:,:,:), wg(:)
    allocate(phi_tmp(nx,ny,nz))
    phi_tmp = phi
    d = -x(1) + x(2)
    sigma = dble(n) * d / 3.d0
    allocate(wg(2*n+1))
    do ix = -n, n
      wg(ix+n+1) = exp(-((ix*d)**2)/(2.d0*sigma**2))
    enddo
    wg = wg / sum(wg)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = 0.d0
          do ix = -n, n
            if (i+ix >= 1 .and. i+ix <= nx) then
              tmp = tmp + phi_tmp(i+ix,j,k)*wg(ix+n+1)
            else
              tmp = tmp + phi_tmp(i,j,k)
            endif
          enddo
          phi(i,j,k) = tmp
    enddo;enddo;enddo
    deallocate(phi_tmp)
  end subroutine calc_Gaussian_filter_x


  subroutine calc_Gaussian_filter_y(nx, ny, nz, n, y, phi)
    integer, intent(in)    :: nx, ny, nz, n
    real(8), intent(in)    :: y(ny)
    real(8), intent(inout) :: phi(nx,ny,nz)
    real(8) sigma, d, tmp
    integer :: i, j, k, iy
    real(8), allocatable :: phi_tmp(:,:,:), wg(:)
    allocate(phi_tmp(nx,ny,nz))
    phi_tmp = phi
    d = -y(ny-1) + y(ny)
    sigma = dble(n) * d / 3.d0
    allocate(wg(2*n+1))
    do iy = -n, n
      wg(iy+n+1) = exp(-((iy*d)**2)/(2.d0*sigma**2))
    enddo
    wg = wg / sum(wg)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = 0.d0
          do iy = -n, n
            if (j+iy >= 1 .and. j+iy <= ny) then
              tmp = tmp + phi_tmp(i,j+iy,k)*wg(iy+n+1)
            else
              tmp = tmp + phi_tmp(i,j,k)
            endif
          enddo
          phi(i,j,k) = tmp
    enddo;enddo;enddo
    deallocate(phi_tmp, wg)
  end subroutine calc_Gaussian_filter_y


  subroutine calc_Gaussian_filter_z(nx, ny, nz, n, z, phi)
    integer, intent(in)    :: nx, ny, nz, n
    real(8), intent(in)    :: z(nz)
    real(8), intent(inout) :: phi(nx,ny,nz)
    real(8) sigma, d, tmp
    integer :: i, j, k, iz
    real(8), allocatable :: phi_tmp(:,:,:), wg(:)
    allocate(phi_tmp(nx,ny,nz))
    phi_tmp = phi
    d = -z(1) + z(2)
    sigma = dble(n) * d / 3.d0
    allocate(wg(2*n+1))
    do iz = -n, n
      wg(iz+n+1) = exp(-((iz*d)**2)/(2.d0*sigma**2))
    enddo
    wg = wg / sum(wg)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = 0.d0
          do iz = -n, n
            if (k+iz >= 1 .and. k+iz <= nz) then
              tmp = tmp + phi_tmp(i,j,k+iz)*wg(iz+n+1)
            else
              tmp = tmp + phi_tmp(i,j,k)
            endif
          enddo
          phi(i,j,k) = tmp
    enddo;enddo;enddo
    deallocate(phi_tmp, wg)
  end subroutine calc_Gaussian_filter_z


  subroutine set_bc_cyclic_x_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(inout) :: ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz)
    integer j, k
    do k = 1, nz
      do j = 1, ny
        ustd(1,j,k)    = ustd(nx-5,j,k)
        ustd(2,j,k)    = ustd(nx-4,j,k)
        ustd(3,j,k)    = ustd(nx-3,j,k)
        ustd(nx-2,j,k) = ustd(4,j,k)
        ustd(nx-1,j,k) = ustd(5,j,k)
        ustd(nx,j,k)   = ustd(6,j,k)
        vstd(1,j,k)    = vstd(nx-5,j,k)
        vstd(2,j,k)    = vstd(nx-4,j,k)
        vstd(3,j,k)    = vstd(nx-3,j,k)
        vstd(nx-2,j,k) = vstd(4,j,k)
        vstd(nx-1,j,k) = vstd(5,j,k)
        vstd(nx,j,k)   = vstd(6,j,k)
        wstd(1,j,k)    = wstd(nx-5,j,k)
        wstd(2,j,k)    = wstd(nx-4,j,k)
        wstd(3,j,k)    = wstd(nx-3,j,k)
        wstd(nx-2,j,k) = wstd(4,j,k)
        wstd(nx-1,j,k) = wstd(5,j,k)
        wstd(nx,j,k)   = wstd(6,j,k)
        Tstd(1,j,k)    = Tstd(nx-5,j,k)
        Tstd(2,j,k)    = Tstd(nx-4,j,k)
        Tstd(3,j,k)    = Tstd(nx-3,j,k)
        Tstd(nx-2,j,k) = Tstd(4,j,k)
        Tstd(nx-1,j,k) = Tstd(5,j,k)
        Tstd(nx,j,k)   = Tstd(6,j,k)
    enddo;enddo
  end subroutine set_bc_cyclic_x_cpu


  subroutine set_bc_cyclic_z_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(inout) :: ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz)
    integer i, j
    do j = 1, ny
      do i = 1, nx
        ustd(i,j,1)    = ustd(i,j,nz-5)
        ustd(i,j,2)    = ustd(i,j,nz-4)
        ustd(i,j,3)    = ustd(i,j,nz-3)
        ustd(i,j,nz-2) = ustd(i,j,4)
        ustd(i,j,nz-1) = ustd(i,j,5)
        ustd(i,j,nz)   = ustd(i,j,6)
        vstd(i,j,1)    = vstd(i,j,nz-5)
        vstd(i,j,2)    = vstd(i,j,nz-4)
        vstd(i,j,3)    = vstd(i,j,nz-3)
        vstd(i,j,nz-2) = vstd(i,j,4)
        vstd(i,j,nz-1) = vstd(i,j,5)
        vstd(i,j,nz)   = vstd(i,j,6)
        wstd(i,j,1)    = wstd(i,j,nz-5)
        wstd(i,j,2)    = wstd(i,j,nz-4)
        wstd(i,j,3)    = wstd(i,j,nz-3)
        wstd(i,j,nz-2) = wstd(i,j,4)
        wstd(i,j,nz-1) = wstd(i,j,5)
        wstd(i,j,nz)   = wstd(i,j,6)
        Tstd(i,j,1)    = Tstd(i,j,nz-5)
        Tstd(i,j,2)    = Tstd(i,j,nz-4)
        Tstd(i,j,3)    = Tstd(i,j,nz-3)
        Tstd(i,j,nz-2) = Tstd(i,j,4)
        Tstd(i,j,nz-1) = Tstd(i,j,5)
        Tstd(i,j,nz)   = Tstd(i,j,6)
    enddo;enddo
  end subroutine set_bc_cyclic_z_cpu


  subroutine calc_rms(nx, ny, nz, phi, rms)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: phi(nx,ny,nz)
    real(8), intent(out) :: rms
    real(8), allocatable :: phi2(:,:,:)
    real(8) phi_mean, phi2_mean
    allocate(phi2(nx,ny,nz))
    phi2 = phi * phi
    phi2_mean = sum(phi2) / size(phi2)
    phi_mean  = sum(phi)  / size(phi)
    rms = sqrt(phi2_mean - phi_mean * phi_mean)
    deallocate(phi2)
  end subroutine calc_rms


  subroutine set_init_tbl(nx, ny, nz, xs, ys, zs, rand, blt0, blt, rf, u0, p0, T0, M0, Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(in)  :: rand, blt0, blt, rf, u0, p0, T0, M0
    real(8), intent(out) :: Q(5,nx,ny,nz)
    integer i, j, k
    real(8) :: eta, rho, u, v, w, T, Tw, Taw, p_wall
    real(8) :: fd, pi = acos(-1.d0)
    ! random
    real(8), allocatable :: randum(:,:,:,:), ustd(:,:,:), vstd(:,:,:), wstd(:,:,:), Tstd(:,:,:)
    integer ir, jr, kr, nxr, nyr, nzr, n
    ! generate randum
    nxr = (nx+9) / 10
    nyr = (ny+1) / 2
    nzr = (nz+1) / 2
    allocate(randum(4,nxr,nyr,nzr), ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz))
    ! generate random number
    do k = 1, nzr
      do j = 1, nyr
        do i = 1, nxr
          call random_number(randum(1,i,j,k))
          call random_number(randum(2,i,j,k))
          call random_number(randum(3,i,j,k))
          call random_number(randum(4,i,j,k))
          randum(1,i,j,k) = 2.d0 * randum(1,i,j,k) - 1.d0
          randum(2,i,j,k) = 2.d0 * randum(2,i,j,k) - 1.d0
          randum(3,i,j,k) = 2.d0 * randum(3,i,j,k) - 1.d0
          randum(4,i,j,k) = 2.d0 * randum(4,i,j,k) - 1.d0
    enddo;enddo;enddo
    ! generate fluctuation
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (ys(j) <= blt) then
            ir = i / 10 + 1
            jr = j / 2  + 1
            kr = k / 2  + 1
            ustd(i,j,k) = 2.d0 * rand * u0 * randum(1,ir,jr,kr)
            vstd(i,j,k) =        rand * u0 * randum(2,ir,jr,kr)
            wstd(i,j,k) =        rand * u0 * randum(3,ir,jr,kr)
            Tstd(i,j,k) = T0 * gamma_1 * M0**2 * randum(4,ir,jr,kr) * rand
          else
            ustd(i,j,k) = 0.d0
            vstd(i,j,k) = 0.d0
            wstd(i,j,k) = 0.d0
            Tstd(i,j,k) = 0.d0
          endif
    enddo;enddo;enddo
    ! damping
    do k = 1, nz
      do j = 1, ny
        if (j >= 10) then
          fd = sin(pi * ys(j) / (2.d0 * blt))
        else
          fd = 0.d0
        endif
        do i = 1, nx
          ustd(i,j,k) = fd * ustd(i,j,k)
          vstd(i,j,k) = fd * vstd(i,j,k)
          wstd(i,j,k) = fd * wstd(i,j,k)
          Tstd(i,j,k) = fd * Tstd(i,j,k)
    enddo;enddo;enddo
    block
      real(8) urms1, vrms1, wrms1, Trms1, urms2, vrms2, wrms2, Trms2
      call calc_rms(nx, ny, nz, ustd, urms1)
      call calc_rms(nx, ny, nz, vstd, vrms1)
      call calc_rms(nx, ny, nz, wstd, wrms1)
      call calc_rms(nx, ny, nz, Tstd, Trms1)
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, ustd)
      call calc_Gaussian_filter_y(nx, ny, nz, 3, ys, ustd)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, ustd)
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, vstd)
      call calc_Gaussian_filter_y(nx, ny, nz, 3, ys, vstd)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, vstd)
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, wstd)
      call calc_Gaussian_filter_y(nx, ny, nz, 3, ys, wstd)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, wstd)
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, Tstd)
      call calc_Gaussian_filter_y(nx, ny, nz, 3, ys, Tstd)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, Tstd)
      call set_bc_cyclic_x_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
      call set_bc_cyclic_z_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
      call calc_rms(nx, ny, nz, ustd, urms2)
      call calc_rms(nx, ny, nz, vstd, vrms2)
      call calc_rms(nx, ny, nz, wstd, wrms2)
      call calc_rms(nx, ny, nz, Tstd, Trms2)
      ustd = ustd * urms1 / urms2
      vstd = vstd * vrms1 / vrms2
      wstd = wstd * wrms1 / wrms2
      Tstd = Tstd * Trms1 / Trms2
    end block
    ! add fluctuation
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          eta = 5.d0 * ys(j) / blt0
          u   = min(u0, u0 * (0.0015d0 * eta**4 - 0.0181d0 * eta**3 + 0.029d0 * eta**2 + 0.3192 * eta + 0.0003d0))
          v   = 0.d0
          Taw = T0 * (1.d0 + rf * 0.5d0 * gamma_1 * M0**2)
          Tw  = Taw
          T   = Tw + (Taw - Tw) * u / u0 - rf * u**2 / (2.d0 * Cp)
          u   = u + ustd(i,j,k)
          v   = v + vstd(i,j,k)
          w   = wstd(i,j,k)
          T   = T + Tstd(i,j,k)
          rho = p0 / (R * T)
          Q(1,i,j,k) = rho
          Q(2,i,j,k) = Q(1,i,j,k) * u
          Q(3,i,j,k) = Q(1,i,j,k) * v
          Q(4,i,j,k) = Q(1,i,j,k) * w
          Q(5,i,j,k) = p0 * over_gamma_1 + 0.5d0 * (Q(2,i,j,k)**2 + Q(3,i,j,k)**2 + Q(4,i,j,k)**2) / Q(1,i,j,k)
    enddo;enddo;enddo
    deallocate(randum, ustd, vstd, wstd, Tstd)
    ! bottom
    Q(1,:,1,:) = Q(1,:,2,:)
    Q(2,:,1,:) = 0.d0
    Q(3,:,1,:) = 0.d0
    Q(4,:,1,:) = 0.d0
    p_wall = gamma_1 * (Q(5,2,2,2) - 0.5d0 * (Q(2,2,2,2)**2 + Q(3,2,2,2)**2 + Q(4,2,2,2)**2) / Q(1,2,2,2))
    Q(5,:,1,:) = p_wall * over_gamma_1
  end subroutine set_init_tbl
end module set_init_common

