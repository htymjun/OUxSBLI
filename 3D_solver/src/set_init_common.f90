module set_init_common
  use mod_globals, only : gamma, R, Taw, rf
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  use set_compressible_bl
  implicit none
contains
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
    deallocate(phi_tmp, wg)
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


  subroutine set_init_tbl(nx, ny, nz, x, y, z, rand, blt0, blt, u0, p0, T0, M0, Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(in)  :: rand, blt0, blt, u0, p0, T0, M0
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    integer(4) sz
    integer(4), allocatable :: seed(:)
    real(8) :: p_wall
    real(8) :: fd, pi = acos(-1.d0)
    ! random
    real(8), allocatable :: rho(:), u(:), v(:), T(:), randum(:,:,:,:), ustd(:,:,:), vstd(:,:,:), wstd(:,:,:), Tstd(:,:,:)
    integer ir, jr, kr, nxr, nyr, nzr, n
    ! generate randum
    nxr = (nx+9) / 10
    nyr = (ny+1) / 2
    nzr = (nz+1) / 2
    allocate(rho(ny), u(ny), v(ny), T(ny))
    call calc_HD_Blasius(ny, y, blt0, u0, T0, p0, M0, rho, u, v, T)
    allocate(randum(4,nxr,nyr,nzr), ustd(nx,ny,nz), vstd(nx,ny,nz), wstd(nx,ny,nz), Tstd(nx,ny,nz))
    call random_seed(size=sz)
    allocate(seed(sz))
    call random_seed(get=seed)
    seed(1) = 48
    call random_seed(put=seed)
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
          if (y(j) <= blt) then
            ir = i / 10 + 1
            jr = j / 2  + 1
            kr = k / 2  + 1
            ustd(i,j,k) = 2.d0 * rand * u0 * randum(1,ir,jr,kr)
            vstd(i,j,k) =        rand * u0 * randum(2,ir,jr,kr)
            wstd(i,j,k) =        rand * u0 * randum(3,ir,jr,kr)
            !Tstd(i,j,k) = T0 * gamma_1 * M0**2 * randum(4,ir,jr,kr) * rand ! old version. This is not wrong.
            Tstd(i,j,k) = T0 * gamma_1 * M0**2 * randum(4,ir,jr,kr) * rand * u(j) / u0 ! fixed 2026/03/30
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
          fd = sin(pi * y(j) / (2.d0 * blt))
        else
          fd = 0.d0
        endif
        do i = 1, nx
          ustd(i,j,k) = fd * ustd(i,j,k)
          vstd(i,j,k) = fd * vstd(i,j,k)
          wstd(i,j,k) = fd * wstd(i,j,k)
          Tstd(i,j,k) = fd * Tstd(i,j,k)
    enddo;enddo;enddo
    call calc_Gaussian_filter_x(nx, ny, nz, 3, x, ustd)
    call calc_Gaussian_filter_y(nx, ny, nz, 3, y, ustd)
    call calc_Gaussian_filter_z(nx, ny, nz, 3, z, ustd)
    call calc_Gaussian_filter_x(nx, ny, nz, 3, x, vstd)
    call calc_Gaussian_filter_y(nx, ny, nz, 3, y, vstd)
    call calc_Gaussian_filter_z(nx, ny, nz, 3, z, vstd)
    call calc_Gaussian_filter_x(nx, ny, nz, 3, x, wstd)
    call calc_Gaussian_filter_y(nx, ny, nz, 3, y, wstd)
    call calc_Gaussian_filter_z(nx, ny, nz, 3, z, wstd)
    call calc_Gaussian_filter_x(nx, ny, nz, 3, x, Tstd)
    call calc_Gaussian_filter_y(nx, ny, nz, 3, y, Tstd)
    call calc_Gaussian_filter_z(nx, ny, nz, 3, z, Tstd)
    call set_bc_cyclic_x_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
    call set_bc_cyclic_z_cpu(nx, ny, nz, ustd, vstd, wstd, Tstd)
    ! add fluctuation
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,j,k,1) = p0 / (R * (T(j) + Tstd(i,j,k)))!rho(j)
          Q(i,j,k,2) = Q(i,j,k,1) * (u(j) + ustd(i,j,k))
          Q(i,j,k,3) = Q(i,j,k,1) * (v(j) + vstd(i,j,k))
          Q(i,j,k,4) = Q(i,j,k,1) * wstd(i,j,k)
          Q(i,j,k,5) = p0 * over_gamma_1 + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo
    deallocate(rho, u, v, T, randum, ustd, vstd, wstd, Tstd, seed)
    ! bottom
    Q(:,1,:,1) = Q(:,2,:,1)
    Q(:,1,:,2) = 0.d0
    Q(:,1,:,3) = 0.d0
    Q(:,1,:,4) = 0.d0
    p_wall = gamma_1 * (Q(2,2,2,5) - 0.5d0 * (Q(2,2,2,2)**2 + Q(2,2,2,3)**2 + Q(2,2,2,4)**2) / Q(2,2,2,1))
    Q(:,1,:,5) = p_wall * over_gamma_1
  end subroutine set_init_tbl
end module set_init_common

