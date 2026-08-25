module set
  use cudafor
  use mod_globals, only : gamma, p0, rho0, rho1, rho2, u1, u2, du, delta_theta0, &
                          density_ratio, disturbance_intensity, peak_wavelengths_x
  use set_bc_common, only : set_bc_mut_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1

    dx1 = Lx / dble(nx - 6)
    dy1 = Ly / dble(ny - 6)
    dz1 = Lz / dble(nz - 6)
    dx(:) = dx1
    dy(:) = dy1
    dz(:) = dz1

    do i = 1, nx
      x(i) = (dble(i) - 3.5d0) * dx1
    enddo
    do j = 1, ny
      y(j) = (dble(j) - 3.5d0) * dy1 - 0.5d0 * Ly
    enddo
    do k = 1, nz
      z(k) = (dble(k) - 3.5d0) * dz1
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k, offset
    real(8) rho, u, p, eta, envelope
    real(8) up, vp, wp, kx, kz, phase1, phase2

    call ghost_width(offset)
    Q(:,:,:,:) = 0.d0

    kx = 2.d0 * acos(-1.d0) * dble(peak_wavelengths_x) / (x(nx-offset) - x(1+offset))
    kz = kx

    do k = 1 + offset, nz - offset
      do j = 1 + offset, ny - offset
        eta = y(j) / delta_theta0
        envelope = exp(-eta**2)
        do i = 1 + offset, nx - offset
          call base_state(y(j), rho, u, p)

          phase1 = kx * x(i)
          phase2 = kz * z(k)
          up = -2.d0 * disturbance_intensity * du * eta * envelope * sin(phase1) * sin(phase2)
          vp = -disturbance_intensity * du * envelope * &
               (delta_theta0 * kx * cos(phase1) + delta_theta0 * kz * sin(phase1)) * sin(phase2)
          wp =  2.d0 * disturbance_intensity * du * eta * envelope * sin(phase1) * cos(phase2)

          Q(i,j,k,1) = rho
          Q(i,j,k,2) = rho * (u + up)
          Q(i,j,k,3) = rho * vp
          Q(i,j,k,4) = rho * wp
          Q(i,j,k,5) = p / (gamma - 1.d0) + 0.5d0 * rho * ((u + up)**2 + vp**2 + wp**2)
    enddo;enddo;enddo

    call set_bc_ps2002_init(nx, ny, nz, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q_1, Q_2, Q_3, Q_4, Q_5, &
                    Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in), value            :: myrank, nx, ny, nz
    real(8), intent(in), device           :: Jacobian(nx,ny)
    real(8), intent(inout), device        :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(in), device, optional :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), &
                                             Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), &
                                             Qre_5(ny*(nz-6))

    call set_bc_ps2002(nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
  end subroutine set_bc


  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)

    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut


  subroutine ghost_width(offset)
    use mod_constant, only : id_accuracy
    integer, intent(out) :: offset

    if (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    else
      offset = 3
    endif
  end subroutine ghost_width


  subroutine base_state(y, rho, u, p)
    real(8), intent(in)  :: y
    real(8), intent(out) :: rho, u, p
    real(8) profile, lambda

    profile = tanh(-y / (2.d0 * delta_theta0))
    lambda = (density_ratio - 1.d0) / (density_ratio + 1.d0)
    rho = rho0 * (1.d0 + lambda * profile)
    u = 0.5d0 * du * profile
    p = p0
  end subroutine base_state


  subroutine far_state(is_top, rho, u, p)
    logical, intent(in)  :: is_top
    real(8), intent(out) :: rho, u, p

    if (is_top) then
      rho = rho1
      u = u1
    else
      rho = rho2
      u = u2
    endif
    p = p0
  end subroutine far_state


  subroutine set_bc_ps2002_init(nx, ny, nz, Q)
    integer, intent(in)     :: nx, ny, nz
    real(8), intent(inout)  :: Q(nx,ny,nz,5)
    integer i, j, k
    real(8) rho, u, p

    do k = 4, nz - 3
      do j = 4, ny - 3
        Q(1:3,j,k,:) = Q(nx-5:nx-3,j,k,:)
        Q(nx-2:nx,j,k,:) = Q(4:6,j,k,:)
    enddo;enddo

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1:3,:) = Q(i,j,nz-5:nz-3,:)
        Q(i,j,nz-2:nz,:) = Q(i,j,4:6,:)
    enddo;enddo

    call far_state(.false., rho, u, p)
    do k = 1, nz
      do j = 1, 3
        do i = 1, nx
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = rho * u
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0
          Q(i,j,k,5) = p / (gamma - 1.d0) + 0.5d0 * rho * u**2
    enddo;enddo;enddo

    call far_state(.true., rho, u, p)
    do k = 1, nz
      do j = ny - 2, ny
        do i = 1, nx
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = rho * u
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0
          Q(i,j,k,5) = p / (gamma - 1.d0) + 0.5d0 * rho * u**2
    enddo;enddo;enddo
  end subroutine set_bc_ps2002_init


  subroutine set_bc_ps2002(nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer i, j, k
    real(8), parameter :: e1 = p0 / (gamma - 1.d0) + 0.5d0 * rho1 * u1**2
    real(8), parameter :: e2 = p0 / (gamma - 1.d0) + 0.5d0 * rho2 * u2**2

    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz - 3
      do j = 4, ny - 3
        Q_1(1,j,k) = Q_1(nx-5,j,k); Q_1(2,j,k) = Q_1(nx-4,j,k); Q_1(3,j,k) = Q_1(nx-3,j,k)
        Q_1(nx-2,j,k) = Q_1(4,j,k); Q_1(nx-1,j,k) = Q_1(5,j,k); Q_1(nx,j,k) = Q_1(6,j,k)
        Q_2(1,j,k) = Q_2(nx-5,j,k); Q_2(2,j,k) = Q_2(nx-4,j,k); Q_2(3,j,k) = Q_2(nx-3,j,k)
        Q_2(nx-2,j,k) = Q_2(4,j,k); Q_2(nx-1,j,k) = Q_2(5,j,k); Q_2(nx,j,k) = Q_2(6,j,k)
        Q_3(1,j,k) = Q_3(nx-5,j,k); Q_3(2,j,k) = Q_3(nx-4,j,k); Q_3(3,j,k) = Q_3(nx-3,j,k)
        Q_3(nx-2,j,k) = Q_3(4,j,k); Q_3(nx-1,j,k) = Q_3(5,j,k); Q_3(nx,j,k) = Q_3(6,j,k)
        Q_4(1,j,k) = Q_4(nx-5,j,k); Q_4(2,j,k) = Q_4(nx-4,j,k); Q_4(3,j,k) = Q_4(nx-3,j,k)
        Q_4(nx-2,j,k) = Q_4(4,j,k); Q_4(nx-1,j,k) = Q_4(5,j,k); Q_4(nx,j,k) = Q_4(6,j,k)
        Q_5(1,j,k) = Q_5(nx-5,j,k); Q_5(2,j,k) = Q_5(nx-4,j,k); Q_5(3,j,k) = Q_5(nx-3,j,k)
        Q_5(nx-2,j,k) = Q_5(4,j,k); Q_5(nx-1,j,k) = Q_5(5,j,k); Q_5(nx,j,k) = Q_5(6,j,k)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        Q_1(i,j,1) = Q_1(i,j,nz-5); Q_1(i,j,2) = Q_1(i,j,nz-4); Q_1(i,j,3) = Q_1(i,j,nz-3)
        Q_1(i,j,nz-2) = Q_1(i,j,4); Q_1(i,j,nz-1) = Q_1(i,j,5); Q_1(i,j,nz) = Q_1(i,j,6)
        Q_2(i,j,1) = Q_2(i,j,nz-5); Q_2(i,j,2) = Q_2(i,j,nz-4); Q_2(i,j,3) = Q_2(i,j,nz-3)
        Q_2(i,j,nz-2) = Q_2(i,j,4); Q_2(i,j,nz-1) = Q_2(i,j,5); Q_2(i,j,nz) = Q_2(i,j,6)
        Q_3(i,j,1) = Q_3(i,j,nz-5); Q_3(i,j,2) = Q_3(i,j,nz-4); Q_3(i,j,3) = Q_3(i,j,nz-3)
        Q_3(i,j,nz-2) = Q_3(i,j,4); Q_3(i,j,nz-1) = Q_3(i,j,5); Q_3(i,j,nz) = Q_3(i,j,6)
        Q_4(i,j,1) = Q_4(i,j,nz-5); Q_4(i,j,2) = Q_4(i,j,nz-4); Q_4(i,j,3) = Q_4(i,j,nz-3)
        Q_4(i,j,nz-2) = Q_4(i,j,4); Q_4(i,j,nz-1) = Q_4(i,j,5); Q_4(i,j,nz) = Q_4(i,j,6)
        Q_5(i,j,1) = Q_5(i,j,nz-5); Q_5(i,j,2) = Q_5(i,j,nz-4); Q_5(i,j,3) = Q_5(i,j,nz-3)
        Q_5(i,j,nz-2) = Q_5(i,j,4); Q_5(i,j,nz-1) = Q_5(i,j,5); Q_5(i,j,nz) = Q_5(i,j,6)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        do j = 1, 3
          Q_1(i,j,k) = rho2
          Q_2(i,j,k) = rho2 * u2
          Q_3(i,j,k) = 0.d0
          Q_4(i,j,k) = 0.d0
          Q_5(i,j,k) = e2
        enddo
        do j = ny - 2, ny
          Q_1(i,j,k) = rho1
          Q_2(i,j,k) = rho1 * u1
          Q_3(i,j,k) = 0.d0
          Q_4(i,j,k) = 0.d0
          Q_5(i,j,k) = e1
        enddo
    enddo;enddo
  end subroutine set_bc_ps2002
end module set
