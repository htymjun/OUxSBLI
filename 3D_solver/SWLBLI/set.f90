module set
  use cudafor
  use mod_globals, only : gamma, rho0, u0, p0, rho2, p2, ux, uy, beta, Lx, Ly, Lz, x_in, Xsh, i_LE, nx, dt
  use mod_constant, only : gamma_1, over_gamma_1
  use set_bc_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dz1
    real(8) tanh_s, yi
    real(8), parameter :: s = 2.4d0 !1.6 ! tanh wall-clustering stretch
    dx1 = Lx / dble(nx-1)
    dz1 = Lz / dble(nz-1)

    x(1) = x_in
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    tanh_s = tanh(s)
    do j = 1, ny
      yi = dble(j-1) / dble(ny-1)
      y(j) = Ly * (1 - tanh(s * (1.d0 - yi)) / tanh_s)
    enddo
    do j = 1, ny-1
      dy(j) = y(j+1) - y(j)
    enddo

    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
    z(:) = z(:) - 0.5d0 * Lz
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,j,k,1) = rho0
          Q(i,j,k,2) = rho0 * u0
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0
          Q(i,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
        enddo
        ! The inlet column i=1 is never updated afterwards (interior-only RK kernels,
        ! no inlet BC in set_bc), so it stays frozen at these values: below the
        ! incident-shock trace it holds the freestream, above it the post-shock
        ! state, acting as the oblique-shock generator.
        ! if (ys(j) > (Xsh - x_in) * dtan(beta)) then
        !   Q(1,j,k,1) = rho2
        !   Q(1,j,k,2) = rho2 * ux
        !   Q(1,j,k,3) = rho2 * uy
        !   Q(1,j,k,4) = 0.d0
        !   Q(1,j,k,5) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
        ! endif
    enddo;enddo
    do k = 1, nz
      do i = 1, nx    
        if (Ly / dtan(beta) + Lx / dble(nx-1) * dble(i-1) + x_in >= Xsh .and. Lx / dble(nx-1) * dble(i-1) + x_in < 0.09d0) then
          Q(i,ny,k,1) = rho2
          Q(i,ny,k,2) = rho2 * ux
          Q(i,ny,k,3) = rho2 * uy
          Q(i,ny,k,4) = 0.d0
          Q(i,ny,k,5) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
        endif
    enddo;enddo
    do k = 1, nz
      do i = i_LE, nx
        Q(i,1,k,1) = rho0
        Q(i,1,k,2) = 0.d0
        Q(i,1,k,3) = 0.d0
        Q(i,1,k,4) = 0.d0
        Q(i,1,k,5) = Q(i,2,k,5)
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, xs, zs, t_now, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), zs(nz)
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8) Jacobian_tmp
    integer i, j, k
    real(8) :: p_wall
    real(8) :: rhoin, uin, vin, pin, cin
    real(8) :: Rp, Rm, vb, cb, sb, ub, rhob, pb
    ! exterior (post-shock) sound speed for the top Riemann state
    real(8), parameter :: c_ext = sqrt(gamma * p2 / rho2)
    real(8), parameter :: c0    = sqrt(gamma * p0 / rho0)
    ! region of blowing (laminar-to-turbulent transition)
    real(8), parameter :: A = 0.02d0
    real(8) :: f_x, g_z, h_t
    real(8), parameter :: x_a = 0.01d0 !beginning of the blowing and suction zone
    real(8), parameter :: x_b = 0.02d0 !end of the blowing and suction zone
    integer :: l = 1, l_max = 10, m = 1, m_max = 10
    real(8), parameter :: beta_force = 75000 !(Hz)
    real(8) :: Z_l(l_max), T_m(m_max)
    real(8) :: phi_l(l_max), phi_m(m_max)
    real(8), device ::phi_l_gpu(l_max), phi_m_gpu(m_max)
    real(8), device :: Z_l_gpu(l_max), T_m_gpu(m_max)
    real(8), device :: xs_gpu(nx), zs_gpu(nz)
    real(8) r, theta, t_now, t_now_local
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real(8), parameter :: ratio = 1.25d0

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 2, ny-1
        ! outlet
        QJ_1(nx,j,k) = QJ_1(nx-1,j,k); QJ_2(nx,j,k) = QJ_2(nx-1,j,k); QJ_3(nx,j,k) = QJ_3(nx-1,j,k)
        QJ_4(nx,j,k) = QJ_4(nx-1,j,k); QJ_5(nx,j,k) = QJ_5(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        if (Ly / dtan(beta) + Lx / dble(nx-1) * dble(i-1) + x_in < Xsh) then
          ! ---- pre-shock state (j = ny) ----
          ! ---- interior state (j = ny-1) ----
          rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
          uin   = QJ_2(i,ny-1,k) / QJ_1(i,ny-1,k)
          vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
          pin   = gamma_1 * ( QJ_5(i,ny-1,k) * Jacobian(i,ny-1) &
                - 0.5d0 * rhoin * (uin**2 + vin**2) )
          cin   = sqrt(gamma * pin / rhoin)

          ! ---- Riemann invariant (normal = y direction) ----
          Rp = vin + 2.d0 * cin * over_gamma_1
          Rm = -2.d0 * c0 * over_gamma_1 !uy - 2.d0 * c_ext * over_gamma_1
          vb = 0.5d0 * (Rp + Rm)
          cb = 0.25d0 * gamma_1 * (Rp - Rm)

          if (vb >= 0.d0) then
            sb = pin / rhoin**gamma
            ub = uin
          else
            sb = p0 / rho0**gamma !sb = p2 / rho2**gamma
            ub = u0 !ub = ux
          endif

          rhob = (cb**2 / (gamma * sb)) ** over_gamma_1
          pb   = sb * rhob**gamma

          Jacobian_tmp = 1.d0 / Jacobian(i,ny)
          QJ_1(i,ny,k) = rhob * Jacobian_tmp
          QJ_2(i,ny,k) = rhob * ub * Jacobian_tmp
          QJ_3(i,ny,k) = rhob * vb * Jacobian_tmp
          QJ_4(i,ny,k) = 0.d0
          QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
        
        else if (Ly / dtan(beta) + Lx / dble(nx-1) * dble(i-1) + x_in >= Xsh .and. Lx / dble(nx-1) * dble(i-1) + x_in < 0.09d0) then
          ! ---- post-shock state (j = ny) ----
          Jacobian_tmp = 1.d0 / Jacobian(i,ny)
          QJ_1(i,ny,k) = rho2 * Jacobian_tmp
          QJ_2(i,ny,k) = rho2 * ux * Jacobian_tmp
          QJ_3(i,ny,k) = rho2 * uy * Jacobian_tmp
          QJ_4(i,ny,k) = 0.d0
          QJ_5(i,ny,k) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
        
        else
          ! ---- Neumann ----
          QJ_1(i,ny,k) = QJ_1(i,ny-1,k); QJ_2(i,ny,k) = QJ_2(i,ny-1,k); QJ_3(i,ny,k) = QJ_3(i,ny-1,k)
          QJ_4(i,ny,k) = QJ_4(i,ny-1,k); QJ_5(i,ny,k) = QJ_5(i,ny-1,k)
        endif
    enddo;enddo
    
    !prepareing for the region of blowing (laminar-to-turbulent transition)
    call random_number(phi_l)
    call random_number(phi_m)

    r = 1.d0 / ratio ! r = 0.8
    Z_l(1) = (1.0d0 - r) / (1.0d0 - r**dble(l_max))
    do l = 2, l_max
      Z_l(l) = Z_l(l-1) * r
    enddo
    T_m(1) = (1.0d0 - r) / (1.0d0 - r**dble(m_max))
    do m = 2, m_max
      T_m(m) = T_m(m-1) * r
    enddo
    t_now_local = dt * dble(t_now)

    phi_l_gpu = phi_l
    phi_m_gpu = phi_m
    Z_l_gpu   = Z_l
    T_m_gpu   = T_m
    xs_gpu    = xs
    zs_gpu    = zs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        if(i < i_LE) then
          !Neumann
          QJ_1(i,1,k) = QJ_1(i,2,k)
          QJ_2(i,1,k) = QJ_2(i,2,k)
          QJ_3(i,1,k) = 0.d0
          QJ_4(i,1,k) = 0.d0
          QJ_5(i,1,k) = QJ_5(i,2,k) - 0.5d0 * QJ_3(i,2,k)**2 / QJ_1(i,2,k)
        else
          !NoSlip
          QJ_1(i,1,k) = QJ_1(i,2,k)
          QJ_2(i,1,k) = 0.d0
          QJ_3(i,1,k) = 0.d0
          QJ_4(i,1,k) = 0.d0
          p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
          QJ_5(i,1,k) = p_wall * over_gamma_1
        endif

        ! region of blowing (laminar-to-turbulent transition)
        if(xs_gpu(i) >= x_a .and. xs_gpu(i) <= x_b) then
          g_z = 0.d0; h_t = 0.d0
          do l = 1, l_max
            g_z = g_z + Z_l_gpu(l) * sin(2.d0 * pi * dble(l) * (zs_gpu(k) / Lz + phi_l_gpu(l)))
          enddo
          
          do m = 1, m_max
            h_t = h_t + T_m_gpu(m) * sin(beta_force * t_now_local + 2.d0 * pi * phi_m_gpu(m))
          enddo

          theta = 2.d0 * pi * (xs_gpu(i) - x_a) / (x_b - x_a)
          f_x = 4.d0 * sin(theta) * (1.d0 -cos(theta)) / sqrt(27.d0)

          QJ_3(i,1,k) = QJ_1(i,1,k) * A * u0 * f_x * g_z * h_t
        endif
    enddo;enddo

    

    call set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
  end subroutine set_bc


  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(inout), device  :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        ! inlet
        mut(1,j,k)  = mut(2,j,k)
        qc2(1,j,k)  = qc2(2,j,k)
        ! outlet
        mut(nx,j,k) = mut(nx-1,j,k)
        qc2(nx,j,k) = qc2(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        ! wall
        mut(i,1,k) = 0.d0
        qc2(i,1,k) = 0.d0
        ! top
        mut(i,ny,k) = mut(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,ny-1,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        ! span
        mut(i,j,1)    = mut(i,j,nz-5)
        mut(i,j,2)    = mut(i,j,nz-4)
        mut(i,j,3)    = mut(i,j,nz-3)
        mut(i,j,nz-2) = mut(i,j,4)
        mut(i,j,nz-1) = mut(i,j,5)
        mut(i,j,nz)   = mut(i,j,6)
        qc2(i,j,1)    = qc2(i,j,nz-5)
        qc2(i,j,2)    = qc2(i,j,nz-4)
        qc2(i,j,3)    = qc2(i,j,nz-3)
        qc2(i,j,nz-2) = qc2(i,j,4)
        qc2(i,j,nz-1) = qc2(i,j,5)
        qc2(i,j,nz)   = qc2(i,j,6)
    enddo;enddo
  end subroutine set_bc_mut
end module set
