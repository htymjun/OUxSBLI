module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, Pr, rho0, u0, p0, T0, M0, blt, rho2, p2, ux, uy, rf, Taw, beta
  use mod_constant, only : id_rescale, Cp, gamma_1, over_gamma, over_gamma_1, mu0_T0_S_over_T0_2_3
  use set_bc_common
  use set_init_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    integer i, j, ny_b
    real(8) dx1, dy1
    real(8) s, tanh_s, yi
    dx1 = Lx / dble(nx-1)
    dy1 = dx1

    x(1) = -0.016d0 !0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    s = 1.6d0
    tanh_s = tanh(s)
    do j = 1, ny
      yi = dble(j-1) / dble(ny-1)
      y(j) = Ly * (1 - tanh(s * (1.d0 - yi)) / tanh_s)
    enddo
    do j = 1, ny-1
      dy(j) = y(j+1) - y(j)
    enddo
    ! y(1) = 0.d0
    ! do j = 1, ny-1
    !   if (y(j) <= 3.d0 * blt) then
    !     dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
    !     ny_b  = j
    !   else
    !     dy(j) = dy1 * (1.d0 + 0.75d0 * dble(j-ny_b) / dble(ny-ny_b))
    !   endif
    !   y(j+1) = y(j) + dy(j)
    ! enddo
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    integer i, j
    
    do j = 1, ny
      do i = 1, nx
        Q(i,1,j) = rho0
        Q(i,2,j) = rho0 * u0
        Q(i,3,j) = 0.d0
        Q(i,4,j) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
      if((ys(j) / 0.08d0) > 1.2d0 * dtan(beta)) then
        Q(1,1,j) = rho2
        Q(1,2,j) = rho2 * ux
        Q(1,3,j) = rho2 * uy
        Q(1,4,j) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
      endif
    enddo
    do i = nx / 11 + 1, nx
      Q(i,1,1) = rho0
      Q(i,2,1) = 0.d0
      Q(i,3,1) = 0.d0
      Q(i,4,1) = Q(i,4,2)
    enddo    
  end subroutine set_init

  subroutine set_bc(myrank, nx, ny, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny) !Q / Jacobian
    real(8) Jacobian_tmp
    integer i, j, l, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall
    real(8) :: rhoin, uin, vin, pin, cin
    real(8) :: c_ext
    real(8) :: Rp, Rm, vb, cb, sb, ub, rhob, pb

    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      do l = 1, 4
        ! outlet
        QJ(nx,l,j) = QJ(nx-1,l,j)
    enddo;enddo

    c_ext = sqrt(gamma * p2 / rho2)

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! ---- interior state (j = ny-1) ----
      rhoin = QJ(i,1,ny-1) * Jacobian(i,ny-1)
      uin   = QJ(i,2,ny-1) / QJ(i,1,ny-1)
      vin   = QJ(i,3,ny-1) / QJ(i,1,ny-1)
      pin   = gamma_1 * ( QJ(i,4,ny-1) * Jacobian(i,ny-1) &
            - 0.5d0 * rhoin * (uin**2 + vin**2) )
      cin   = sqrt(gamma * pin / rhoin)

      ! ---- Riemann invariant（normal=y direction） ----
      Rp = vin   + 2.d0 * cin * over_gamma_1
      Rm = uy - 2.d0 * c_ext * over_gamma_1
      vb = 0.5d0 * (Rp + Rm)
      cb = 0.25d0 * gamma_1 * (Rp - Rm)

      if (vb >= 0.d0) then
        sb = pin / rhoin**gamma
        ub = uin
      else
        sb = p2 / rho2**gamma
        ub = ux
      endif

      rhob = (cb**2 / (gamma * sb)) ** over_gamma_1
      pb   = sb * rhob**gamma

      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      QJ(i,1,ny) = rhob * Jacobian_tmp
      QJ(i,2,ny) = rhob * ub * Jacobian_tmp
      QJ(i,3,ny) = rhob * vb * Jacobian_tmp
      QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      if(i <= nx / 11) then
        !Neumann
        QJ(i,1,1) = QJ(i,1,2)
        QJ(i,2,1) = QJ(i,2,2)
        QJ(i,3,1) = 0.d0 !QJ(i,3,2)
        QJ(i,4,1) = QJ(i,4,2) - 0.5d0 * QJ(i,3,2)**2 / QJ(i,1,2)
      else
        !NoSlip
        QJ(i,1,1) = QJ(i,1,2)
        QJ(i,2,1) = 0.d0
        QJ(i,3,1) = 0.d0
        p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
        QJ(i,4,1) = p_wall * over_gamma_1
      endif
    enddo
  end subroutine set_bc
end module set