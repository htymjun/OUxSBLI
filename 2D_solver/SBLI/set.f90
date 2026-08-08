module set
  use cudafor
  use mod_globals, only : gamma, rho0, u0, p0, rho2, p2, ux, uy, beta, x_in, Xsh, i_LE
  use mod_constant, only : gamma_1, over_gamma_1
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    integer i, j
    real(8) dx1
    real(8) tanh_s, yi
    real(8), parameter :: s = 1.6d0 ! tanh wall-clustering stretch
    dx1 = Lx / dble(nx-1)

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
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,ny,4)
    integer i, j

    do j = 1, ny
      do i = 1, nx
        Q(i,j,1) = rho0
        Q(i,j,2) = rho0 * u0
        Q(i,j,3) = 0.d0
        Q(i,j,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
      ! The inlet column i=1 is never updated afterwards (interior-only RK kernels,
      ! no inlet BC in set_bc), so it stays frozen at these values: below the
      ! incident-shock trace it holds the freestream, above it the post-shock
      ! state, acting as the oblique-shock generator.
      if (ys(j) > (Xsh - x_in) * dtan(beta)) then
        Q(1,j,1) = rho2
        Q(1,j,2) = rho2 * ux
        Q(1,j,3) = rho2 * uy
        Q(1,j,4) = p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)
      endif
    enddo
    do i = i_LE, nx
      Q(i,1,1) = rho0
      Q(i,1,2) = 0.d0
      Q(i,1,3) = 0.d0
      Q(i,1,4) = Q(i,2,4)
    enddo
  end subroutine set_init

  subroutine set_bc(myrank, nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    real(8) Jacobian_tmp
    integer i, j
    real(8) :: p_wall
    real(8) :: rhoin, uin, vin, pin, cin
    real(8) :: Rp, Rm, vb, cb, sb, ub, rhob, pb
    ! exterior (post-shock) sound speed for the top Riemann state
    real(8), parameter :: c_ext = sqrt(gamma * p2 / rho2)

    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      ! outlet
      QJ_1(nx,j) = QJ_1(nx-1,j); QJ_2(nx,j) = QJ_2(nx-1,j); QJ_3(nx,j) = QJ_3(nx-1,j); QJ_4(nx,j) = QJ_4(nx-1,j)
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! ---- interior state (j = ny-1) ----
      rhoin = QJ_1(i,ny-1) * Jacobian(i,ny-1)
      uin   = QJ_2(i,ny-1) / QJ_1(i,ny-1)
      vin   = QJ_3(i,ny-1) / QJ_1(i,ny-1)
      pin   = gamma_1 * ( QJ_4(i,ny-1) * Jacobian(i,ny-1) &
            - 0.5d0 * rhoin * (uin**2 + vin**2) )
      cin   = sqrt(gamma * pin / rhoin)

      ! ---- Riemann invariant (normal = y direction) ----
      Rp = vin + 2.d0 * cin * over_gamma_1
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
      QJ_1(i,ny) = rhob * Jacobian_tmp
      QJ_2(i,ny) = rhob * ub * Jacobian_tmp
      QJ_3(i,ny) = rhob * vb * Jacobian_tmp
      QJ_4(i,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      if(i < i_LE) then
        !Neumann
        QJ_1(i,1) = QJ_1(i,2)
        QJ_2(i,1) = QJ_2(i,2)
        QJ_3(i,1) = 0.d0
        QJ_4(i,1) = QJ_4(i,2) - 0.5d0 * QJ_3(i,2)**2 / QJ_1(i,2)
      else
        !NoSlip
        QJ_1(i,1) = QJ_1(i,2)
        QJ_2(i,1) = 0.d0
        QJ_3(i,1) = 0.d0
        p_wall = gamma_1 * (QJ_4(i,2) - 0.5d0 * (QJ_2(i,2)**2 + QJ_3(i,2)**2) / QJ_1(i,2))
        QJ_4(i,1) = p_wall * over_gamma_1
      endif
    enddo
  end subroutine set_bc
end module set
