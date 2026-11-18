module set
#ifdef PREC_DF
  use mod_prec
#endif
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
    REAL_W, intent(in), device     :: Jacobian(nx,ny)
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    call bc_1<<<(ny+127)/128, 128>>>(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    call bc_2<<<(nx+127)/128, 128>>>(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    call bc_3<<<(nx+127)/128, 128>>>(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
  end subroutine set_bc


  !> loop j = 2 .. ny-1 of the former CUF kernel 1
  attributes(global) subroutine bc_1(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: nx, ny
    REAL_W, intent(in), device     :: Jacobian(nx,ny)
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    integer j
    j = (blockIdx%x - 1) * blockDim%x + threadIdx%x + 1
    if (j > ny-1) return
    ! outlet
    QJ_1(nx,j) = QJ_1(nx-1,j); QJ_2(nx,j) = QJ_2(nx-1,j); QJ_3(nx,j) = QJ_3(nx-1,j); QJ_4(nx,j) = QJ_4(nx-1,j)
  end subroutine bc_1


  !> loop i = 1 .. nx of the former CUF kernel 2
  attributes(global) subroutine bc_2(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    use mod_globals, only : g_gamma => gamma, g_rho2 => rho2, g_p2 => p2, g_ux => ux, g_uy => uy
    integer, intent(in), value     :: nx, ny
    REAL_W, intent(in), device     :: Jacobian(nx,ny)
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    ! constants in the working kind (no FP64 arithmetic when PREC = 'FP32')
    REAL_W, parameter :: gamma = g_gamma, rho2 = g_rho2, p2 = g_p2, ux = g_ux, uy = g_uy
    REAL_W, parameter :: half = 0.5d0, quarter = 0.25d0, two = 2.d0, one = 1.d0, zero = 0.d0
    REAL_W, parameter :: c_ext = sqrt(gamma * p2 / rho2)
    REAL_W Jacobian_tmp
    REAL_T rhoin, uin, vin, pin, cin, Rp, Rm, vb, cb, sb, ub, rhob, pb
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx) return
    ! ---- interior state (j = ny-1) ----
    rhoin = QJ_1(i,ny-1) * Jacobian(i,ny-1)
    uin   = QJ_2(i,ny-1) / QJ_1(i,ny-1)
    vin   = QJ_3(i,ny-1) / QJ_1(i,ny-1)
    pin   = gamma_1 * ( QJ_4(i,ny-1) * Jacobian(i,ny-1) &
          - half * rhoin * (uin**2 + vin**2) )
    cin   = sqrt(gamma * pin / rhoin)

    ! ---- Riemann invariant (normal = y direction) ----
    Rp = vin + two * cin * over_gamma_1
    Rm = uy - two * c_ext * over_gamma_1
    vb = half * (Rp + Rm)
    cb = quarter * gamma_1 * (Rp - Rm)

    if (vb >= zero) then
      sb = pin / rhoin**gamma
      ub = uin
    else
      sb = p2 / rho2**gamma
      ub = ux
    endif

    rhob = (cb**2 / (gamma * sb)) ** over_gamma_1
    pb   = sb * rhob**gamma

    Jacobian_tmp = one / Jacobian(i,ny)
    QJ_1(i,ny) = rhob * Jacobian_tmp
    QJ_2(i,ny) = rhob * ub * Jacobian_tmp
    QJ_3(i,ny) = rhob * vb * Jacobian_tmp
    QJ_4(i,ny) = (pb * over_gamma_1 + half * rhob * (ub**2 + vb**2)) * Jacobian_tmp
  end subroutine bc_2


  !> loop i = 1 .. nx of the former CUF kernel 3
  attributes(global) subroutine bc_3(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: nx, ny
    REAL_W, intent(in), device     :: Jacobian(nx,ny)
    REAL_T, intent(inout), device  :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny)
    ! constants in the working kind (no FP64 arithmetic when PREC = 'FP32')
    REAL_W, parameter :: half = 0.5d0, zero = 0.d0
    REAL_T p_wall
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (i > nx) return
    if(i < i_LE) then
      !Neumann
      QJ_1(i,1) = QJ_1(i,2)
      QJ_2(i,1) = QJ_2(i,2)
      QJ_3(i,1) = zero
      QJ_4(i,1) = QJ_4(i,2) - half * QJ_3(i,2)**2 / QJ_1(i,2)
    else
      !NoSlip
      QJ_1(i,1) = QJ_1(i,2)
      QJ_2(i,1) = zero
      QJ_3(i,1) = zero
      p_wall = gamma_1 * (QJ_4(i,2) - half * (QJ_2(i,2)**2 + QJ_3(i,2)**2) / QJ_1(i,2))
      QJ_4(i,1) = p_wall * over_gamma_1
    endif
  end subroutine bc_3
end module set
