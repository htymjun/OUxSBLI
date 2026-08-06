module set
  use cudafor
  use mod_globals, only : gamma, rho0, u0, p0, i_LE
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1
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

    ! place the leading edge (first no-slip wall point, i = i_LE) exactly at x = 0
    x(1) = -dble(i_LE - 1) * dx1
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

    ! Uniform freestream everywhere. The inlet column i=1 is never updated
    ! afterwards (interior-only RK kernels, no inlet BC in set_bc), so it stays
    ! frozen at the freestream state and acts as a Dirichlet inflow.
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1) = rho0
        Q(i,j,2) = rho0 * u0
        Q(i,j,3) = 0.d0
        Q(i,j,4) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo
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
    integer i, j
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, vb, cb, pb
    ! cache
    real(8) Jacobian_tmp
    ! freestream state at the top boundary
    real(8), parameter :: v0      = 0.d0
    real(8), parameter :: c0      = sqrt(gamma * p0 / rho0)
    real(8), parameter :: over_c0 = 1.d0 / c0

    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      ! outlet
      QJ_1(nx,j) = QJ_1(nx-1,j); QJ_2(nx,j) = QJ_2(nx-1,j); QJ_3(nx,j) = QJ_3(nx-1,j); QJ_4(nx,j) = QJ_4(nx-1,j)
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! top
      ! Riemann invariants
      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      pin   = gamma_1 * (QJ_4(i,ny-1) - 0.5d0 * (QJ_2(i,ny-1)**2 + QJ_3(i,ny-1)**2) &
              / QJ_1(i,ny-1)) * Jacobian(i,ny-1)
      rhoin = QJ_1(i,ny-1) * Jacobian(i,ny-1)
      cin   = sqrt(gamma * pin / rhoin)
      vin   = QJ_3(i,ny-1) / QJ_1(i,ny-1)
      Rp   = vin + 2.d0 * cin * over_gamma_1
      Rm   = v0  - 2.d0 * c0  * over_gamma_1
      vb   = 0.5d0 * (Rp + Rm)
      cb   = 0.25d0 * gamma_1 * (Rp - Rm)
      rhob = (cb * over_c0)**(2.d0 * over_gamma_1) * rho0
      pb   = (rhob * cb**2) * over_gamma
      QJ_1(i,ny) = rhob * Jacobian_tmp
      QJ_2(i,ny) = rhob * u0 * Jacobian_tmp
      QJ_3(i,ny) = rhob * vb * Jacobian_tmp
      QJ_4(i,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
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
