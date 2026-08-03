module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, Cp, Pr, u0, rho0, p0, T0, M0, blt, beta, &
                          rho2, p2, ux, uy, rf, Taw, rho3, p3, ux3, uy3
  use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  use set_init_common
  implicit none
  integer No
contains
  subroutine set_grid(myrank, nx, ny, Lx, Ly, x, y, dx, dy)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: Lx, Ly
    real(8), intent(out) :: x(nx), y(ny), dx(nx-1), dy(ny-1)
    integer i, j, ny_b
    real(8) dx1, dy1, ximp, Lx_s
    dx1  = Lx / dble(nx-1)
    dy1  = dx1
    ximp = 70.d0 * blt

    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    x(:) = x(:) - ximp

    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
        ny_b  = j
      else
        dy(j) = dy1 * (1.d0 + 0.75d0 * dble(j-ny_b) / dble(ny-ny_b))
      endif
      y(j+1) = y(j) + dy(j)
    enddo

    Lx_s = y(ny) / dble(beta) / blt
    do i = 1, nx
      if (x(i) / blt + Lx_s >= 0.d0) then
        No = i
        exit
      endif
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,ny,4)
    call set_init_tbl(nx, ny, xs, ys, blt, blt, u0, p0, T0, M0, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny) ! Q / Jacobian
    integer i, j
    real(8) :: p_wall, pre, rho, rhou, rhov, rhow, p, e
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb, v0 = 0.d0
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    real(8), parameter :: T    = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0 = p0 / (R * T)
    real(8), parameter :: c0   = sqrt(gamma * p0 / rho0)
    real(8), parameter :: c3   = sqrt(gamma * p3 / rho3)
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      ! cyclic
      QJ_1(1,j)    = QJ_1(nx-5,j); QJ_2(1,j)    = QJ_2(nx-5,j); QJ_3(1,j)    = QJ_3(nx-5,j); QJ_4(1,j)    = QJ_4(nx-5,j)
      QJ_1(2,j)    = QJ_1(nx-4,j); QJ_2(2,j)    = QJ_2(nx-4,j); QJ_3(2,j)    = QJ_3(nx-4,j); QJ_4(2,j)    = QJ_4(nx-4,j)
      QJ_1(3,j)    = QJ_1(nx-3,j); QJ_2(3,j)    = QJ_2(nx-3,j); QJ_3(3,j)    = QJ_3(nx-3,j); QJ_4(3,j)    = QJ_4(nx-3,j)
      QJ_1(nx-2,j) = QJ_1(4,j);    QJ_2(nx-2,j) = QJ_2(4,j);    QJ_3(nx-2,j) = QJ_3(4,j);    QJ_4(nx-2,j) = QJ_4(4,j)
      QJ_1(nx-1,j) = QJ_1(5,j);    QJ_2(nx-1,j) = QJ_2(5,j);    QJ_3(nx-1,j) = QJ_3(5,j);    QJ_4(nx-1,j) = QJ_4(5,j)
      QJ_1(nx,j)   = QJ_1(6,j);    QJ_2(nx,j)   = QJ_2(6,j);    QJ_3(nx,j)   = QJ_3(6,j);    QJ_4(nx,j)   = QJ_4(6,j)
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! top
      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      QJ_1(i,ny) = rho0 * Jacobian_tmp
      QJ_2(i,ny) = rho0 * u0 * Jacobian_tmp
      QJ_3(i,ny) = 0.d0
      QJ_4(i,ny) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      ! NoSlip
      QJ_1(i,1) = QJ_1(i,2)
      QJ_2(i,1) = 0.d0
      QJ_3(i,1) = 0.d0
      p_wall = gamma_1 * (QJ_4(i,2) - 0.5d0 * (QJ_2(i,2)**2 + QJ_3(i,2)**2) / QJ_1(i,2))
      QJ_4(i,1) = p_wall * over_gamma_1
    enddo
  end subroutine set_bc
end module set

