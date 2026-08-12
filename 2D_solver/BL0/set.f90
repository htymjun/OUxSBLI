module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy, rf, Taw
  use mod_constant, only : id_rescale, Cp, gamma_1, over_gamma, over_gamma_1
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

    x(1) = 0.d0
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    s = 2.2d0
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
    real(8) :: rho0 = p0 / (R * T0)
    integer i, j
    !call set_init_tbl(nx, ny, xs, ys, 0.75d0*blt, blt, u0, p0, T0, M0, Q)
    do j = 1, ny
      do i = 1, nx
        Q(i,1,j) = rho0
        Q(i,2,j) = rho0 * u0
        Q(i,3,j) = 0.d0
        Q(i,4,j) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
    enddo;enddo
    do i = nx / 10 + 1, nx
      Q(i,1,1) = rho0
      Q(i,2,1) = 0.d0
      Q(i,3,1) = 0.d0
      Q(i,4,1) = Q(i,4,2)
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny) ! Q / Jacobian
    integer i, j, l, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb, v0 = 0.d0
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    real(8), parameter :: T       = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0    = p0 / (R * T)
    real(8), parameter :: c0      = sqrt(gamma * p0 / rho0)
    real(8), parameter :: over_c0 = 1.d0 / c0
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      Jacobian_tmp = 1.d0 / Jacobian(1,j)
      QJ(1,1,j) = rho0 * Jacobian_tmp
      QJ(1,2,j) = rho0 * u0 * Jacobian_tmp
      QJ(1,3,j) = 0.d0
      QJ(1,4,j) = (p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      do l = 1, 4
        ! outlet
        QJ(nx,l,j) = QJ(nx-1,l,j)
    enddo;enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! top
      ! Riemann invariants
      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      pin   = gamma_1 * (QJ(i,4,ny-1) - 0.5d0 * (QJ(i,2,ny-1)**2 + QJ(i,3,ny-1)**2) &
              / QJ(i,1,ny-1)) * Jacobian(i,ny-1)
      rhoin = QJ(i,1,ny-1) * Jacobian(i,ny-1)
      cin   = sqrt(gamma * pin / rhoin)
      vin   = QJ(i,3,ny-1) / QJ(i,1,ny-1)
      Rp   = vin + 2.d0 * cin * over_gamma_1
      Rm   = v0  - 2.d0 * c0  * over_gamma_1
      vb   = 0.5d0 * (Rp + Rm)
      cb   = 0.25d0 * gamma_1 * (Rp - Rm)
      rhob = (cb * over_c0)**(2.d0 * over_gamma_1) * rho0
      pb   = (rhob * cb**2) * over_gamma
      QJ(i,1,ny) = rhob * Jacobian_tmp
      QJ(i,2,ny) = rhob * u0 * Jacobian_tmp
      QJ(i,3,ny) = rhob * vb * Jacobian_tmp
      QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      if(i <= int(nx / 10)) then      
        ! Neumann
        QJ(i,1,1) = QJ(i,1,2)
        QJ(i,2,1) = QJ(i,2,2)
        QJ(i,3,1) = QJ(i,3,2)
        QJ(i,4,1) = QJ(i,4,2)
      else
        ! NoSlip
        QJ(i,1,1) = QJ(i,1,2)
        QJ(i,2,1) = 0.d0
        QJ(i,3,1) = 0.d0
        p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
        QJ(i,4,1) = p_wall * over_gamma_1
      endif
    enddo
  end subroutine set_bc
end module set

