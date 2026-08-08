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
    real(8), parameter :: s = 2.4d0 ! tanh wall-clustering stretch
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
    ! top far field: interior (j = ny-1) state and boundary state
    real(8) :: rhoin, uin, vin, pin, rhob, ub, vb, pb, sb
    ! cache
    real(8) Jacobian_tmp
    ! freestream entropy at the top boundary
    real(8), parameter :: s0 = p0 / rho0**gamma

    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      ! outlet: subsonic pressure outflow.
      !
      ! Extrapolating all four conservatives imposes nothing where subsonic
      ! outflow needs exactly one condition, so the exit pressure floated: it sat
      ! 0.008*q_inf below p0 and dragged a favorable gradient ~25 mm back up the
      ! plate, inflating Cf there by 6%. A zero-incidence flat plate has dp/dx = 0,
      ! so p = p0 is the physically right single condition and matches what the
      ! top boundary already imposes. Density and both momenta still come from the
      ! interior; only the energy is rebuilt around p0.
      QJ_1(nx,j) = QJ_1(nx-1,j)
      QJ_2(nx,j) = QJ_2(nx-1,j)
      QJ_3(nx,j) = QJ_3(nx-1,j)
      QJ_4(nx,j) = p0 * over_gamma_1 / Jacobian(nx,j) &
                 + 0.5d0 * (QJ_2(nx,j)**2 + QJ_3(nx,j)**2) / QJ_1(nx,j)
    enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! top: constant-pressure far field.
      !
      ! The Riemann-invariant far field this used to be (and that SBLI still uses)
      ! holds Rm = v_ext - 2*c_ext/(gamma-1) fixed, which in steady state forces
      ! p_b - p0 = rho0*c0*v_b, i.e. (p_b - p0)/q_inf = 2*v_b/(u0*M0). That is a
      ! transient non-reflection property, not p -> p0, and the 1/M0 factor makes it
      ! unusable at M0 = 0.1: the displacement-induced v_b (0.11 m/s at the leading
      ! edge, decaying downstream) became a 0.07*q_inf favorable pressure gradient
      ! along the plate, accelerating the edge flow 3% and inflating Cf by >10% at
      ! the trailing end. At SBLI's M0 = 2.15 there is no 1/M0 amplification, which
      ! is why the same form is fine there.
      !
      ! For a steady zero-incidence plate the far field is p = p0, so impose that and
      ! take everything else from the outgoing characteristics: v is always
      ! extrapolated (this is also the zero-shear far-field condition, replacing the
      ! old u = u0 clamp that put a spurious shear layer across the outer half of the
      ! domain), while entropy and tangential momentum are upwinded on the sign of v.
      ! Subsonic outflow takes 1 condition from outside, subsonic inflow takes 3.
      rhoin = QJ_1(i,ny-1) * Jacobian(i,ny-1)
      uin   = QJ_2(i,ny-1) / QJ_1(i,ny-1)
      vin   = QJ_3(i,ny-1) / QJ_1(i,ny-1)
      pin   = gamma_1 * (QJ_4(i,ny-1) * Jacobian(i,ny-1) &
              - 0.5d0 * rhoin * (uin**2 + vin**2))

      pb = p0
      vb = vin
      if (vb >= 0.d0) then
        ! outflow: entropy and tangential momentum leave the domain
        sb = pin / rhoin**gamma
        ub = uin
      else
        ! inflow: both come from the freestream (rhob then reduces to rho0)
        sb = s0
        ub = u0
      endif
      rhob = (pb / sb)**over_gamma

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
