!> SWTBLI: grid, initial condition and boundary conditions of the two coupled
!> domains driven by src/sbli.f90 (ranks 0/1 = TBL generator, ranks 2/3 = SBLI).
!>
!> Design (see mod_globals.f90 for the numbers):
!>  * Both domains share dx0, the wall-normal grid function and the spanwise grid,
!>    so the TBL planes nre2-5..nre2-3 can be copied verbatim (Q/J with identical
!>    Jacobians) into the SBLI inlet columns i=1..3: the 6-point stencils see a
!>    continuous field across the junction.
!>  * Above the TBL domain (j >= ny1) the SBLI inlet continues a running time mean
!>    (tau = 10 delta/u0) of the spanwise-mean TBL row ny1-1, so there is no step in
!>    the mean at y = Ly1 and no acoustic noise is injected into the outer flow.
!>  * The incident shock enters through the SBLI top boundary at
!>    x_top = x_imp - Ly2/tan(beta1). The top boundary is a 1D characteristic
!>    (Riemann-invariant) condition everywhere: R+ is extrapolated, R-, entropy and
!>    the tangential velocity come from an exterior state that is state 1 upstream,
!>    blends smoothly (tanh over ~4 cells) into state 2 at x_top (state 2 rescaled
!>    to the actual upstream level, see set_bc) and stays state 2 from x_top + 3 delta
!>    on while the flow under the top still comes down through the incident shock.
!>    Fixing all five variables would reflect outgoing acoustic waves, because the
!>    normal Mach number behind the incident shock is only ~0.3.
!>    Behind the reflected (separation) shock, detected from the running time mean
!>    of the local spanwise-mean v, the top row is instead copied from row ny-1 at
!>    the foot of the right-running characteristic, so that shock and its waves
!>    leave the domain.
!>  * The top 2 delta of the SBLI domain (gy_top per cell) and the last 9.4 delta in
!>    x (3 %/cell) are stretched to damp outgoing waves before the boundaries.
!>  * The inter-domain transfer is staged through host memory (1.6 MB per stage):
!>    a device-pointer exchange can corrupt silently on some MPI stacks
!>    (3D_solver/CLAUDE.md, COMMZ section).
!>
!> Calibration of blt/Cf_init/dnu0 (Re_theta = 500): Guarini et al.'s composite
!> profile (Reichardt + Finley wake, Van Driest II inversion, Crocco-Busemann
!> temperature) for M0 = 2.28, T_tot = 290 K, p_tot = 50 kPa gives
!> Cf = 3.783e-3, delta99 = 1.0724 mm, delta_nu = 9.952 um, Re_tau = 107.8.
module set
  use cudafor
  use mpi
  use mod_globals, only : ny1, nz1, nbuf, nre2, nin, gamma, R, u0, p0, T0, rho0, M0, blt, dt, &
                          rf, Taw, Cf_init, dx0, dyw0, y_str, gy_top, theta_deg, ximp_in, Lz1
  use mod_constant, only : gamma_1, over_gamma_1, mu0_T0_S_over_T0_2_3, id_rescale
  use set_bc_common
  use set_init_common, only : calc_Gaussian_filter_x, calc_Gaussian_filter_z
  use calc_para, only : flatten_rescale, reconstruct_sbli_inlet
  implicit none
  private
  public set_grid, set_init, set_bc, set_bc_mut

  real(8), parameter :: pi = 3.14159265358979323846d0
  ! state behind the incident shock (2), filled by init_shock
  real(8), save :: rho2s = 0.d0, u2s = 0.d0, v2s = 0.d0, p2s = 0.d0
  real(8), save :: beta1 = 0.d0, x_imp = 0.d0, x_top = 0.d0, x_dyn = 0.d0
  ! host copy of this rank's x (set_grid), turned into per-column exterior states
  real(8), allocatable, save :: x_host(:)
  ! device data built on the first set_bc call (the device is selected later than set_grid)
  logical, save :: dev_ready = .false.
  real(8), allocatable, device, save :: s_d(:)
  ! running mean (tau_mavg) of the spanwise-mean primitive state (rho, u, v, p) next to
  ! the top, per column; downstream of x_dyn it decides between the Riemann condition
  ! on state 2 and the characteristic copy, and gives the characteristic's angle
  real(8), allocatable, device, save :: etop_d(:,:), itop_d(:,:)
  logical, save :: etop_init = .false.
  ! x on the device and the spacing of the top cell, for the characteristic extrapolation
  real(8), allocatable, device, save :: x_d(:)
  real(8), save :: dy_top = 0.d0
  integer, allocatable, device, save :: idyn_d(:)
  real(8), allocatable, device, save :: q1d_d(:), mrow_d(:,:), mavg_d(:,:)
  logical, save :: mavg_init = .false.
  ! running mean of the SBLI top pressure upstream of x_top, over p0 (see set_bc)
  real(8), save :: f_top = 1.d0
  logical, save :: ftop_init = .false.
  integer, save :: i_up0 = 0, i_up1 = 0, ncall_top = 0
  real(8), allocatable, save :: ptop_h(:)   ! host copy of the top pressure of those columns
  ! time constant of the running means (inlet rows above the TBL domain, top boundary)
  real(8), parameter :: tau_mavg = 10.d0   ! [blt/u0]
  ! their relaxation factor per set_bc call: set_bc runs once per RK stage, 3 times per dt
  real(8), parameter :: eps_mavg = dt / (3.d0 * tau_mavg * blt / u0)
  real(8), allocatable, pinned, save :: q1d_h(:)
  integer, parameter :: tag_inlet = 4242
contains
  ! ======================================================================
  ! grid
  ! ======================================================================
  !> Wall-normal grid shared by both domains: y(1) = 0, first cell dyw0, +3.5 %/cell
  !> up to dy = dx0, uniform to y_str, gy_top per cell above
  !> (ouxsbli.analysis.swtbli_stats.make_y is the Python twin). The TBL domain uses the
  !> first ny1 points of the same sequence, so the junction rows coincide exactly.
  subroutine make_y(ny, y)
    integer, intent(in)  :: ny
    real(8), intent(out) :: y(ny)
    real(8) dy
    integer j
    y(1) = 0.d0
    dy = dyw0
    do j = 2, ny
      y(j) = y(j-1) + dy
      if (y(j) < y_str) then
        dy = min(dy * 1.035d0, dx0)
      else
        dy = dy * gy_top
      endif
    enddo
  end subroutine make_y


  !> rank_domain = 0 for the TBL generator, 2 for the SBLI domain (sbli.f90 passes
  !> myrank for even ranks and myrank-1 for the I/O ranks).
  subroutine set_grid(rank_domain, nx, ny, nz, Lx, Ly, Lz, Lx1_in, x, y, z, dx, dy, dz)
    integer, intent(in)  :: rank_domain, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz, Lx1_in
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k, nu
    real(8) d, dz1

    if (rank_domain < 2) then
      do i = 1, nx
        x(i) = dble(i-1) * dx0
      enddo
    else
      ! inlet column i=1 coincides with TBL plane nin = nre2-5
      nu = nx - 1 - nbuf
      x(1) = dble(nin-1) * dx0
      do i = 2, nu + 1
        x(i) = x(1) + dble(i-1) * dx0
      enddo
      d = dx0
      do i = nu + 2, nx
        d = min(d * 1.03d0, 4.d0 * dx0)
        x(i) = x(i-1) + d
      enddo
    endif
    do i = 1, nx-1
      dx(i) = x(i+1) - x(i)
    enddo

    call make_y(ny, y)
    do j = 1, ny-1
      dy(j) = y(j+1) - y(j)
    enddo

    ! periodic span Lz over nz-6 unique planes; plane k=4 at z=-Lz/2
    dz1 = Lz / dble(nz-6)
    do k = 1, nz
      z(k) = dble(k-4) * dz1 - 0.5d0 * Lz
    enddo
    do k = 1, nz-1
      dz(k) = dz1
    enddo

    if (allocated(x_host)) deallocate(x_host)
    allocate(x_host(nx))
    x_host = x
    dy_top = y(ny) - y(ny-1)
    call init_shock(rank_domain, nx, ny, x, y)
  end subroutine set_grid


  ! ======================================================================
  ! oblique-shock states
  ! ======================================================================
  !> weak-solution shock angle for deflection th [rad] at Mach M (bisection)
  function beta_weak(M, th) result(b)
    real(8), intent(in) :: M, th
    real(8) b, lo, hi, mid, f_lo, f_mid
    integer it
    lo = asin(1.d0 / M) + 1.d-10
    hi = 64.d0 * pi / 180.d0
    f_lo = tan(th) - 2.d0 / tan(lo) * (M**2 * sin(lo)**2 - 1.d0) / (M**2 * (gamma + cos(2.d0*lo)) + 2.d0)
    do it = 1, 200
      mid = 0.5d0 * (lo + hi)
      f_mid = tan(th) - 2.d0 / tan(mid) * (M**2 * sin(mid)**2 - 1.d0) / (M**2 * (gamma + cos(2.d0*mid)) + 2.d0)
      if (f_lo * f_mid <= 0.d0) then
        hi = mid
      else
        lo = mid; f_lo = f_mid
      endif
    enddo
    b = 0.5d0 * (lo + hi)
  end function beta_weak


  !> Rankine-Hugoniot across an oblique shock of angle b (relative to the
  !> upstream flow) for upstream (rho, p, Mach M); returns rho, p and Mach behind.
  subroutine oblique_rh(M, b, th, rho_a, p_a, rho_b, p_b, M_b)
    real(8), intent(in)  :: M, b, th, rho_a, p_a
    real(8), intent(out) :: rho_b, p_b, M_b
    real(8) Mn, Mn2, Mnb
    Mn  = M * sin(b)
    Mn2 = Mn**2
    p_b   = p_a * (1.d0 + 2.d0 * gamma / (gamma + 1.d0) * (Mn2 - 1.d0))
    rho_b = rho_a * (gamma + 1.d0) * Mn2 / ((gamma - 1.d0) * Mn2 + 2.d0)
    Mnb = sqrt((1.d0 + 0.5d0 * (gamma - 1.d0) * Mn2) / (gamma * Mn2 - 0.5d0 * (gamma - 1.d0)))
    M_b = Mnb / sin(b - th)
  end subroutine oblique_rh


  subroutine init_shock(rank_domain, nx, ny, x, y)
    integer, intent(in) :: rank_domain, nx, ny
    real(8), intent(in) :: x(nx), y(ny)
    real(8) th, M2, M3, br, c2, V2, rho3, p3
    th = theta_deg * pi / 180.d0
    if (theta_deg <= 0.d0) then
      rho2s = rho0; u2s = u0; v2s = 0.d0; p2s = p0
      p3 = p0
      beta1 = asin(1.d0 / M0)
    else
      beta1 = beta_weak(M0, th)
      call oblique_rh(M0, beta1, th, rho0, p0, rho2s, p2s, M2)
      c2 = sqrt(gamma * p2s / rho2s)
      V2 = M2 * c2
      u2s =  V2 * cos(th)
      v2s = -V2 * sin(th)
      br = beta_weak(M2, th)
      ! state 3 (behind the reflected shock) is only reported, not imposed
      call oblique_rh(M2, br, th, rho2s, p2s, rho3, p3, M3)
    endif
    if (rank_domain >= 2) then
      x_imp = x(1) + ximp_in * blt
      x_top = x_imp - y(ny) / tan(beta1)
      x_dyn = x_top + 3.d0 * blt
      if (rank_domain == 2) then
        print '(a,f6.2,a,f7.3,a,f7.4,a,f7.4)', " SBLI: theta = ", theta_deg, " deg, beta1 = ", beta1*180.d0/pi, &
              " deg, p2/p1 = ", p2s/p0, ", p3/p1 = ", p3/p0
        print '(a,f8.3,a,f8.3,a,f8.3,a)', " SBLI: x(1) = ", x(1)/blt, " delta, x_top = ", x_top/blt, &
              " delta, x_imp = ", x_imp/blt, " delta"
      endif
    endif
  end subroutine init_shock


  ! ======================================================================
  ! initial condition: turbulent composite profile + coherent perturbations
  ! ======================================================================
  function mu_suth(T) result(mu)
    real(8), intent(in) :: T
    real(8) mu
    mu = mu0_T0_S_over_T0_2_3 * T**1.5d0 / (T + 110.4d0)
  end function mu_suth


  !> Reichardt inner law + Finley wake in Van Driest-transformed units
  function ucp(yp, dp) result(u)
    real(8), intent(in) :: yp, dp
    real(8) u, e, C1
    real(8), parameter :: kappa = 0.41d0, Cl = 5.2d0, Pi_w = 0.25d0
    C1 = -log(kappa) / kappa + Cl
    e = min(yp / dp, 1.d0)
    u = log(1.d0 + kappa * yp) / kappa &
        + C1 * (1.d0 - exp(-yp / 11.d0) - (yp / 11.d0) * exp(-0.33d0 * yp)) &
        + (e**2 - e**3 + 6.d0 * Pi_w * e**2 - 4.d0 * Pi_w * e**3) / kappa
  end function ucp


  !> Composite mean profile (Guarini et al. 2000, eq. 3.3) for skin friction
  !> Cf_init; ouxsbli/analysis/tbl_stats.reference_profile is the Python twin.
  subroutine composite_profile(ny, y, u, T, dnu)
    integer, intent(in)  :: ny
    real(8), intent(in)  :: y(ny)
    real(8), intent(out) :: u(ny), T(ny), dnu
    real(8) rhow, A, muw, uip, ut, lo, hi, mid, target
    integer j, it
    rhow = p0 / (R * Taw)
    A    = sqrt(rf * 0.5d0 * (gamma - 1.d0) * M0**2 / (Taw / T0))
    muw  = mu_suth(Taw)
    uip  = 1.d0 / sqrt(0.5d0 * Cf_init * (rho0 / rhow))
    ut   = u0 / uip
    dnu  = muw / rhow / ut
    target = asin(A) / A * uip
    lo = 2.d0; hi = 1.d6
    do it = 1, 200
      mid = 0.5d0 * (lo + hi)
      if (ucp(mid, mid) > target) then
        hi = mid
      else
        lo = mid
      endif
    enddo
    do j = 1, ny
      u(j) = min(max(sin(A * ucp(y(j) / dnu, mid) / uip) / A * u0, 0.d0), u0)
    enddo
    u(1) = 0.d0
    do j = 1, ny
      T(j) = Taw - rf * 0.5d0 * (gamma - 1.d0) * M0**2 * T0 * (u(j) / u0)**2
    enddo
  end subroutine composite_profile


  !> Mean: composite turbulent profile with delta99 = blt (not Blasius).
  !> Perturbations chosen to by-pass the slow linear transition:
  !>  (1) near-wall streaks, 0.12 u0, peak at y+ = 15, spanwise wavelength
  !>      Lz1/2 (1 delta, ~108 wall units, in the 2 delta span the runs were started
  !>      with; 2 delta for Lz1 = 4 delta), modulated in x over 4 delta;
  !>  (2) divergence-free (in y-z) quasi-streamwise vortex pairs, |v'|,|w'| ~ 0.05 u0,
  !>      centred at y+ = 30, same spanwise wavelength, phase-shifted in x;
  !>  (3) Gaussian-filtered random noise, 5 % of u0, inside the boundary layer.
  !> The recycling then keeps the disturbances in the domain.
  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    real(8), allocatable :: um(:), Tm(:), up(:,:,:), vp(:,:,:), wp(:,:,:), rnd(:,:,:,:)
    real(8) dnu, yp, fs, g, gp, yv, lz_w, lx_w, psi0, fd, ax, bx, rho, u, v, w, T
    integer i, j, k, ir, jr, kr, nxr, nyr, nzr, sz
    integer, allocatable :: seed(:)

    allocate(um(ny), Tm(ny), up(nx,ny,nz), vp(nx,ny,nz), wp(nx,ny,nz))
    call composite_profile(ny, ys, um, Tm, dnu)
    if (myrank == 0) print '(a,f8.3,a,f8.2)', " IC: wall unit [um] = ", dnu*1.d6, ",  Re_tau = ", blt/dnu

    lz_w = 0.5d0 * Lz1          ! 1 delta for the 2 delta span of the first runs, 2 delta for Lz1 = 4 delta
    lx_w = 4.d0 * blt
    yv   = 30.d0 * dnu
    psi0 = 0.05d0 * u0 / (2.d0 * pi / lz_w)   ! max |v'| = 0.05 u0 (g peaks at 1)

    ! (1) + (2) coherent part
    do k = 1, nz
      do j = 1, ny
        yp = ys(j) / dnu
        fs = (yp / 15.d0) * exp(1.d0 - yp / 15.d0)
        g  = (ys(j) / yv)**2 * exp(1.d0 - (ys(j) / yv)**2)
        gp = 2.d0 * ys(j) / yv**2 * exp(1.d0 - (ys(j) / yv)**2) * (1.d0 - (ys(j) / yv)**2)
        do i = 1, nx
          ax = 0.5d0 + 0.5d0 * sin(2.d0 * pi * xs(i) / lx_w)
          bx = cos(2.d0 * pi * xs(i) / lx_w + 0.25d0 * pi)
          up(i,j,k) = 0.12d0 * u0 * fs * cos(2.d0 * pi * zs(k) / lz_w) * ax
          vp(i,j,k) =  psi0 * g  * (2.d0 * pi / lz_w) * cos(2.d0 * pi * zs(k) / lz_w + 0.5d0 * pi) * bx
          wp(i,j,k) = -psi0 * gp * sin(2.d0 * pi * zs(k) / lz_w + 0.5d0 * pi) * bx
    enddo;enddo;enddo

    ! (3) random part on a coarse lattice (10 x 2 x 2 cells), filtered
    nxr = nx / 10 + 2; nyr = ny / 2 + 2; nzr = nz / 2 + 2
    allocate(rnd(3,nxr,nyr,nzr))
    call random_seed(size=sz)
    allocate(seed(sz))
    seed = 48 + 7 * myrank
    call random_seed(put=seed)
    call random_number(rnd)
    rnd = 2.d0 * rnd - 1.d0
    block
      real(8), allocatable :: nu_(:,:,:), nv_(:,:,:), nw_(:,:,:)
      allocate(nu_(nx,ny,nz), nv_(nx,ny,nz), nw_(nx,ny,nz))
      do k = 1, nz
        do j = 1, ny
          if (ys(j) <= blt) then
            fd = sin(0.5d0 * pi * min(ys(j) / (20.d0 * dnu), 1.d0)) * (1.d0 - ys(j) / blt)
          else
            fd = 0.d0
          endif
          do i = 1, nx
            ir = i / 10 + 1; jr = j / 2 + 1; kr = k / 2 + 1
            nu_(i,j,k) = 0.05d0 * u0 * fd * rnd(1,ir,jr,kr) * 2.d0
            nv_(i,j,k) = 0.05d0 * u0 * fd * rnd(2,ir,jr,kr)
            nw_(i,j,k) = 0.05d0 * u0 * fd * rnd(3,ir,jr,kr)
      enddo;enddo;enddo
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, nu_)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, nu_)
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, nv_)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, nv_)
      call calc_Gaussian_filter_x(nx, ny, nz, 3, xs, nw_)
      call calc_Gaussian_filter_z(nx, ny, nz, 3, zs, nw_)
      up = up + nu_; vp = vp + nv_; wp = wp + nw_
      deallocate(nu_, nv_, nw_)
    end block

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          u = um(j) + up(i,j,k)
          v = vp(i,j,k)
          w = wp(i,j,k)
          if (j == 1) then
            u = 0.d0; v = 0.d0; w = 0.d0
          endif
          T   = Tm(j)
          rho = p0 / (R * T)
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = rho * u
          Q(i,j,k,3) = rho * v
          Q(i,j,k,4) = rho * w
          Q(i,j,k,5) = p0 * over_gamma_1 + 0.5d0 * rho * (u**2 + v**2 + w**2)
    enddo;enddo;enddo

    ! periodic ghost planes (the random part is not periodic by construction)
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1,:)    = Q(i,j,nz-5,:)
        Q(i,j,2,:)    = Q(i,j,nz-4,:)
        Q(i,j,3,:)    = Q(i,j,nz-3,:)
        Q(i,j,nz-2,:) = Q(i,j,4,:)
        Q(i,j,nz-1,:) = Q(i,j,5,:)
        Q(i,j,nz,:)   = Q(i,j,6,:)
    enddo;enddo
    deallocate(um, Tm, up, vp, wp, rnd, seed)
  end subroutine set_init


  ! ======================================================================
  ! boundary conditions
  ! ======================================================================
  !> One-time device setup: per-column exterior state of the top boundary and the
  !> inter-domain buffers. Done here, not in set_grid, because the GPU is selected
  !> inside RungeKutta (check_gpu), after set_grid has run.
  subroutine setup_device(myrank, nx, ny, nz)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), allocatable :: sw(:)
    integer, allocatable :: idyn(:)
    integer i, n
    allocate(sw(nx), idyn(nx))
    do i = 1, nx
      ! weight of state 2 in the exterior state: 0 upstream, 1 behind the incident
      ! shock, tanh over ~4 cells at x_top
      sw(i) = 0.d0
      if (myrank >= 2) sw(i) = 0.5d0 * (1.d0 + tanh((x_host(i) - x_top) / (2.d0 * dx0)))
      idyn(i) = 0
      if (myrank >= 2 .and. x_host(i) >= x_dyn) idyn(i) = 1
    enddo
    allocate(s_d(nx), idyn_d(nx), etop_d(nx,4), itop_d(nx,4), x_d(nx))
    x_d = x_host
    s_d = sw; idyn_d = idyn
    n = 3 * (ny1 - 2) * (nz - 6) * 5
    allocate(q1d_d(n), q1d_h(n), mrow_d(3,4), mavg_d(3,4))
    deallocate(sw, idyn)
    ! columns used for the upstream top-pressure level: x_in + 1 delta .. x_top - 1 delta
    if (myrank >= 2) then
      do i = 1, nx
        if (i_up0 == 0 .and. x_host(i) >= x_host(1) + blt) i_up0 = i
        if (x_host(i) <= x_top - blt) i_up1 = i
      enddo
      if (i_up1 < i_up0 + 5) i_up1 = i_up0 + 5
      allocate(ptop_h(i_up1 - i_up0 + 1))
    endif
    dev_ready = .true.
  end subroutine setup_device


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(in), device, optional :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    integer i, j, k, offset, n, ierr
    integer mstat(MPI_STATUS_SIZE)
    real(8) p_wall, Jt, rhoin, uin, vin, win, pin, cin, rhox, ux, vx, px, cx, Rp, Rm, vb, cb, sb, ub, wb, rhob, pb
    ! exterior states as plain scalars for the CUF kernels
    real(8) v2half, psum, fp, fr, r2, u2, v2, p2, sx
    real(8) dyt, alpha, xsrc, wa, ra, ua, va, wwa, pa, rb2, ub2, vb2, wb2, pb2
    integer ii

    if (.not. dev_ready) call setup_device(myrank, nx, ny, nz)
    v2half = 0.5d0 * v2s
    r2 = rho2s; u2 = u2s; v2 = v2s; p2 = p2s
    dyt = dy_top
    n = 3 * (ny1 - 2) * (nz - 6) * 5

    if (myrank == 0) then
      ! ---------------- TBL generator ----------------
      if (kind(id_rescale) == 4 .and. present(Qre_1)) then
        !$cuf kernel do(2)<<<*,*>>>
        do k = 1, nz-6
          do j = 2, ny-1
            offset = ny*(k-1) + j
            QJ_1(1,j,k+3) = Qre_1(offset)
            QJ_2(1,j,k+3) = Qre_2(offset)
            QJ_3(1,j,k+3) = Qre_3(offset)
            QJ_4(1,j,k+3) = Qre_4(offset)
            QJ_5(1,j,k+3) = Qre_5(offset)
        enddo;enddo
      endif
    else
      ! ---------------- SBLI domain: inlet from the TBL planes ----------------
      call MPI_RECV(q1d_h, n, MPI_REAL8, 0, tag_inlet, MPI_COMM_WORLD, mstat, ierr)
      q1d_d = q1d_h
      call reconstruct_sbli_inlet(nx, ny1, ny, nz, 3, q1d_d, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
      ! spanwise mean of the TBL row ny1-1 (primitive), per inlet column
      call span_mean_row(nx, ny, nz, 3, ny1-1, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, mrow_d)
      ! Above y = Ly1 impose a running (exponential) time mean of that row, not the
      ! instantaneous spanwise mean: the latter carries the TBL domain's acoustic
      ! fluctuations as column-uniform density/pressure pulses, which convect down the
      ! SBLI domain as vertical fronts (seen in the first A1 schlieren).
      call running_mean(3, mrow_d, mavg_d, mavg_init)
      !$cuf kernel do(3)<<<*,*>>>
      do k = 4, nz-3
        do j = ny1, ny-1
          do i = 1, 3
            Jt = 1.d0 / Jacobian(i,j)
            rhob = mavg_d(i,1); ub = mavg_d(i,2); vb = mavg_d(i,3); pb = mavg_d(i,4)
            QJ_1(i,j,k) = rhob * Jt
            QJ_2(i,j,k) = rhob * ub * Jt
            QJ_3(i,j,k) = rhob * vb * Jt
            QJ_4(i,j,k) = 0.d0
            QJ_5(i,j,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jt
      enddo;enddo;enddo
    endif

    ! ---------------- outlet (both domains): zero gradient ----------------
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        QJ_1(nx,j,k) = QJ_1(nx-1,j,k)
        QJ_2(nx,j,k) = QJ_2(nx-1,j,k)
        QJ_3(nx,j,k) = QJ_3(nx-1,j,k)
        QJ_4(nx,j,k) = QJ_4(nx-1,j,k)
        QJ_5(nx,j,k) = QJ_5(nx-1,j,k)
    enddo;enddo

    ! Level of the exterior state behind the incident shock at the SBLI top. Upstream
    ! of the shock the outer flow is not at p0 but ~2 % above it (compression by the
    ! growing boundary layer's displacement), so the nominal state 2 fights the actual
    ! flow and radiates weak waves behind the incident shock (seen in A2a). State 2 is
    ! rescaled isentropically by the running mean (tau = 10 delta/u0) of the pressure
    ! next to the top upstream of x_top: p*f, rho*f**(1/gamma), same velocities.
    ! State 1 (upstream exterior) stays nominal: rescaling it too would feed the
    ! measured level back into itself (f drifted upward in the first A2c attempt).
    fp = 1.d0; fr = 1.d0
    if (myrank >= 2) then
      ! spanwise-mean state next to the top, per column (rho, u, v, p)
      call span_mean_row(nx, ny, nz, nx, ny-1, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, itop_d)
      ! Mean pressure of the upstream columns, summed on the host in a fixed order: a
      ! `cuf kernel` sum reduction over (i, k) returned run-to-run different last bits,
      ! which made the whole SBLI domain irreproducible from the first step.
      ptop_h = itop_d(i_up0:i_up1,4)
      psum = sum(ptop_h) / dble(i_up1 - i_up0 + 1) / p0
      if (.not. ftop_init) then
        f_top = psum
        ftop_init = .true.
      else
        f_top = f_top + eps_mavg * (psum - f_top)
      endif
      fp = f_top
      fr = f_top**(1.d0 / gamma)
      ncall_top = ncall_top + 1
      if (mod(ncall_top, 15000) == 1) print '(a,i8,a,f9.5,a,f9.5)', " SBLI top level: call ", ncall_top, &
                                                ", instantaneous p/p0 = ", psum, ", running f = ", f_top
      r2 = r2 * fr; p2 = p2 * fp
      ! Downstream of x_dyn the flow under the top boundary has crossed the reflected
      ! (separation) shock and is not a known uniform state: with a large bubble it is
      ! neither state 2 nor state 3. The former SWTBLI switch between the two, decided
      ! cell by cell from the local v, kept injecting waves (A2c segment 1 schlieren:
      ! a fan of parallel waves from the top where the reflected shock leaves). The
      ! running mean of the local state now makes that decision per column (and gives
      ! the characteristic direction), see the top-boundary kernel below; imposing
      ! that mean itself through the Riemann condition was tried and rejected.
      call running_mean(nx, itop_d, etop_d, etop_init)
    endif

    ! ---------------- top (Riemann invariants) and adiabatic no-slip wall ----------------
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
        uin   = QJ_2(i,ny-1,k) / QJ_1(i,ny-1,k)
        vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
        win   = QJ_4(i,ny-1,k) / QJ_1(i,ny-1,k)
        pin   = gamma_1 * (QJ_5(i,ny-1,k) * Jacobian(i,ny-1) - 0.5d0 * rhoin * (uin**2 + vin**2 + win**2))
        cin   = sqrt(gamma * pin / rhoin)
        if (idyn_d(i) == 1 .and. etop_d(i,3) >= v2half) then
          ! Behind the reflected (separation) shock the flow is supersonic and nearly
          ! parallel to the top: copy the state from row ny-1 at the foot of the
          ! right-running characteristic through (x_i, y_ny), x_src = x_i - dy/tan(mu+phi)
          ! (mean local Mach angle mu and flow angle phi). The reflected shock and the
          ! Mach waves of the same family, which run at about that angle, then leave
          ! the domain instead of being reflected by a normal-direction condition.
          alpha = asin(min(1.d0, sqrt(gamma * etop_d(i,4) / etop_d(i,1)) &
                  / sqrt(etop_d(i,2)**2 + etop_d(i,3)**2))) + atan(etop_d(i,3) / etop_d(i,2))
          xsrc = x_d(i) - dyt / tan(alpha)
          ii = min(i, nx-1)   ! keeps ii+1 inside the row if the characteristic runs downstream
          do while (ii > 1 .and. x_d(ii) > xsrc)
            ii = ii - 1
          enddo
          wa = (xsrc - x_d(ii)) / (x_d(ii+1) - x_d(ii))
          wa = min(max(wa, 0.d0), 1.d0)
          ra  = QJ_1(ii,ny-1,k) * Jacobian(ii,ny-1)
          ua  = QJ_2(ii,ny-1,k) / QJ_1(ii,ny-1,k)
          va  = QJ_3(ii,ny-1,k) / QJ_1(ii,ny-1,k)
          wwa = QJ_4(ii,ny-1,k) / QJ_1(ii,ny-1,k)
          pa  = gamma_1 * (QJ_5(ii,ny-1,k) * Jacobian(ii,ny-1) - 0.5d0 * ra * (ua**2 + va**2 + wwa**2))
          rb2 = QJ_1(ii+1,ny-1,k) * Jacobian(ii+1,ny-1)
          ub2 = QJ_2(ii+1,ny-1,k) / QJ_1(ii+1,ny-1,k)
          vb2 = QJ_3(ii+1,ny-1,k) / QJ_1(ii+1,ny-1,k)
          wb2 = QJ_4(ii+1,ny-1,k) / QJ_1(ii+1,ny-1,k)
          pb2 = gamma_1 * (QJ_5(ii+1,ny-1,k) * Jacobian(ii+1,ny-1) - 0.5d0 * rb2 * (ub2**2 + vb2**2 + wb2**2))
          rhob = (1.d0 - wa) * ra  + wa * rb2
          ub   = (1.d0 - wa) * ua  + wa * ub2
          vb   = (1.d0 - wa) * va  + wa * vb2
          wb   = (1.d0 - wa) * wwa + wa * wb2
          pb   = (1.d0 - wa) * pa  + wa * pb2
        else
        if (idyn_d(i) == 1) then
          ! still under the incident-shock flow (state 2 entering through the top)
          rhox = r2; ux = u2; vx = v2; px = p2
        else
          sx = s_d(i)
          rhox = (1.d0 - sx) * rho0 + sx * r2
          ux   = (1.d0 - sx) * u0   + sx * u2
          vx   =                      sx * v2
          px   = (1.d0 - sx) * p0   + sx * p2
        endif
        cx = sqrt(gamma * px / rhox)
        Rp = vin + 2.d0 * cin * over_gamma_1
        Rm = vx  - 2.d0 * cx  * over_gamma_1
        vb = 0.5d0 * (Rp + Rm)
        cb = 0.25d0 * gamma_1 * (Rp - Rm)
        if (vb >= 0.d0) then
          sb = pin / rhoin**gamma
          ub = uin
          wb = win
        else
          sb = px / rhox**gamma
          ub = ux
          wb = 0.d0
        endif
        rhob = (cb**2 / (gamma * sb))**over_gamma_1
        pb   = sb * rhob**gamma
        endif
        Jt = 1.d0 / Jacobian(i,ny)
        QJ_1(i,ny,k) = rhob * Jt
        QJ_2(i,ny,k) = rhob * ub * Jt
        QJ_3(i,ny,k) = rhob * vb * Jt
        QJ_4(i,ny,k) = rhob * wb * Jt
        QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2 + wb**2)) * Jt
        ! wall: u = 0, dp/dy = 0, drho/dy = 0 (adiabatic)
        QJ_1(i,1,k) = QJ_1(i,2,k)
        QJ_2(i,1,k) = 0.d0
        QJ_3(i,1,k) = 0.d0
        QJ_4(i,1,k) = 0.d0
        p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
        QJ_5(i,1,k) = p_wall * over_gamma_1
    enddo;enddo

    call set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)

    if (myrank == 0) then
      ! hand the planes nre2-5..nre2-3 to the SBLI inlet (through host memory)
      call flatten_rescale(nx, ny, nz, nre2, 3, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, q1d_d)
      q1d_h = q1d_d
      call MPI_SEND(q1d_h, n, MPI_REAL8, 2, tag_inlet, MPI_COMM_WORLD, ierr)
    endif
  end subroutine set_bc


  !> Spanwise mean of the primitive state (rho, u, v, p) on row j, columns 1..n
  subroutine span_mean_row(nx, ny, nz, n, j, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, m)
    integer, intent(in), value   :: nx, ny, nz, n, j
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(out), device :: m(n,4)
    real(8) m1, m2, m3, m4, rr
    integer i, k
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, n
      m1 = 0.d0; m2 = 0.d0; m3 = 0.d0; m4 = 0.d0
      do k = 4, nz-3
        rr = QJ_1(i,j,k) * Jacobian(i,j)
        m1 = m1 + rr
        m2 = m2 + QJ_2(i,j,k) / QJ_1(i,j,k)
        m3 = m3 + QJ_3(i,j,k) / QJ_1(i,j,k)
        m4 = m4 + gamma_1 * Jacobian(i,j) * (QJ_5(i,j,k) - 0.5d0 * (QJ_2(i,j,k)**2 &
             + QJ_3(i,j,k)**2 + QJ_4(i,j,k)**2) / QJ_1(i,j,k))
      enddo
      m(i,1) = m1 / dble(nz-6)
      m(i,2) = m2 / dble(nz-6)
      m(i,3) = m3 / dble(nz-6)
      m(i,4) = m4 / dble(nz-6)
    enddo
  end subroutine span_mean_row


  !> Exponential running mean (time constant tau_mavg) of inst; the first call copies it
  subroutine running_mean(n, inst, avg, init)
    integer, intent(in), value     :: n
    real(8), intent(in), device    :: inst(n,4)
    real(8), intent(inout), device :: avg(n,4)
    logical, intent(inout)         :: init
    integer i, k
    if (.not. init) then
      avg = inst
      init = .true.
    else
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, 4
        do i = 1, n
          avg(i,k) = (1.d0 - eps_mavg) * avg(i,k) + eps_mavg * inst(i,k)
      enddo;enddo
    endif
  end subroutine running_mean


  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(inout), device  :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        mut(1,j,k)  = mut(2,j,k)
        qc2(1,j,k)  = qc2(2,j,k)
        mut(nx,j,k) = mut(nx-1,j,k)
        qc2(nx,j,k) = qc2(nx-1,j,k)
    enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        mut(i,1,k)  = 0.d0
        qc2(i,1,k)  = 0.d0
        mut(i,ny,k) = mut(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,ny-1,k)
    enddo;enddo
    call set_bc_cyclic_z_scalar(nx, ny, nz, mut)
    call set_bc_cyclic_z_scalar(nx, ny, nz, qc2)
  end subroutine set_bc_mut


  subroutine set_bc_cyclic_z_scalar(nx, ny, nz, f)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: f(nx,ny,nz)
    integer i, j
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        f(i,j,1)    = f(i,j,nz-5)
        f(i,j,2)    = f(i,j,nz-4)
        f(i,j,3)    = f(i,j,nz-3)
        f(i,j,nz-2) = f(i,j,4)
        f(i,j,nz-1) = f(i,j,5)
        f(i,j,nz)   = f(i,j,6)
    enddo;enddo
  end subroutine set_bc_cyclic_z_scalar
end module set
