module set
  use cudafor
  use mod_globals, only : nx, ny, nz, Lz, chord, aoa, far_r, Ma_inf, gamma, R, Pr
  use set_coordinate, only : set_metrics_curv
  use set_bc_common
  implicit none
contains
  ! Generate O-grid around NACA 0012.  xi (i) wraps around airfoil (periodic,
  ! i=1 and i=nx are ghost cells), eta (j) is wall-normal (j=1 wall, j=ny far-field),
  ! zeta (k) is uniform spanwise.
  ! Outer boundary: circle of radius far_r_in centered at (chord_in/2, 0).
  ! Wall-normal distribution uses geometric stretching (ratio 1.05).
  subroutine set_grid_c_wing(nx_in, ny_in, nz_in, chord_in, aoa_in, far_r_in, Lz_in, &
                              x_phys, y_phys, zc, dz)
    integer, intent(in)  :: nx_in, ny_in, nz_in
    real(8), intent(in)  :: chord_in, aoa_in, far_r_in, Lz_in
    real(8), intent(out) :: x_phys(nx_in,ny_in), y_phys(nx_in,ny_in)
    real(8), intent(out) :: zc(nz_in), dz(nz_in-1)
    ! stretch=1.05: first wall-normal cell ≈0.011 chord (Euler; no need for viscous clustering)
    real(8), parameter   :: stretch = 1.05d0
    real(8) :: pi, phi, xn, yt, xi(nx_in), yi(nx_in), xo(nx_in), yo(nx_in)
    real(8) :: tj, xp, yp, xr, yr
    integer :: i, j, k
    pi = 4.d0 * atan(1.d0)
    ! Interior cells i=2..nx_in-1 span a full 2*pi so that ghost cells can be
    ! exact periodic copies: x_phys(1,j)=x_phys(nx_in-1,j), x_phys(nx_in,j)=x_phys(2,j).
    do i = 2, nx_in-1
      phi  = 2.d0*pi*dble(i-2)/dble(nx_in-2)
      xn   = 0.5d0*(1.d0 + cos(phi))
      yt   = (0.12d0/0.2d0) * (0.2969d0*sqrt(max(xn, 1.d-14)) &
             - 0.1260d0*xn - 0.3516d0*xn**2 + 0.2843d0*xn**3 - 0.1015d0*xn**4)
      xi(i) = xn*chord_in
      yi(i) = merge(-yt*chord_in, yt*chord_in, phi <= pi)
    enddo
    xi(1) = xi(nx_in-1);  yi(1) = yi(nx_in-1)
    xi(nx_in) = xi(2);    yi(nx_in) = yi(2)
    do i = 2, nx_in-1
      phi   = 2.d0*pi*dble(i-2)/dble(nx_in-2)
      xo(i) = 0.5d0*chord_in + far_r_in*cos(phi)
      yo(i) = -far_r_in*sin(phi)
    enddo
    xo(1) = xo(nx_in-1);  yo(1) = yo(nx_in-1)
    xo(nx_in) = xo(2);    yo(nx_in) = yo(2)
    do j = 1, ny_in
      tj = (stretch**(j-1) - 1.d0) / (stretch**(ny_in-1) - 1.d0)
      do i = 1, nx_in
        x_phys(i,j) = xi(i) + tj*(xo(i) - xi(i))
        y_phys(i,j) = yi(i) + tj*(yo(i) - yi(i))
      enddo
    enddo
    if (abs(aoa_in) > 1.d-14) then
      xp = 0.25d0*chord_in;  yp = 0.d0
      do j = 1, ny_in
        do i = 1, nx_in
          xr = cos(aoa_in)*(x_phys(i,j)-xp) - sin(aoa_in)*(y_phys(i,j)-yp) + xp
          yr = sin(aoa_in)*(x_phys(i,j)-xp) + cos(aoa_in)*(y_phys(i,j)-yp) + yp
          x_phys(i,j) = xr;  y_phys(i,j) = yr
        enddo
      enddo
    endif
    do k = 1, nz_in
      zc(k) = Lz_in*dble(k-1)/dble(nz_in-1)
    enddo
    do k = 1, nz_in-1
      dz(k) = zc(k+1) - zc(k)
    enddo
  end subroutine set_grid_c_wing


  !> Set up O-grid around NACA 0012.
  !> Lx and Ly are dummies passed by main_curv; grid geometry comes from mod_globals
  !> (chord, aoa, far_r).  Lz is used for spanwise extent.
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz, &
                      x_phys_g, y_phys_g)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: Lx, Ly, Lz
    real(8), intent(out), allocatable :: xc(:), yc(:), zc(:)
    real(8), intent(out), allocatable :: dx(:), dy(:), dz(:)
    real(8), intent(out), allocatable :: x_phys_g(:,:), y_phys_g(:,:)
    integer :: i, j
    allocate(xc(nx), yc(ny), zc(nz))
    allocate(dx(nx-1), dy(ny-1), dz(nz-1))
    allocate(x_phys_g(nx,ny), y_phys_g(nx,ny))
    call set_grid_c_wing(nx, ny, nz, chord, aoa, far_r, Lz, &
                         x_phys_g, y_phys_g, zc, dz)
    xc = (/ (dble(i), i=1,nx) /)
    yc = (/ (dble(j), j=1,ny) /)
    dx = 1.d0
    dy = 1.d0
    if (myrank == 0) then
      print *, "NACA 0012 O-grid: nx=", nx, " ny=", ny, " nz=", nz
      print *, "  chord=", chord, " far_r=", far_r, " aoa=", aoa*180.d0/acos(-1.d0), " deg"
      print *, "  Ma_inf=", Ma_inf
    endif
  end subroutine set_grid


  !> Compute 2D metrics for NACA O-grid.
  !> set_metrics_curv fills ALL face normals (n_xi, n_eta) and cell-centre metrics
  !> for interior m=2..nx-1, n=2..ny-1.
  !> This routine then:
  !>   (1) fills eta boundary cells j=1 (airfoil wall) and j=ny (far-field) using
  !>       1-sided stencils for m=2..nx-1, and
  !>   (2) fills xi ghost cells i=1 and i=nx by periodic copy from i=nx-1 and i=2.
  subroutine set_metrics(nx, ny, x_phys_g, y_phys_g, &
                         n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                         xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Jac_cpu)
    integer, intent(in)  :: nx, ny
    real(8), intent(in)  :: x_phys_g(nx,ny), y_phys_g(nx,ny)
    real(8), intent(out), allocatable :: n_xi_x_cpu(:,:), n_xi_y_cpu(:,:)
    real(8), intent(out), allocatable :: n_eta_x_cpu(:,:), n_eta_y_cpu(:,:)
    real(8), intent(out), allocatable :: xi_x_cpu(:,:), xi_y_cpu(:,:)
    real(8), intent(out), allocatable :: eta_x_cpu(:,:), eta_y_cpu(:,:)
    real(8), intent(out), allocatable :: Jac_cpu(:,:)
    real(8) :: xxi, yxi, xeta, yeta, J2
    integer :: m, n
    allocate(n_xi_x_cpu(nx-1,ny-2), n_xi_y_cpu(nx-1,ny-2))
    allocate(n_eta_x_cpu(nx-2,ny-1), n_eta_y_cpu(nx-2,ny-1))
    allocate(xi_x_cpu(nx,ny), xi_y_cpu(nx,ny))
    allocate(eta_x_cpu(nx,ny), eta_y_cpu(nx,ny))
    allocate(Jac_cpu(nx,ny))
    ! Fill ALL face normals and interior cell-centre metrics
    call set_metrics_curv(nx, ny, x_phys_g, y_phys_g, &
                          n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                          xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Jac_cpu)
    ! (1) eta boundary cells: j=1 (airfoil wall, 1-sided forward) and
    !     j=ny (far-field, 1-sided backward), for interior xi (m=2..nx-1)
    do m = 2, nx-1
      ! j=1: forward difference in eta
      xxi  = 0.5d0*(x_phys_g(m+1,1)-x_phys_g(m-1,1))
      yxi  = 0.5d0*(y_phys_g(m+1,1)-y_phys_g(m-1,1))
      xeta = x_phys_g(m,2)-x_phys_g(m,1);  yeta = y_phys_g(m,2)-y_phys_g(m,1)
      J2             = xxi*yeta - xeta*yxi
      Jac_cpu(m,1)   = J2
      xi_x_cpu(m,1)  =  yeta/J2;  xi_y_cpu(m,1)  = -xeta/J2
      eta_x_cpu(m,1) = -yxi /J2;  eta_y_cpu(m,1) =  xxi /J2
      ! j=ny: backward difference in eta
      xxi  = 0.5d0*(x_phys_g(m+1,ny)-x_phys_g(m-1,ny))
      yxi  = 0.5d0*(y_phys_g(m+1,ny)-y_phys_g(m-1,ny))
      xeta = x_phys_g(m,ny)-x_phys_g(m,ny-1);  yeta = y_phys_g(m,ny)-y_phys_g(m,ny-1)
      J2              = xxi*yeta - xeta*yxi
      Jac_cpu(m,ny)   = J2
      xi_x_cpu(m,ny)  =  yeta/J2;  xi_y_cpu(m,ny)  = -xeta/J2
      eta_x_cpu(m,ny) = -yxi /J2;  eta_y_cpu(m,ny) =  xxi /J2
    enddo
    ! (2) xi ghost cells: periodic copy.
    !     x_phys_g(1,:) = x_phys_g(nx-1,:) and x_phys_g(nx,:) = x_phys_g(2,:),
    !     so metrics at ghost cells equal metrics at the corresponding interior cells.
    !     Step (1) has already filled j=1 and j=ny for m=2..nx-1, so corners are
    !     covered correctly here.
    do n = 1, ny
      Jac_cpu(1,n)    = Jac_cpu(nx-1,n)
      xi_x_cpu(1,n)   = xi_x_cpu(nx-1,n)
      xi_y_cpu(1,n)   = xi_y_cpu(nx-1,n)
      eta_x_cpu(1,n)  = eta_x_cpu(nx-1,n)
      eta_y_cpu(1,n)  = eta_y_cpu(nx-1,n)
      Jac_cpu(nx,n)   = Jac_cpu(2,n)
      xi_x_cpu(nx,n)  = xi_x_cpu(2,n)
      xi_y_cpu(nx,n)  = xi_y_cpu(2,n)
      eta_x_cpu(nx,n) = eta_x_cpu(2,n)
      eta_y_cpu(nx,n) = eta_y_cpu(2,n)
    enddo
  end subroutine set_metrics


  !> Initialize flow to uniform subsonic free-stream.
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : Ma_inf, rho_inf, p_inf, T_inf
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    real(8) u_inf, v_inf, E_inf
    u_inf   = Ma_inf * sqrt(gamma * R * T_inf)
    v_inf   = 0.d0
    E_inf   = p_inf / (gamma - 1.d0) + 0.5d0 * rho_inf * (u_inf**2 + v_inf**2)
    Q(:,:,:,1) = rho_inf
    Q(:,:,:,2) = rho_inf * u_inf
    Q(:,:,:,3) = rho_inf * v_inf
    Q(:,:,:,4) = 0.d0
    Q(:,:,:,5) = E_inf
    if (myrank == 0) then
      print *, "Flow: M_inf =", Ma_inf, " aoa =", aoa*180.d0/acos(-1.d0), " deg"
      print *, "  rho_inf =", rho_inf, " p_inf =", p_inf, " E_inf =", E_inf
    endif
  end subroutine set_init


  !> Apply boundary conditions on the GPU.
  !> (a) xi periodic:  i=1  ghost ← i=nx-1 interior; i=nx ghost ← i=2 interior
  !> (b) eta j=1:      Euler slip wall on airfoil surface
  !> (c) eta j=ny:     Dirichlet far-field free-stream
  !> (d) z-periodic:   k=1 and k=nz ghost cells
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, eta_x, eta_y, Q_1, Q_2, Q_3, Q_4, Q_5)
    use mod_globals, only : Ma_inf, rho_inf, p_inf, T_inf
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in),    device :: Jacobian(nx,ny)
    real(8), intent(in),    device :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer :: i, j, k
    real(8) :: u_inf, v_inf, E_inf
    real(8) :: nxw, nyw, nmag, u_int, v_int, u_n, Jratio
    u_inf   = Ma_inf * sqrt(gamma * R * T_inf)! * cos(aoa)
    v_inf   = 0.d0!Ma_inf * sqrt(gamma * R * T_inf) * sin(aoa)
    E_inf   = p_inf / (gamma - 1.d0) + 0.5d0 * rho_inf * (u_inf**2 + v_inf**2)
    ! (a) xi periodic: O-grid seam at trailing edge
    !     Jacobian ratio = 1 exactly since grid is periodic, so direct copy.
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do j = 1, ny
        Q_1(1,j,k) = Q_1(nx-1,j,k)
        Q_2(1,j,k) = Q_2(nx-1,j,k)
        Q_3(1,j,k) = Q_3(nx-1,j,k)
        Q_4(1,j,k) = Q_4(nx-1,j,k)
        Q_5(1,j,k) = Q_5(nx-1,j,k)
        Q_1(nx,j,k) = Q_1(2,j,k)
        Q_2(nx,j,k) = Q_2(2,j,k)
        Q_3(nx,j,k) = Q_3(2,j,k)
        Q_4(nx,j,k) = Q_4(2,j,k)
        Q_5(nx,j,k) = Q_5(2,j,k)
      enddo
    enddo
    ! (b) eta j=1: wall on airfoil
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do i = 1, nx
        ! Euler slip wall
        !nxw   = eta_x(i,1);  nyw = eta_y(i,1)
        !nmag  = sqrt(nxw*nxw + nyw*nyw)
        !nxw   = nxw / nmag;  nyw = nyw / nmag
        !u_int = Q_2(i,2,k) / Q_1(i,2,k)
        !v_int = Q_3(i,2,k) / Q_1(i,2,k)
        !u_n   = u_int*nxw + v_int*nyw
        Jratio = Jacobian(i,2) / Jacobian(i,1)
        !Q_1(i,1,k) = Q_1(i,2,k) * Jratio
        !Q_2(i,1,k) = (Q_2(i,2,k) - 2.d0*u_n*nxw*Q_1(i,2,k)) * Jratio
        !Q_3(i,1,k) = (Q_3(i,2,k) - 2.d0*u_n*nyw*Q_1(i,2,k)) * Jratio
        !Q_4(i,1,k) = Q_4(i,2,k) * Jratio
        !Q_5(i,1,k) = Q_5(i,2,k) * Jratio
        ! NS no-slip wall
        Q_1(i,1,k) = Q_1(i,2,k) * Jratio
        Q_2(i,1,k) =-Q_2(i,2,k) * Jratio
        Q_3(i,1,k) =-Q_3(i,2,k) * Jratio
        Q_4(i,1,k) =-Q_4(i,2,k) * Jratio
        Q_5(i,1,k) = Q_5(i,2,k) * Jratio
      enddo
    enddo
    ! (c) eta j=ny: Dirichlet far-field free-stream
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do i = 1, nx
        Q_1(i,ny,k) = rho_inf / Jacobian(i,ny)
        Q_2(i,ny,k) = rho_inf * u_inf / Jacobian(i,ny)
        Q_3(i,ny,k) = rho_inf * v_inf / Jacobian(i,ny)
        Q_4(i,ny,k) = 0.d0
        Q_5(i,ny,k) = E_inf / Jacobian(i,ny)
      enddo
    enddo
    ! (d) z-periodic
    Q_1(:,:,1)  = Q_1(:,:,nz-1)
    Q_2(:,:,1)  = Q_2(:,:,nz-1)
    Q_3(:,:,1)  = Q_3(:,:,nz-1)
    Q_4(:,:,1)  = Q_4(:,:,nz-1)
    Q_5(:,:,1)  = Q_5(:,:,nz-1)
    Q_1(:,:,nz) = Q_1(:,:,2)
    Q_2(:,:,nz) = Q_2(:,:,2)
    Q_3(:,:,nz) = Q_3(:,:,2)
    Q_4(:,:,nz) = Q_4(:,:,2)
    Q_5(:,:,nz) = Q_5(:,:,2)
  end subroutine set_bc
end module set
