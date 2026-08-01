module set
  use cudafor
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, theta, x_corner, Ma_inf, gamma, R, Pr, &
                          rho_inf, u_inf, v_inf, p_inf
  use set_coordinate, only : set_grid_c_corner, set_metrics_curv
  use set_bc_common
  implicit none
contains
  ! Compression corner reflection test case. Flat plate (y=0) for x <= x_corner,
  ! then ramp at angle theta for x > x_corner. Upper wall flat at y=Ly.
  ! xi-direction: uniform from x=0 (inlet) to x=Lx (outlet), non-periodic.
  ! eta-direction: linear blending from lower to upper wall.
  ! zeta-direction: uniform spanwise.
  subroutine set_grid_c_corner(nx, ny, nz, Lx, Ly, Lz, theta, x_corner, &
                                x_phys, y_phys, zc, dz)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz, theta, x_corner
    real(8), intent(out) :: x_phys(nx,ny), y_phys(nx,ny)
    real(8), intent(out) :: zc(nz), dz(nz-1)
    real(8) :: dxi, tj, xi(nx), yi(nx), xo(nx), yo(nx)
    integer :: i, j, k
    ! Interior cells i=2..nx-1 span xi from 0 to Lx
    dxi = Lx / dble(nx-2)
    do i = 2, nx-1
      xi(i) = dxi * dble(i-2)
      if (xi(i) <= x_corner) then
        yi(i) = 0.d0
      else
        yi(i) = (xi(i) - x_corner) * tan(theta)
      endif
      xo(i) = xi(i)
      yo(i) = Ly
    enddo
    xi(1)  = xi(2)    - dxi;  yi(1)  = 0.d0
    xi(nx) = xi(nx-1) + dxi;  yi(nx) = yi(nx-1) + dxi * tan(theta)
    xo(1)  = xi(1);   yo(1)  = Ly
    xo(nx) = xi(nx);  yo(nx) = Ly
    do j = 1, ny
      tj = dble(j-1) / dble(ny-1)
      do i = 1, nx
        x_phys(i,j) = xi(i) + tj * (xo(i) - xi(i))
        y_phys(i,j) = yi(i) + tj * (yo(i) - yi(i))
      enddo
    enddo
    do k = 1, nz
      zc(k) = Lz * dble(k-1) / dble(nz-1)
    enddo
    do k = 1, nz-1
      dz(k) = zc(k+1) - zc(k)
    enddo
  end subroutine set_grid_c_corner


  !> Set up body-fitted grid for compression corner reflection test case
  !> Lower boundary: flat plate for x <= x_corner, then ramp at angle theta
  !> Upper boundary: flat wall at y = Ly
  !> Populates x_phys_g, y_phys_g with physical coordinates
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz, &
                      x_phys_g, y_phys_g)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: Lx, Ly, Lz
    real(8), intent(out), allocatable :: xc(:), yc(:), zc(:)
    real(8), intent(out), allocatable :: dx(:), dy(:), dz(:)
    real(8), intent(out), allocatable :: x_phys_g(:,:), y_phys_g(:,:)
    integer i, j
    allocate(xc(nx), yc(ny), zc(nz))
    allocate(dx(nx-1), dy(ny-1), dz(nz-1))
    allocate(x_phys_g(nx,ny), y_phys_g(nx,ny))
    ! Generate body-fitted grid: flat plate + compression corner ramp
    call set_grid_c_corner(nx, ny, nz, Lx, Ly, Lz, theta, x_corner, &
                           x_phys_g, y_phys_g, zc, dz)
    ! Dummy 1D coordinate arrays (computational spacing is uniform: Δξ = Δη = 1)
    xc = (/ (dble(i), i=1,nx) /)
    yc = (/ (dble(j), j=1,ny) /)
    ! zc already filled by set_grid_c_corner
    ! Uniform computational spacing
    dx = 1.d0
    dy = 1.d0
    ! dz already filled
    if (myrank == 0) then
      print *, "CORN grid: nx=", nx, " ny=", ny, " nz=", nz
      print *, "  Domain: Lx=", Lx, " Ly=", Ly, " Lz=", Lz
      print *, "  Corner: x_corner=", x_corner, " theta=", theta*180.d0/acos(-1.d0), " deg"
    endif
  end subroutine set_grid


  !> Compute 2D metrics and face normals for curvilinear body-fitted grid.
  !> After set_metrics_curv fills the interior, this routine fills boundary and
  !> corner cells using 1-sided stencils appropriate for the CORN case:
  !>   eta j=1 (lower wall) and j=ny (upper wall): 1-sided forward/backward
  !>   xi  i=1 (inlet)      and i=nx (outlet):     1-sided forward/backward
  !> No periodic ghost-cell enforcement: xi boundaries are Dirichlet/zero-gradient.
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
    call set_metrics_curv(nx, ny, x_phys_g, y_phys_g, &
                          n_xi_x_cpu, n_xi_y_cpu, n_eta_x_cpu, n_eta_y_cpu, &
                          xi_x_cpu, xi_y_cpu, eta_x_cpu, eta_y_cpu, Jac_cpu)
    ! eta-boundary cells: j=1 (lower wall, 1-sided forward) and j=ny (upper wall, backward)
    do m = 2, nx-1
      xxi  = 0.5d0*(x_phys_g(m+1,1)-x_phys_g(m-1,1))
      yxi  = 0.5d0*(y_phys_g(m+1,1)-y_phys_g(m-1,1))
      xeta = x_phys_g(m,2)-x_phys_g(m,1);  yeta = y_phys_g(m,2)-y_phys_g(m,1)
      J2             = xxi*yeta - xeta*yxi
      Jac_cpu(m,1)   = J2
      xi_x_cpu(m,1)  =  yeta/J2;  xi_y_cpu(m,1)  = -xeta/J2
      eta_x_cpu(m,1) = -yxi /J2;  eta_y_cpu(m,1) =  xxi /J2
      xxi  = 0.5d0*(x_phys_g(m+1,ny)-x_phys_g(m-1,ny))
      yxi  = 0.5d0*(y_phys_g(m+1,ny)-y_phys_g(m-1,ny))
      xeta = x_phys_g(m,ny)-x_phys_g(m,ny-1);  yeta = y_phys_g(m,ny)-y_phys_g(m,ny-1)
      J2              = xxi*yeta - xeta*yxi
      Jac_cpu(m,ny)   = J2
      xi_x_cpu(m,ny)  =  yeta/J2;  xi_y_cpu(m,ny)  = -xeta/J2
      eta_x_cpu(m,ny) = -yxi /J2;  eta_y_cpu(m,ny) =  xxi /J2
    enddo
    ! xi-boundary cells: i=1 (inlet, 1-sided fwd) and i=nx (outlet, backward), all j
    do n = 1, ny
      xxi = x_phys_g(2,n)-x_phys_g(1,n);  yxi = y_phys_g(2,n)-y_phys_g(1,n)
      if (n == 1) then
        xeta = x_phys_g(1,2)-x_phys_g(1,1);  yeta = y_phys_g(1,2)-y_phys_g(1,1)
      elseif (n == ny) then
        xeta = x_phys_g(1,ny)-x_phys_g(1,ny-1);  yeta = y_phys_g(1,ny)-y_phys_g(1,ny-1)
      else
        xeta = 0.5d0*(x_phys_g(1,n+1)-x_phys_g(1,n-1))
        yeta = 0.5d0*(y_phys_g(1,n+1)-y_phys_g(1,n-1))
      endif
      J2            = xxi*yeta - xeta*yxi
      Jac_cpu(1,n)  = J2
      xi_x_cpu(1,n)  =  yeta/J2;  xi_y_cpu(1,n)  = -xeta/J2
      eta_x_cpu(1,n) = -yxi /J2;  eta_y_cpu(1,n) =  xxi /J2
      xxi = x_phys_g(nx,n)-x_phys_g(nx-1,n);  yxi = y_phys_g(nx,n)-y_phys_g(nx-1,n)
      if (n == 1) then
        xeta = x_phys_g(nx,2)-x_phys_g(nx,1);  yeta = y_phys_g(nx,2)-y_phys_g(nx,1)
      elseif (n == ny) then
        xeta = x_phys_g(nx,ny)-x_phys_g(nx,ny-1);  yeta = y_phys_g(nx,ny)-y_phys_g(nx,ny-1)
      else
        xeta = 0.5d0*(x_phys_g(nx,n+1)-x_phys_g(nx,n-1))
        yeta = 0.5d0*(y_phys_g(nx,n+1)-y_phys_g(nx,n-1))
      endif
      J2             = xxi*yeta - xeta*yxi
      Jac_cpu(nx,n)  = J2
      xi_x_cpu(nx,n)  =  yeta/J2;  xi_y_cpu(nx,n)  = -xeta/J2
      eta_x_cpu(nx,n) = -yxi /J2;  eta_y_cpu(nx,n) =  xxi /J2
    enddo
  end subroutine set_metrics


  !> Initialize flow field with supersonic free-stream conditions
  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    real(8) :: E_inf = p_inf / (gamma - 1.d0) + 0.5d0 * rho_inf * (u_inf**2 + v_inf**2)
    ! Initialize all cells to free-stream
    Q = 0.d0
    Q(:,:,:,1) = rho_inf
    Q(:,:,:,2) = rho_inf * u_inf
    Q(:,:,:,3) = rho_inf * v_inf
    Q(:,:,:,4) = 0.d0  ! w = 0 (quasi-2D in z)
    Q(:,:,:,5) = E_inf
    if (myrank == 0) then
      print *, "Flow: M_inf =", Ma_inf, " (supersonic)"
      print *, "  rho_inf =", rho_inf, " p_inf =", p_inf, " E_inf =", E_inf
    endif
  end subroutine set_init


  !> Apply boundary conditions
  !> (a) xi inlet (i=1): Dirichlet free-stream
  !> (b) xi outlet (i=nx): zero-gradient extrapolation
  !> (c) eta lower wall (j=1): Euler slip wall (flat plate)
  !> (d) eta upper wall (j=ny): Euler slip wall (flat reflector for shock reflection)
  !> (e) z-periodic
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, eta_x, eta_y, Q_1, Q_2, Q_3, Q_4, Q_5)
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in), device :: Jacobian(nx,ny)
    real(8), intent(in), device :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer :: i, j, k
    real(8) :: nxw, nyw, nmag, u_int, v_int, u_n, Jratio
    real(8) :: E_inf = p_inf / (gamma - 1.d0) + 0.5d0 * rho_inf * (u_inf**2 + v_inf**2)
    ! (a) xi inlet ghost (i=1): Dirichlet free-stream
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do j = 1, ny
        Q_1(1,j,k) = rho_inf / Jacobian(1,j)
        Q_2(1,j,k) = rho_inf * u_inf / Jacobian(1,j)
        Q_3(1,j,k) = rho_inf * v_inf / Jacobian(1,j)
        Q_4(1,j,k) = 0.d0
        Q_5(1,j,k) = E_inf / Jacobian(1,j)
      enddo
    enddo
    ! (b) xi outlet ghost (i=nx): zero-gradient (physical) extrapolation
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do j = 1, ny
        Jratio = Jacobian(nx-1,j) / Jacobian(nx,j)
        Q_1(nx,j,k) = Q_1(nx-1,j,k) * Jratio
        Q_2(nx,j,k) = Q_2(nx-1,j,k) * Jratio
        Q_3(nx,j,k) = Q_3(nx-1,j,k) * Jratio
        Q_4(nx,j,k) = Q_4(nx-1,j,k) * Jratio
        Q_5(nx,j,k) = Q_5(nx-1,j,k) * Jratio
      enddo
    enddo
    ! (c) eta lower wall (j=1): Euler slip wall (flat plate, normal = eta direction)
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do i = 1, nx
        nxw   = eta_x(i,1);  nyw = eta_y(i,1)
        nmag  = sqrt(nxw*nxw + nyw*nyw)
        nxw   = nxw / nmag;  nyw = nyw / nmag
        u_int = Q_2(i,2,k) / Q_1(i,2,k)
        v_int = Q_3(i,2,k) / Q_1(i,2,k)
        u_n   = u_int*nxw + v_int*nyw
        Jratio = Jacobian(i,2) / Jacobian(i,1)
        Q_1(i,1,k) = Q_1(i,2,k) * Jratio
        Q_2(i,1,k) = (Q_2(i,2,k) - 2.d0*u_n*nxw*Q_1(i,2,k)) * Jratio
        Q_3(i,1,k) = (Q_3(i,2,k) - 2.d0*u_n*nyw*Q_1(i,2,k)) * Jratio
        Q_4(i,1,k) = Q_4(i,2,k) * Jratio
        Q_5(i,1,k) = Q_5(i,2,k) * Jratio
      enddo
    enddo
    ! (d) eta upper wall (j=ny): Euler slip wall (flat reflector, normal = eta direction)
    !$cuf kernel do(2) <<<*,(16,16)>>>
    do k = 1, nz
      do i = 1, nx
        nxw   = eta_x(i,ny);  nyw = eta_y(i,ny)
        nmag  = sqrt(nxw*nxw + nyw*nyw)
        nxw   = nxw / nmag;  nyw = nyw / nmag
        u_int = Q_2(i,ny-1,k) / Q_1(i,ny-1,k)
        v_int = Q_3(i,ny-1,k) / Q_1(i,ny-1,k)
        u_n   = u_int*nxw + v_int*nyw
        Jratio = Jacobian(i,ny-1) / Jacobian(i,ny)
        Q_1(i,ny,k) = Q_1(i,ny-1,k) * Jratio
        Q_2(i,ny,k) = (Q_2(i,ny-1,k) - 2.d0*u_n*nxw*Q_1(i,ny-1,k)) * Jratio
        Q_3(i,ny,k) = (Q_3(i,ny-1,k) - 2.d0*u_n*nyw*Q_1(i,ny-1,k)) * Jratio
        Q_4(i,ny,k) = Q_4(i,ny-1,k) * Jratio
        Q_5(i,ny,k) = Q_5(i,ny-1,k) * Jratio
      enddo
    enddo
    ! (e) z-periodic: ghost cells k=1 and k=nz wrap around interior k=2..nz-1
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
