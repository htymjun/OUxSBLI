module set
  use cudafor
  use mod_globals, only : gamma, rho0, u0, p0, i_LE
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1
  use set_bc_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dz1
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

    ! z: uniform and periodic, matching set_grid_cyclic6_3D's own convention
    ! (dz1 = Lz/(nz-6), z(k) = dz1*(k-4)); no absolute origin to respect here
    ! (unlike x), so the same formula is valid for every k, ghost or interior.
    dz1 = Lz / dble(nz-6)
    do k = 1, nz
      z(k) = dz1 * dble(k-4)
    enddo
    dz(:) = dz1
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k

    ! Uniform freestream everywhere. The inlet column i=1 is never updated
    ! afterwards (interior-only RK kernels, no inlet BC in set_bc), so it stays
    ! frozen at the freestream state and acts as a Dirichlet inflow.
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          Q(i,j,k,1) = rho0
          Q(i,j,k,2) = rho0 * u0
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0 ! no spanwise velocity: the BL is uniform in z by construction
          Q(i,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2
      enddo;enddo;enddo
    do k = 1, nz
      do i = i_LE, nx
        Q(i,1,k,1) = rho0
        Q(i,1,k,2) = 0.d0
        Q(i,1,k,3) = 0.d0
        Q(i,1,k,4) = 0.d0
        Q(i,1,k,5) = Q(i,2,k,5)
    enddo;enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    integer i, j, k
    real(8) :: p_wall
    ! top far field: interior (j = ny-1) state and boundary state
    real(8) :: rhoin, uin, vin, win, pin, rhob, ub, vb, wb, pb, sb
    ! cache
    real(8) Jacobian_tmp
    ! freestream entropy at the top boundary
    real(8), parameter :: s0 = p0 / rho0**gamma

    ! inlet: active Dirichlet freestream, unlike 2D_solver/BL's "leave i=1
    ! frozen, no inlet BC" idiom. Confirmed by direct testing: with BC_X=True
    ! and i=1 left untouched, the boundary-aware x convective kernel corrupts
    ! columns i=2..7 (exactly ORDER columns in from the boundary) into NaN
    ! within a single RK step -- reproduced with both SLAU and KEEP, with and
    ! without viscosity, and at every nx tried, so it is not scheme, physics,
    ! or resolution dependent. No existing 3D case has ever combined BC_X=True
    ! with an untouched i=1: SBLI and TBL (the only other BC_X=True 3D cases)
    ! both actively rewrite their leftmost columns every step via rescaling.
    ! Re-asserting the same freestream state here every step is numerically a
    ! no-op relative to leaving it frozen -- and empirically avoids whatever
    ! that untested combination trips in the boundary-aware kernel.
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        Jacobian_tmp = 1.d0 / Jacobian(1,j)
        QJ_1(1,j,k) = rho0 * Jacobian_tmp
        QJ_2(1,j,k) = rho0 * u0 * Jacobian_tmp
        QJ_3(1,j,k) = 0.d0
        QJ_4(1,j,k) = 0.d0
        QJ_5(1,j,k) = (p0 / (gamma - 1.d0) + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 2, ny-1
        ! outlet: subsonic pressure outflow (see 2D_solver/BL/set.f90 for the
        ! full rationale -- p = p0 is the correct single condition for a
        ! zero-incidence flat plate). w extrapolates like the other momenta.
        QJ_1(nx,j,k) = QJ_1(nx-1,j,k)
        QJ_2(nx,j,k) = QJ_2(nx-1,j,k)
        QJ_3(nx,j,k) = QJ_3(nx-1,j,k)
        QJ_4(nx,j,k) = QJ_4(nx-1,j,k)
        QJ_5(nx,j,k) = p0 * over_gamma_1 / Jacobian(nx,j) &
                     + 0.5d0 * (QJ_2(nx,j,k)**2 + QJ_3(nx,j,k)**2 + QJ_4(nx,j,k)**2) / QJ_1(nx,j,k)
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        ! top: constant-pressure far field (see 2D_solver/BL/set.f90). v is
        ! always extrapolated (wall-normal, zero-shear far field); entropy and
        ! tangential momentum (u, and here also w) are upwinded on sign(v).
        rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
        uin   = QJ_2(i,ny-1,k) / QJ_1(i,ny-1,k)
        vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
        win   = QJ_4(i,ny-1,k) / QJ_1(i,ny-1,k)
        pin   = gamma_1 * (QJ_5(i,ny-1,k) * Jacobian(i,ny-1) &
                - 0.5d0 * rhoin * (uin**2 + vin**2 + win**2))

        pb = p0
        vb = vin
        if (vb >= 0.d0) then
          ! outflow: entropy and tangential momentum leave the domain
          sb = pin / rhoin**gamma
          ub = uin
          wb = win
        else
          ! inflow: both come from the freestream (rhob then reduces to rho0)
          sb = s0
          ub = u0
          wb = 0.d0
        endif
        rhob = (pb / sb)**over_gamma

        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        QJ_1(i,ny,k) = rhob * Jacobian_tmp
        QJ_2(i,ny,k) = rhob * ub * Jacobian_tmp
        QJ_3(i,ny,k) = rhob * vb * Jacobian_tmp
        QJ_4(i,ny,k) = rhob * wb * Jacobian_tmp
        QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2 + wb**2)) * Jacobian_tmp
    enddo;enddo

    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        if (i < i_LE) then
          ! symmetry (slip) wall upstream of the leading edge: v (wall-normal)
          ! forced to zero and its kinetic energy removed from the interior
          ! state; u and w (both tangential to this plane) are Neumann-copied.
          QJ_1(i,1,k) = QJ_1(i,2,k)
          QJ_2(i,1,k) = QJ_2(i,2,k)
          QJ_3(i,1,k) = 0.d0
          QJ_4(i,1,k) = QJ_4(i,2,k)
          QJ_5(i,1,k) = QJ_5(i,2,k) - 0.5d0 * QJ_3(i,2,k)**2 / QJ_1(i,2,k)
        else
          ! no-slip wall: all three velocity components vanish; wall pressure
          ! is recovered by removing the full (u,v,w) kinetic energy from the
          ! interior state.
          QJ_1(i,1,k) = QJ_1(i,2,k)
          QJ_2(i,1,k) = 0.d0
          QJ_3(i,1,k) = 0.d0
          QJ_4(i,1,k) = 0.d0
          p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
          QJ_5(i,1,k) = p_wall * over_gamma_1
        endif
    enddo;enddo

    call set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
  end subroutine set_bc


  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut
end module set
