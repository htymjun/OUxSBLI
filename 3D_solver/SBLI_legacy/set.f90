module set
  use mod_globals, only : id_rescale, nx, ny, nz, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-1)
    dy1 = dx1
    dz1 = Lz / dble(nz-1)
    
    x(1) = Lx * 0.5d0 * dble(myrank)
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo

    print *, "myrank is ", myrank, " x(1) = ", x(1)

    y(1) = 0.d0
    do j = 1, ny
      if (y(j) <= 7.5e-3) then
        ! DNS 0.05: yp=0.5
        dy(j)  = min(1.d0, max(0.075d0, dble(j)/dble(128))) * dy1
        y(j+1) = y(j) + dy(j)
      else
        ! buffer
        dy(j)  = 2.d0 * dy1
        y(j+1) = y(j) + dy(j)
      endif
    enddo

    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid

  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer i, j, k
    real(8) :: blt0 = 0.5d0 * blt
    real(8) :: Cp   = gamma * R / (gamma - 1.d0), rf = 0.89d0
    real(8) :: eta, rho, u, v, w, T, Tw, Taw, p_wall
    ! random
    real(8) :: std, ustd, vstd, wstd, Tstd
    real(8), allocatable :: randum(:,:,:,:)
    integer ir, jr, kr, nxr, nyr, nzr

    ! generate randum
    nxr = (nx+9) / 10
    nyr = (ny+1) / 2
    nzr = (nz+1) / 2
    allocate(randum(nxr,nyr,nzr,4))

    do k = 1, nzr
      do j = 1, nyr
        do i = 1, nxr
          call random_number(randum(i,j,k,1))
          call random_number(randum(i,j,k,2))
          call random_number(randum(i,j,k,3))
          call random_number(randum(i,j,k,4))
          randum(i,j,k,1) = 2.d0 * randum(i,j,k,1) - 1.d0
          randum(i,j,k,2) = 2.d0 * randum(i,j,k,2) - 1.d0
          randum(i,j,k,3) = 2.d0 * randum(i,j,k,3) - 1.d0
          randum(i,j,k,4) = 2.d0 * randum(i,j,k,4) - 1.d0
    enddo;enddo;enddo

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          eta = 5.d0 * ys(j) / blt0
          u   = min(u0, u0 * (0.0015d0 * eta**4 - 0.0181d0 * eta**3 + 0.029d0 * eta**2 + 0.3192 * eta + 0.0003d0))
          v   = 0.d0
          Taw = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
          Tw  = Taw
          T   = Tw + (Taw - Tw) * u / u0 - rf * u**2 / (2.d0 * (gamma * R / (gamma - 1.d0)))
          !call calc_Blasius(eta,d,u,v)
          if (10 < j .and. ys(j) <= blt) then
            ir = i / 10 + 1
            jr = j / 2  + 1
            kr = k / 2  + 1
            ustd = 0.2d0 * u0 * randum(ir,jr,kr,1)
            vstd = 0.1d0 * u0 * randum(ir,jr,kr,2)
            wstd = 0.1d0 * u0 * randum(ir,jr,kr,3)
            Tstd = T0 * (gamma - 1.d0) * M0**2 * 0.2d0 * randum(ir,jr,kr,4)
          else
            ustd = 0.d0
            vstd = 0.d0
            wstd = 0.d0
            Tstd = 0.d0
          endif
          u   = u + ustd
          v   = v + vstd
          w   = wstd
          T   = T + Tstd
          rho = p0 / (R * T)
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = Q(i,j,k,1) * u
          Q(i,j,k,3) = Q(i,j,k,1) * v
          Q(i,j,k,4) = Q(i,j,k,1) * w
          Q(i,j,k,5) = p0 / (gamma - 1.d0) + 0.5d0 * (Q(i,j,k,2)**2 + Q(i,j,k,3)**2 + Q(i,j,k,4)**2) / Q(i,j,k,1)
    enddo;enddo;enddo

    deallocate(randum)

    ! bottom
    Q(:,1,:,1) = Q(:,2,:,1)
    Q(:,1,:,2) = 0.d0
    Q(:,1,:,3) = 0.d0
    Q(:,1,:,4) = 0.d0
    p_wall = (gamma - 1.d0) * (Q(2,2,2,5) - 0.5d0 * (Q(2,2,2,2)**2 + Q(2,2,2,3)**2 + Q(2,2,2,4)**2) / Q(2,2,2,1))
    Q(:,1,:,5) = p_wall / (gamma - 1.d0)
  end subroutine set_init
  
  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(ny)
    real(8), intent(inout), device :: QJ(nx,ny,nz,5) ! Q / Jacobian
    real(8), intent(in), device    :: Qre(ny*(nz-6)*5)
    integer i, j, k, l, No
    real(8) :: Cp = gamma * R / (gamma - 1.d0), rf = 0.89d0
    real(8) :: p_wall
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb
    real(8) :: rho0, c0, v0 = 0.d0, Taw, T

    if (kind(id_rescale) == 4 .and. myrank == 0) then
      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 1, nz-6
          do j = 1, ny
            QJ(1,j,k+3,l)  = Qre(ny*(nz-6)*(l-1)+ny*(k-1)+j)
      enddo;enddo;enddo
    elseif (kind(id_rescale) == 4 .and. myrank == 14) then
      !$cuf kernel do(3)<<<*,*>>>
      do l = 1, 5
        do k = 4, nz-3
          do j = 2, ny-1
            QJ(nx,j,k,l) = QJ(nx-1,j,k,l)
      enddo;enddo;enddo
    endif

    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        ! top
        ! Riemann invariants
        pin   = (gamma - 1.d0) * (QJ(i,ny-1,k,5) - 0.5d0 * (QJ(i,ny-1,k,2)**2 + QJ(i,ny-1,k,3)**2 + QJ(i,ny-1,k,4)**2) &
                / QJ(i,ny-1,k,1)) * Jacobian(ny-1)
        rhoin = QJ(i,ny-1,k,1) * Jacobian(ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ(i,ny-1,k,3) / QJ(i,ny-1,k,1)
        ! temperature and density at top
        Taw   = T0 * (1.d0 + rf * 0.5d0 * (gamma - 1.d0) * M0**2)
        T     = Taw - rf * u0**2 / (2.d0 * (gamma * R / (gamma - 1.d0)))
        rho0  = p0 / (R * T)
        c0    = sqrt(gamma * p0 / rho0)

        Rp    = vin + 2.d0 * cin / (gamma - 1.d0)
        Rm    = v0  - 2.d0 * c0  / (gamma - 1.d0)
        vb    = 0.5d0 * (Rp + Rm)
        cb    = 0.25d0 * (gamma - 1.d0) * (Rp - Rm)
        rhob  = (cb / c0)**(2.d0 / (gamma - 1.d0)) * rho0
        pb    = (rhob * cb**2) / gamma

        QJ(i,ny,k,1) = rhob / Jacobian(ny)
        QJ(i,ny,k,2) = rhob * u0 / Jacobian(ny)
        QJ(i,ny,k,3) = rhob * vb / Jacobian(ny)
        QJ(i,ny,k,4) = 0.d0
        QJ(i,ny,k,5) = (pb / (gamma - 1.d0) + 0.5d0 * rhob * (u0**2 + vb**2)) / Jacobian(ny)

        ! NoSlip
        QJ(i,1,k,1) = QJ(i,2,k,1)
        QJ(i,1,k,2) = 0.d0
        QJ(i,1,k,3) = 0.d0
        QJ(i,1,k,4) = 0.d0
        p_wall = (gamma - 1.d0) * (QJ(i,2,k,5) - 0.5d0 * (QJ(i,2,k,2)**2 + QJ(i,2,k,3)**2 + QJ(i,2,k,4)**2) / QJ(i,2,k,1))
        QJ(i,1,k,5) = p_wall / (gamma - 1.d0)
    enddo;enddo

    if (6 == myrank) then
      No = int(0.2d0 * nx)
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz
        do i = No, nx
          QJ(i,ny,k,1) = rho2 / Jacobian(ny)
          QJ(i,ny,k,2) = rho2 * ux / Jacobian(ny)
          QJ(i,ny,k,3) = rho2 * uy / Jacobian(ny)
          QJ(i,ny,k,4) = 0.d0
          QJ(i,ny,k,5) = (p2 / (gamma - 1.d0) + 0.5d0 * rho2 * (ux**2 + uy**2)) / Jacobian(ny)
      enddo;enddo
    elseif (8 <= myrank) then
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz
        do i = 1, nx
          pin   = (gamma - 1.d0) * (QJ(i,ny-1,k,5) - 0.5d0 * (QJ(i,ny-1,k,2)**2 + QJ(i,ny-1,k,3)**2 + QJ(i,ny-1,k,4)**2) &
                  / QJ(i,ny-1,k,1)) * Jacobian(ny-1)
          rhoin = QJ(i,ny-1,k,1) * Jacobian(ny-1)
          cin   = sqrt(gamma * pin / rhoin)
          vin   = QJ(i,ny-1,k,3) / QJ(i,ny-1,k,1)
          c0    = sqrt(gamma * p2 / rho2)
          Rp    = vin + 2.d0 * cin / (gamma - 1.d0)
          Rm    = uy  - 2.d0 * c0  / (gamma - 1.d0)
          vb    = 0.5d0 * (Rp + Rm)
          cb    = 0.25d0 * (gamma - 1.d0) * (Rp - Rm)
          rhob  = cin * rhoin / cb
          pb    = (rhob * cb**2) / gamma
          ub    = sqrt(2.d0 * gamma * (p2 / rho2 - pb / rhob) / (gamma - 1.d0) + ux**2 + uy**2 - vb**2)
          QJ(i,ny,k,1) = rhob / Jacobian(ny)
          QJ(i,ny,k,2) = rhob * ub / Jacobian(ny)
          QJ(i,ny,k,3) = rhob * vb / Jacobian(ny)
          QJ(i,ny,k,4) = 0.d0
          QJ(i,ny,k,5) = (pb / (gamma - 1.d0) + 0.5d0 * rhob * (ub**2 + vb**2)) / Jacobian(ny)
      enddo;enddo
    endif

    ! cyclic
    !$cuf kernel do(3)<<<*,*>>>
    do l = 1, 5
      do j = 1, ny
        do i = 1, nx
          QJ(i,j,1,l) = QJ(i,j,nz-5,l)
          QJ(i,j,2,l) = QJ(i,j,nz-4,l)
          QJ(i,j,3,l) = QJ(i,j,nz-3,l)
          QJ(i,j,nz-2,l) = QJ(i,j,4,l)
          QJ(i,j,nz-1,l) = QJ(i,j,5,l)
          QJ(i,j,nz,l)   = QJ(i,j,6,l)
    enddo;enddo;enddo
  end subroutine set_bc

  subroutine set_bc_mut(nx,ny,nz,mut,qc2)
    integer, intent(in), value      :: nx, ny, nz
    real(8), intent(inout), device  :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        ! inlet
        mut(1,j,k)  = mut(nre2,j,k)
        qc2(1,j,k)  = qc2(nre2,j,k)
        ! outlet
        mut(nx,j,k) = mut(nx-1,j,k)
        qc2(nx,j,k) = qc2(nx-1,j,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        ! wall
        mut(i,1,k) = 0.d0
        qc2(i,1,k) = 0.d0
        ! top
        mut(i,ny,k) = mut(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,ny-1,k)
    enddo;enddo

    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        ! span
        mut(i,j,1)    = mut(i,j,nz-5)
        mut(i,j,2)    = mut(i,j,nz-4)
        mut(i,j,3)    = mut(i,j,nz-3)
        mut(i,j,nz-2) = mut(i,j,4)
        mut(i,j,nz-1) = mut(i,j,5)
        mut(i,j,nz)   = mut(i,j,6)
        qc2(i,j,1)    = qc2(i,j,nz-5)
        qc2(i,j,2)    = qc2(i,j,nz-4)
        qc2(i,j,3)    = qc2(i,j,nz-3)
        qc2(i,j,nz-2) = qc2(i,j,4)
        qc2(i,j,nz-1) = qc2(i,j,5)
        qc2(i,j,nz)   = qc2(i,j,6)
    enddo;enddo
  end subroutine set_bc_mut
end module set

