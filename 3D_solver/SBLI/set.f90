module set
  use cudafor
  use mpi
  use mod_globals, only : ny1, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, beta, &
                          ny2, rho2, p2, ux, uy, rf, Taw, rho3, p3, ux3, uy3
  use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1, id_rescale
  use set_bc_common
  use set_bc_tbl_sbli
  use set_init_common
  use calc_para
  implicit none
  integer No
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, Lx1, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz, Lx1
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k, nx1, ny_b
    real(8) dx1, dy1, dz1, ximp, Lx_s
    dx1  = 20.d0 * blt / dble(512)
    dy1  = dx1
    dz1  = Lz / dble(nz-1)
    ximp = 0.9d0 * Lx1 + 30.d0 * blt

    if (myrank == 0) then
      x(1) = 0.d0
      do i = 1, nx-1
        dx(i) = dx1
        x(i+1) = x(i) + dx(i)
      enddo
    else
      x(1) = Lx1
      nx1  = int(0.9d0 * dble(nx-1))
      ! computational region
      do i = 1, nx1
        dx(i)  = dx1
        x(i+1) = x(i) + dx(i)
      enddo
      ! buffer region
      do i = nx1 + 1, nx-1
        dx(i)  = dx1 * (1.d0 + 3.d0 * dble(i-nx1) / dble(nx-nx1))
        x(i+1) = x(i) + dx(i)
      enddo
    endif
    x(:) = x(:) - ximp

    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
        ny_b  = j
      else
        dy(j) = dy1 * (1.d0 + 0.75d0 * dble(j-ny_b) / dble(ny2-ny_b))
      endif
      y(j+1) = y(j) + dy(j)
    enddo

    if (myrank == 2) then
      Lx_s = y(ny) / dble(beta) / blt
      do i = 1, nx
        if (x(i) / blt + Lx_s >= 0.d0) then
          No = i
          exit
        endif
      enddo
    endif

    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
    z(:) = z(:) - 0.5d0 * Lz
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    call set_init_tbl(nx, ny, nz, xs, ys, zs, 0.1d0, 0.75d0*blt, blt, u0, p0, T0, M0, Q)
  end subroutine set_init


  subroutine set_bc_Gaussian(nx, ny, nz, nxg, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: nx, ny, nz, nxg
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), device :: tmp(4)
    integer i, k
    ! y direction one-sided
    !$cuf kernel do(1)<<<*,*>>>
    do k = 1, nz
      do i = nxg, nx
        tmp(:) = QJ_1(i,ny-3:ny,k) * Jacobian(i,ny-3:ny)
        QJ_1(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_2(i,ny-3:ny,k) * Jacobian(i,ny-3:ny)
        QJ_2(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_3(i,ny-3:ny,k) * Jacobian(i,ny-3:ny)
        QJ_3(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_4(i,ny-3:ny,k) * Jacobian(i,ny-3:ny)
        QJ_4(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_5(i,ny-3:ny,k) * Jacobian(i,ny-3:ny)
        QJ_5(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
    enddo;enddo
    ! x direction one-sided
    !$cuf kernel do(1)<<<*,*>>>
    do k = 1, nz
      do i = nxg, nx
        tmp(:) = QJ_1(i-3:i,ny,k) * Jacobian(i-3:i,ny)
        QJ_1(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_2(i-3:i,ny,k) * Jacobian(i-3:i,ny)
        QJ_2(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_3(i-3:i,ny,k) * Jacobian(i-3:i,ny)
        QJ_3(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_4(i-3:i,ny,k) * Jacobian(i-3:i,ny)
        QJ_4(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
        tmp(:) = QJ_5(i-3:i,ny,k) * Jacobian(i-3:i,ny)
        QJ_5(i,ny,k) = (0.05d0 * tmp(1) + 0.15d0 * tmp(2) + 0.3d0 * tmp(3) + 0.5d0 * tmp(4)) / Jacobian(i,ny)
    enddo;enddo
  end subroutine set_bc_Gaussian


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, optional :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    integer i, j, k, offset, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall, pre, rho, rhou, rhov, rhow, p, e
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb, v0 = 0.d0
    ! parallel
    real(8), device :: Q1d(3*(ny1-2)*(nz-6)*5)
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    real(8), parameter :: T    = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0 = p0 / (R * T)
    real(8), parameter :: c0   = sqrt(gamma * p0 / rho0)
    real(8), parameter :: c3   = sqrt(gamma * p3 / rho3)
    if (myrank == 0) then
      if (kind(id_rescale) == 4) then
        !$cuf kernel do(2)<<<*,*>>>
        do k = 1, nz-6
          do j = 2, ny-1
            offset = ny*(k-1) + j
            QJ_1(1,j,k+3)  = Qre_1(offset)
            QJ_2(1,j,k+3)  = Qre_2(offset)
            QJ_3(1,j,k+3)  = Qre_3(offset)
            QJ_4(1,j,k+3)  = Qre_4(offset)
            QJ_5(1,j,k+3)  = Qre_5(offset)
            QJ_1(nx,j,k+3) = QJ_1(nx-1,j,k+3)
            QJ_2(nx,j,k+3) = QJ_2(nx-1,j,k+3)
            QJ_3(nx,j,k+3) = QJ_3(nx-1,j,k+3)
            QJ_4(nx,j,k+3) = QJ_4(nx-1,j,k+3)
            QJ_5(nx,j,k+3) = QJ_5(nx-1,j,k+3)
        enddo;enddo
      else
        !$cuf kernel do(2)<<<*,*>>>
        do k = 4, nz-3
          do j = 2, ny-1
            ! inlet
            QJ_1(1,j,k) = QJ_1(nx-5,j,k)
            QJ_2(1,j,k) = QJ_2(nx-5,j,k)
            QJ_3(1,j,k) = QJ_3(nx-5,j,k)
            QJ_4(1,j,k) = QJ_4(nx-5,j,k)
            QJ_5(1,j,k) = QJ_5(nx-5,j,k)
            QJ_1(2,j,k) = QJ_1(nx-4,j,k)
            QJ_2(2,j,k) = QJ_2(nx-4,j,k)
            QJ_3(2,j,k) = QJ_3(nx-4,j,k)
            QJ_4(2,j,k) = QJ_4(nx-4,j,k)
            QJ_5(2,j,k) = QJ_5(nx-4,j,k)
            QJ_1(3,j,k) = QJ_1(nx-3,j,k)
            QJ_2(3,j,k) = QJ_2(nx-3,j,k)
            QJ_3(3,j,k) = QJ_3(nx-3,j,k)
            QJ_4(3,j,k) = QJ_4(nx-3,j,k)
            QJ_5(3,j,k) = QJ_5(nx-3,j,k)
            ! outlet
            QJ_1(nx-2,j,k) = QJ_1(4,j,k)
            QJ_2(nx-2,j,k) = QJ_2(4,j,k)
            QJ_3(nx-2,j,k) = QJ_3(4,j,k)
            QJ_4(nx-2,j,k) = QJ_4(4,j,k)
            QJ_5(nx-2,j,k) = QJ_5(4,j,k)
            QJ_1(nx-1,j,k) = QJ_1(5,j,k)
            QJ_2(nx-1,j,k) = QJ_2(5,j,k)
            QJ_3(nx-1,j,k) = QJ_3(5,j,k)
            QJ_4(nx-1,j,k) = QJ_4(5,j,k)
            QJ_5(nx-1,j,k) = QJ_5(5,j,k)
            QJ_1(nx,j,k)   = QJ_1(6,j,k)
            QJ_2(nx,j,k)   = QJ_2(6,j,k)
            QJ_3(nx,j,k)   = QJ_3(6,j,k)
            QJ_4(nx,j,k)   = QJ_4(6,j,k)
            QJ_5(nx,j,k)   = QJ_5(6,j,k)
        enddo;enddo
      endif
      call flatten_rescale(nx, ny1, nz, nre2, 3, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Q1d)
      !Q_cpu = Q1d ! This is safe but very slow
      call MPI_ISEND(Q1d, 5*3*(ny1-2)*(nz-6), MPI_REAL8, myrank+2, 0, MPI_COMM_WORLD, ireq, ierr)
    else
      call MPI_IRECV(Q1d, 5*3*(ny1-2)*(nz-6), MPI_REAL8, myrank-2, 0, MPI_COMM_WORLD, ireq, ierr)
      call MPI_WAIT(ireq, istat, ierr)
      !Q1d = Q_cpu ! This is safe but very slow
      ! inlet boundary layer
      call reconstruct_sbli_inlet(nx, ny1, ny, nz, 3, Q1d, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
      !$cuf kernel do(2)<<<*,*>>>
      do k = 4, nz-3
        do j = ny1-1, ny-1
          ! inlet free stream flow
          Jacobian_tmp = 1.d0 / Jacobian(1,j)
          QJ_1(1,j,k) = rho0 * Jacobian_tmp     !rho  * Jacobian_tmp
          QJ_2(1,j,k) = rho0 * u0 * Jacobian_tmp!rhou * Jacobian_tmp
          QJ_3(1,j,k) = 0.d0                    !rhov * Jacobian_tmp
          QJ_4(1,j,k) = 0.d0                    !rhow * Jacobian_tmp
          QJ_5(1,j,k) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp!e * Jacobian_tmp
      enddo;enddo
      !$cuf kernel do(2)<<<*,*>>>
      do k = 4, nz-3
        do j = 2, ny-1
          ! outlet
          QJ_1(nx,j,k) = QJ_1(nx-1,j,k)
          QJ_2(nx,j,k) = QJ_2(nx-1,j,k)
          QJ_3(nx,j,k) = QJ_3(nx-1,j,k)
          QJ_4(nx,j,k) = QJ_4(nx-1,j,k)
          QJ_5(nx,j,k) = QJ_5(nx-1,j,k)
      enddo;enddo
    endif

    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        ! top
        ! Riemann invariants
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        pin   = gamma_1 * (QJ_5(i,ny-1,k) - 0.5d0 * (QJ_2(i,ny-1,k)**2 + QJ_3(i,ny-1,k)**2 + QJ_4(i,ny-1,k)**2) &
                / QJ_1(i,ny-1,k)) * Jacobian(i,ny-1)
        rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
        Rp    = vin + 2.d0 * cin * over_gamma_1
        Rm    = v0  - 2.d0 * c0  * over_gamma_1
        vb    = 0.5d0 * (Rp + Rm)
        cb    = 0.25d0 * gamma_1 * (Rp - Rm)
        rhob  = (cb / c0)**(2.d0 * over_gamma_1) * rho0
        pb    = (rhob * cb**2) / gamma
        QJ_1(i,ny,k) = rhob * Jacobian_tmp
        QJ_2(i,ny,k) = rhob * u0 * Jacobian_tmp
        QJ_3(i,ny,k) = rhob * vb * Jacobian_tmp
        QJ_4(i,ny,k) = 0.d0
        QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
        ! NoSlip
        QJ_1(i,1,k) = QJ_1(i,2,k)
        QJ_2(i,1,k) = 0.d0
        QJ_3(i,1,k) = 0.d0
        QJ_4(i,1,k) = 0.d0
        p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
        QJ_5(i,1,k) = p_wall * over_gamma_1
    enddo;enddo
    !call set_bc_Riemann_tbl_top_down(nx, ny, nz, 3, 1, nx, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    !call set_bc_Neumann_tbl_top_down(nx, ny, nz, 3, 1, nx, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)

    if (myrank == 2) then
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz
        do i = No, nx
          Jacobian_tmp = 1.d0 / Jacobian(i,ny)
          vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
          if (0.5d0 * uy > vin) then
            QJ_1(i,ny,k) = rho2 * Jacobian_tmp
            QJ_2(i,ny,k) = rho2 * ux * Jacobian_tmp
            QJ_3(i,ny,k) = rho2 * uy * Jacobian_tmp
            QJ_4(i,ny,k) = 0.d0
            QJ_5(i,ny,k) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
          else
            pin   = gamma_1 * (QJ_5(i,ny-1,k) - 0.5d0 * (QJ_2(i,ny-1,k)**2 + QJ_3(i,ny-1,k)**2 + QJ_4(i,ny-1,k)**2) &
                    / QJ_1(i,ny-1,k)) * Jacobian(i,ny-1)
            rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
            cin   = sqrt(gamma * pin / rhoin)
            Rp    = vin + 2.d0 * cin * over_gamma_1
            Rm    = uy3 - 2.d0 * c3  * over_gamma_1
            vb    = 0.5d0 * (Rp + Rm)
            cb    = 0.25d0 * gamma_1 * (Rp - Rm)
            rhob  = cin * rhoin / cb
            pb    = (rhob * cb**2) * over_gamma
            ub    = sqrt(2.d0 * gamma * (p3 / rho3 - pb / rhob) * over_gamma_1 + ux3**2 + uy3**2 - vb**2)
            QJ_1(i,ny,k) = rhob * Jacobian_tmp
            QJ_2(i,ny,k) = rhob * ub * Jacobian_tmp
            QJ_3(i,ny,k) = rhob * vb * Jacobian_tmp
            QJ_4(i,ny,k) = 0.d0
            QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
          endif
      enddo;enddo
    endif

    call set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    call set_bc_Gaussian(nx, ny, nz, int(0.5d0 * nx), Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)

    if (myrank == 0) then
      call MPI_WAIT(ireq, istat, ierr)
    endif
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

