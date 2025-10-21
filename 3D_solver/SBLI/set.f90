module set
  use cudafor
  use mpi
  use mod_globals, only : id_accuracy, id_rescale, ny1, nre2, gamma, R, rf, Cp, Pr, u0, p0, T0, M0, blt, rho2, p2, ux, uy
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  use set_bc_common
  use set_bc_tbl_sbli
  use set_init_common
  use calc_para
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, Lx1, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz, Lx1
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1
    dx1 = Lx / dble(nx-1)
    dy1 = dx1
    dz1 = Lz / dble(nz-1)
    if (myrank == 0) then
      x(1) = 0.d0
    else
      x(1) = Lx1
    endif
    do i = 1, nx-1
      dx(i) = dx1
      x(i+1) = x(i) + dx(i)
    enddo
    y(1) = 0.d0
    do j = 1, ny-1
      if (y(j) <= 3.d0 * blt) then
        dy(j) = min(1.d0, max(0.07d0, dble(j)/dble(128))) * dy1
      elseif (3.d0 * blt <= y(j) .and. y(j) <= 8.d0 * blt) then
        dy(j) = 1.5d0 * dy1
      else
        dy(j) = 1.75d0 * dy1
      endif
      y(j+1) = y(j) + dy(j)
    enddo
    print *, y(ny)
    z(1) = 0.d0
    do k = 1, nz-1
      dz(k) = dz1
      z(k+1) = z(k) + dz(k)
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, ny, nz, xs, ys, zs, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(out) :: Q(5,nx,ny,nz)
    call set_init_tbl(nx, ny, nz, xs, ys, zs, 0.05d0, 0.75d0*blt, blt, rf, u0, p0, T0, M0, Q)
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
    integer i, j, k, l, No, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall, rf = 0.89d0
    ! Riemann invariants
    real(8) :: rhoin, pin, cin, vin, Rp, Rm, rhob, ub, vb, cb, pb
    real(8) :: rho0, c0, v0 = 0.d0, Taw, T
    ! parallel
    real(8), device :: Q1d(3*(ny1-2)*(nz-6)*5)
    !real(8) Q_cpu(3*(ny1-2)*(nz-6)*5)
    ! cache
    real(8) Jacobian_tmp
    ! temperature and density at top
    Taw  = T0 * (1.d0 + rf * 0.5d0 * gamma_1 * M0**2)
    T    = Taw - rf * u0**2 / (2.d0 * (gamma * R * over_gamma_1))
    rho0 = p0 / (R * T)
    c0   = sqrt(gamma * p0 / rho0)
    if (myrank == 0) then
      ! inlet & outlet !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (kind(id_rescale) == 4) then
        !$cuf kernel do(2)<<<*,*>>>
        do k = 1, nz-6
          do j = 2, ny-1
            do l = 1, 5
              ! inlet
              QJ(l,1,j,k+3)  = Qre(ny*5*(k-1)+5*(j-1)+l)
              ! outlet
              QJ(l,nx,j,k+3) = QJ(l,nx-1,j,k+3)
        enddo;enddo;enddo
      else
        !$cuf kernel do(2)<<<*,*>>>
        do k = 4, nz-3
          do j = 2, ny-1
            do l = 1, 5
              ! inlet
              QJ(l,1,j,k) = QJ(l,nx-5,j,k)
              QJ(l,2,j,k) = QJ(l,nx-4,j,k)
              QJ(l,3,j,k) = QJ(l,nx-3,j,k)
              ! outlet
              QJ(l,nx-2,j,k) = QJ(l,4,j,k)
              QJ(l,nx-1,j,k) = QJ(l,5,j,k)
              QJ(l,nx,j,k)   = QJ(l,6,j,k)
        enddo;enddo;enddo
      endif
      call flatten_rescale(nx, ny1, nz, nre2, 3, QJ, Q1d)
      !Q_cpu = Q1d ! This is safe but very slow
      call MPI_ISEND(Q1d, 5*3*(ny1-2)*(nz-6), MPI_REAL8, myrank+2, 0, MPI_COMM_WORLD, ireq, ierr)
    else
      call MPI_IRECV(Q1d, 5*3*(ny1-2)*(nz-6), MPI_REAL8, myrank-2, 0, MPI_COMM_WORLD, ireq, ierr)
      call MPI_WAIT(ireq, istat, ierr)
      !Q1d = Q_cpu ! This is safe but very slow
      ! inlet boundary layer
      call reconstruct_sbli_inlet(nx, ny1, ny, nz, 3, Q1d, QJ)
      !$cuf kernel do(2)<<<*,*>>>
      do k = 4, nz-3
        do j = ny1-1, ny-1
          ! inlet free stream flow
          Jacobian_tmp = 1.d0 / Jacobian(1,j)
          QJ(1,1,j,k) = rho0 * Jacobian_tmp
          QJ(2,1,j,k) = rho0 * u0 * Jacobian_tmp
          QJ(3,1,j,k) = 0.d0
          QJ(4,1,j,k) = 0.d0
          QJ(5,1,j,k) = (p0 * over_gamma_1 + 0.5d0 * rho0 * u0**2) * Jacobian_tmp
      enddo;enddo
      !$cuf kernel do(2)<<<*,*>>>
      do k = 4, nz-3
        do j = 2, ny-1
          do l = 1, 5
            ! outlet
            QJ(l,nx,j,k) = QJ(l,nx-1,j,k)
      enddo;enddo;enddo
    endif

    !call set_bc_Riemann_tbl_top_down(nx, ny, nz, 3, 1, nx, Jacobian, QJ)
    call set_bc_Neumann_tbl_top_down(nx, ny, nz, 3, 1, nx, Jacobian, QJ)
    
    if (myrank == 2) then
      No = int(dble(nx)*0.1d0)!int(dble(nx) * 0.33d0 / 35.d0)
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz
        do i = No, nx
          Jacobian_tmp = 1.d0 / Jacobian(i,ny)
          vin   = QJ(3,i,ny-1,k) / QJ(1,i,ny-1,k)
          if (0.5d0 * uy > vin) then 
            QJ(1,i,ny,k) = rho2 * Jacobian_tmp
            QJ(2,i,ny,k) = rho2 * ux * Jacobian_tmp
            QJ(3,i,ny,k) = rho2 * uy * Jacobian_tmp
            QJ(4,i,ny,k) = 0.d0
            QJ(5,i,ny,k) = (p2 * over_gamma_1 + 0.5d0 * rho2 * (ux**2 + uy**2)) * Jacobian_tmp
          else
            pin   = gamma_1 * (QJ(5,i,ny-1,k) - 0.5d0 * (QJ(2,i,ny-1,k)**2 + QJ(3,i,ny-1,k)**2 + QJ(4,i,ny-1,k)**2) &
                    / QJ(1,i,ny-1,k)) * Jacobian(i,ny-1)
            rhoin = QJ(1,i,ny-1,k) * Jacobian(i,ny-1)
            cin   = sqrt(gamma * pin / rhoin)
            c0    = sqrt(gamma * p2 / rho2)
            Rp    = vin + 2.d0 * cin * over_gamma_1
            Rm    = v0  - 2.d0 * c0  * over_gamma_1
            vb    = 0.5d0 * (Rp + Rm)
            cb    = 0.25d0 * gamma_1 * (Rp - Rm)
            rhob  = cin * rhoin / cb
            pb    = (rhob * cb**2) / gamma
            ub    = sqrt(2.d0 * gamma * (p2 / rho2 - pb / rhob) * over_gamma_1 + ux**2 + uy**2 - vb**2)
            QJ(1,i,ny,k) = rhob * Jacobian_tmp
            QJ(2,i,ny,k) = rhob * ub * Jacobian_tmp
            QJ(3,i,ny,k) = rhob * vb * Jacobian_tmp
            QJ(4,i,ny,k) = 0.d0
            QJ(5,i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp
          endif
      enddo;enddo
    endif
    call set_bc_cyclic_z(nx, ny, nz, QJ)
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

