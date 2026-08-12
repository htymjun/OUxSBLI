module set
  use cudafor
  use mpi
  use mod_globals, only : ny1, nre2, gamma, R, Cp, Pr, u0, p0, T0, M0, blt, rf, Taw
  use mod_constant, only : Cp, gamma_1, over_gamma, over_gamma_1, id_rescale
  use set_bc_common
  use set_init_common
  implicit none
contains
  subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: Lx, Ly, Lz
    real(8), intent(out) :: x(nx), y(ny), z(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    integer i, j, k
    real(8) dx1, dy1, dz1, tanh_s, yi
    real(8), parameter :: s = 2.4d0
    dx1 = Lx / dble(nx-1)
    dz1 = Lz / dble(nz-1)

    x(1) = 0.d0
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


  subroutine set_bc(myrank, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, optional :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    integer i, j, k, offset 
    real(8) p_wall, Jacobian_tmp
    real(8) rhoin, uin, vin, win, pin, rhob, ub, vb, wb, pb, sb
    real(8), parameter :: Tinf = Taw - rf * u0**2 / (2.0d0 * Cp) ! Crocco-Busemann
    real(8), parameter :: s0   = p0 / (p0 / (R * Tinf))**gamma

    if (kind(id_rescale) == 4 .and. present(Qre_1)) then
      !$cuf kernel do(2)<<<*,*>>>
      do k = 1, nz-6
        do j = 2, ny-1
          offset = ny*(k-1) + j
          ! inlet
          QJ_1(1,j,k+3)  = Qre_1(offset)
          QJ_2(1,j,k+3)  = Qre_2(offset)
          QJ_3(1,j,k+3)  = Qre_3(offset)
          QJ_4(1,j,k+3)  = Qre_4(offset)
          QJ_5(1,j,k+3)  = Qre_5(offset)
          ! outlet
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

    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 1, nx
        ! top
        ! Riemann invariants
        rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
        uin   = QJ_2(i,ny-1,k) / QJ_1(i,ny-1,k)
        vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
        win   = QJ_4(i,ny-1,k) / QJ_1(i,ny-1,k)
        pin   = gamma_1 * (QJ_5(i,ny-1,k) * Jacobian(i,ny-1) &
                - 0.5d0 * rhoin * (uin**2 + vin**2 + win**2))
        pb = p0
        vb = vin
        if (vb >= 0.d0) then
          sb = pin / rhoin**gamma
          ub = uin
          wb = win
        else
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
        ! NoSlip
        QJ_1(i,1,k) = QJ_1(i,2,k)
        QJ_2(i,1,k) = 0.d0
        QJ_3(i,1,k) = 0.d0
        QJ_4(i,1,k) = 0.d0
        p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
        QJ_5(i,1,k) = p_wall * over_gamma_1
    enddo;enddo

    call set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
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
