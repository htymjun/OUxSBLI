module set
  use cudafor
  use mpi
  use mod_globals, only : gamma, R, Pr, rho0, u0, p0, T0, M0, blt, rho2, p2, ux, uy, rf, Taw
  use mod_constant, only : id_rescale, Cp, gamma_1, over_gamma, over_gamma_1, mu0_T0_S_over_T0_2_3
  use set_bc_common
  use set_init_common
  use mod_shock !!!!!!!!!!!!!!!!!!!!!!!oblique shock
  implicit none
  real(8) ptbl, rho0_ptbl
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

    s = 1.6d0
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


  subroutine interpolate_Q(nx, nyi, ny, y, Qi, Qp)
    integer, intent(in)  :: nx, nyi, ny
    real(8), intent(in)  :: y(ny)
    real(8), intent(in)  :: Qi(5,nyi)
    real(8), intent(out) :: Qp(4,ny)
    integer :: i, j
    real(8) :: y1, y2, alpha
    do j = 1, ny
      if (y(j) <= Qi(1,1)) then
        Qp(:,j) = Qi(2:5,1)
      elseif (y(j) >= Qi(1,nyi)) then
        Qp(:,j) = Qi(2:5,nyi)
      else
        do i = 1, nyi-1
          y1 = Qi(1,i)
          y2 = Qi(1,i+1)
          if (y1 <= y(j) .and. y(j) <= y2) then
            alpha = (y(j) - y1) / (y2 - y1)
            Qp(:,j) = (1.0d0 - alpha) * Qi(2:5,i) + alpha * Qi(2:5,i+1)
            exit
          endif
        enddo
      endif
    enddo
  end subroutine interpolate_Q


  subroutine set_init(myrank, nx, ny, xs, ys, Q)
    integer, intent(in)  :: myrank, nx, ny
    real(8), intent(in)  :: xs(nx), ys(ny)
    real(8), intent(out) :: Q(nx,4,ny)
    character(len=40) filename
    integer i, j, nyi, filesize, ios
    real(8) mu0, nu0, disp_thic, Re_disp, deta !!!!displacement thickness
    real(8) rho, u, v, p
    real(8), allocatable :: Qi(:,:), Qp(:,:)
    write(filename, "(a)") "./Qin.dat"
    open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file Qin.dat"
    endif
    inquire(10, size=filesize)
    nyi = filesize / (8 * 5)
    allocate(Qi(5,nyi))
    read(10) Qi
    close(10)
    allocate(Qp(4,ny))
    call interpolate_Q(nx, nyi, ny, ys, Qi, Qp)
    ! ptbl = 0.d0
    ! do j = 1, nyi
    !   ptbl = ptbl + Qp(4,j)
    ! enddo
    ! ptbl = ptbl / dble(nyi)
    !ptbl = Qp(4,ny) !!!!!!!!!!!!!!!!!!!!!!!
    do j = 1, ny
      rho = Qp(1,j)
      u   = Qp(2,j)
      v   = Qp(3,j)
      Q(:,1,j) = rho
      Q(:,2,j) = rho * u
      Q(:,3,j) = rho * v
      !Q(:,4,j) = ptbl * over_gamma_1 + 0.5d0 * rho * (u**2 + v**2)
      Q(:,4,j) = Qp(4,j) * over_gamma_1 + 0.5d0 * rho * (u**2 + v**2)
    enddo

    !call calc_p_rho_init(ptbl)
    call Qin_init(Qp(:,ny))
    
    ptbl = 0.d0
    do j = 1, nyi
      ptbl = ptbl + Qp(4,j)
    enddo
    ptbl = ptbl / dble(nyi)
    rho0_ptbl = ptbl / (R * T0)    

    mu0 = mu0_T0_S_over_T0_2_3 / (T0 + 110.4d0) * T0**1.5d0 !!!!!!!!!!!!!!!!!!!!!displacement thickness
    nu0 = mu0 * R * T0 / ptbl
    disp_thic = 0.d0
    do j = 2, ny
        deta = -ys(j-1) + ys(j)
        disp_thic = disp_thic + 0.5d0 * ((1 - (Qp(1,j-1) * Qp(2,j-1)) / (rho0_ptbl * u0)) + (1 - (Qp(1,j) * Qp(2,j)) / (rho0_ptbl * u0))) * deta
    enddo
    Re_disp = u0 * disp_thic / nu0
    print *, "Re_disp = ", Re_disp

    deallocate(Qi, Qp)

    do i = int(0.65d0 * nx), nx
      Q(i,1,ny) = rho2_init
      Q(i,2,ny) = rho2_init * ux_init
      Q(i,3,ny) = rho2_init * uy_init
      Q(i,4,ny) = p2_init * over_gamma_1 + 0.5d0 * rho2_init * (ux_init**2 + uy_init**2)
    enddo
  end subroutine set_init


  subroutine set_bc(myrank, nx, ny, Jacobian, QJ)
    integer, intent(in), value     :: myrank, nx, ny
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ(nx,4,ny) ! Q / Jacobian
    real(8) Jacobian_tmp
    integer i, j, l, ireq, ierr, istat(MPI_STATUS_SIZE)
    real(8) :: p_wall
    real(8) :: rho_ext, u_ext, v_ext, p_ext, c_ext
    real(8) :: rhoin, uin, vin, pin, cin
    real(8) :: Rp, Rm, vb, cb, sb, ub, rhob, pb
    integer :: i_switch

    i_switch = int(0.65d0 * nx)

    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      do l = 1, 4
        ! outlet
        QJ(nx,l,j) = QJ(nx-1,l,j)
    enddo;enddo

    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      ! Jacobian_tmp = 1.d0 / Jacobian(1,ny-1)
      ! ! Neumann
      ! if (i < int(0.7d0 * nx)) then !0.1nx
      !   QJ(i,1,ny) = QJ(i,1,ny-1) 
      !   QJ(i,2,ny) = QJ(i,2,ny-1)
      !   QJ(i,3,ny) = QJ(i,3,ny-1)
      !   QJ(i,4,ny) = QJ(i,4,ny-1)
      ! else
      !   QJ(i,1,ny) = rho2_init * Jacobian_tmp
      !   QJ(i,2,ny) = rho2_init * ux * Jacobian_tmp
      !   QJ(i,3,ny) = rho2_init * uy * Jacobian_tmp
      !   QJ(i,4,ny) = (p2_init * over_gamma_1 + 0.5d0 * rho2_init * (ux**2 + uy**2)) * Jacobian_tmp
      ! endif

      ! ---- 外部基準状態（x依存） ----
      if (i < i_switch) then
        rho_ext = rho0_init
        u_ext   = u0_init
        v_ext   = v0_init !0.d0
        p_ext   = p0_init
      else
        rho_ext = rho2_init
        u_ext   = ux_init
        v_ext   = uy_init
        p_ext   = p2_init
      endif
      c_ext = sqrt(gamma * p_ext / rho_ext)

      ! ---- 内部状態 (j = ny-1) ----
      rhoin = QJ(i,1,ny-1) * Jacobian(i,ny-1)
      uin   = QJ(i,2,ny-1) / QJ(i,1,ny-1)
      vin   = QJ(i,3,ny-1) / QJ(i,1,ny-1)
      pin   = gamma_1 * ( QJ(i,4,ny-1) * Jacobian(i,ny-1) &
            - 0.5d0 * rhoin * (uin**2 + vin**2) )
      cin   = sqrt(gamma * pin / rhoin)

      ! ---- リーマン不変量（法線=y方向） ----
      Rp = vin   + 2.d0 * cin   * over_gamma_1   ! 常に流出（内部から）
      Rm = v_ext - 2.d0 * c_ext * over_gamma_1   ! 常に流入（外部基準から）
      vb = 0.5d0 * (Rp + Rm)
      cb = 0.25d0 * gamma_1 * (Rp - Rm)

      ! ---- エントロピー・接線速度：vbの符号で切替 ----
      if (vb >= 0.d0) then
        sb = pin / rhoin**gamma      ! 流出：内部から
        ub = uin
      else
        sb = p_ext / rho_ext**gamma  ! 流入：外部基準から
        ub = u_ext
      endif

      rhob = (cb**2 / (gamma * sb)) ** over_gamma_1
      pb   = sb * rhob**gamma

      Jacobian_tmp = 1.d0 / Jacobian(i,ny)
      QJ(i,1,ny) = rhob * Jacobian_tmp
      QJ(i,2,ny) = rhob * ub * Jacobian_tmp
      QJ(i,3,ny) = rhob * vb * Jacobian_tmp
      QJ(i,4,ny) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2)) * Jacobian_tmp

      ! NoSlip
      QJ(i,1,1) = QJ(i,1,2)
      QJ(i,2,1) = 0.d0
      QJ(i,3,1) = 0.d0
      p_wall = gamma_1 * (QJ(i,4,2) - 0.5d0 * (QJ(i,2,2)**2 + QJ(i,3,2)**2) / QJ(i,1,2))
      QJ(i,4,1) = p_wall * over_gamma_1
    enddo
  end subroutine set_bc
end module set
