module set
  use mod_globals , only : nx, ny, nz, Lx, Ly, Lz, gamma, R,dt,nt, &
                           Lx_main, Lx_buf,Ly_main, Ly_buf,Lz_main,Lz_buf,nx_main, nx_buf,ny_main, ny_buf, nz_main, nz_buf
  use set_bc_common
  use set_coordinate
  implicit none

  ! グリッド生成用の内部ルーチン(他モジュールから見せない)
  private :: partitioned_axis, stretched_widths

  real(8), device :: rho_target(ny,nz), u_target(ny,nz), inlet_envelope(ny,nz)
  real(8), device :: sigma_x_1d(nx), sigma_y_1d(ny), sigma_z_1d(nz)
  real(8), device :: rho_target_x(ny,nz), u_target_x(ny,nz)
  real(8), device :: v_target_x(ny,nz), w_target_x(ny,nz), p_target_x(ny,nz)
  real(8), device :: over_jacobian_tab(nx,ny)
  ! QJ retains a reference z width; residuals are scaled by dz_ref/dz_cell.
  real(8), device :: z_residual_scale(nz)
  integer, save :: inlet_bc_calls = 0
  integer, parameter :: nsp_x = nx_buf, nsp_y = 7
  integer, parameter :: nsp_z = max(1, nint(dble(nz_buf*nsp_y)/dble(ny_buf)))
  real(8), parameter :: sigma_max_x = 0.06d0, sigma_max_y = 0.01d0

contains


!=====================================================================
! Deterministic, stateless pseudo-random numbers for the GPU inlet.
! The SplitMix64 hash makes the result independent of CUDA scheduling.
!=====================================================================
pure attributes(device) function inlet_uniform(j, k, frame, stream, seed) result(r)
  integer, intent(in), value :: j, k, frame, stream, seed
  integer(8) :: x, bits
  real(8) :: r

  x = int(seed,8) + 104729_8*int(j,8) + 130363_8*int(k,8) &
      + 15485863_8*int(frame,8) + 32452843_8*int(stream,8)
  x = x + int(z'9E3779B97F4A7C15',8)
  x = ieor(x, ishft(x,-30)) * int(z'BF58476D1CE4E5B9',8)
  x = ieor(x, ishft(x,-27)) * int(z'94D049BB133111EB',8)
  x = ieor(x, ishft(x,-31))
  bits = iand(ishft(x,-11), int(z'001FFFFFFFFFFFFF',8))
  r = (real(bits,8) + 0.5d0) * 1.1102230246251565d-16
end function inlet_uniform


pure attributes(device) function inlet_normal(j, k, frame, component, seed) result(g)
  integer, intent(in), value :: j, k, frame, component, seed
  real(8) :: g, r1, r2
  real(8), parameter :: two_pi = 6.2831853071795864769d0

  r1 = inlet_uniform(j, k, frame, 2*component-1, seed)
  r2 = inlet_uniform(j, k, frame, 2*component,   seed)
  g = sqrt(-2.d0*log(max(r1,1.d-15))) * cos(two_pi*r2)
end function inlet_normal


subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  use mod_constant, only : id_accuracy
  use mod_globals, only : dt, CFL, u1, u2, gamma, R, T1, T2, endT, np, Pr
  integer, intent(in)  :: myrank, nx, ny, nz
  real(8), intent(in)  :: Lx, Ly, Lz
  real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)

  real(8) :: dx_min, dy_min, dz_min, dmin, c1, c2, umax, cmax, dt_suggest
  real(8) :: vmax_est, wmax_est, dt_conv, dt_diff, mu_max, rho_min
  integer :: nt_suggest

  call set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
                            nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)
  dx_min = minval(dx)
  dy_min = minval(dy)
  dz_min = minval(dz)

  c1 = sqrt(gamma*R*T1)
  c2 = sqrt(gamma*R*T2)
  umax = max(abs(u1), abs(u2))
  cmax = max(c1, c2)

  ! --- 渦発達を見込んだ速度の安全係数（経験的に1.5〜2倍） ---
  vmax_est = 1.5d0 * umax    ! 渦のピーク速度を見込む
  wmax_est = 1.5d0 * umax

  ! --- 3方向を合算した対流CFL（より保守的） ---
  dt_conv = CFL / ( (umax+cmax)/dx_min + (vmax_est+cmax)/dy_min + (wmax_est+cmax)/dz_min )

  ! --- 拡散(粘性)のCFL制約も確認 ---
  mu_max  = 1.716d-5 * (T2/273.15d0)**1.5d0 * (273.15d0+110.4d0)/(T2+110.4d0)
  rho_min = 8.d4 / (R * max(T1,T2))

  dt_diff = CFL * 0.5d0 * min(dx_min,dy_min,dz_min)**2 * rho_min / mu_max

  dt_suggest = min(dt_conv, dt_diff)
  nt_suggest = int(endT / (dble(np) * dt_suggest))

  print *, "===================================================="
  print *, "dx_min, dy_min, dz_min =", dx_min, dy_min, dz_min
  print *, "dt_conv (対流CFL, 3D合算) =", dt_conv
  print *, "dt_diff (拡散/粘性制約)   =", dt_diff
  print *, "Suggested dt (min)        =", dt_suggest
  print *, "Current  dt (mod_globals) =", dt
  print *, "Ratio current/suggested   =", dt / dt_suggest
  print *, "Suggested nt =", nt_suggest
  print *, "===================================================="

end subroutine set_grid

!=====================================================================
! 等比数列の公比rを求める（一定伸び率版：現在は使用しないが互換のため残す）
! 条件: dr_0 * (r^n - 1)/(r-1) = L_target
!=====================================================================
real(8) function calc_stretch_ratio(dr_0, n, L_target) result(r)
  real(8), intent(in) :: dr_0, L_target
  integer, intent(in) :: n
  real(8) :: r_lo, r_hi, r_mid, f_mid
  integer :: iter
  real(8), parameter :: tol = 1.d-12

  if (abs(dr_0*dble(n) - L_target) < tol) then
    r = 1.0d0
    return
  endif

  if (dr_0*dble(n) < L_target) then
    r_lo = 1.0d0 + tol
    r_hi = 10.0d0
  else
    r_lo = tol
    r_hi = 1.0d0 - tol
  endif

  do iter = 1, 200
    r_mid = 0.5d0*(r_lo + r_hi)
    f_mid = dr_0*(r_mid**dble(n) - 1.0d0)/(r_mid - 1.0d0) - L_target
    if (f_mid > 0.d0) then
      r_hi = r_mid
    else
      r_lo = r_mid
    endif
    if (abs(r_hi - r_lo) < tol) exit
  enddo
  r = 0.5d0*(r_lo + r_hi)

end function calc_stretch_ratio


!=====================================================================
! 滑らかランプ版: セル幅の伸び率を接続点の局所伸び率 r0 から
! 終端伸び率 r_target までコサインランプで滑らかに変化させ、
! 合計長さが L_target になる r_target を二分法で求める。
! ジャンクション(主計算領域との接続点)での伸び率の急変を避けるための関数。
!=====================================================================
real(8) function calc_smooth_stretch_ratio(dr_0, r0, n, L_target) result(r_target)
  real(8), intent(in) :: dr_0    ! 接続点でのセル幅(主計算領域の最後のセル幅)
  real(8), intent(in) :: r0      ! 接続点での局所伸び率(主計算領域側から引き継ぐ、通常1に近い)
  integer, intent(in) :: n       ! バッファのセル数
  real(8), intent(in) :: L_target
  real(8) :: r_lo, r_hi, r_mid, total_mid
  integer :: iter
  real(8), parameter :: tol = 1.d-12

  if (n < 2 .or. dr_0 <= 0.d0 .or. L_target <= 0.d0 .or. r0 < 1.d0-1.d-12) &
    error stop 'Invalid stretched buffer parameters'
  r_lo = 1.d0
  if (calc_ramp_total(dr_0, r0, r_lo, n) > L_target+tol) &
    error stop 'Buffer too short for main-domain edge spacing'
  r_hi = max(1.01d0, r0)
  do while (calc_ramp_total(dr_0, r0, r_hi, n) < L_target)
    r_hi = 2.d0*r_hi
  enddo

  do iter = 1, 200
    r_mid = 0.5d0*(r_lo + r_hi)
    total_mid = calc_ramp_total(dr_0, r0, r_mid, n)
    if (total_mid < L_target) then
      r_lo = r_mid
    else
      r_hi = r_mid
    endif
    if (abs(r_hi - r_lo) < tol) exit
  enddo
  r_target = 0.5d0*(r_lo + r_hi)

contains
  real(8) function calc_ramp_total(dr0_in, r0_in, r1_in, n_in) result(total)
    real(8), intent(in) :: dr0_in, r0_in, r1_in
    integer, intent(in) :: n_in
    real(8) :: dr_local, eta, r_local, pi_local
    integer :: k
    pi_local = acos(-1.d0)
    dr_local = dr0_in
    total = 0.d0
    do k = 1, n_in
      eta = dble(k-1) / dble(max(n_in-1,1))
      r_local = r0_in + (r1_in - r0_in) * (1.0d0 - cos(pi_local*eta)) / 2.0d0
      total = total + dr_local     ! 先に加算(代入してからdrを伸ばす)
      if (total > L_target) return ! bracket comparison only; avoid overflow
      dr_local = dr_local * r_local
    enddo
  end function calc_ramp_total
end function calc_smooth_stretch_ratio


!=====================================================================
! Host-only main-domain mesh generation, independent of buffer parameters.
! 中心一様区間 + 外側伸長区間(5次smoothstepで接続点の格子幅を滑らかに接続)
!=====================================================================
subroutine partitioned_axis(n, length, uniform_length, n_uniform, symmetric, coord)
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  integer, intent(in) :: n, n_uniform
  real(8), intent(in) :: length, uniform_length
  logical, intent(in) :: symmetric
  real(8), intent(out) :: coord(n)
  real(8) :: h, width(n-1), side_length
  integer :: m, n_side, first, j

  if (n < 2 .or. n_uniform < 2 .or. n_uniform > n) error stop 'Invalid uniform point count'
  if (.not.ieee_is_finite(length) .or. .not.ieee_is_finite(uniform_length)) &
    error stop 'Nonfinite mesh length'
  if (uniform_length <= 0.d0 .or. uniform_length > length) error stop 'Invalid uniform length'
  h = uniform_length / dble(n_uniform-1)

  if (.not.symmetric) then
    coord(1) = 0.d0
    do j = 2, n_uniform
      coord(j) = dble(j-1)*h
    enddo
    m = n-n_uniform
    call stretched_widths(m, length-uniform_length, h, width(1:m))
    do j = 1, m
      coord(n_uniform+j) = coord(n_uniform+j-1)+width(j)
    enddo
  else
    if (mod(n-n_uniform,2) /= 0) error stop 'Symmetric mesh requires even n-n_uniform'
    n_side = (n-n_uniform)/2
    first = n_side+1
    side_length = 0.5d0*(length-uniform_length)
    do j = 0, n_uniform-1
      coord(first+j) = side_length+dble(j)*h
    enddo
    call stretched_widths(n_side, side_length, h, width(1:n_side))
    do j = 1, n_side
      coord(first-j) = coord(first-j+1)-width(j)
      coord(first+n_uniform-1+j) = coord(first+n_uniform-2+j)+width(j)
    enddo
  endif

  if (abs(coord(1)) > 1.d-11*length .or. abs(coord(n)-length) > 1.d-11*length) &
    error stop 'Main mesh length mismatch'
end subroutine partitioned_axis


subroutine stretched_widths(n, length, h, width)
  integer, intent(in) :: n
  real(8), intent(in) :: length, h
  real(8), intent(out) :: width(n)
  real(8) :: t, excess, weight(n), amplitude, tol
  integer :: j

  tol = 1.d-12*max(h, length)
  excess = length-dble(n)*h
  if (n == 0) then
    if (abs(length) > tol) error stop 'Nonzero length without stretched intervals'
    return
  endif
  if (excess < -tol) &
    error stop 'Stretch section too short: increase uniform point count or reduce uniform length'
  if (n == 1) then
    if (abs(excess) > tol) error stop 'Need at least two stretched intervals for smooth matching'
    width = h
    return
  endif
  ! Width equals h at the junction; quintic ramp has zero slope at both ends.
  do j = 1, n
    t = dble(j-1)/dble(n-1)
    weight(j) = t**3*(10.d0-15.d0*t+6.d0*t*t)
  enddo
  amplitude = max(0.d0, excess)/sum(weight)
  width = h+amplitude*weight
end subroutine stretched_widths


subroutine set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
  nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)

  use mod_globals, only : Ly_uniform, ny_uniform, Lz_uniform, nz_uniform
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: Lx_main, Lx_buf, Ly_main, Ly_buf, Lz
  integer, intent(in) :: nx_main, nx_buf, ny_main, ny_buf
  real(8), intent(out) :: xc(nx), yc(ny), zc(nz)
  real(8), intent(out) :: dx(nx-1), dy(ny-1), dz(nz-1)

  real(8) :: x(nx+1), y(ny+1), z(nz+1)
  real(8) :: dr, eta, pi
  integer :: i, j, k

  ! 格子伸長パラメータ
  real(8)  :: r0_x, r1_x               ! x方向バッファ: 開始/終端伸び率
  real(8)  :: r_local

  pi = acos(-1.d0)

  !=================================================================
  ! x方向: 主計算領域は全域一様 + 流出側バッファ(滑らかランプ伸長)
  !=================================================================
  if (nx_main < 3 .or. Lx_main <= 0.d0) error stop 'Invalid uniform x mesh'
  do i = 2, nx_main+1
    x(i) = Lx_main*dble(i-2)/dble(nx_main-1)
  enddo

  ! 接続点の格子幅。主計算領域は一様なので局所伸び率は1
  dr   = x(nx_main+1) - x(nx_main)
  r0_x = 1.d0

  r1_x = calc_smooth_stretch_ratio(dr, r0_x, nx_buf, Lx_buf)

  do i = nx_main+2, nx-1
    eta = dble(i - (nx_main+2)) / dble(max(nx_buf-1,1))
    r_local = r0_x + (r1_x - r0_x) * (1.0d0 - cos(pi*eta)) / 2.0d0
    x(i) = x(i-1) + dr
    dr = dr * r_local
  enddo

  ! 左ゴースト(i=1): 主計算領域流入端の格子幅を複製
  x(1) = x(2) - (x(3)-x(2))
  ! 右ゴースト
  x(nx)   = x(nx-1) + dr
  x(nx+1) = x(nx)   + dr

  print *, "x(nx-1) - x(nx_main+1) =", x(nx-1) - x(nx_main+1), " (目標:", Lx_buf, ")"
  print *, "xバッファ 開始伸び率=", r0_x, " 終端伸び率=", r1_x

  !=================================================================
  ! y/z: 中心一様 + 外側伸長の主領域と、滑らかランプのバッファ(y,z共通)
  !=================================================================
  call transverse_faces(ny, ny_main, ny_buf, Ly_main, Ly_buf, Ly_uniform, ny_uniform, y)
  call transverse_faces(nz, nz_main, nz_buf, Lz_main, Lz_buf, Lz_uniform, nz_uniform, z)

  !=================================================================
  ! セル中心・格子間隔
  !=================================================================
  do i = 1, nx; xc(i) = 0.5d0*(x(i)+x(i+1)); enddo
  do j = 1, ny; yc(j) = 0.5d0*(y(j)+y(j+1)); enddo
  do k = 1, nz; zc(k) = 0.5d0*(z(k)+z(k+1)); enddo

  do i = 1, nx-1; dx(i) = xc(i+1)-xc(i); enddo
  do j = 1, ny-1; dy(j) = yc(j+1)-yc(j); enddo
  do k = 1, nz-1; dz(k) = zc(k+1)-zc(k); enddo
!=================================================================
! Cache the same cell-center cross profile used by set_init.
!=================================================================
block
  real(8) :: rho_h(ny,nz), u_h(ny,nz), env_h(ny,nz)
  real(8) :: sx_h(nx), sy_h(ny), sz_h(nz), temp_h
  integer :: jj, kk, ii
  do kk = 1, nz
    do jj = 1, ny
      call cross_profile(yc(jj), zc(kk), u_h(jj,kk), temp_h, rho_h(jj,kk), env_h(jj,kk))
    enddo
  enddo
  do ii = 1, nx
    sx_h(ii) = sponge_weight(dble(ii-(nx-1-nsp_x))/dble(nsp_x), sigma_max_x)
  enddo
  do jj = 1, ny
    sy_h(jj) = sponge_weight(max(dble(1+nsp_y-jj), dble(jj-(ny-nsp_y)))/dble(nsp_y), sigma_max_y)
  enddo
  do kk = 1, nz
    sz_h(kk) = sponge_weight(max(dble(1+nsp_z-kk), dble(kk-(nz-nsp_z)))/dble(nsp_z), sigma_max_y)
  enddo
  rho_target = rho_h
  u_target = u_h
  inlet_envelope = env_h
  sigma_x_1d = sx_h
  sigma_y_1d = sy_h
  sigma_z_1d = sz_h
end block
  !=================================================================
  ! 確認出力
  !=================================================================
  print *, "=== 格子生成確認 ==="
  print *, "x: 主計算領域 ", x(2), "~", x(nx_main+1), " バッファ ~", x(nx-1)
  print *, "dx_main(一様) =", x(3)-x(2), " dx_buf_max =", x(nx-1)-x(nx-2)
  print *, "y: バッファ下 ", y(2), "~", y(ny_buf+2)
  print *, "y: 主計算領域 ", y(ny_buf+2), "~", y(ny_buf+ny_main+1)
  print *, "y: バッファ上 ", y(ny_buf+ny_main+1), "~", y(ny-1)
  print *, "dy_center(最小格子幅) =", y(ny_buf+2+(ny_main/2)+1)-y(ny_buf+2+(ny_main/2))
  print *, "dy_edge  (主計算端)   =", y(ny_buf+3)-y(ny_buf+2)
  print *, "dy_buf_max(バッファ端) =", y(ny-1)-y(ny-2)
  print *, "dy_edge/dy_center     =", (y(ny_buf+3)-y(ny_buf+2)) / &
                                       (y(ny_buf+2+(ny_main/2)+1)-y(ny_buf+2+(ny_main/2)))
  print *, "dy_buf_max/dy_edge    =", (y(ny-1)-y(ny-2))/(y(ny_buf+3)-y(ny_buf+2))

  print *, "z: main faces =", z(nz_buf+2), z(nz_buf+nz_main+1)
  print *, "z: buffer extents =", z(nz_buf+2)-z(2), z(nz-1)-z(nz_buf+nz_main+1)
  print *, "dz min/max =", minval(dz), maxval(dz)
end subroutine set_grid_main_buffer

  !=====================================================================
  ! y/z共通の面座標生成:
  !   主領域 = 中心一様区間 + 外側伸長(partitioned_axis, 対称)
  !   バッファ = 接続点の伸び率を引き継ぐコサインランプ伸長
  !=====================================================================
  subroutine transverse_faces(n, nmain, nbuf, lmain, lbuf, luni, nuni, faces)
    integer, intent(in) :: n, nmain, nbuf, nuni
    real(8), intent(in) :: lmain, lbuf, luni
    real(8), intent(out) :: faces(n+1)
    real(8) :: dr, eta, pi, r_local, r0_lo, r1_lo, r0_hi, r1_hi
    integer :: m
    if (nmain < 4 .or. nbuf < 2 .or. n /= nmain+2*nbuf+2) error stop 'Invalid transverse mesh'
    pi = acos(-1.d0)

    ! 主領域: 0..lmain で生成し、下バッファ長さ分だけ平行移動
    call partitioned_axis(nmain, lmain, luni, nuni, .true., &
                          faces(nbuf+2:nbuf+nmain+1))
    faces(nbuf+2:nbuf+nmain+1) = faces(nbuf+2:nbuf+nmain+1) + lbuf

    ! --- 下バッファ: 接続点の格子幅と、その外側への局所伸び率を引き継ぐ ---
    dr = faces(nbuf+3) - faces(nbuf+2)
    r0_lo = dr / (faces(nbuf+4) - faces(nbuf+3))
    r1_lo = calc_smooth_stretch_ratio(dr, r0_lo, nbuf, lbuf)

    do m = nbuf+1, 2, -1
      eta = dble((nbuf+2-1) - m) / dble(max(nbuf-1,1))
      r_local = r0_lo + (r1_lo - r0_lo) * (1.0d0 - cos(pi*eta)) / 2.0d0
      faces(m) = faces(m+1) - dr
      dr = dr * r_local
    enddo

    ! --- 上バッファ ---
    dr = faces(nbuf+nmain+1) - faces(nbuf+nmain)
    r0_hi = dr / (faces(nbuf+nmain) - faces(nbuf+nmain-1))
    r1_hi = calc_smooth_stretch_ratio(dr, r0_hi, nbuf, lbuf)

    do m = nbuf+nmain+2, n-1
      eta = dble(m - (nbuf+nmain+2)) / dble(max(nbuf-1,1))
      r_local = r0_hi + (r1_hi - r0_hi) * (1.0d0 - cos(pi*eta)) / 2.0d0
      faces(m) = faces(m-1) + dr
      dr = dr * r_local
    enddo

    ! 下ゴースト(m=1)
    faces(1)   = faces(2) - (faces(3)-faces(2))
    ! 上ゴースト
    faces(n)   = faces(n-1) + dr
    faces(n+1) = faces(n)   + dr
  end subroutine transverse_faces

  !=====================================================================
  ! ヤコビアン: 内部は中央差分、ゴーストセル(1点)は隣接内部セルの値を複製
  !=====================================================================
  subroutine set_Jacobian_xy3_stretch(nx, ny, nz, xc, yc, zc, Jacobian)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: xc(nx), yc(ny), zc(nz)
    real(8), intent(out) :: Jacobian(nx,ny)
    integer :: i, j, k
    real(8) :: dx(nx-1), dy(ny-1), dzdzeta, z_scale_h(nz)

    dx = xc(2:nx) - xc(1:nx-1)
    dy = yc(2:ny) - yc(1:ny-1)

    ! Reference volume Vref=dx_cell*dy_cell*dz_ref, not the physical volume.
    ! calc_R multiplies the area-weighted residual by Vref/Vcell.
    dzdzeta = 0.5d0*(zc(nz/2+1)-zc(nz/2-1))
    do k = 2, nz-1
      z_scale_h(k) = dzdzeta/(0.5d0*(zc(k+1)-zc(k-1)))
    enddo
    z_scale_h(1) = z_scale_h(2)
    z_scale_h(nz) = z_scale_h(nz-1)
    z_residual_scale = z_scale_h

    do j = 2, ny-1
      do i = 2, nx-1
        ! Same spacing-based form as set_Jacobian_xy3, retaining dz_ref.
        Jacobian(i,j) = 8.d0 / &
          ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * 2.d0*dzdzeta)
      enddo;enddo

    Jacobian(1,:)  = Jacobian(2,:)
    Jacobian(nx,:) = Jacobian(nx-1,:)
    Jacobian(:,1)  = Jacobian(:,2)
    Jacobian(:,ny) = Jacobian(:,ny-1)
    Jacobian(1,1)   = Jacobian(2,2)
    Jacobian(1,ny)  = Jacobian(2,ny-1)
    Jacobian(nx,1)  = Jacobian(nx-1,2)
    Jacobian(nx,ny) = Jacobian(nx-1,ny-1)

    ! 1/Jacobian も同時に1回だけ計算してGPUへ転送しておく(set_bc内での毎回の割り算を排除)
    over_jacobian_tab = 1.0d0 / Jacobian

  end subroutine set_Jacobian_xy3_stretch


  pure subroutine cross_profile(y, z, velocity, temperature, density, envelope)
    use mod_globals, only : u1, u2, T1, T2, p, delta_bl
    real(8), intent(in) :: y, z
    real(8), intent(out) :: velocity, temperature, density, envelope
    real(8) :: ty, tz, shape
    ty = tanh(2.d0*(y-0.5d0*Ly)/delta_bl)
    tz = tanh(2.d0*(z-0.5d0*Lz)/delta_bl)
    shape = ty*tz
    velocity = 0.5d0*(u1+u2) + 0.5d0*(u1-u2)*shape
    temperature = 0.5d0*(T1+T2) + 0.5d0*(T1-T2)*shape
    density = p/(R*temperature)
    envelope = max(0.d0, 1.d0-shape**2)
  end subroutine cross_profile

  pure function sponge_weight(eta, strength) result(weight)
    real(8), intent(in) :: eta, strength
    real(8) :: weight, t, ramp
    t = min(1.d0, max(0.d0, eta))
    ramp = t**3*(10.d0-15.d0*t+6.d0*t*t)
    weight = 1.d0-(1.d0-strength)**ramp
  end function sponge_weight

  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : p
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,5,ny,nz)
    integer :: i, j, k
    real(8) :: vel, temp, rho, env
    do k = 1, nz
      do j = 1, ny
        call cross_profile(y(j), z(k), vel, temp, rho, env)
        do i = 1, nx
          Q(i,1,j,k) = rho
          Q(i,2,j,k) = rho*vel
          Q(i,3:4,j,k) = 0.d0
          Q(i,5,j,k) = p/(gamma-1.d0)+0.5d0*rho*vel**2
        enddo
      enddo
    enddo
  end subroutine set_init

subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
  use mod_constant, only : id_accuracy
  use mod_globals, only : u1, rho1, u2, rho2, p, amp, dt, gamma, Ly, Lz, Lx,T1,T2,delta_bl, step_offset, &
                          inlet_fluctuation_rms, inlet_random_seed, inlet_rk_stages
  use set_coordinate
  integer, intent(in), value            :: myrank, nx, ny, nz
  real(8), intent(in), device           :: jacobian(nx,ny)
  real(8), intent(inout), device        :: Q(nx,5,ny,nz)
  real(8), intent(in), device, optional :: Qre(:)
  real(8) :: u, v, w, u_init, rho_init
  integer :: i, j, k, l, i_target, inlet_frame

  real(8) :: sigma_x, sigma_y
  real(8) :: fluctuation_scale
  ! スポンジ層での基本変数ブレンド用(block構文はCUFカーネル内でエラーになりうるため
  ! サブルーチンレベルのスカラーとして宣言する)
  real(8) :: rho_now, u_now, v_now, w_now, p_now
  real(8) :: rho_blend, u_blend, v_blend, w_blend, p_blend

    !===================================================================
    ! x = 1 flow inlet: shear-layer-localized Gaussian fluctuations in
    ! u, v and w.
    ! One independent Gaussian random field is retained for all RK stages
    ! and renewed once per physical time step.  At every physical (y,z)
    ! point, N(0,1) samples are multiplied by
    ! inlet_fluctuation_rms*abs(u_init)*(1-tanh_y**2*tanh_z**2).
    !===================================================================
    inlet_frame = inlet_bc_calls / inlet_rk_stages
    inlet_bc_calls = inlet_bc_calls + 1

    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        i = 1
        rho_init = rho_target(j,k)
        u_init   = u_target(j,k)
        fluctuation_scale = inlet_fluctuation_rms*abs(u_init)*inlet_envelope(j,k)
        u = u_init + fluctuation_scale*inlet_normal(j, k, inlet_frame, 1, inlet_random_seed)
        v =fluctuation_scale*inlet_normal(j, k, inlet_frame, 2, inlet_random_seed)
        w =fluctuation_scale*inlet_normal(j, k, inlet_frame, 3, inlet_random_seed)

        Q(i,1,j,k) = rho_init*over_jacobian_tab(i,j)
        Q(i,2,j,k) = rho_init*u*over_jacobian_tab(i,j)
        Q(i,3,j,k) = rho_init*v*over_jacobian_tab(i,j)
        Q(i,4,j,k) = rho_init*w*over_jacobian_tab(i,j)
        Q(i,5,j,k) = (p/(gamma-1.d0) + 0.5d0*rho_init*(u**2+v**2+w**2))*over_jacobian_tab(i,j)
      enddo;enddo

    !===================================================================
    ! Capture the full yz plane before damping.  A z average would erase
    ! the reversed high/low streams.  Keep each (j,k) target independently.
    !===================================================================
    i_target = nx - 1 - nsp_x
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        rho_now = Q(i_target,1,j,k)*jacobian(i_target,j)
        u_now = Q(i_target,2,j,k)/Q(i_target,1,j,k)
        v_now = Q(i_target,3,j,k)/Q(i_target,1,j,k)
        w_now = Q(i_target,4,j,k)/Q(i_target,1,j,k)
        p_now = (gamma-1.d0)*(Q(i_target,5,j,k)*jacobian(i_target,j) &
                -0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
        rho_target_x(j,k) = rho_now
        u_target_x(j,k) = u_now
        v_target_x(j,k) = v_now
        w_target_x(j,k) = w_now
        p_target_x(j,k) = p_now
      enddo
    enddo

    !===================================================================
    ! xスポンジだけを処理する。全領域を走査せず、対象となる末端nsp_x点に
    ! カーネル範囲を限定してGPU上の不要なメモリアクセスを避ける。
    !===================================================================
    !$cuf kernel do(3)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do i = nx-nsp_x, nx-1
          sigma_x = sigma_x_1d(i)
          if (sigma_x > 0.d0) then
            rho_now = Q(i,1,j,k) * jacobian(i,j)
            u_now   = Q(i,2,j,k) / Q(i,1,j,k)
            v_now   = Q(i,3,j,k) / Q(i,1,j,k)
            w_now   = Q(i,4,j,k) / Q(i,1,j,k)
            p_now   = (gamma-1.d0) * (Q(i,5,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))

            rho_blend = (1.d0-sigma_x)*rho_now + sigma_x*rho_target_x(j,k)
            u_blend   = (1.d0-sigma_x)*u_now   + sigma_x*u_target_x(j,k)
            v_blend   = (1.d0-sigma_x)*v_now   + sigma_x*v_target_x(j,k)
            w_blend   = (1.d0-sigma_x)*w_now   + sigma_x*w_target_x(j,k)
            p_blend   = (1.d0-sigma_x)*p_now   + sigma_x*p_target_x(j,k)
            Q(i,1,j,k) = rho_blend * over_jacobian_tab(i,j)
            Q(i,2,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
            Q(i,3,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
            Q(i,4,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
            Q(i,5,j,k) = (p_blend/(gamma-1.d0) &
                         + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)

          endif
        enddo;enddo;enddo

    !===================================================================
    ! Apply y/z sponge weights together to preserve symmetry at corners.
    ! Both sides relax to the same yz-dependent cross profile as the inlet.
    !===================================================================
    !$cuf kernel do(3)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          sigma_y = 1.d0-(1.d0-sigma_y_1d(j))*(1.d0-sigma_z_1d(k))
          if (sigma_y > 0.d0) then
            rho_now = Q(i,1,j,k) * jacobian(i,j)
            u_now   = Q(i,2,j,k) / Q(i,1,j,k)
            v_now   = Q(i,3,j,k) / Q(i,1,j,k)
            w_now   = Q(i,4,j,k) / Q(i,1,j,k)
            p_now   = (gamma-1.d0) * (Q(i,5,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
            rho_init = rho_target(j,k)
            u_init   = u_target(j,k)
            rho_blend = (1.d0-sigma_y)*rho_now + sigma_y*rho_init
            u_blend   = (1.d0-sigma_y)*u_now   + sigma_y*u_init
            v_blend   = (1.d0-sigma_y)*v_now
            w_blend   = (1.d0-sigma_y)*w_now
            p_blend   = (1.d0-sigma_y)*p_now + sigma_y*p
            Q(i,1,j,k) = rho_blend * over_jacobian_tab(i,j)
            Q(i,2,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
            Q(i,3,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
            Q(i,4,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
            Q(i,5,j,k) = (p_blend/(gamma-1.d0) &
                         + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)
          endif
        enddo;enddo;enddo

    !===================================================================
    ! x = nx 流出境界(1点ゴースト): 物理量で1次外挿してJ補正
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do l = 1, 5
          Q(nx,l,j,k) = ( 2.0d0*(Q(nx-1,l,j,k)*jacobian(nx-1,j)) &
                             - (Q(nx-2,l,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)
        enddo;enddo;enddo 
    !===================================================================
    ! y方向境界(j=1, j=ny): 0次外挿。物理保存量をコピーしてJ補正。
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
        do l = 1, 5
          Q(i,l,ny,k) = Q(i,l,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)

          Q(i,l,1,k) = Q(i,l,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
        enddo;enddo;enddo 
    !===================================================================
    ! z: zeroth-order extrapolation copies the adjacent interior state.
    ! Jacobian is independent of z, so its factors cancel here.
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(i,l,j,1) = Q(i,l,j,2)
          Q(i,l,j,nz) = Q(i,l,j,nz-1)
        enddo
      enddo
    enddo

  end subroutine set_bc

  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer :: i, j, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        mut(1,j,k) = mut(2,j,k); mut(nx,j,k) = mut(nx-1,j,k)
        qc2(1,j,k) = qc2(2,j,k); qc2(nx,j,k) = qc2(nx-1,j,k)
      enddo
    enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 1, nx
        mut(i,1,k) = mut(i,2,k); mut(i,ny,k) = mut(i,ny-1,k)
        qc2(i,1,k) = qc2(i,2,k); qc2(i,ny,k) = qc2(i,ny-1,k)
      enddo
    enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        mut(i,j,1) = mut(i,j,2); mut(i,j,nz) = mut(i,j,nz-1)
        qc2(i,j,1) = qc2(i,j,2); qc2(i,j,nz) = qc2(i,j,nz-1)
      enddo
    enddo
  end subroutine set_bc_mut

end module set

! module set
!   use mod_globals , only : nx, ny, nz, Lx, Ly, Lz, gamma, R,dt,nt, &
!                            Lx_main, Lx_buf,Ly_main, Ly_buf,Lz_main,Lz_buf,nx_main, nx_buf,ny_main, ny_buf, nz_main, nz_buf
!   use set_bc_common
!   use set_coordinate
!   implicit none


!   real(8), device :: rho_target(ny,nz), u_target(ny,nz), inlet_envelope(ny,nz)
!   real(8), device :: sigma_x_1d(nx), sigma_y_1d(ny), sigma_z_1d(nz)
!   real(8), device :: rho_target_x(ny,nz), u_target_x(ny,nz)
!   real(8), device :: v_target_x(ny,nz), w_target_x(ny,nz), p_target_x(ny,nz)
!   real(8), device :: over_jacobian_tab(nx,ny)
!   ! QJ retains a reference z width; residuals are scaled by dz_ref/dz_cell.
!   real(8), device :: z_residual_scale(nz)
!   integer, save :: inlet_bc_calls = 0
!   integer, parameter :: nsp_x = nx_buf, nsp_y = 7
!   integer, parameter :: nsp_z = max(1, nint(dble(nz_buf*nsp_y)/dble(ny_buf)))
!   real(8), parameter :: sigma_max_x = 0.06d0, sigma_max_y = 0.01d0

! contains


! !=====================================================================
! ! Deterministic, stateless pseudo-random numbers for the GPU inlet.
! ! The SplitMix64 hash makes the result independent of CUDA scheduling.
! !=====================================================================
! pure attributes(device) function inlet_uniform(j, k, frame, stream, seed) result(r)
!   integer, intent(in), value :: j, k, frame, stream, seed
!   integer(8) :: x, bits
!   real(8) :: r

!   x = int(seed,8) + 104729_8*int(j,8) + 130363_8*int(k,8) &
!       + 15485863_8*int(frame,8) + 32452843_8*int(stream,8)
!   x = x + int(z'9E3779B97F4A7C15',8)
!   x = ieor(x, ishft(x,-30)) * int(z'BF58476D1CE4E5B9',8)
!   x = ieor(x, ishft(x,-27)) * int(z'94D049BB133111EB',8)
!   x = ieor(x, ishft(x,-31))
!   bits = iand(ishft(x,-11), int(z'001FFFFFFFFFFFFF',8))
!   r = (real(bits,8) + 0.5d0) * 1.1102230246251565d-16
! end function inlet_uniform


! pure attributes(device) function inlet_normal(j, k, frame, component, seed) result(g)
!   integer, intent(in), value :: j, k, frame, component, seed
!   real(8) :: g, r1, r2
!   real(8), parameter :: two_pi = 6.2831853071795864769d0

!   r1 = inlet_uniform(j, k, frame, 2*component-1, seed)
!   r2 = inlet_uniform(j, k, frame, 2*component,   seed)
!   g = sqrt(-2.d0*log(max(r1,1.d-15))) * cos(two_pi*r2)
! end function inlet_normal


! subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
!   use mod_constant, only : id_accuracy
!   use mod_globals, only : dt, CFL, u1, u2, gamma, R, T1, T2, endT, np, Pr
!   integer, intent(in)  :: myrank, nx, ny, nz
!   real(8), intent(in)  :: Lx, Ly, Lz
!   real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)

!   real(8) :: dx_min, dy_min, dz_min, dmin, c1, c2, umax, cmax, dt_suggest
!   real(8) :: vmax_est, wmax_est, dt_conv, dt_diff, mu_max, rho_min
!   integer :: nt_suggest

!   call set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
!                             nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)
!   dx_min = minval(dx)
!   dy_min = minval(dy)
!   dz_min = minval(dz)

!   c1 = sqrt(gamma*R*T1)
!   c2 = sqrt(gamma*R*T2)
!   umax = max(abs(u1), abs(u2))
!   cmax = max(c1, c2)

!   ! --- 渦発達を見込んだ速度の安全係数（経験的に1.5〜2倍） ---
!   vmax_est = 1.5d0 * umax    ! 渦のピーク速度を見込む
!   wmax_est = 1.5d0 * umax

!   ! --- 3方向を合算した対流CFL（より保守的） ---
!   dt_conv = CFL / ( (umax+cmax)/dx_min + (vmax_est+cmax)/dy_min + (wmax_est+cmax)/dz_min )

!   ! --- 拡散(粘性)のCFL制約も確認 ---
!   mu_max  = 1.716d-5 * (T2/273.15d0)**1.5d0 * (273.15d0+110.4d0)/(T2+110.4d0)
!   rho_min = 8.d4 / (R * max(T1,T2))

!   dt_diff = CFL * 0.5d0 * min(dx_min,dy_min,dz_min)**2 * rho_min / mu_max

!   dt_suggest = min(dt_conv, dt_diff)
!   nt_suggest = int(endT / (dble(np) * dt_suggest))

!   print *, "===================================================="
!   print *, "dx_min, dy_min, dz_min =", dx_min, dy_min, dz_min
!   print *, "dt_conv (対流CFL, 3D合算) =", dt_conv
!   print *, "dt_diff (拡散/粘性制約)   =", dt_diff
!   print *, "Suggested dt (min)        =", dt_suggest
!   print *, "Current  dt (mod_globals) =", dt
!   print *, "Ratio current/suggested   =", dt / dt_suggest
!   print *, "Suggested nt =", nt_suggest
!   print *, "===================================================="

! end subroutine set_grid

! !=====================================================================
! ! 等比数列の公比rを求める（一定伸び率版：現在は使用しないが互換のため残す）
! ! 条件: dr_0 * (r^n - 1)/(r-1) = L_target
! !=====================================================================
! real(8) function calc_stretch_ratio(dr_0, n, L_target) result(r)
!   real(8), intent(in) :: dr_0, L_target
!   integer, intent(in) :: n
!   real(8) :: r_lo, r_hi, r_mid, f_mid
!   integer :: iter
!   real(8), parameter :: tol = 1.d-12

!   if (abs(dr_0*dble(n) - L_target) < tol) then
!     r = 1.0d0
!     return
!   endif

!   if (dr_0*dble(n) < L_target) then
!     r_lo = 1.0d0 + tol
!     r_hi = 10.0d0
!   else
!     r_lo = tol
!     r_hi = 1.0d0 - tol
!   endif

!   do iter = 1, 200
!     r_mid = 0.5d0*(r_lo + r_hi)
!     f_mid = dr_0*(r_mid**dble(n) - 1.0d0)/(r_mid - 1.0d0) - L_target
!     if (f_mid > 0.d0) then
!       r_hi = r_mid
!     else
!       r_lo = r_mid
!     endif
!     if (abs(r_hi - r_lo) < tol) exit
!   enddo
!   r = 0.5d0*(r_lo + r_hi)

! end function calc_stretch_ratio


! !=====================================================================
! ! 滑らかランプ版: セル幅の伸び率を接続点の局所伸び率 r0 から
! ! 終端伸び率 r_target までコサインランプで滑らかに変化させ、
! ! 合計長さが L_target になる r_target を二分法で求める。
! ! ジャンクション(主計算領域との接続点)での伸び率の急変を避けるための関数。
! !=====================================================================
! real(8) function calc_smooth_stretch_ratio(dr_0, r0, n, L_target) result(r_target)
!   real(8), intent(in) :: dr_0    ! 接続点でのセル幅(主計算領域の最後のセル幅)
!   real(8), intent(in) :: r0      ! 接続点での局所伸び率(主計算領域側から引き継ぐ、通常1に近い)
!   integer, intent(in) :: n       ! バッファのセル数
!   real(8), intent(in) :: L_target
!   real(8) :: r_lo, r_hi, r_mid, total_mid
!   integer :: iter
!   real(8), parameter :: tol = 1.d-12

!   r_lo = r0
!   r_hi = 1000.0d0

!   do iter = 1, 200
!     r_mid = 0.5d0*(r_lo + r_hi)
!     total_mid = calc_ramp_total(dr_0, r0, r_mid, n)
!     if (total_mid < L_target) then
!       r_lo = r_mid
!     else
!       r_hi = r_mid
!     endif
!     if (abs(r_hi - r_lo) < tol) exit
!   enddo
!   r_target = 0.5d0*(r_lo + r_hi)

! contains
!   real(8) function calc_ramp_total(dr0_in, r0_in, r1_in, n_in) result(total)
!     real(8), intent(in) :: dr0_in, r0_in, r1_in
!     integer, intent(in) :: n_in
!     real(8) :: dr_local, eta, r_local, pi_local
!     integer :: k
!     pi_local = acos(-1.d0)
!     dr_local = dr0_in
!     total = 0.d0
!     do k = 1, n_in
!       eta = dble(k-1) / dble(max(n_in-1,1))
!       r_local = r0_in + (r1_in - r0_in) * (1.0d0 - cos(pi_local*eta)) / 2.0d0
!       total = total + dr_local     ! 先に加算(現行と同じ順序: 代入してからdrを伸ばす)
!       dr_local = dr_local * r_local
!     enddo
!   end function calc_ramp_total
! end function calc_smooth_stretch_ratio


! subroutine set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
!   nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)

!   integer, intent(in) :: nx, ny, nz
!   real(8), intent(in) :: Lx_main, Lx_buf, Ly_main, Ly_buf, Lz
!   integer, intent(in) :: nx_main, nx_buf, ny_main, ny_buf
!   real(8), intent(out) :: xc(nx), yc(ny), zc(nz)
!   real(8), intent(out) :: dx(nx-1), dy(ny-1), dz(nz-1)

!   real(8) :: x(nx+1), y(ny+1), z(nz+1)
!   real(8) :: dx_main
!   real(8) :: s, dr, eta, pi
!   integer :: i, j, k

!   ! 格子伸長パラメータ
!   real(8), parameter :: hx_pow   = 1.1d0 !1.0d0  ! x方向主計算領域(流入側集中)冪乗
!   real(8)  :: r0_x, r1_x               ! x方向バッファ: 開始/終端伸び率
!   real(8)  :: r_local

!   pi = acos(-1.d0)

!   dx_main = Lx_main * hx_pow * (1.0d0)**(hx_pow-1.0d0) / dble(nx_main-1)

!   !=================================================================
!   ! x方向: 主計算領域(冪乗、流入側集中) + 流出側バッファ(滑らかランプ伸長)
!   !=================================================================
!   do i = 2, nx_main+1
!     s = dble(i-2) / dble(nx_main-1)
!     x(i) = Lx_main * s**hx_pow
!   enddo

!   ! 接続点の格子幅とその直前の局所伸び率(主計算領域側から引き継ぐ)
!   dr  = x(nx_main+1) - x(nx_main)
!   r0_x = dr / (x(nx_main) - x(nx_main-1))

!   r1_x = calc_smooth_stretch_ratio(dr, r0_x, nx_buf, Lx_buf)

!   do i = nx_main+2, nx-1
!     eta = dble(i - (nx_main+2)) / dble(max(nx_buf-1,1))
!     r_local = r0_x + (r1_x - r0_x) * (1.0d0 - cos(pi*eta)) / 2.0d0
!     x(i) = x(i-1) + dr
!     dr = dr * r_local
!   enddo

!   ! 左ゴースト(i=1): 主計算領域流入端の格子幅を複製
!   x(1) = x(2) - (x(3)-x(2))
!   ! 右ゴースト
!   x(nx)   = x(nx-1) + dr
!   x(nx+1) = x(nx)   + dr

!   print *, "x(nx-1) - x(nx_main+1) =", x(nx-1) - x(nx_main+1), " (目標:", Lx_buf, ")"
!   print *, "xバッファ 開始伸び率=", r0_x, " 終端伸び率=", r1_x

!   !=================================================================
!   ! y/z: identical center-clustered core and cosine-ramp buffer stretching.
!   !=================================================================
!   call transverse_faces(ny, ny_main, ny_buf, Ly_main, Ly_buf, y)
!   call transverse_faces(nz, nz_main, nz_buf, Lz_main, Lz_buf, z)

!   !=================================================================
!   ! セル中心・格子間隔
!   !=================================================================
!   do i = 1, nx; xc(i) = 0.5d0*(x(i)+x(i+1)); enddo
!   do j = 1, ny; yc(j) = 0.5d0*(y(j)+y(j+1)); enddo
!   do k = 1, nz; zc(k) = 0.5d0*(z(k)+z(k+1)); enddo

!   do i = 1, nx-1; dx(i) = xc(i+1)-xc(i); enddo
!   do j = 1, ny-1; dy(j) = yc(j+1)-yc(j); enddo
!   do k = 1, nz-1; dz(k) = zc(k+1)-zc(k); enddo
! !=================================================================
! ! Cache the same cell-center cross profile used by set_init.
! !=================================================================
! block
!   real(8) :: rho_h(ny,nz), u_h(ny,nz), env_h(ny,nz)
!   real(8) :: sx_h(nx), sy_h(ny), sz_h(nz), temp_h
!   integer :: jj, kk, ii
!   do kk = 1, nz
!     do jj = 1, ny
!       call cross_profile(yc(jj), zc(kk), u_h(jj,kk), temp_h, rho_h(jj,kk), env_h(jj,kk))
!     enddo
!   enddo
!   do ii = 1, nx
!     sx_h(ii) = sponge_weight(dble(ii-(nx-1-nsp_x))/dble(nsp_x), sigma_max_x)
!   enddo
!   do jj = 1, ny
!     sy_h(jj) = sponge_weight(max(dble(1+nsp_y-jj), dble(jj-(ny-nsp_y)))/dble(nsp_y), sigma_max_y)
!   enddo
!   do kk = 1, nz
!     sz_h(kk) = sponge_weight(max(dble(1+nsp_z-kk), dble(kk-(nz-nsp_z)))/dble(nsp_z), sigma_max_y)
!   enddo
!   rho_target = rho_h
!   u_target = u_h
!   inlet_envelope = env_h
!   sigma_x_1d = sx_h
!   sigma_y_1d = sy_h
!   sigma_z_1d = sz_h
! end block
!   !=================================================================
!   ! 確認出力
!   !=================================================================
!   print *, "=== 格子生成確認 ==="
!   print *, "x: 主計算領域 ", x(2), "~", x(nx_main+1), " バッファ ~", x(nx-1)
!   print *, "y: バッファ下 ", y(2), "~", y(ny_buf+2)
!   print *, "y: 主計算領域 ", y(ny_buf+2), "~", y(ny_buf+ny_main+1)
!   print *, "y: バッファ上 ", y(ny_buf+ny_main+1), "~", y(ny-1)
!   print *, "dy_center(最小格子幅) =", y(ny_buf+2+(ny_main/2)+1)-y(ny_buf+2+(ny_main/2))
!   print *, "dy_edge  (主計算端)   =", y(ny_buf+3)-y(ny_buf+2)
!   print *, "dy_buf_max(バッファ端) =", y(ny-1)-y(ny-2)
!   print *, "dy_edge/dy_center     =", (y(ny_buf+3)-y(ny_buf+2)) / &
!                                        (y(ny_buf+2+(ny_main/2)+1)-y(ny_buf+2+(ny_main/2)))
!   print *, "dy_buf_max/dy_edge    =", (y(ny-1)-y(ny-2))/(y(ny_buf+3)-y(ny_buf+2))

!   print *, "z: main faces =", z(nz_buf+2), z(nz_buf+nz_main+1)
!   print *, "z: buffer extents =", z(nz_buf+2)-z(2), z(nz-1)-z(nz_buf+nz_main+1)
!   print *, "dz min/max =", minval(dz), maxval(dz)
! end subroutine set_grid_main_buffer

!   subroutine transverse_faces(n, nmain, nbuf, lmain, lbuf, faces)
!     integer, intent(in) :: n, nmain, nbuf
!     real(8), intent(in) :: lmain, lbuf
!     real(8), intent(out) :: faces(n+1)
!     real(8), parameter :: beta = 1.7d0
!     real(8) :: s, dr, eta, pi, r_local, r0_lo, r1_lo, r0_hi, r1_hi
!     integer :: m
!     if (nmain < 4 .or. nbuf < 2 .or. n /= nmain+2*nbuf+2) error stop 'Invalid transverse mesh'
!     pi = acos(-1.d0)
!   do m = nbuf+2, nbuf+nmain+1
!     s = 2.0d0*dble(m-(nbuf+2))/dble(nmain-1) - 1.0d0
!     faces(m) = lbuf + 0.5d0*lmain*(1.0d0 + atanh(s*tanh(beta))/beta)
!   enddo

!   ! --- 下バッファ: 接続点の格子幅と、その外側への局所伸び率を引き継ぐ ---
!   dr = faces(nbuf+3) - faces(nbuf+2)
!   r0_lo = dr / (faces(nbuf+4) - faces(nbuf+3))

!   r1_lo = calc_smooth_stretch_ratio(dr, r0_lo, nbuf, lbuf)


!   do m = nbuf+1, 2, -1
!     eta = dble((nbuf+2-1) - m) / dble(max(nbuf-1,1))
!     r_local = r0_lo + (r1_lo - r0_lo) * (1.0d0 - cos(pi*eta)) / 2.0d0
!     faces(m) = faces(m+1) - dr
!     dr = dr * r_local
!   enddo


!   ! --- 上バッファ ---
!   dr = faces(nbuf+nmain+1) - faces(nbuf+nmain)
!   r0_hi = dr / (faces(nbuf+nmain) - faces(nbuf+nmain-1))

!   r1_hi = calc_smooth_stretch_ratio(dr, r0_hi, nbuf, lbuf)


!   do m = nbuf+nmain+2, n-1
!     eta = dble(m - (nbuf+nmain+2)) / dble(max(nbuf-1,1))
!     r_local = r0_hi + (r1_hi - r0_hi) * (1.0d0 - cos(pi*eta)) / 2.0d0
!     faces(m) = faces(m-1) + dr
!     dr = dr * r_local
!   enddo


!   ! 下ゴースト(m=1)
!   faces(1) = faces(2) - (faces(3)-faces(2))
!   ! 上ゴースト
!   faces(n)   = faces(n-1) + dr
!   faces(n+1) = faces(n)   + dr

!   end subroutine transverse_faces

!   !=====================================================================
!   ! ヤコビアン: 内部は中央差分、ゴーストセル(1点)は隣接内部セルの値を複製
!   !=====================================================================
!   subroutine set_Jacobian_xy3_stretch(nx, ny, nz, xc, yc, zc, Jacobian)
!     integer, intent(in)  :: nx, ny, nz
!     real(8), intent(in)  :: xc(nx), yc(ny), zc(nz)
!     real(8), intent(out) :: Jacobian(nx,ny)
!     integer :: i, j, k
!     real(8) :: dx(nx-1), dy(ny-1), dzdzeta, z_scale_h(nz)

!     dx = xc(2:nx) - xc(1:nx-1)
!     dy = yc(2:ny) - yc(1:ny-1)

!     ! Reference volume Vref=dx_cell*dy_cell*dz_ref, not the physical volume.
!     ! calc_R multiplies the area-weighted residual by Vref/Vcell.
!     dzdzeta = 0.5d0*(zc(nz/2+1)-zc(nz/2-1))
!     do k = 2, nz-1
!       z_scale_h(k) = dzdzeta/(0.5d0*(zc(k+1)-zc(k-1)))
!     enddo
!     z_scale_h(1) = z_scale_h(2)
!     z_scale_h(nz) = z_scale_h(nz-1)
!     z_residual_scale = z_scale_h

!     do j = 2, ny-1
!       do i = 2, nx-1
!         ! Same spacing-based form as set_Jacobian_xy3, retaining dz_ref.
!         Jacobian(i,j) = 8.d0 / &
!           ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * 2.d0*dzdzeta)
!       enddo;enddo

!     Jacobian(1,:)  = Jacobian(2,:)
!     Jacobian(nx,:) = Jacobian(nx-1,:)
!     Jacobian(:,1)  = Jacobian(:,2)
!     Jacobian(:,ny) = Jacobian(:,ny-1)
!     Jacobian(1,1)   = Jacobian(2,2)
!     Jacobian(1,ny)  = Jacobian(2,ny-1)
!     Jacobian(nx,1)  = Jacobian(nx-1,2)
!     Jacobian(nx,ny) = Jacobian(nx-1,ny-1)

!     ! 1/Jacobian も同時に1回だけ計算してGPUへ転送しておく(set_bc内での毎回の割り算を排除)
!     over_jacobian_tab = 1.0d0 / Jacobian

!   end subroutine set_Jacobian_xy3_stretch


!   pure subroutine cross_profile(y, z, velocity, temperature, density, envelope)
!     use mod_globals, only : u1, u2, T1, T2, p, delta_bl
!     real(8), intent(in) :: y, z
!     real(8), intent(out) :: velocity, temperature, density, envelope
!     real(8) :: ty, tz, shape
!     ty = tanh(2.d0*(y-0.5d0*Ly)/delta_bl)
!     tz = tanh(2.d0*(z-0.5d0*Lz)/delta_bl)
!     shape = ty*tz
!     velocity = 0.5d0*(u1+u2) + 0.5d0*(u1-u2)*shape
!     temperature = 0.5d0*(T1+T2) + 0.5d0*(T1-T2)*shape
!     density = p/(R*temperature)
!     envelope = max(0.d0, 1.d0-shape**2)
!   end subroutine cross_profile

!   pure function sponge_weight(eta, strength) result(weight)
!     real(8), intent(in) :: eta, strength
!     real(8) :: weight, t, ramp
!     t = min(1.d0, max(0.d0, eta))
!     ramp = t**3*(10.d0-15.d0*t+6.d0*t*t)
!     weight = 1.d0-(1.d0-strength)**ramp
!   end function sponge_weight

!   subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
!     use mod_globals, only : p
!     integer, intent(in) :: myrank, nx, ny, nz
!     real(8), intent(in) :: x(nx), y(ny), z(nz)
!     real(8), intent(out) :: Q(nx,5,ny,nz)
!     integer :: i, j, k
!     real(8) :: vel, temp, rho, env
!     do k = 1, nz
!       do j = 1, ny
!         call cross_profile(y(j), z(k), vel, temp, rho, env)
!         do i = 1, nx
!           Q(i,1,j,k) = rho
!           Q(i,2,j,k) = rho*vel
!           Q(i,3:4,j,k) = 0.d0
!           Q(i,5,j,k) = p/(gamma-1.d0)+0.5d0*rho*vel**2
!         enddo
!       enddo
!     enddo
!   end subroutine set_init

! subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
!   use mod_constant, only : id_accuracy
!   use mod_globals, only : u1, rho1, u2, rho2, p, amp, dt, gamma, Ly, Lz, Lx,T1,T2,delta_bl, step_offset, &
!                           inlet_fluctuation_rms, inlet_random_seed, inlet_rk_stages
!   use set_coordinate
!   integer, intent(in), value            :: myrank, nx, ny, nz
!   real(8), intent(in), device           :: jacobian(nx,ny)
!   real(8), intent(inout), device        :: Q(nx,5,ny,nz)
!   real(8), intent(in), device, optional :: Qre(:)
!   real(8) :: u, v, w, u_init, rho_init
!   integer :: i, j, k, l, i_target, inlet_frame

!   real(8) :: sigma_x, sigma_y
!   real(8) :: fluctuation_scale
!   ! スポンジ層での基本変数ブレンド用(block構文はCUFカーネル内でエラーになりうるため
!   ! サブルーチンレベルのスカラーとして宣言する)
!   real(8) :: rho_now, u_now, v_now, w_now, p_now
!   real(8) :: rho_blend, u_blend, v_blend, w_blend, p_blend

!     !===================================================================
!     ! x = 1 flow inlet: shear-layer-localized Gaussian fluctuations in
!     ! u, v and w.
!     ! One independent Gaussian random field is retained for all RK stages
!     ! and renewed once per physical time step.  At every physical (y,z)
!     ! point, N(0,1) samples are multiplied by
!     ! inlet_fluctuation_rms*abs(u_init)*(1-tanh_y**2*tanh_z**2).
!     !===================================================================
!     inlet_frame = inlet_bc_calls / inlet_rk_stages
!     inlet_bc_calls = inlet_bc_calls + 1

!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 2, nz-1
!       do j = 2, ny-1
!         i = 1
!         rho_init = rho_target(j,k)
!         u_init   = u_target(j,k)
!         fluctuation_scale = inlet_fluctuation_rms*abs(u_init)*inlet_envelope(j,k)
!         u = u_init + fluctuation_scale*inlet_normal(j, k, inlet_frame, 1, inlet_random_seed)
!         v =fluctuation_scale*inlet_normal(j, k, inlet_frame, 2, inlet_random_seed)
!         w =fluctuation_scale*inlet_normal(j, k, inlet_frame, 3, inlet_random_seed)

!         Q(i,1,j,k) = rho_init*over_jacobian_tab(i,j)
!         Q(i,2,j,k) = rho_init*u*over_jacobian_tab(i,j)
!         Q(i,3,j,k) = rho_init*v*over_jacobian_tab(i,j)
!         Q(i,4,j,k) = rho_init*w*over_jacobian_tab(i,j)
!         Q(i,5,j,k) = (p/(gamma-1.d0) + 0.5d0*rho_init*(u**2+v**2+w**2))*over_jacobian_tab(i,j)
!       enddo;enddo

!     !===================================================================
!     ! Capture the full yz plane before damping.  A z average would erase
!     ! the reversed high/low streams.  Keep each (j,k) target independently.
!     !===================================================================
!     i_target = nx - 1 - nsp_x
!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 2, nz-1
!       do j = 2, ny-1
!         rho_now = Q(i_target,1,j,k)*jacobian(i_target,j)
!         u_now = Q(i_target,2,j,k)/Q(i_target,1,j,k)
!         v_now = Q(i_target,3,j,k)/Q(i_target,1,j,k)
!         w_now = Q(i_target,4,j,k)/Q(i_target,1,j,k)
!         p_now = (gamma-1.d0)*(Q(i_target,5,j,k)*jacobian(i_target,j) &
!                 -0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
!         rho_target_x(j,k) = rho_now
!         u_target_x(j,k) = u_now
!         v_target_x(j,k) = v_now
!         w_target_x(j,k) = w_now
!         p_target_x(j,k) = p_now
!       enddo
!     enddo

!     !===================================================================
!     ! xスポンジだけを処理する。全領域を走査せず、対象となる末端nsp_x点に
!     ! カーネル範囲を限定してGPU上の不要なメモリアクセスを避ける。
!     !===================================================================
!     !$cuf kernel do(3)<<<*,*>>>
!     do k = 2, nz-1
!       do j = 2, ny-1
!         do i = nx-nsp_x, nx-1
!           sigma_x = sigma_x_1d(i)
!           if (sigma_x > 0.d0) then
!             rho_now = Q(i,1,j,k) * jacobian(i,j)
!             u_now   = Q(i,2,j,k) / Q(i,1,j,k)
!             v_now   = Q(i,3,j,k) / Q(i,1,j,k)
!             w_now   = Q(i,4,j,k) / Q(i,1,j,k)
!             p_now   = (gamma-1.d0) * (Q(i,5,j,k)*jacobian(i,j) &
!                       - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))

!             rho_blend = (1.d0-sigma_x)*rho_now + sigma_x*rho_target_x(j,k)
!             u_blend   = (1.d0-sigma_x)*u_now   + sigma_x*u_target_x(j,k)
!             v_blend   = (1.d0-sigma_x)*v_now   + sigma_x*v_target_x(j,k)
!             w_blend   = (1.d0-sigma_x)*w_now   + sigma_x*w_target_x(j,k)
!             p_blend   = (1.d0-sigma_x)*p_now   + sigma_x*p_target_x(j,k)
!             Q(i,1,j,k) = rho_blend * over_jacobian_tab(i,j)
!             Q(i,2,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
!             Q(i,3,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
!             Q(i,4,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
!             Q(i,5,j,k) = (p_blend/(gamma-1.d0) &
!                          + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)

!           endif
!         enddo;enddo;enddo

!     !===================================================================
!     ! Apply y/z sponge weights together to preserve symmetry at corners.
!     ! Both sides relax to the same yz-dependent cross profile as the inlet.
!     !===================================================================
!     !$cuf kernel do(3)<<<*,*>>>
!     do k = 2, nz-1
!       do j = 2, ny-1
!         do i = 2, nx-1
!           sigma_y = 1.d0-(1.d0-sigma_y_1d(j))*(1.d0-sigma_z_1d(k))
!           if (sigma_y > 0.d0) then
!             rho_now = Q(i,1,j,k) * jacobian(i,j)
!             u_now   = Q(i,2,j,k) / Q(i,1,j,k)
!             v_now   = Q(i,3,j,k) / Q(i,1,j,k)
!             w_now   = Q(i,4,j,k) / Q(i,1,j,k)
!             p_now   = (gamma-1.d0) * (Q(i,5,j,k)*jacobian(i,j) &
!                       - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
!             rho_init = rho_target(j,k)
!             u_init   = u_target(j,k)
!             rho_blend = (1.d0-sigma_y)*rho_now + sigma_y*rho_init
!             u_blend   = (1.d0-sigma_y)*u_now   + sigma_y*u_init
!             v_blend   = (1.d0-sigma_y)*v_now
!             w_blend   = (1.d0-sigma_y)*w_now
!             p_blend   = (1.d0-sigma_y)*p_now + sigma_y*p
!             Q(i,1,j,k) = rho_blend * over_jacobian_tab(i,j)
!             Q(i,2,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
!             Q(i,3,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
!             Q(i,4,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
!             Q(i,5,j,k) = (p_blend/(gamma-1.d0) &
!                          + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)
!           endif
!         enddo;enddo;enddo

!     !===================================================================
!     ! x = nx 流出境界(1点ゴースト): 物理量で1次外挿してJ補正
!     !===================================================================
!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 1, nz
!       do j = 1, ny
!         do l = 1, 5
!           Q(nx,l,j,k) = ( 2.0d0*(Q(nx-1,l,j,k)*jacobian(nx-1,j)) &
!                              - (Q(nx-2,l,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)
!         enddo;enddo;enddo 
!     !===================================================================
!     ! y方向境界(j=1, j=ny): 0次外挿。物理保存量をコピーしてJ補正。
!     !===================================================================
!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 1, nz
!       do i = 1, nx
!         do l = 1, 5
!           Q(i,l,ny,k) = Q(i,l,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)

!           Q(i,l,1,k) = Q(i,l,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
!         enddo;enddo;enddo 
!     !===================================================================
!     ! z: zeroth-order extrapolation copies the adjacent interior state.
!     ! Jacobian is independent of z, so its factors cancel here.
!     !===================================================================
!     !$cuf kernel do(2)<<<*,*>>>
!     do j = 1, ny
!       do i = 1, nx
!         do l = 1, 5
!           Q(i,l,j,1) = Q(i,l,j,2)
!           Q(i,l,j,nz) = Q(i,l,j,nz-1)
!         enddo
!       enddo
!     enddo

!   end subroutine set_bc

!   subroutine set_bc_mut(nx, ny, nz, mut, qc2)
!     integer, intent(in), value     :: nx, ny, nz
!     real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
!     integer :: i, j, k
!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 2, nz-1
!       do j = 2, ny-1
!         mut(1,j,k) = mut(2,j,k); mut(nx,j,k) = mut(nx-1,j,k)
!         qc2(1,j,k) = qc2(2,j,k); qc2(nx,j,k) = qc2(nx-1,j,k)
!       enddo
!     enddo
!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 2, nz-1
!       do i = 1, nx
!         mut(i,1,k) = mut(i,2,k); mut(i,ny,k) = mut(i,ny-1,k)
!         qc2(i,1,k) = qc2(i,2,k); qc2(i,ny,k) = qc2(i,ny-1,k)
!       enddo
!     enddo
!     !$cuf kernel do(2)<<<*,*>>>
!     do j = 1, ny
!       do i = 1, nx
!         mut(i,j,1) = mut(i,j,2); mut(i,j,nz) = mut(i,j,nz-1)
!         qc2(i,j,1) = qc2(i,j,2); qc2(i,j,nz) = qc2(i,j,nz-1)
!       enddo
!     enddo
!   end subroutine set_bc_mut

! end module set
