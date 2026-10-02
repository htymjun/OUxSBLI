program convert_2d_to_3d
  implicit none
  real(8), parameter :: blt = 1.d-3
  ! mesh
  real(8), parameter :: Lx = 165.d0 * blt
  real(8), parameter :: Ly = 56.d0 * blt
  real(8), parameter :: Lz = 4.d0 * blt   ! ★ スパン方向長さ（3D版mod_globalsのLzと一致させること）
  ! DNS
  integer, parameter :: nx = 257
  integer, parameter :: ny = 480
  integer, parameter :: nz_new = 54
  integer, parameter :: nz_int = nz_new - 6   ! 1ランクあたりの内部平面数（ゴースト3枚ずつ除く）
  integer, parameter :: n_compute_ranks = 1   ! ★ 3D計算で使う計算ランク数（GPU数）に合わせる
  ! flat-plate geometry
  real(8), parameter :: x_in = -5.d0 * blt
  real(8), parameter :: Xsh  = 100.d0 * blt
  integer, parameter :: i_LE = nint(-x_in * dble(nx-1) / Lx) + 1

  ! --- 擾乱パラメータ ---
  real(8), parameter :: eps      = 0.0001d0         ! 擾乱振幅（u0に対する比率などお好みで調整）
  integer, parameter :: n_beta   = 3               ! スパン方向に入れる周期数（波数 = 2*pi*n_beta/Lz）
  real(8), parameter :: pi       = 4.d0 * atan(1.d0)
  real(8), parameter :: delta99  = 1.d0 * blt      ! Xshでのδ99 = 1mm
  real(8), parameter :: y_center = 0.1d0 * delta99 ! 境界層中央付近に擾乱を集中
  real(8), parameter :: y_width  = 0.05d0 * delta99 ! 境界層厚みの広がり
  real(8), parameter :: s        = 2.8d0           ! set_gridと同じtanhクラスタリング係数

  real(8), parameter :: R     = 287.15d0
  real(8), parameter :: gamma = 1.4d0
  real(8), parameter :: M0    = 2.15d0
  real(8), parameter :: p_tot = 100.d3 !25.d3 Re_x = 1.25 * 10**5 (Xsh = 100.d0 * blt)
  real(8), parameter :: p0    = p_tot / ((1.d0 + 0.5d0 * (gamma - 1.d0) * M0**2)**(gamma/(gamma-1.d0)))
  real(8), parameter :: T0    = 288.15d0
  real(8), parameter :: rho0  = p0 / (R * T0)
  real(8), parameter :: u0    = M0 * sqrt(gamma * R * T0)   ! 3D版mod_globalsのu0と同じ計算式

  real(8), allocatable :: Q2d(:,:,:), Q3d(:,:,:,:)
  real(8), allocatable :: zcoord(:), ycoord(:)
  integer i, j, k, irank, iz_offset
  real(8) dz1, yi, tanh_s, shape_y
  real(8) rho_local, rhov_old, rhov_new, delta
  character(len=40) filename

  allocate(Q2d(nx,ny,4))
  allocate(Q3d(nx,ny,nz_new,5))
  allocate(zcoord(nz_new))
  allocate(ycoord(ny))

  ! --- 2次元のリスタートファイルを読み込み（全ランク共通の元データ） ---
  open(10, file="recal/Q00001_nx257_np170_restart.dat", form="unformatted", access="stream", status="old")
  read(10) Q2d
  close(10)

  ! --- yグリッドをset_gridと同じ規約で再構築 ---
  tanh_s = tanh(s)
  do j = 1, ny
    yi = dble(j-1) / dble(ny-1)
    ycoord(j) = Ly * (1.d0 - tanh(s * (1.d0 - yi)) / tanh_s)
  enddo

  ! --- グローバルなdz（set_gridと同じ規約：全ランク分の内部平面で1周期） ---
  dz1 = Lz / dble(n_compute_ranks * nz_int)

  ! --- 計算ランクごとにファイルを生成 ---
  do irank = 0, n_compute_ranks - 1
    iz_offset = irank * nz_int

    ! 全z平面（ゴースト含む）に同じ値を複製し、rho*w=0を挿入
    do k = 1, nz_new
      Q3d(:,:,k,1) = Q2d(:,:,1)   ! rho
      Q3d(:,:,k,2) = Q2d(:,:,2)   ! rho*u
      Q3d(:,:,k,3) = Q2d(:,:,3)   ! rho*v
      Q3d(:,:,k,4) = 0.d0         ! rho*w = 0
      Q3d(:,:,k,5) = Q2d(:,:,4)   ! E
      ! set_gridと同じ規約でこのランクが担当するグローバルz座標を計算
      zcoord(k) = dble(iz_offset + k - 4) * dz1 - 0.5d0 * Lz
    enddo

    ! --- 壁法線方向速度(v)にスパン正弦波擾乱を追加。yはガウシアンで滑らかに分配 ---
    ! 壁(j=1)はNoSlipで上書きされるので j=2 から対象にする
    do k = 1, nz_new
      do j = 2, ny
        shape_y = exp(-((ycoord(j) - y_center) / y_width)**2)
        do i = i_LE, nx
          rho_local = Q3d(i,j,k,1)
          rhov_old  = Q3d(i,j,k,3)
          delta     = rho_local * eps * u0 * shape_y &
                    * sin(2.d0 * pi * dble(n_beta) * zcoord(k) / Lz)
          rhov_new  = rhov_old + delta

          ! 圧力を保ったまま、運動エネルギー変化分だけEを補正
          Q3d(i,j,k,5) = Q3d(i,j,k,5) + 0.5d0 * (rhov_new**2 - rhov_old**2) / rho_local
          Q3d(i,j,k,3) = rhov_new
        enddo
      enddo
    enddo

    ! --- このランク用のリスタートファイルとして保存 ---
    write(filename, "(a, i5.5, a)") "recal/Q_3d", irank+1, ".dat"
    open(20, file=filename, form="unformatted", access="stream", status="replace")
    write(20) Q3d
    close(20)
  enddo
end program convert_2d_to_3d