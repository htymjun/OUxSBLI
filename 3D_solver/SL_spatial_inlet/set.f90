module set
  use mod_globals , only : nx, ny, nz, Lx, Ly, Lz, gamma, R,dt,nt, &
                           Lx_main, Lx_buf,Ly_main, Ly_buf,nx_main, nx_buf,ny_main, ny_buf, &
                           npx, nx_global
  use set_bc_common
  use set_coordinate
  implicit none
!   call set_Jacobian_xy3_stretch(nx, ny, nz, x, y, z, Jacobian)
  ! グリッド生成用の内部ルーチン(他モジュールから見せない)
  private :: partitioned_axis, stretched_widths

  !共通変数
  real(8), device  :: xi(nx), yj(ny), zk(nz)
  real(8), device  :: rho_target_1d(ny), u_target_1d(ny)
  real(8), device  :: inlet_envelope_1d(ny)           ! shear-layer-localized inlet noise envelope
  real(8), device  :: amp_u_env_1d(ny)                ! amp*u_init(j)*env(j) を事前計算(追加最適化)
  real(8), device  :: sigma_x_1d(nx), sigma_y_1d(ny)  ! スポンジ係数(x,y)を事前計算
  real(8), device  :: rho_target_x_1d(ny), u_target_x_1d(ny)
  real(8), device  :: v_target_x_1d(ny), w_target_x_1d(ny), p_target_x_1d(ny)
  real(8), device  :: over_jacobian_tab(nx,ny)        ! 1/Jacobian を事前計算(格子は時間不変)
  integer, save    :: inlet_bc_calls = 0

  ! x方向MPI分割(mod_globalsのnpx、config.fyppはCOMMZ=False): 偶数(計算)ランクごとに1スラブ、
  ! 奇数(I/O)ランクは相方と同じスラブ。global i = local i + i_offset。接続側のゴーストは3面で、
  ! set_bcの最後(exchange_x)で隣スラブの内部3面を受け取る。共有カーネルは配列端から3面目以降の
  ! 面/セルを高次で計算するので、内部セルは分割なしと同じステンシルになる。
  ! zは分割しない: Lzは周期スパン、nzは全体の点数(ゴースト3+3点を含む)で、周期は自スラブ内コピー。
  integer, save    :: my_slab   = 0     ! このランクのxスラブ番号(0..npx-1)
  integer, save    :: i_offset  = 0     ! ローカル→全体のx番号オフセット
  integer, save    :: comm_compute      ! 計算ランクだけのコミュニケータ(x袖交換用)
  integer, parameter :: n_compute = 1   ! zスラブ数(z方向は分割しない)
  integer, parameter :: k_offset  = 0   ! zのグローバル番号オフセット(同上)
  ! 配列端のヤコビアンは隣の列のコピー(set_Jacobian_xy3)なので、端のゴースト面だけ
  ! QJ = Q/J を真のJとの比で換算して受け取る(一様部では1)。
  real(8), save    :: x_fac_lo = 1.d0, x_fac_hi = 1.d0
  real(8), allocatable, device :: xs_lo_d(:), xs_hi_d(:), xr_lo_d(:), xr_hi_d(:)
  ! ホスト側の送受信バッファはピン留めメモリ(ページ可能メモリより転送が速い)
  real(8), allocatable, pinned :: xs_lo_h(:), xs_hi_h(:), xr_lo_h(:), xr_hi_h(:)

  ! スポンジ・強制振動パラメータ(set_grid_main_bufferとset_bc双方から使うためモジュールスコープへ)
  ! ---- パッシブスカラー(混合分率 xi: 流れ1で1、流れ2で0)。mod_globalsのscalar_onで有効化 ----
  ! 流れ場の1ステップごと(最終RK段のset_bc)に、更新後の速度場でSSP-RK3により1ステップ進める。
  real(8), allocatable, device :: sc_phi(:,:,:)       ! xi (時刻n)
  real(8), allocatable, device :: sc_phis(:,:,:)      ! RK中間段
  real(8), allocatable, device :: sc_phib(:,:,:)      ! RK 段2の結果
  real(8), allocatable, device :: sc_gam(:,:,:)       ! rho*D = mu/Sc
  real(8), device  :: sc_mx(nx), sc_my(ny)            ! セル中心のメトリック di/dx, dj/dy
  real(8), device  :: sc_idx(nx-1), sc_idy(ny-1)      ! 1/(xc(i+1)-xc(i)), 1/(yc(j+1)-yc(j))
  real(8), device  :: sc_in(ny)                       ! 流入分布(速度のtanh分布と同じ形)
  real(8), save    :: sc_mz = 0.d0                    ! 1/dz
  real(4), save    :: sc_xc(nx), sc_yc(ny), sc_zc(nz)  ! VTK出力用のセル中心座標(流れ場のVTKと同じ)

  integer, parameter :: nsp_x = nx_buf
  integer, parameter :: nsp_y = ny_buf/2
  real(8), parameter :: sigma_max_x = 0.06d0
  real(8), parameter :: sigma_max_y = 0.01d0

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


!=====================================================================
! パッシブスカラー: 風上側から再構成したセル界面値(Korenリミッタ、TVD)。
! a=風上のさらに風上, b=風上, c=風下。値は必ずbとcの間に収まる。
! 勾配比 r=(c-b)/(b-a) の割り算を避け、psi(r)*(b-a) を直接評価する(割り算版と数学的に同じ):
!   s=sign(b-a), psi*(b-a) = s*max(0, min(2*s*(c-b), (|b-a|+2*s*(c-b))/3, 2*|b-a|))
!=====================================================================
pure attributes(device) function sc_face(a, b, c) result(f)
  real(8), intent(in), value :: a, b, c
  real(8) :: f, da, dc, s, m
  da = b - a
  dc = c - b
  s  = sign(1.d0, da)
  m  = max(0.d0, min(2.d0*s*dc, min((abs(da) + 2.d0*s*dc)/3.d0, 2.d0*abs(da))))
  f  = b + 0.5d0*s*m
end function sc_face


!=====================================================================
! パッシブスカラー: z方向ゴースト3面(周期)。zは分割しないので自スラブ内コピー。
!=====================================================================
subroutine sc_exchange_z(myrank, nx, ny, nz, f)
  integer, intent(in), value     :: myrank, nx, ny, nz
  real(8), intent(inout), device :: f(nx,ny,nz)
  integer :: i, j, k

  !$cuf kernel do(2)<<<*,*>>>
  do j = 1, ny
    do i = 1, nx
      do k = 1, 3
        f(i,j,k)      = f(i,j,nz-6+k)
        f(i,j,nz-3+k) = f(i,j,3+k)
  enddo;enddo;enddo
end subroutine sc_exchange_z


!=====================================================================
! パッシブスカラー: 境界条件。流入は固定分布、流出と上下は勾配ゼロ、zは周期。
!=====================================================================
subroutine sc_bc(myrank, nx, ny, nz, f)
  integer, intent(in), value     :: myrank, nx, ny, nz
  real(8), intent(inout), device :: f(nx,ny,nz)
  integer :: i, j, k

  if (my_slab == 0) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        f(1,j,k)  = sc_in(j)
    enddo;enddo
  endif
  if (my_slab == npx-1) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        f(nx,j,k) = f(nx-1,j,k)
    enddo;enddo
  endif

  !$cuf kernel do(2)<<<*,*>>>
  do k = 4, nz-3
    do i = 1, nx
      f(i,1,k)  = f(i,2,k)
      f(i,ny,k) = f(i,ny-1,k)
  enddo;enddo

  call sc_exchange_z(myrank, nx, ny, nz, f)
  call sc_exchange_x(nx, ny, nz, f)
end subroutine sc_bc


!=====================================================================
! パッシブスカラー: x方向の袖交換(非周期)。流れ場のexchange_xと同じ面を交換する。
! 送受信バッファは流れ場用(xs_*, xr_*)の先頭 3*ny*nz 要素を借りる。
!=====================================================================
subroutine sc_exchange_x(nx, ny, nz, f)
  use mpi
  integer, intent(in), value     :: nx, ny, nz
  real(8), intent(inout), device :: f(nx,ny,nz)
  integer :: left, right, n, ierr
  integer :: istat(MPI_STATUS_SIZE)

  if (npx == 1) return
  call alloc_exchange_x(ny, nz)
  n = 3*ny*nz
  left  = MPI_PROC_NULL;  if (my_slab > 0)     left  = my_slab - 1
  right = MPI_PROC_NULL;  if (my_slab < npx-1) right = my_slab + 1
  if (right /= MPI_PROC_NULL) then
    call sc_pack_x(nx, ny, nz, nx-5, f, xs_hi_d)
    xs_hi_h(1:n) = xs_hi_d(1:n)
  endif
  if (left /= MPI_PROC_NULL) then
    call sc_pack_x(nx, ny, nz, 4, f, xs_lo_d)
    xs_lo_h(1:n) = xs_lo_d(1:n)
  endif
  call MPI_SENDRECV(xs_hi_h, n, MPI_REAL8, right, 43, xr_lo_h, n, MPI_REAL8, left,  43, comm_compute, istat, ierr)
  call MPI_SENDRECV(xs_lo_h, n, MPI_REAL8, left,  44, xr_hi_h, n, MPI_REAL8, right, 44, comm_compute, istat, ierr)
  if (left /= MPI_PROC_NULL) then
    xr_lo_d(1:n) = xr_lo_h(1:n)
    call sc_unpack_x(nx, ny, nz, 1, xr_lo_d, f)
  endif
  if (right /= MPI_PROC_NULL) then
    xr_hi_d(1:n) = xr_hi_h(1:n)
    call sc_unpack_x(nx, ny, nz, nx-2, xr_hi_d, f)
  endif
end subroutine sc_exchange_x

subroutine sc_pack_x(nx, ny, nz, i0, f, buf)
  integer, intent(in), value     :: nx, ny, nz, i0
  real(8), intent(in), device    :: f(nx,ny,nz)
  real(8), intent(inout), device :: buf(3*ny*nz*5)
  integer :: j, k, m
  !$cuf kernel do(2)<<<*,*>>>
  do k = 1, nz
    do j = 1, ny
      do m = 1, 3
        buf(m + 3*((j-1) + ny*(k-1))) = f(i0+m-1,j,k)
  enddo;enddo;enddo
end subroutine sc_pack_x

subroutine sc_unpack_x(nx, ny, nz, i0, buf, f)
  integer, intent(in), value     :: nx, ny, nz, i0
  real(8), intent(in), device    :: buf(3*ny*nz*5)
  real(8), intent(inout), device :: f(nx,ny,nz)
  integer :: j, k, m
  !$cuf kernel do(2)<<<*,*>>>
  do k = 1, nz
    do j = 1, ny
      do m = 1, 3
        f(i0+m-1,j,k) = buf(m + 3*((j-1) + ny*(k-1)))
  enddo;enddo;enddo
end subroutine sc_unpack_x


!=====================================================================
! x方向の袖交換バッファ(流れ場5変数x3面。スカラーは先頭の1変数分を使う)
!=====================================================================
subroutine alloc_exchange_x(ny, nz)
  integer, intent(in) :: ny, nz
  integer :: n
  if (allocated(xs_lo_d)) return
  n = 3*ny*nz*5
  allocate(xs_lo_d(n), xs_hi_d(n), xr_lo_d(n), xr_hi_d(n))
  allocate(xs_lo_h(n), xs_hi_h(n), xr_lo_h(n), xr_hi_h(n))
end subroutine alloc_exchange_x


!=====================================================================
! x方向の袖交換(非周期: 両端のスラブは外側に相手を持たない)。
! 右隣へ内部の最後の3面(nx-5..nx-3)、左隣へ最初の3面(4..6)を送り、
! ゴースト(1..3, nx-2..nx)に受け取る。y/zのゴーストも含む全(j,k)。
! GPU上で1次元バッファに詰め、ホスト経由でMPI_SENDRECVする。
!=====================================================================
subroutine exchange_x(nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
  use mpi
  integer, intent(in), value     :: nx, ny, nz
  real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  integer :: left, right, n, ierr
  integer :: istat(MPI_STATUS_SIZE)

  if (npx == 1) return
  call alloc_exchange_x(ny, nz)
  n = 3*ny*nz*5
  left  = MPI_PROC_NULL;  if (my_slab > 0)     left  = my_slab - 1
  right = MPI_PROC_NULL;  if (my_slab < npx-1) right = my_slab + 1

  if (right /= MPI_PROC_NULL) then
    call pack_x(nx, ny, nz, nx-5, Q_1, Q_2, Q_3, Q_4, Q_5, xs_hi_d)
    xs_hi_h = xs_hi_d
  endif
  if (left /= MPI_PROC_NULL) then
    call pack_x(nx, ny, nz, 4, Q_1, Q_2, Q_3, Q_4, Q_5, xs_lo_d)
    xs_lo_h = xs_lo_d
  endif
  call MPI_SENDRECV(xs_hi_h, n, MPI_REAL8, right, 41, xr_lo_h, n, MPI_REAL8, left,  41, comm_compute, istat, ierr)
  call MPI_SENDRECV(xs_lo_h, n, MPI_REAL8, left,  42, xr_hi_h, n, MPI_REAL8, right, 42, comm_compute, istat, ierr)
  if (left /= MPI_PROC_NULL) then
    xr_lo_d = xr_lo_h
    call unpack_x(nx, ny, nz, 1, 1, x_fac_lo, xr_lo_d, Q_1, Q_2, Q_3, Q_4, Q_5)
  endif
  if (right /= MPI_PROC_NULL) then
    xr_hi_d = xr_hi_h
    call unpack_x(nx, ny, nz, nx-2, nx, x_fac_hi, xr_hi_d, Q_1, Q_2, Q_3, Q_4, Q_5)
  endif
end subroutine exchange_x

!> x面 i0..i0+2 を1次元バッファへ詰める。
subroutine pack_x(nx, ny, nz, i0, Q_1, Q_2, Q_3, Q_4, Q_5, buf)
  integer, intent(in), value     :: nx, ny, nz, i0
  real(8), intent(in), device    :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  real(8), intent(inout), device :: buf(3*ny*nz*5)
  integer :: j, k, m, idx, n3
  n3 = 3*ny*nz
  !$cuf kernel do(2)<<<*,*>>>
  do k = 1, nz
    do j = 1, ny
      do m = 1, 3
        idx = m + 3*((j-1) + ny*(k-1))
        buf(idx)        = Q_1(i0+m-1,j,k)
        buf(idx+n3)     = Q_2(i0+m-1,j,k)
        buf(idx+2*n3)   = Q_3(i0+m-1,j,k)
        buf(idx+3*n3)   = Q_4(i0+m-1,j,k)
        buf(idx+4*n3)   = Q_5(i0+m-1,j,k)
  enddo;enddo;enddo
end subroutine pack_x

!> 1次元バッファをx面 i0..i0+2 へ展開する。配列端の面 i_edge だけ、
!> Q/J を真のJ/ローカルのJ の比 fac で換算する(モジュール先頭の説明参照)。
subroutine unpack_x(nx, ny, nz, i0, i_edge, fac, buf, Q_1, Q_2, Q_3, Q_4, Q_5)
  integer, intent(in), value     :: nx, ny, nz, i0, i_edge
  real(8), intent(in), value     :: fac
  real(8), intent(in), device    :: buf(3*ny*nz*5)
  real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  integer :: j, k, m, idx, n3
  real(8) :: f
  n3 = 3*ny*nz
  !$cuf kernel do(2)<<<*,*>>>
  do k = 1, nz
    do j = 1, ny
      do m = 1, 3
        idx = m + 3*((j-1) + ny*(k-1))
        f = 1.d0
        if (i0+m-1 == i_edge) f = fac
        Q_1(i0+m-1,j,k) = buf(idx)*f
        Q_2(i0+m-1,j,k) = buf(idx+n3)*f
        Q_3(i0+m-1,j,k) = buf(idx+2*n3)*f
        Q_4(i0+m-1,j,k) = buf(idx+3*n3)*f
        Q_5(i0+m-1,j,k) = buf(idx+4*n3)*f
  enddo;enddo;enddo
end subroutine unpack_x


!=====================================================================
! パッシブスカラー: 拡散係数 rho*D = mu/Sc (ソルバーと同じサザーランド則)
!=====================================================================
subroutine sc_calc_gam(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5, gam)
  use mod_globals, only : Sc
  use mod_constant, only : mu0_T0_S_over_T0_2_3
  integer, intent(in), value  :: myrank, nx, ny, nz
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  real(8), intent(inout), device :: gam(nx,ny,nz)
  integer :: i, j, k
  real(8) :: rho, u, v, w, pres, temp, over_Sc

  over_Sc = 1.d0/Sc
  !$cuf kernel do(3)<<<*,*>>>
  do k = 4, nz-3
    do j = 1, ny
      do i = 1, nx
        rho  = Q_1(i,j,k)*jacobian(i,j)
        u    = Q_2(i,j,k)/Q_1(i,j,k)
        v    = Q_3(i,j,k)/Q_1(i,j,k)
        w    = Q_4(i,j,k)/Q_1(i,j,k)
        pres = (gamma-1.d0)*(Q_5(i,j,k)*jacobian(i,j) - 0.5d0*rho*(u*u + v*v + w*w))
        temp = pres/(R*rho)
        gam(i,j,k) = mu0_T0_S_over_T0_2_3/(temp + 110.4d0)*(temp*sqrt(temp))*over_Sc
  enddo;enddo;enddo

  call sc_exchange_z(myrank, nx, ny, nz, gam)
end subroutine sc_calc_gam


!=====================================================================
! パッシブスカラー: 1セルの右辺  d(xi)/dt = -u.grad(xi) + (1/rho) div(rho*D grad(xi))
! 移流はKorenリミッタ付き3次風上(端は添字クランプで1次風上に落ちる)、拡散は2次中心。
! 近傍のxiは最初にまとめて読み込む(互いに独立な読み出しをまとめて発行してレイテンシを重ねる)。
! 読み込んだ値は界面値と拡散項の計算で使い回す。全セルを計算する(省略はしない)。
!=====================================================================
attributes(device) function sc_rhs_cell(nx, ny, nz, i, j, k, mz, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam) result(r)
  integer, intent(in), value  :: nx, ny, nz, i, j, k
  real(8), intent(in), value  :: mz
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
  real(8), intent(in), device :: f(nx,ny,nz), gam(nx,ny,nz)
  real(8) :: r
  integer :: im2, ip2, jm2, jp2
  real(8) :: fc, fim1, fip1, fim2, fip2, fjm1, fjp1, fjm2, fjp2, fkm1, fkp1, fkm2, fkp2
  real(8) :: rho, u, v, w, adv, dif, d, dz_adv

  im2 = max(i-2, 1);  ip2 = min(i+2, nx)
  jm2 = max(j-2, 1);  jp2 = min(j+2, ny)

  fc   = f(i,j,k)
  fim1 = f(i-1,j,k);  fip1 = f(i+1,j,k);  fim2 = f(im2,j,k);  fip2 = f(ip2,j,k)
  fjm1 = f(i,j-1,k);  fjp1 = f(i,j+1,k);  fjm2 = f(i,jm2,k);  fjp2 = f(i,jp2,k)
  fkm1 = f(i,j,k-1);  fkp1 = f(i,j,k+1);  fkm2 = f(i,j,k-2);  fkp2 = f(i,j,k+2)

  rho = Q_1(i,j,k)*jacobian(i,j)
  u   = Q_2(i,j,k)/Q_1(i,j,k)
  v   = Q_3(i,j,k)/Q_1(i,j,k)
  w   = Q_4(i,j,k)/Q_1(i,j,k)

  ! --- x ---
  if (u >= 0.d0) then
    d = sc_face(fim1, fc, fip1) - sc_face(fim2, fim1, fc)
  else
    d = sc_face(fip2, fip1, fc) - sc_face(fip1, fc, fim1)
  endif
  adv = u*d*sc_mx(i)
  dif = sc_mx(i)*( 0.5d0*(gam(i,j,k)+gam(i+1,j,k))*(fip1-fc)*sc_idx(i) &
                 - 0.5d0*(gam(i,j,k)+gam(i-1,j,k))*(fc-fim1)*sc_idx(i-1) )

  ! --- y ---
  if (v >= 0.d0) then
    d = sc_face(fjm1, fc, fjp1) - sc_face(fjm2, fjm1, fc)
  else
    d = sc_face(fjp2, fjp1, fc) - sc_face(fjp1, fc, fjm1)
  endif
  adv = adv + v*d*sc_my(j)
  dif = dif + sc_my(j)*( 0.5d0*(gam(i,j,k)+gam(i,j+1,k))*(fjp1-fc)*sc_idy(j) &
                       - 0.5d0*(gam(i,j,k)+gam(i,j-1,k))*(fc-fjm1)*sc_idy(j-1) )

  ! --- z の風上差分(面値の差) ---
  if (w >= 0.d0) then
    dz_adv = sc_face(fkm1, fc, fkp1) - sc_face(fkm2, fkm1, fc)
  else
    dz_adv = sc_face(fkp2, fkp1, fc) - sc_face(fkp1, fc, fkm1)
  endif
  ! --- z (一様格子、ゴースト3面) ---
  adv = adv + w*dz_adv*mz
  dif = dif + mz*mz*( 0.5d0*(gam(i,j,k)+gam(i,j,k+1))*(fkp1-fc) &
                    - 0.5d0*(gam(i,j,k)+gam(i,j,k-1))*(fc-fkm1) )

  r = -adv + dif/rho
end function sc_rhs_cell


!=====================================================================
! パッシブスカラー: SSP-RK3の1段 (右辺の評価と更新を1つのカーネルにまとめたもの)
!   sc_stage : out = ca*phi0 + cb*( f + dt*L(f) )     (段1: ca=0,cb=1  段2: ca=3/4,cb=1/4)
!   sc_stage3: phi = ca*phi  + cb*( f + dt*L(f) )     (段3: ca=1/3,cb=2/3, phiをその場で更新)
! outとfは別の配列(隣のセルを読むので同じ配列には書かない)。段3は読むのがfの近傍とphi自身の
! 同じ点だけなので、phiにその場で書いても競合しない。
!=====================================================================
subroutine sc_stage(nx, ny, nz, ca, cb, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam, phi0, out)
  integer, intent(in), value  :: nx, ny, nz
  real(8), intent(in), value  :: ca, cb
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
  real(8), intent(in), device :: f(nx,ny,nz), gam(nx,ny,nz), phi0(nx,ny,nz)
  real(8), intent(inout), device :: out(nx,ny,nz)
  integer :: i, j, k
  real(8) :: mz
  mz = sc_mz

  !$cuf kernel do(3)<<<*,(32,4,1)>>>
  do k = 4, nz-3
    do j = 2, ny-1
      do i = 2, nx-1
        out(i,j,k) = ca*phi0(i,j,k) + cb*( f(i,j,k) + dt*sc_rhs_cell(nx, ny, nz, i, j, k, mz, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam) )
  enddo;enddo;enddo
end subroutine sc_stage


subroutine sc_stage3(nx, ny, nz, ca, cb, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam, phi)
  integer, intent(in), value  :: nx, ny, nz
  real(8), intent(in), value  :: ca, cb
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
  real(8), intent(in), device :: f(nx,ny,nz), gam(nx,ny,nz)
  real(8), intent(inout), device :: phi(nx,ny,nz)
  integer :: i, j, k
  real(8) :: mz
  mz = sc_mz

  !$cuf kernel do(3)<<<*,(32,4,1)>>>
  do k = 4, nz-3
    do j = 2, ny-1
      do i = 2, nx-1
        phi(i,j,k) = ca*phi(i,j,k) + cb*( f(i,j,k) + dt*sc_rhs_cell(nx, ny, nz, i, j, k, mz, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam) )
  enddo;enddo;enddo
end subroutine sc_stage3


!=====================================================================
! パッシブスカラー: 流入分布を全域に与える(初期値)
!=====================================================================
subroutine sc_fill_inlet_profile(nx, ny, nz, f)
  integer, intent(in), value     :: nx, ny, nz
  real(8), intent(inout), device :: f(nx,ny,nz)
  integer :: i, j, k

  !$cuf kernel do(3)<<<*,*>>>
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        f(i,j,k) = sc_in(j)
  enddo;enddo;enddo
end subroutine sc_fill_inlet_profile


!=====================================================================
! パッシブスカラー: 1ステップ進める(SSP-RK3、速度場は更新後の値で固定)
!   段1: phis = phi + dt L(phi)
!   段2: phib = 3/4 phi + 1/4 ( phis + dt L(phis) )
!   段3: phi  = 1/3 phi + 2/3 ( phib + dt L(phib) )
!=====================================================================
subroutine sc_advance(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)
  integer, intent(in), value  :: myrank, nx, ny, nz
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)

  call sc_calc_gam(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5, sc_gam)

  call sc_stage(nx, ny, nz, 0.d0, 1.d0, jacobian, Q_1, Q_2, Q_3, Q_4, sc_phi, sc_gam, sc_phi, sc_phis)
  call sc_bc(myrank, nx, ny, nz, sc_phis)

  call sc_stage(nx, ny, nz, 0.75d0, 0.25d0, jacobian, Q_1, Q_2, Q_3, Q_4, sc_phis, sc_gam, sc_phi, sc_phib)
  call sc_bc(myrank, nx, ny, nz, sc_phib)

  call sc_stage3(nx, ny, nz, 1.d0/3.d0, 2.d0/3.d0, jacobian, Q_1, Q_2, Q_3, Q_4, sc_phib, sc_gam, sc_phi)
  call sc_bc(myrank, nx, ny, nz, sc_phi)
end subroutine sc_advance


!=====================================================================
! パッシブスカラー: 出力。data/(<I/Oランク>/)xiNNNNN.vtr
! 流れ場のQNNNNN.vtrと同じ格子・同じ番号のVTK(RectilinearGrid、単精度、点データ"xi")。
! ParaViewで直接開ける。zゴースト面も流れ場と同じく含む。
!=====================================================================
subroutine sc_write(myrank, nx, ny, nz, idx)
  use mod_globals, only : step_offset
  integer, intent(in), value :: myrank, nx, ny, nz, idx
  real(8), allocatable :: h(:,:,:)
  real(4), allocatable :: h4(:,:,:)
  character(len=64) :: filename
  character(len=12) :: off1, off2, off3
  character(len=6)  :: e1, e2, e3
  character(len=1)  :: lf
  integer    :: u
  integer(4) :: nb_x, nb_y, nb_z, nb_f
  real(8)    :: fmin, fmax
  logical    :: out_of_range

  allocate(h(nx,ny,nz), h4(nx,ny,nz))
  h  = sc_phi
  h4 = real(h, 4)
  ! 値域の検査(NaNも「範囲外」になる)。正常なときは何も表示しない。
  out_of_range = .not. all(h(2:nx-1,2:ny-1,4:nz-3) >= -1.d-6 .and. h(2:nx-1,2:ny-1,4:nz-3) <= 1.d0+1.d-6)
  fmin = minval(h(2:nx-1,2:ny-1,4:nz-3))
  fmax = maxval(h(2:nx-1,2:ny-1,4:nz-3))
  if (npx >= 2) then
    write(filename, "(a, i0, a, i5.5, a)") "data/", myrank+1, "/xi", idx+step_offset, ".vtr"
  else
    write(filename, "(a, i5.5, a)") "data/xi", idx+step_offset, ".vtr"
  endif

  lf = char(10)
  nb_x = 4*nx;  nb_y = 4*ny;  nb_z = 4*nz;  nb_f = 4*nx*ny*nz   ! 各ブロックのデータ長[byte]
  write(e1, '(i6)') nx-1;  write(e2, '(i6)') ny-1;  write(e3, '(i6)') nz-1
  write(off1, '(i12)') 4_8 + nb_x
  write(off2, '(i12)') 8_8 + nb_x + nb_y
  write(off3, '(i12)') 12_8 + nb_x + nb_y + nb_z
  open(newunit=u, file=filename, status="replace", action="write", form="unformatted", access="stream", convert="little_endian")
  write(u) '<?xml version="1.0"?>'//lf
  write(u) '<VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian" header_type="UInt32">'//lf
  write(u) '  <RectilinearGrid WholeExtent="0 '//e1//' 0 '//e2//' 0 '//e3//'">'//lf
  write(u) '    <Piece Extent="0 '//e1//' 0 '//e2//' 0 '//e3//'">'//lf
  write(u) '      <Coordinates>'//lf
  write(u) '        <DataArray type="Float32" Name="x" format="appended" offset="0"/>'//lf
  write(u) '        <DataArray type="Float32" Name="y" format="appended" offset="'//off1//'"/>'//lf
  write(u) '        <DataArray type="Float32" Name="z" format="appended" offset="'//off2//'"/>'//lf
  write(u) '      </Coordinates>'//lf
  write(u) '      <PointData Scalars="xi">'//lf
  write(u) '        <DataArray type="Float32" Name="xi" NumberOfComponents="1" format="appended" offset="'//off3//'"/>'//lf
  write(u) '      </PointData>'//lf
  write(u) '    </Piece>'//lf
  write(u) '  </RectilinearGrid>'//lf
  write(u) '  <AppendedData encoding="raw">'//lf
  write(u) '  _', nb_x, sc_xc, nb_y, sc_yc, nb_z, sc_zc, nb_f, h4, lf
  write(u) '  </AppendedData>'//lf
  write(u) '</VTKFile>'//lf
  close(u)
  deallocate(h, h4)
  if (out_of_range) &
    print '(1x,a,i0,a,a,a,es13.6,a,es13.6)', "myrank is ", myrank, "  WARNING scalar out of [0,1]: ", trim(filename), &
          "  min=", fmin, "  max=", fmax
end subroutine sc_write


!=====================================================================
! パッシブスカラー: リスタートファイル recal/xiNNNNN.dat (スラブごと、real(8)、ゴースト込み)。
! 流れ場のrecal/QNNNNN.datと同じく計算の最後に1回だけ書く。
!=====================================================================
subroutine sc_write_restart(myrank, nx, ny, nz)
  integer, intent(in), value :: myrank, nx, ny, nz
  real(8), allocatable :: h(:,:,:)
  character(len=64) :: filename
  integer :: u

  allocate(h(nx,ny,nz))
  h = sc_phi
  write(filename, "(a, i5.5, a)") "recal/xi", myrank/2+1, ".dat"
  open(newunit=u, file=filename, status="replace", action="write", form="unformatted", access="stream")
  write(u) h
  close(u)
  deallocate(h)
  print *, "myrank is ", myrank, " scalar restart file written: ", trim(filename)
end subroutine sc_write_restart


!=====================================================================
! パッシブスカラー: 初期化。RESTART=Trueでrecal/xiNNNNN.datがあれば読み、
! なければ流入分布(tanh)を全域に与える。
!=====================================================================
subroutine sc_init(myrank, nx, ny, nz)
  use mod_constant, only : id_recal
  integer, intent(in), value :: myrank, nx, ny, nz
  real(8), allocatable :: h(:,:,:)
  character(len=64) :: filename
  integer :: u, ios, ierr
  logical :: have_file

  allocate(sc_phi(nx,ny,nz), sc_phis(nx,ny,nz), sc_phib(nx,ny,nz), sc_gam(nx,ny,nz), stat=ierr)
  if (ierr /= 0) then
    print *, "myrank is ", myrank, " scalar: GPU memory allocation failed", ierr
    error stop 'scalar allocation failed'
  endif
  if (npx == 1) call execute_command_line("mkdir -p data", wait=.true., exitstat=ierr)
  sc_phib = 0.d0
  sc_gam = 0.d0

  have_file = .false.
  if (id_recal) then
    write(filename, "(a, i5.5, a)") "recal/xi", myrank/2+1, ".dat"
    inquire(file=filename, exist=have_file)
  endif
  if (have_file) then
    allocate(h(nx,ny,nz))
    open(newunit=u, file=filename, status="old", action="read", form="unformatted", access="stream")
    read(u, iostat=ios) h
    close(u)
    if (ios /= 0) error stop 'scalar restart file has the wrong size'
    sc_phi = h
    deallocate(h)
    print *, "myrank is ", myrank, " scalar restarted from ", trim(filename)
  else
    if (id_recal) print *, "myrank is ", myrank, " scalar: no restart file, starting from the inlet profile"
    call sc_fill_inlet_profile(nx, ny, nz, sc_phi)
  endif
  call sc_bc(myrank, nx, ny, nz, sc_phi)
  sc_phis = sc_phi
  call sc_write(myrank, nx, ny, nz, 0)
end subroutine sc_init


!=====================================================================
! パッシブスカラー: set_bcの最後から呼ばれる入口。最初の呼び出しで初期化し、
! 各ステップの最終RK段(流れ場がQ^(n+1)になった時点)でxiを1ステップ進める。
!=====================================================================
subroutine sc_step(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)
  use mod_globals, only : np, inlet_rk_stages, scalar_output_every
  integer, intent(in), value  :: myrank, nx, ny, nz
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  integer :: n

  if (.not. allocated(sc_phi)) call sc_init(myrank, nx, ny, nz)
  if (mod(inlet_bc_calls, inlet_rk_stages) /= 0) return
  n = inlet_bc_calls/inlet_rk_stages                 ! 完了したステップ数
  call sc_advance(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)
  if (mod(n, nt*scalar_output_every) == 0) call sc_write(myrank, nx, ny, nz, n/nt)
  if (n == np*nt) call sc_write_restart(myrank, nx, ny, nz)
end subroutine sc_step


!=====================================================================
! このランクのGPUを選ぶ。main.f90のmygpuと同じ式(ノード内ランク/2 を GPU数で割った余り)。
! モジュール変数のdevice配列(over_jacobian_tab, sigma_*, sc_*など)は GPU ごとに別のコピーが
! あり、書き込んだ時点で選ばれている GPU のコピーに入る。本来の選択(check_gpu)は set_grid より
! 後のRungeKutta内なので、複数GPUでは GPU 0 に表を書いて別のGPUで読むことになる。
! そのため、表に書き込む前(set_grid の最初)に同じ GPU を選んでおく。GPU 1枚なら常に 0。
!=====================================================================
subroutine select_my_gpu()
  use mpi
  use cudafor
  integer :: comm_node, rank_node, ndev, ierr, stat

  call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_node, ierr)
  call MPI_COMM_RANK(comm_node, rank_node, ierr)
  call MPI_COMM_FREE(comm_node, ierr)
  stat = cudaGetDeviceCount(ndev)
  if (stat == 0 .and. ndev >= 1) stat = cudaSetDevice(mod(rank_node/2, ndev))
end subroutine select_my_gpu


!=====================================================================
! 診断(最初の数回のset_bcだけ): CUDAエラー(カーネル起動失敗など)の報告、このランクの
! GPU(型番とcc)、device表の中身(空でないか)を表示する。通常時は1行(GPUの情報)だけ出る。
!=====================================================================
subroutine diag_first_calls(myrank, nx, ny)
  use cudafor
  integer, intent(in), value :: myrank, nx, ny
  integer :: stat, dev
  type(cudaDeviceProp) :: prop
  real(8) :: t

  if (inlet_bc_calls > 7) return
  stat = cudaDeviceSynchronize()
  if (stat /= 0) print '(1x,a,i0,a,i0,a,a)', "myrank is ", myrank, "  CUDA error before set_bc call ", &
                       inlet_bc_calls+1, ": ", trim(cudaGetErrorString(stat))
  stat = cudaGetLastError()
  if (stat /= 0) print '(1x,a,i0,a,i0,a,a)', "myrank is ", myrank, "  CUDA error before set_bc call ", &
                       inlet_bc_calls+1, ": ", trim(cudaGetErrorString(stat))
  if (inlet_bc_calls == 0) then
    stat = cudaGetDevice(dev)
    stat = cudaGetDeviceProperties(prop, dev)
    t = over_jacobian_tab(nx/2, ny/2)
    print '(1x,a,i0,a,i0,a,a,a,i0,a,i0,a,es10.3)', "myrank is ", myrank, "  set_bc runs on GPU ", dev, " (", &
          trim(prop%name), ", cc ", prop%major, ".", prop%minor, ")  over_jacobian_tab sample =", t
    if (.not. (t > 0.d0 .and. t < huge(1.d0))) print '(1x,a,i0,a)', "myrank is ", myrank, &
          "  ERROR: the module device tables are empty on this GPU (they were written on another device)"
  endif
end subroutine diag_first_calls


!=====================================================================
! 配列中のNaN/Infの個数と、最初の(最小の線形番号の)位置を数える。
!=====================================================================
subroutine nan_scan(nx, ny, nz, Q, cnt, first)
  integer, intent(in), value  :: nx, ny, nz
  real(8), intent(in), device :: Q(nx,ny,nz)
  integer, intent(out)        :: cnt, first
  integer :: i, j, k, c, f, lin

  c = 0
  f = huge(1)
  !$cuf kernel do(3)<<<*,*>>>
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if (.not. (abs(Q(i,j,k)) <= 1.d300)) then
          c = c + 1
          lin = i + nx*((j-1) + ny*(k-1))
          f = min(f, lin)
        endif
  enddo;enddo;enddo
  cnt = c
  first = f
end subroutine nan_scan


!=====================================================================
! NaN/Infの検査。見つけたら場所を表示して全ランクを強制終了する(MPI_ABORT)。
! 最初の8回のset_bc(= 最初の2ステップ)は毎回、その後は nan_check_every ステップごとの
! 4段目(Q^(n+1)になった時点)で、流れ場の保存変数5つと濃度xiを検査する。
! nan_check_every = 0 なら検査しない。set_bcの入口で呼ぶので、直前のRK段の結果を見る。
!=====================================================================
subroutine check_nan_abort(myrank, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
  use mpi
  use mod_globals, only : nan_check_every, inlet_rk_stages, dt
  integer, intent(in), value  :: myrank, nx, ny, nz
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  integer :: m, cnt, first, ii, jj, kk, ierr, stage, step
  logical :: bad
  character(len=8) :: name

  if (nan_check_every <= 0) return
  stage = mod(inlet_bc_calls, inlet_rk_stages) + 1          ! この呼び出しが見るRK段の結果(1..4)
  step  = inlet_bc_calls/inlet_rk_stages + 1                ! その段が属するステップ番号
  if (inlet_bc_calls > 7) then
    if (stage /= inlet_rk_stages) return
    if (mod(step, max(nan_check_every, 1)) /= 0) return
  endif

  bad = .false.
  do m = 1, 6
    select case (m)
    case (1); name = "Q_1";  call nan_scan(nx, ny, nz, Q_1, cnt, first)
    case (2); name = "Q_2";  call nan_scan(nx, ny, nz, Q_2, cnt, first)
    case (3); name = "Q_3";  call nan_scan(nx, ny, nz, Q_3, cnt, first)
    case (4); name = "Q_4";  call nan_scan(nx, ny, nz, Q_4, cnt, first)
    case (5); name = "Q_5";  call nan_scan(nx, ny, nz, Q_5, cnt, first)
    case (6)
      name = "xi"
      cnt = 0;  first = huge(1)
      if (allocated(sc_phi)) call nan_scan(nx, ny, nz, sc_phi, cnt, first)
    end select
    if (cnt > 0) then
      bad = .true.
      ii = mod(first-1, nx) + 1
      jj = mod((first-1)/nx, ny) + 1
      kk = (first-1)/(nx*ny) + 1
      print '(1x,a,i0,a,i0,a,i0,a,i0,a,es10.3,a)', "myrank is ", myrank, "  NaN/Inf detected: step ", step, &
            ", RK stage ", stage, " (set_bc call ", inlet_bc_calls+1, ", t = ", dble(step-1)*dt, " s)"
      print '(1x,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)', "   array ", trim(name), ": ", cnt, " cells; first at local (i,j,k) = ", &
            ii, ",", jj, ",", kk, "   (global i = ", ii + i_offset, ")"
    endif
  enddo
  if (bad) then
    print '(1x,a,i0,a)', "myrank is ", myrank, "  aborting all MPI ranks because of NaN/Inf"
    call flush(6)
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    stop 1
  endif
end subroutine check_nan_abort


subroutine set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
  use mpi
  use mod_constant, only : id_accuracy
  use mod_globals, only : dt, CFL, u1, u2, gamma, R, T1, T2, endT, np, Pr
  integer, intent(in)  :: myrank, nx, ny, nz
  real(8), intent(in)  :: Lx, Ly, Lz
  real(8), intent(out) :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)

  real(8) :: dx_min, dy_min, dz_min, dmin, c1, c2, umax, cmax, dt_suggest
  real(8) :: vmax_est, wmax_est, dt_conv, dt_diff, mu_max, rho_min
  integer :: nt_suggest
  integer :: nranks, ierr
  character(len=64) :: rank_dir
  real(8) :: jac_h(nx,ny)

  ! 表(device配列)に書き込む前に、このランクのGPUを選ぶ(select_my_gpuの説明を参照)。
  call select_my_gpu()

  ! set_gridは全ランクから呼ばれるので、ここでx分割の情報とコミュニケータを用意する。
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  if (nranks /= 2*npx) then
    print *, "x decomposition: npx =", npx, " needs", 2*npx, " MPI ranks, but got", nranks
    error stop 'MPI rank count must be 2*npx'
  endif
  if (npx < 1 .or. mod(nx_global-6, npx) /= 0) then
    print *, "x decomposition: nx_global-6 =", nx_global-6, " is not divisible by npx =", npx
    error stop '(nx_main+nx_buf-4) must be divisible by npx'
  endif
  my_slab  = myrank/2
  i_offset = my_slab*(nx-6)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, mod(myrank,2), myrank, comm_compute, ierr)
  if (npx > 1) then
    ! 流出スポンジとその目標面(i = nx-1-nsp_x)は最終スラブの内部に収まること。
    if (nx-1-nsp_x < 4) error stop 'x decomposition: the outflow sponge must fit inside the last slab (reduce npx or nx_buf)'
    ! 奇数(I/O)ランクは data/<myrank>/ にスラブごとの出力を書く。
    ! ディレクトリが無いと最初のopenで落ちるので、ここで作っておく。
    if (mod(myrank,2) == 1) then
      write(rank_dir, "(a, i0)") "mkdir -p data/", myrank
      call execute_command_line(trim(rank_dir), wait=.true., exitstat=ierr)
    endif
  endif

  call set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
                            nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)

  ! set_bcが使う 1/Jacobian テーブル。main.f90と同じset_Jacobian_xy3で作る
  ! (set_Jacobian_xy3_stretchはこのソルバーのmainからは呼ばれない)。
  call set_Jacobian_xy3(nx, ny, nz, dx, dy, dz, jac_h)
  over_jacobian_tab = 1.0d0 / jac_h

  ! パッシブスカラー用のメトリック(セル中心間隔から)
  block
    real(8) :: mx_h(nx), my_h(ny), idx_h(nx-1), idy_h(ny-1)
    integer :: ii
    do ii = 2, nx-1; mx_h(ii) = 2.d0/(dx(ii-1) + dx(ii)); enddo
    mx_h(1) = mx_h(2);  mx_h(nx) = mx_h(nx-1)
    do ii = 2, ny-1; my_h(ii) = 2.d0/(dy(ii-1) + dy(ii)); enddo
    my_h(1) = my_h(2);  my_h(ny) = my_h(ny-1)
    idx_h = 1.d0/dx
    idy_h = 1.d0/dy
    sc_mx  = mx_h;   sc_my  = my_h
    sc_idx = idx_h;  sc_idy = idy_h
    sc_mz  = dble(n_compute*(nz-6))/Lz
    sc_xc  = real(xc, 4);  sc_yc = real(yc, 4);  sc_zc = real(zc, 4)
  end block

  dx_min = minval(dx)
  dy_min = minval(dy)
  dz_min = Lz / dble(n_compute*(nz-6))

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
      total = total + dr_local     ! 先に加算(現行と同じ順序: 代入してからdrを伸ばす)
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

  use mod_globals, only : Ly_uniform, ny_uniform
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: Lx_main, Lx_buf, Ly_main, Ly_buf, Lz
  integer, intent(in) :: nx_main, nx_buf, ny_main, ny_buf
  real(8), intent(out) :: xc(nx), yc(ny), zc(nz)
  real(8), intent(out) :: dx(nx-1), dy(ny-1), dz(nz-1)

  ! x は全体格子(nxg点)で生成し、このスラブの範囲を i_offset で切り出す。
  integer, parameter :: nxg = nx_global
  real(8) :: x(nxg+1), xcg(nxg), dxg(nxg-1), y(ny+1), z(nz+1)
  real(8) :: dz1
  real(8) :: s, dr, eta, pi
  integer :: i, j, k, ig

  ! 格子伸長パラメータ
  real(8)  :: r0_x, r1_x               ! x方向バッファ: 開始/終端伸び率
  real(8)  :: r0_y_lo, r1_y_lo         ! y方向下バッファ: 開始/終端伸び率
  real(8)  :: r0_y_hi, r1_y_hi         ! y方向上バッファ: 開始/終端伸び率
  real(8)  :: r_local

  pi = acos(-1.d0)

  ! x main domain is uniform; only the separate downstream buffer stretches.
  if (nx_main < 3 .or. Lx_main <= 0.d0) error stop 'Invalid uniform x mesh'
  if (i_offset + nx > nxg) error stop 'x slab exceeds the global grid'
  do i = 2, nx_main+1
    x(i) = Lx_main*dble(i-2)/dble(nx_main-1)
  enddo

  ! 接続点の格子幅とその直前の局所伸び率(主計算領域側から引き継ぐ)
  dr  = x(nx_main+1) - x(nx_main)
  r0_x = 1.d0

  r1_x = calc_smooth_stretch_ratio(dr, r0_x, nx_buf, Lx_buf)

  do i = nx_main+2, nxg-1
    eta = dble(i - (nx_main+2)) / dble(max(nx_buf-1,1))
    r_local = r0_x + (r1_x - r0_x) * (1.0d0 - cos(pi*eta)) / 2.0d0
    x(i) = x(i-1) + dr
    dr = dr * r_local
  enddo

  ! 左ゴースト(i=1): 主計算領域流入端の格子幅を複製
  x(1) = x(2) - (x(3)-x(2))
  ! 右ゴースト
  x(nxg)   = x(nxg-1) + dr
  x(nxg+1) = x(nxg)   + dr

  print *, "x(nx-1) - x(nx_main+1) =", x(nxg-1) - x(nx_main+1), " (目標:", Lx_buf, ")"
  print *, "xバッファ 開始伸び率=", r0_x, " 終端伸び率=", r1_x

  !=================================================================
  ! y方向: 下バッファ(滑らかランプ) + 主計算領域(中心一様・外側伸長) + 上バッファ(滑らかランプ)
  !=================================================================
  ! y main domain: symmetric stretching outside the central uniform section.
  call partitioned_axis(ny_main, Ly_main, Ly_uniform, ny_uniform, .true., &
                        y(ny_buf+2:ny_buf+ny_main+1))
  y(ny_buf+2:ny_buf+ny_main+1) = y(ny_buf+2:ny_buf+ny_main+1)+Ly_buf

  ! --- 下バッファ: 接続点の格子幅と、その外側への局所伸び率を引き継ぐ ---
  dr = y(ny_buf+3) - y(ny_buf+2)
  r0_y_lo = dr / (y(ny_buf+4) - y(ny_buf+3))

  r1_y_lo = calc_smooth_stretch_ratio(dr, r0_y_lo, ny_buf, Ly_buf)

  print *, "下バッファ 開始伸び率=", r0_y_lo, " 終端伸び率=", r1_y_lo
  print *, "下バッファ合計長さ確認 ="

  do j = ny_buf+1, 2, -1
    eta = dble((ny_buf+2-1) - j) / dble(max(ny_buf-1,1))
    r_local = r0_y_lo + (r1_y_lo - r0_y_lo) * (1.0d0 - cos(pi*eta)) / 2.0d0
    y(j) = y(j+1) - dr
    dr = dr * r_local
  enddo

  print *, "y(ny_buf+2) - y(2) =", y(ny_buf+2) - y(2), " (目標:", Ly_buf, ")"

  ! --- 上バッファ ---
  dr = y(ny_buf+ny_main+1) - y(ny_buf+ny_main)
  r0_y_hi = dr / (y(ny_buf+ny_main) - y(ny_buf+ny_main-1))

  r1_y_hi = calc_smooth_stretch_ratio(dr, r0_y_hi, ny_buf, Ly_buf)

  print *, "上バッファ 開始伸び率=", r0_y_hi, " 終端伸び率=", r1_y_hi

  do j = ny_buf+ny_main+2, ny-1
    eta = dble(j - (ny_buf+ny_main+2)) / dble(max(ny_buf-1,1))
    r_local = r0_y_hi + (r1_y_hi - r0_y_hi) * (1.0d0 - cos(pi*eta)) / 2.0d0
    y(j) = y(j-1) + dr
    dr = dr * r_local
  enddo

  print *, "y(ny-1) - y(ny_buf+ny_main+1) =", y(ny-1)-y(ny_buf+ny_main+1), &
           " (目標:", Ly_buf, ")"

  ! 下ゴースト(j=1)
  y(1) = y(2) - (y(3)-y(2))
  ! 上ゴースト
  y(ny)   = y(ny-1) + dr
  y(ny+1) = y(ny)   + dr

  !=================================================================
  ! z方向: 一様(周期境界、3点ゴースト)
  ! Lzは周期スパン、nz-6が全体の点数(zは分割しない: n_compute=1, k_offset=0)。
  !=================================================================
  dz1 = Lz / dble(n_compute*(nz-6))
  do k = 4, nz-2
    z(k) = dz1 * dble(k_offset+k-4)
  enddo
  z(3) = z(4) - dz1
  z(2) = z(3) - dz1
  z(1) = z(2) - dz1
  z(nz-1) = z(nz-2) + dz1
  z(nz)   = z(nz-1) + dz1
  z(nz+1) = z(nz)   + dz1

  !共通変数に代入 CPU --- GPU
  xi = x(i_offset+1:i_offset+nx)
  yj = y(1:ny)
  zk = z(1:nz)

  !=================================================================
  ! セル中心・格子間隔
  !=================================================================
  do i = 1, nxg; xcg(i) = 0.5d0*(x(i)+x(i+1)); enddo
  do i = 1, nxg-1; dxg(i) = xcg(i+1)-xcg(i); enddo
  do j = 1, ny; yc(j) = 0.5d0*(y(j)+y(j+1)); enddo
  do k = 1, nz; zc(k) = 0.5d0*(z(k)+z(k+1)); enddo

  ! このスラブの切り出し(値は全体格子と同一)
  do i = 1, nx;   xc(i) = xcg(i+i_offset); enddo
  do i = 1, nx-1; dx(i) = dxg(i+i_offset); enddo
  ! 接続側の配列端(local 1, nx)で、真のJ / ローカルのJ(隣の列のコピー)。
  x_fac_lo = 1.d0;  x_fac_hi = 1.d0
  if (my_slab > 0) then
    ig = i_offset + 1
    x_fac_lo = (dxg(ig) + dxg(ig+1)) / (dxg(ig-1) + dxg(ig))
  endif
  if (my_slab < npx-1) then
    ig = i_offset + nx
    x_fac_hi = (dxg(ig-2) + dxg(ig-1)) / (dxg(ig-1) + dxg(ig))
  endif
  if (npx > 1) print *, "x slab", my_slab, " of", npx, ": global i =", i_offset+1, "~", i_offset+nx, &
                        " x =", xc(1), "~", xc(nx)
  do j = 1, ny-1; dy(j) = yc(j+1)-yc(j); enddo
  do k = 1, nz-1; dz(k) = zc(k+1)-zc(k); enddo
!=================================================================
! スポンジ・流入境界用の目標値を、事前に1回だけ計算してテーブル化
! (これらはすべて i, k, 時間に依存しない量なので、set_bc内で毎回
!  計算し直す必要がない。格子は時間不変(RESCALE=False)なので安全)
!=================================================================
block
  use mod_globals, only : u1, u2, T1, T2, p, delta_bl, amp
  real(8) :: rho_target_1d_h(ny), u_target_1d_h(ny), amp_u_env_1d_h(ny)
  real(8) :: inlet_envelope_1d_h(ny)
  real(8) :: sigma_x_1d_h(nx), sigma_y_1d_h(ny)
  real(8) :: T_init_local, env_local, over_delta_bl_local, over_Lz_local
  real(8) :: over_nsp_x_local, over_nsp_y_local, eta_local, ramp_local
  real(8) :: sigma_y_lo_local, sigma_y_hi_local
  integer :: jj, kk, ii

  over_delta_bl_local = 1.0d0/delta_bl
  over_Lz_local = 1.0d0/Lz
  over_nsp_x_local = 1.0d0/dble(nsp_x)
  over_nsp_y_local = 1.0d0/dble(nsp_y)

  ! --- y方向テーブル: rho_target, u_target, amp*u_init*env, sigma_y ---
  do jj = 1, ny
    T_init_local = 0.5d0*(T1+T2) + 0.5d0*(T1-T2)*tanh(2.0d0*(y(jj)-0.5d0*Ly)/delta_bl)
    rho_target_1d_h(jj) = p / (R * T_init_local)
    u_target_1d_h(jj)   = 0.5d0*(u1+u2) + 0.5d0*(u1-u2)*tanh(2.0d0*(y(jj)-0.5d0*Ly)/delta_bl)

    ! Localize the 5% Gaussian inlet forcing to the region where the
    ! tanh mean profile has shear.  This is the normalized mean-velocity
    ! gradient: sech^2(2*(y-Ly/2)/delta_bl), evaluated stably as 1-tanh^2.
    env_local = tanh(2.0d0*(y(jj)-0.5d0*Ly)*over_delta_bl_local)
    inlet_envelope_1d_h(jj) = max(0.d0, 1.d0-env_local**2)

    env_local = exp(-((y(jj)-0.5d0*Ly)*over_delta_bl_local)**2)
    amp_u_env_1d_h(jj) = amp * u_target_1d_h(jj) * env_local

    sigma_y_lo_local = 0.d0
    sigma_y_hi_local = 0.d0
    if (jj <= 1+nsp_y) then
      eta_local = dble((1+nsp_y) - jj) * over_nsp_y_local
      if (eta_local < 0.d0) eta_local = 0.d0
      if (eta_local > 1.d0) eta_local = 1.d0
      ! 5次smoothstepを使い、終端の減衰率がsigma_max_yとなるよう
      ! (1-sigma)=(1-sigma_max_y)**ramp から逆算する。
      ramp_local = eta_local**3 * (10.d0 - 15.d0*eta_local + 6.d0*eta_local**2)
      sigma_y_lo_local = 1.d0 - (1.d0-sigma_max_y)**ramp_local
    endif
    if (jj >= ny-nsp_y) then
      eta_local = dble(jj - (ny-nsp_y)) * over_nsp_y_local
      if (eta_local < 0.d0) eta_local = 0.d0
      if (eta_local > 1.d0) eta_local = 1.d0
      ramp_local = eta_local**3 * (10.d0 - 15.d0*eta_local + 6.d0*eta_local**2)
      sigma_y_hi_local = 1.d0 - (1.d0-sigma_max_y)**ramp_local
    endif
    sigma_y_1d_h(jj) = max(sigma_y_lo_local, sigma_y_hi_local)
  enddo

  ! --- x方向テーブル: sigma_x ---
  do ii = 1, nx
    sigma_x_1d_h(ii) = 0.d0
    if (ii+i_offset >= nx_global-nsp_x) then      ! 全体番号で評価(最終スラブだけが非ゼロ)
      eta_local = dble(ii+i_offset - (nx_global-1-nsp_x)) * over_nsp_x_local
      if (eta_local < 0.d0) eta_local = 0.d0
      if (eta_local > 1.d0) eta_local = 1.d0
      ramp_local = eta_local**3 * (10.d0 - 15.d0*eta_local + 6.d0*eta_local**2)
      sigma_x_1d_h(ii) = 1.d0 - (1.d0-sigma_max_x)**ramp_local
    endif
  enddo

  ! CPU -> GPU 転送(いずれも1回のみ、サイズも小さく無視できるコスト)
  rho_target_1d = rho_target_1d_h
  u_target_1d   = u_target_1d_h
  inlet_envelope_1d = inlet_envelope_1d_h
  amp_u_env_1d  = amp_u_env_1d_h
  sigma_x_1d    = sigma_x_1d_h
  sigma_y_1d    = sigma_y_1d_h
  ! パッシブスカラーの流入分布: 速度分布と同じtanh。流れ1(u1側)で1、流れ2(u2側)で0。
  sc_in = (u_target_1d_h - u2)/(u1 - u2)
end block
  !=================================================================
  ! 確認出力
  !=================================================================
  print *, "=== 格子生成確認 ==="
  print *, "x: 主計算領域 ", x(2), "~", x(nx_main+1), " バッファ ~", x(nxg-1)
  print *, "y: バッファ下 ", y(2), "~", y(ny_buf+2)
  print *, "y: 主計算領域 ", y(ny_buf+2), "~", y(ny_buf+ny_main+1)
  print *, "y: バッファ上 ", y(ny_buf+ny_main+1), "~", y(ny-1)
  print *, "dy_center(最小格子幅) =", y(ny_buf+2+(ny_main/2)+1)-y(ny_buf+2+(ny_main/2))
  print *, "dy_edge  (主計算端)   =", y(ny_buf+3)-y(ny_buf+2)
  print *, "dy_buf_max(バッファ端) =", y(ny-1)-y(ny-2)
  print *, "dy_edge/dy_center     =", (y(ny_buf+3)-y(ny_buf+2)) / &
                                       (y(ny_buf+2+(ny_main/2)+1)-y(ny_buf+2+(ny_main/2)))
  print *, "dy_buf_max/dy_edge    =", (y(ny-1)-y(ny-2))/(y(ny_buf+3)-y(ny_buf+2))

end subroutine set_grid_main_buffer




  !=====================================================================
  ! ヤコビアン: 内部は中央差分、ゴーストセル(1点)は隣接内部セルの値を複製
  !=====================================================================
  subroutine set_Jacobian_xy3_stretch(nx, ny, nz, xc, yc, zc, Jacobian)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: xc(nx), yc(ny), zc(nz)
    real(8), intent(out) :: Jacobian(nx,ny)
    integer :: i, j
    real(8) :: dx(nx-1), dy(ny-1), dz1

    dx = xc(2:nx) - xc(1:nx-1)
    dy = yc(2:ny) - yc(1:ny-1)
    dz1 = zc(2) - zc(1)

    do j = 2, ny-1
      do i = 2, nx-1
        ! Same spacing-based expression as set_Jacobian_xy3.
        Jacobian(i,j) = 8.d0 / &
          ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * 2.d0*dz1)
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


  subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
    use mod_globals, only : u1, u2, p, Ly, T1, T2, delta_bl
    use mod_constant, only : id_accuracy
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(in)  :: x(nx), y(ny), z(nz)
    real(8), intent(out) :: Q(nx,ny,nz,5)   ! このソルバーの並びは (i,j,k,成分)
    integer :: i, j, k
    real(8) :: u_init, rho_init, T_init

    ! Start from the unperturbed mean mixing layer.  All velocity
    ! fluctuations enter through the x-inlet boundary in set_bc.
    do k = 1, nz
      do j = 1, ny
        T_init   = 0.5d0*(T1+T2) + 0.5d0*(T1-T2)*tanh(2.0d0*(y(j)-0.5d0*Ly)/delta_bl)
        rho_init = p / (R * T_init)
        u_init = 0.5*(u1+u2) + 0.5*(u1-u2)*tanh(2.0d0*(y(j)-0.5*Ly)/(delta_bl))

        do i = 1, nx
          Q(i,j,k,1) = rho_init
          Q(i,j,k,2) = rho_init * u_init
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0
          Q(i,j,k,5) = p/(gamma-1.d0) + 0.5d0*rho_init*u_init**2
        enddo;enddo;enddo

  end subroutine set_init



subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q_1, Q_2, Q_3, Q_4, Q_5 )
  use mod_constant, only : id_accuracy
  use mod_globals, only : u1, rho1, u2, rho2, p, amp, dt, gamma, Ly, Lz, Lx,T1,T2,delta_bl, step_offset, &
                          inlet_fluctuation_rms, inlet_random_seed, inlet_rk_stages, scalar_on
  use set_coordinate
  use mpi
  integer, intent(in), value            :: myrank, nx, ny, nz
  real(8), intent(in), device           :: jacobian(nx,ny)
  real(8), intent(inout), device        :: Q_1(nx,ny,nz), Q_2(nx,ny,nz),Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  !real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
  real(8) :: u, v, w, u_init, rho_init
  integer :: i, j, k, l, i_target, inlet_frame, k_off, ierr

  real(8) :: sigma_x, sigma_y, over_nz_physical
  real(8) :: rho_sum, u_sum, v_sum, w_sum, p_sum
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
    ! 0.05*abs(u_init)*sech^2(2*(y-Ly/2)/delta_bl).
    !===================================================================
    call diag_first_calls(myrank, nx, ny)
    call check_nan_abort(myrank, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
    inlet_frame = inlet_bc_calls / inlet_rk_stages
    inlet_bc_calls = inlet_bc_calls + 1
    k_off = k_offset

    ! 流入条件は先頭のxスラブだけが持つ。
    if (my_slab == 0) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 2, ny-1
        i = 1
        rho_init = rho_target_1d(j)
        u_init   = u_target_1d(j)
        fluctuation_scale = inlet_fluctuation_rms*abs(u_init)*inlet_envelope_1d(j)
        u = u_init + fluctuation_scale*inlet_normal(j, k+k_off, inlet_frame, 1, inlet_random_seed)
        v =          fluctuation_scale*inlet_normal(j, k+k_off, inlet_frame, 2, inlet_random_seed)
        w =          fluctuation_scale*inlet_normal(j, k+k_off, inlet_frame, 3, inlet_random_seed)

        Q_1(i,j,k) = rho_init*over_jacobian_tab(i,j)
        Q_2(i,j,k) = rho_init*u*over_jacobian_tab(i,j)
        Q_3(i,j,k) = rho_init*v*over_jacobian_tab(i,j)
        Q_4(i,j,k) = rho_init*w*over_jacobian_tab(i,j)
        Q_5(i,j,k) = (p/(gamma-1.d0) + 0.5d0*rho_init*(u**2+v**2+w**2))*over_jacobian_tab(i,j)
      enddo;enddo
    endif

    ! x流出スポンジ(z平均の目標値を含む)は最終のxスラブだけが持つ。
    if (my_slab == npx-1) then
    !===================================================================
    ! xスポンジ入口面の瞬時z平均を計算する。
    ! 発達した混合層の平均厚さを保ったまま、xスポンジ内の変動だけを
    ! 滑らかに減衰させるため、この平均場をx方向の緩和目標とする。
    !===================================================================
    i_target = nx - 1 - nsp_x
    over_nz_physical = 1.d0/dble(n_compute*(nz-6))
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      rho_sum = 0.d0
      u_sum   = 0.d0
      v_sum   = 0.d0
      w_sum   = 0.d0
      p_sum   = 0.d0
      do k = 4, nz-3
        rho_now = Q_1(i_target,j,k) * jacobian(i_target,j)
        u_now   = Q_2(i_target,j,k) / Q_1(i_target,j,k)
        v_now   = Q_3(i_target,j,k) / Q_1(i_target,j,k)
        w_now   = Q_4(i_target,j,k) / Q_1(i_target,j,k)
        p_now   = (gamma-1.d0) * (Q_5(i_target,j,k)*jacobian(i_target,j) &
                  - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
        rho_sum = rho_sum + rho_now
        u_sum   = u_sum   + u_now
        v_sum   = v_sum   + v_now
        w_sum   = w_sum   + w_now
        p_sum   = p_sum   + p_now
      enddo
      rho_target_x_1d(j) = rho_sum*over_nz_physical
      u_target_x_1d(j)   = u_sum*over_nz_physical
      v_target_x_1d(j)   = v_sum*over_nz_physical
      w_target_x_1d(j)   = w_sum*over_nz_physical
      p_target_x_1d(j)   = p_sum*over_nz_physical
    enddo

    !===================================================================
    ! xスポンジだけを処理する。全領域を走査せず、対象となる末端nsp_x点に
    ! カーネル範囲を限定してGPU上の不要なメモリアクセスを避ける。
    !===================================================================
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 2, ny-1
        do i = nx-nsp_x, nx-1
          sigma_x = sigma_x_1d(i)
          if (sigma_x > 0.d0) then
            rho_now = Q_1(i,j,k) * jacobian(i,j)
            u_now   = Q_2(i,j,k) / Q_1(i,j,k)
            v_now   = Q_3(i,j,k) / Q_1(i,j,k)
            w_now   = Q_4(i,j,k) / Q_1(i,j,k)
            p_now   = (gamma-1.d0) * (Q_5(i,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))

            rho_blend = (1.d0-sigma_x)*rho_now + sigma_x*rho_target_x_1d(j)
            u_blend   = (1.d0-sigma_x)*u_now   + sigma_x*u_target_x_1d(j)
            v_blend   = (1.d0-sigma_x)*v_now   + sigma_x*v_target_x_1d(j)
            w_blend   = (1.d0-sigma_x)*w_now   + sigma_x*w_target_x_1d(j)
            p_blend   = (1.d0-sigma_x)*p_now   + sigma_x*p_target_x_1d(j)
            Q_1(i,j,k) = rho_blend * over_jacobian_tab(i,j)
            Q_2(i,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
            Q_3(i,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
            Q_4(i,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
            Q_5(i,j,k) = (p_blend/(gamma-1.d0) &
                         + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)

          endif
        enddo;enddo;enddo
    endif

    !===================================================================
    ! 下側yスポンジ。xスポンジの後に適用して角部を外部一様流へ戻す。
    !===================================================================
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = 2, 1+nsp_y
        do i = 2, nx-1
          sigma_y = sigma_y_1d(j)
          if (sigma_y > 0.d0) then
            rho_now = Q_1(i,j,k) * jacobian(i,j)
            u_now   = Q_2(i,j,k) / Q_1(i,j,k)
            v_now   = Q_3(i,j,k) / Q_1(i,j,k)
            w_now   = Q_4(i,j,k) / Q_1(i,j,k)
            p_now   = (gamma-1.d0) * (Q_5(i,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
            rho_init = rho_target_1d(j)
            u_init   = u_target_1d(j)
            rho_blend = (1.d0-sigma_y)*rho_now + sigma_y*rho_init
            u_blend   = (1.d0-sigma_y)*u_now   + sigma_y*u_init
            v_blend   = (1.d0-sigma_y)*v_now
            w_blend   = (1.d0-sigma_y)*w_now
            p_blend   = (1.d0-sigma_y)*p_now + sigma_y*p
            Q_1(i,j,k) = rho_blend * over_jacobian_tab(i,j)
            Q_2(i,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
            Q_3(i,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
            Q_4(i,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
            Q_5(i,j,k) = (p_blend/(gamma-1.d0) &
                         + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)
          endif
        enddo;enddo;enddo

    !===================================================================
    ! 上側yスポンジ。
    !===================================================================
    !$cuf kernel do(3)<<<*,*>>>
    do k = 1, nz
      do j = ny-nsp_y, ny-1
        do i = 2, nx-1
          sigma_y = sigma_y_1d(j)
          if (sigma_y > 0.d0) then
            rho_now = Q_1(i,j,k) * jacobian(i,j)
            u_now   = Q_2(i,j,k) / Q_1(i,j,k)
            v_now   = Q_3(i,j,k) / Q_1(i,j,k)
            w_now   = Q_4(i,j,k) / Q_1(i,j,k)
            p_now   = (gamma-1.d0) * (Q_5(i,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
            rho_init = rho_target_1d(j)
            u_init   = u_target_1d(j)
            rho_blend = (1.d0-sigma_y)*rho_now + sigma_y*rho_init
            u_blend   = (1.d0-sigma_y)*u_now   + sigma_y*u_init
            v_blend   = (1.d0-sigma_y)*v_now
            w_blend   = (1.d0-sigma_y)*w_now
            p_blend   = (1.d0-sigma_y)*p_now + sigma_y*p
            Q_1(i,j,k) = rho_blend * over_jacobian_tab(i,j)
            Q_2(i,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
            Q_3(i,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
            Q_4(i,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
            Q_5(i,j,k) = (p_blend/(gamma-1.d0) &
                         + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)
          endif
        enddo;enddo;enddo

    !===================================================================
    ! x = nx 流出境界(1点ゴースト): 物理量で1次外挿してJ補正(最終のxスラブのみ)
    !===================================================================
    if (my_slab == npx-1) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
          Q_1(nx,j,k) = ( 2.0d0*(Q_1(nx-1,j,k)*jacobian(nx-1,j)) &
                             - (Q_1(nx-2,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)
          Q_2(nx,j,k) = ( 2.0d0*(Q_2(nx-1,j,k)*jacobian(nx-1,j)) &
                             - (Q_2(nx-2,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)    
          Q_3(nx,j,k) = ( 2.0d0*(Q_3(nx-1,j,k)*jacobian(nx-1,j)) &
                             - (Q_3(nx-2,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)
          Q_4(nx,j,k) = ( 2.0d0*(Q_4(nx-1,j,k)*jacobian(nx-1,j)) &
                             - (Q_4(nx-2,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)
          Q_5(nx,j,k) = ( 2.0d0*(Q_5(nx-1,j,k)*jacobian(nx-1,j)) &
                             - (Q_5(nx-2,j,k)*jacobian(nx-2,j)) ) * over_jacobian_tab(nx,j)                                                                           
        enddo;enddo
    endif
    !===================================================================
    ! y方向境界(j=1, j=ny): 1点ゴースト、物理量で1次外挿してJ補正
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
       
          Q_1(i,ny,k) = ( 2.0d0*(Q_1(i,ny-1,k)*jacobian(i,ny-1)) &
                             - (Q_1(i,ny-2,k)*jacobian(i,ny-2)) ) * over_jacobian_tab(i,ny)
          Q_2(i,ny,k) = ( 2.0d0*(Q_2(i,ny-1,k)*jacobian(i,ny-1)) &
                             - (Q_2(i,ny-2,k)*jacobian(i,ny-2)) ) * over_jacobian_tab(i,ny)
          Q_3(i,ny,k) = ( 2.0d0*(Q_3(i,ny-1,k)*jacobian(i,ny-1)) &
                             - (Q_3(i,ny-2,k)*jacobian(i,ny-2)) ) * over_jacobian_tab(i,ny)
          Q_4(i,ny,k) = ( 2.0d0*(Q_4(i,ny-1,k)*jacobian(i,ny-1)) &
                             - (Q_4(i,ny-2,k)*jacobian(i,ny-2)) ) * over_jacobian_tab(i,ny)
          Q_5(i,ny,k) = ( 2.0d0*(Q_5(i,ny-1,k)*jacobian(i,ny-1)) &
                             - (Q_5(i,ny-2,k)*jacobian(i,ny-2)) ) * over_jacobian_tab(i,ny)
                             
                             
          Q_1(i,1,k) = ( 2.0d0*(Q_1(i,2,k)*jacobian(i,2)) &
                            - (Q_1(i,3,k)*jacobian(i,3)) ) * over_jacobian_tab(i,1)
          Q_2(i,1,k) = ( 2.0d0*(Q_2(i,2,k)*jacobian(i,2)) &
                            - (Q_2(i,3,k)*jacobian(i,3)) ) * over_jacobian_tab(i,1)
          Q_3(i,1,k) = ( 2.0d0*(Q_3(i,2,k)*jacobian(i,2)) &
                            - (Q_3(i,3,k)*jacobian(i,3)) ) * over_jacobian_tab(i,1)
          Q_4(i,1,k) = ( 2.0d0*(Q_4(i,2,k)*jacobian(i,2)) &
                            - (Q_4(i,3,k)*jacobian(i,3)) ) * over_jacobian_tab(i,1)
          Q_5(i,1,k) = ( 2.0d0*(Q_5(i,2,k)*jacobian(i,2)) &
                            - (Q_5(i,3,k)*jacobian(i,3)) ) * over_jacobian_tab(i,1)                                                                                 
        enddo;enddo
    !===================================================================
    ! z方向境界: 周期境界(3点ゴースト)。zは分割しないので自スラブ内コピー。
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do k = 1, 3
          Q_1(i,j,k) = Q_1(i,j,nz-6+k)
          Q_2(i,j,k) = Q_2(i,j,nz-6+k)
          Q_3(i,j,k) = Q_3(i,j,nz-6+k)
          Q_4(i,j,k) = Q_4(i,j,nz-6+k)
          Q_5(i,j,k) = Q_5(i,j,nz-6+k)

          Q_1(i,j,nz-3+k) = Q_1(i,j,3+k)
          Q_2(i,j,nz-3+k) = Q_2(i,j,3+k)
          Q_3(i,j,nz-3+k) = Q_3(i,j,3+k)
          Q_4(i,j,nz-3+k) = Q_4(i,j,3+k)
          Q_5(i,j,nz-3+k) = Q_5(i,j,3+k)
    enddo;enddo;enddo 

    ! x方向の接続面: 隣スラブの内部3面をゴーストに受け取る(npx=1なら何もしない)。
    call exchange_x(nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)

    ! パッシブスカラー(混合分率)。流れ場には一切影響しない。
    if (scalar_on) call sc_step(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)

  end subroutine set_bc

  subroutine set_bc_mut(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    call set_bc_mut_common(nx, ny, nz, mut, qc2)
  end subroutine set_bc_mut

end module set



! module set
!   use mod_globals , only : nx, ny, nz, Lx, Ly, Lz, gamma, R,dt,nt, &
!                            Lx_main, Lx_buf,Ly_main, Ly_buf,nx_main, nx_buf,ny_main, ny_buf
!   use set_bc_common
!   use set_coordinate
!   implicit none


!   !共通変数
!   real(8), device  :: xi(nx), yj(ny), zk(nz)
!   real(8), device  :: rho_target_1d(ny), u_target_1d(ny)
!   real(8), device  :: inlet_envelope_1d(ny)           ! shear-layer-localized inlet noise envelope
!   real(8), device  :: amp_u_env_1d(ny)                ! amp*u_init(j)*env(j) を事前計算(追加最適化)
!   real(8), device  :: sigma_x_1d(nx), sigma_y_1d(ny)  ! スポンジ係数(x,y)を事前計算
!   real(8), device  :: rho_target_x_1d(ny), u_target_x_1d(ny)
!   real(8), device  :: v_target_x_1d(ny), w_target_x_1d(ny), p_target_x_1d(ny)
!   real(8), device  :: over_jacobian_tab(nx,ny)        ! 1/Jacobian を事前計算(格子は時間不変)
!   integer, save    :: inlet_bc_calls = 0

!   ! スポンジ・強制振動パラメータ(set_grid_main_bufferとset_bc双方から使うためモジュールスコープへ)
!   integer, parameter :: nsp_x = nx_buf
!   integer, parameter :: nsp_y = ny_buf/2
!   real(8), parameter :: sigma_max_x = 0.06d0
!   real(8), parameter :: sigma_max_y = 0.01d0

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
!   dz_min = Lz / dble(nz-6)

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

!   if (n < 2 .or. dr_0 <= 0.d0 .or. L_target <= 0.d0 .or. r0 < 1.d0-1.d-12) &
!     error stop 'Invalid stretched buffer parameters'
!   r_lo = 1.d0
!   if (calc_ramp_total(dr_0, r0, r_lo, n) > L_target+tol) &
!     error stop 'Buffer too short for main-domain edge spacing'
!   r_hi = max(1.01d0, r0)
!   do while (calc_ramp_total(dr_0, r0, r_hi, n) < L_target)
!     r_hi = 2.d0*r_hi
!   enddo

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
!       if (total > L_target) return ! bracket comparison only; avoid overflow
!       dr_local = dr_local * r_local
!     enddo
!   end function calc_ramp_total
! end function calc_smooth_stretch_ratio


! subroutine set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
!   nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)

!   use grid_partitioned, only : partitioned_axis
!   use mod_globals, only : Ly_uniform, ny_uniform
!   integer, intent(in) :: nx, ny, nz
!   real(8), intent(in) :: Lx_main, Lx_buf, Ly_main, Ly_buf, Lz
!   integer, intent(in) :: nx_main, nx_buf, ny_main, ny_buf
!   real(8), intent(out) :: xc(nx), yc(ny), zc(nz)
!   real(8), intent(out) :: dx(nx-1), dy(ny-1), dz(nz-1)

!   real(8) :: x(nx+1), y(ny+1), z(nz+1)
!   real(8) :: dz1
!   real(8) :: s, dr, eta, pi
!   integer :: i, j, k

!   ! 格子伸長パラメータ
!   real(8)  :: r0_x, r1_x               ! x方向バッファ: 開始/終端伸び率
!   real(8)  :: r0_y_lo, r1_y_lo         ! y方向下バッファ: 開始/終端伸び率
!   real(8)  :: r0_y_hi, r1_y_hi         ! y方向上バッファ: 開始/終端伸び率
!   real(8)  :: r_local

!   pi = acos(-1.d0)

!   ! x main domain is uniform; only the separate downstream buffer stretches.
!   if (nx_main < 3 .or. Lx_main <= 0.d0) error stop 'Invalid uniform x mesh'
!   do i = 2, nx_main+1
!     x(i) = Lx_main*dble(i-2)/dble(nx_main-1)
!   enddo

!   ! 接続点の格子幅とその直前の局所伸び率(主計算領域側から引き継ぐ)
!   dr  = x(nx_main+1) - x(nx_main)
!   r0_x = 1.d0

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
!   ! y方向: 下バッファ(滑らかランプ) + 主計算領域(中心一様・外側伸長) + 上バッファ(滑らかランプ)
!   !=================================================================
!   ! y main domain: symmetric stretching outside the central uniform section.
!   call partitioned_axis(ny_main, Ly_main, Ly_uniform, ny_uniform, .true., &
!                         y(ny_buf+2:ny_buf+ny_main+1))
!   y(ny_buf+2:ny_buf+ny_main+1) = y(ny_buf+2:ny_buf+ny_main+1)+Ly_buf

!   ! --- 下バッファ: 接続点の格子幅と、その外側への局所伸び率を引き継ぐ ---
!   dr = y(ny_buf+3) - y(ny_buf+2)
!   r0_y_lo = dr / (y(ny_buf+4) - y(ny_buf+3))

!   r1_y_lo = calc_smooth_stretch_ratio(dr, r0_y_lo, ny_buf, Ly_buf)

!   print *, "下バッファ 開始伸び率=", r0_y_lo, " 終端伸び率=", r1_y_lo
!   print *, "下バッファ合計長さ確認 ="

!   do j = ny_buf+1, 2, -1
!     eta = dble((ny_buf+2-1) - j) / dble(max(ny_buf-1,1))
!     r_local = r0_y_lo + (r1_y_lo - r0_y_lo) * (1.0d0 - cos(pi*eta)) / 2.0d0
!     y(j) = y(j+1) - dr
!     dr = dr * r_local
!   enddo

!   print *, "y(ny_buf+2) - y(2) =", y(ny_buf+2) - y(2), " (目標:", Ly_buf, ")"

!   ! --- 上バッファ ---
!   dr = y(ny_buf+ny_main+1) - y(ny_buf+ny_main)
!   r0_y_hi = dr / (y(ny_buf+ny_main) - y(ny_buf+ny_main-1))

!   r1_y_hi = calc_smooth_stretch_ratio(dr, r0_y_hi, ny_buf, Ly_buf)

!   print *, "上バッファ 開始伸び率=", r0_y_hi, " 終端伸び率=", r1_y_hi

!   do j = ny_buf+ny_main+2, ny-1
!     eta = dble(j - (ny_buf+ny_main+2)) / dble(max(ny_buf-1,1))
!     r_local = r0_y_hi + (r1_y_hi - r0_y_hi) * (1.0d0 - cos(pi*eta)) / 2.0d0
!     y(j) = y(j-1) + dr
!     dr = dr * r_local
!   enddo

!   print *, "y(ny-1) - y(ny_buf+ny_main+1) =", y(ny-1)-y(ny_buf+ny_main+1), &
!            " (目標:", Ly_buf, ")"

!   ! 下ゴースト(j=1)
!   y(1) = y(2) - (y(3)-y(2))
!   ! 上ゴースト
!   y(ny)   = y(ny-1) + dr
!   y(ny+1) = y(ny)   + dr

!   !=================================================================
!   ! z方向: 一様(周期境界、3点ゴースト)
!   !=================================================================
!   dz1 = Lz / dble(nz-6)
!   do k = 4, nz-2
!     z(k) = dz1 * dble(k-4)
!   enddo
!   z(3) = z(4) - dz1
!   z(2) = z(3) - dz1
!   z(1) = z(2) - dz1
!   z(nz-1) = z(nz-2) + dz1
!   z(nz)   = z(nz-1) + dz1
!   z(nz+1) = z(nz)   + dz1

!   !共通変数に代入 CPU --- GPU
!   xi = x(1:nx)
!   yj = y(1:ny)
!   zk = z(1:nz)

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
! ! スポンジ・流入境界用の目標値を、事前に1回だけ計算してテーブル化
! ! (これらはすべて i, k, 時間に依存しない量なので、set_bc内で毎回
! !  計算し直す必要がない。格子は時間不変(RESCALE=False)なので安全)
! !=================================================================
! block
!   use mod_globals, only : u1, u2, T1, T2, p, delta_bl, amp
!   real(8) :: rho_target_1d_h(ny), u_target_1d_h(ny), amp_u_env_1d_h(ny)
!   real(8) :: inlet_envelope_1d_h(ny)
!   real(8) :: sigma_x_1d_h(nx), sigma_y_1d_h(ny)
!   real(8) :: T_init_local, env_local, over_delta_bl_local, over_Lz_local
!   real(8) :: over_nsp_x_local, over_nsp_y_local, eta_local, ramp_local
!   real(8) :: sigma_y_lo_local, sigma_y_hi_local
!   integer :: jj, kk, ii

!   over_delta_bl_local = 1.0d0/delta_bl
!   over_Lz_local = 1.0d0/Lz
!   over_nsp_x_local = 1.0d0/dble(nsp_x)
!   over_nsp_y_local = 1.0d0/dble(nsp_y)

!   ! --- y方向テーブル: rho_target, u_target, amp*u_init*env, sigma_y ---
!   do jj = 1, ny
!     T_init_local = 0.5d0*(T1+T2) + 0.5d0*(T1-T2)*tanh(2.0d0*(y(jj)-0.5d0*Ly)/delta_bl)
!     rho_target_1d_h(jj) = p / (R * T_init_local)
!     u_target_1d_h(jj)   = 0.5d0*(u1+u2) + 0.5d0*(u1-u2)*tanh(2.0d0*(y(jj)-0.5d0*Ly)/delta_bl)

!     ! Localize the 5% Gaussian inlet forcing to the region where the
!     ! tanh mean profile has shear.  This is the normalized mean-velocity
!     ! gradient: sech^2(2*(y-Ly/2)/delta_bl), evaluated stably as 1-tanh^2.
!     env_local = tanh(2.0d0*(y(jj)-0.5d0*Ly)*over_delta_bl_local)
!     inlet_envelope_1d_h(jj) = max(0.d0, 1.d0-env_local**2)

!     env_local = exp(-((y(jj)-0.5d0*Ly)*over_delta_bl_local)**2)
!     amp_u_env_1d_h(jj) = amp * u_target_1d_h(jj) * env_local

!     sigma_y_lo_local = 0.d0
!     sigma_y_hi_local = 0.d0
!     if (jj <= 1+nsp_y) then
!       eta_local = dble((1+nsp_y) - jj) * over_nsp_y_local
!       if (eta_local < 0.d0) eta_local = 0.d0
!       if (eta_local > 1.d0) eta_local = 1.d0
!       ! 5次smoothstepを使い、終端の減衰率がsigma_max_yとなるよう
!       ! (1-sigma)=(1-sigma_max_y)**ramp から逆算する。
!       ramp_local = eta_local**3 * (10.d0 - 15.d0*eta_local + 6.d0*eta_local**2)
!       sigma_y_lo_local = 1.d0 - (1.d0-sigma_max_y)**ramp_local
!     endif
!     if (jj >= ny-nsp_y) then
!       eta_local = dble(jj - (ny-nsp_y)) * over_nsp_y_local
!       if (eta_local < 0.d0) eta_local = 0.d0
!       if (eta_local > 1.d0) eta_local = 1.d0
!       ramp_local = eta_local**3 * (10.d0 - 15.d0*eta_local + 6.d0*eta_local**2)
!       sigma_y_hi_local = 1.d0 - (1.d0-sigma_max_y)**ramp_local
!     endif
!     sigma_y_1d_h(jj) = max(sigma_y_lo_local, sigma_y_hi_local)
!   enddo

!   ! --- x方向テーブル: sigma_x ---
!   do ii = 1, nx
!     sigma_x_1d_h(ii) = 0.d0
!     if (ii >= nx-nsp_x) then
!       eta_local = dble(ii - (nx-1-nsp_x)) * over_nsp_x_local
!       if (eta_local < 0.d0) eta_local = 0.d0
!       if (eta_local > 1.d0) eta_local = 1.d0
!       ramp_local = eta_local**3 * (10.d0 - 15.d0*eta_local + 6.d0*eta_local**2)
!       sigma_x_1d_h(ii) = 1.d0 - (1.d0-sigma_max_x)**ramp_local
!     endif
!   enddo

!   ! CPU -> GPU 転送(いずれも1回のみ、サイズも小さく無視できるコスト)
!   rho_target_1d = rho_target_1d_h
!   u_target_1d   = u_target_1d_h
!   inlet_envelope_1d = inlet_envelope_1d_h
!   amp_u_env_1d  = amp_u_env_1d_h
!   sigma_x_1d    = sigma_x_1d_h
!   sigma_y_1d    = sigma_y_1d_h
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

! end subroutine set_grid_main_buffer




!   !=====================================================================
!   ! ヤコビアン: 内部は中央差分、ゴーストセル(1点)は隣接内部セルの値を複製
!   !=====================================================================
!   subroutine set_Jacobian_xy3_stretch(nx, ny, nz, xc, yc, zc, Jacobian)
!     integer, intent(in)  :: nx, ny, nz
!     real(8), intent(in)  :: xc(nx), yc(ny), zc(nz)
!     real(8), intent(out) :: Jacobian(nx,ny)
!     integer :: i, j
!     real(8) :: dx(nx-1), dy(ny-1), dz1

!     dx = xc(2:nx) - xc(1:nx-1)
!     dy = yc(2:ny) - yc(1:ny-1)
!     dz1 = zc(2) - zc(1)

!     do j = 2, ny-1
!       do i = 2, nx-1
!         ! Same spacing-based expression as set_Jacobian_xy3.
!         Jacobian(i,j) = 8.d0 / &
!           ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * 2.d0*dz1)
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


!   subroutine set_init(myrank, nx, ny, nz, x, y, z, Q)
!     use mod_globals, only : u1, u2, p, Ly, T1, T2, delta_bl
!     use mod_constant, only : id_accuracy
!     integer, intent(in)  :: myrank, nx, ny, nz
!     real(8), intent(in)  :: x(nx), y(ny), z(nz)
!     real(8), intent(out) :: Q(nx,5,ny,nz)
!     integer :: i, j, k
!     real(8) :: u_init, rho_init, T_init

!     ! Start from the unperturbed mean mixing layer.  All velocity
!     ! fluctuations enter through the x-inlet boundary in set_bc.
!     do k = 1, nz
!       do j = 1, ny
!         T_init   = 0.5d0*(T1+T2) + 0.5d0*(T1-T2)*tanh(2.0d0*(y(j)-0.5d0*Ly)/delta_bl)
!         rho_init = p / (R * T_init)
!         u_init = 0.5*(u1+u2) + 0.5*(u1-u2)*tanh(2.0d0*(y(j)-0.5*Ly)/(delta_bl))

!         do i = 1, nx
!           Q(i,1,j,k) = rho_init
!           Q(i,2,j,k) = rho_init * u_init
!           Q(i,3,j,k) = 0.d0
!           Q(i,4,j,k) = 0.d0
!           Q(i,5,j,k) = p/(gamma-1.d0) + 0.5d0*rho_init*u_init**2
!         enddo;enddo;enddo

!   end subroutine set_init



! subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q, Qre)
!   use mod_constant, only : id_accuracy
!   use mod_globals, only : u1, rho1, u2, rho2, p, amp, dt, gamma, Ly, Lz, Lx,T1,T2,delta_bl, step_offset, &
!                           inlet_fluctuation_rms, inlet_random_seed, inlet_rk_stages
!   use set_coordinate
!   integer, intent(in), value            :: myrank, nx, ny, nz
!   real(8), intent(in), device           :: jacobian(nx,ny)
!   real(8), intent(inout), device        :: Q(nx,5,ny,nz)
!   real(8), intent(in), device, optional :: Qre(ny*(nz-6)*5)
!   real(8) :: u, v, w, u_init, rho_init
!   integer :: i, j, k, l, i_target, inlet_frame

!   real(8) :: sigma_x, sigma_y, over_nz_physical
!   real(8) :: rho_sum, u_sum, v_sum, w_sum, p_sum
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
!     ! 0.05*abs(u_init)*sech^2(2*(y-Ly/2)/delta_bl).
!     !===================================================================
!     inlet_frame = inlet_bc_calls / inlet_rk_stages
!     inlet_bc_calls = inlet_bc_calls + 1

!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 4, nz-3
!       do j = 2, ny-1
!         i = 1
!         rho_init = rho_target_1d(j)
!         u_init   = u_target_1d(j)
!         fluctuation_scale = inlet_fluctuation_rms*abs(u_init)*inlet_envelope_1d(j)
!         u = u_init + fluctuation_scale*inlet_normal(j, k, inlet_frame, 1, inlet_random_seed)
!         v =          fluctuation_scale*inlet_normal(j, k, inlet_frame, 2, inlet_random_seed)
!         w =          fluctuation_scale*inlet_normal(j, k, inlet_frame, 3, inlet_random_seed)

!         Q(i,1,j,k) = rho_init*over_jacobian_tab(i,j)
!         Q(i,2,j,k) = rho_init*u*over_jacobian_tab(i,j)
!         Q(i,3,j,k) = rho_init*v*over_jacobian_tab(i,j)
!         Q(i,4,j,k) = rho_init*w*over_jacobian_tab(i,j)
!         Q(i,5,j,k) = (p/(gamma-1.d0) + 0.5d0*rho_init*(u**2+v**2+w**2))*over_jacobian_tab(i,j)
!       enddo;enddo

!     !===================================================================
!     ! xスポンジ入口面の瞬時z平均を計算する。
!     ! 発達した混合層の平均厚さを保ったまま、xスポンジ内の変動だけを
!     ! 滑らかに減衰させるため、この平均場をx方向の緩和目標とする。
!     !===================================================================
!     i_target = nx - 1 - nsp_x
!     over_nz_physical = 1.d0/dble(nz-6)
!     !$cuf kernel do(1)<<<*,*>>>
!     do j = 2, ny-1
!       rho_sum = 0.d0
!       u_sum   = 0.d0
!       v_sum   = 0.d0
!       w_sum   = 0.d0
!       p_sum   = 0.d0
!       do k = 4, nz-3
!         rho_now = Q(i_target,1,j,k) * jacobian(i_target,j)
!         u_now   = Q(i_target,2,j,k) / Q(i_target,1,j,k)
!         v_now   = Q(i_target,3,j,k) / Q(i_target,1,j,k)
!         w_now   = Q(i_target,4,j,k) / Q(i_target,1,j,k)
!         p_now   = (gamma-1.d0) * (Q(i_target,5,j,k)*jacobian(i_target,j) &
!                   - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
!         rho_sum = rho_sum + rho_now
!         u_sum   = u_sum   + u_now
!         v_sum   = v_sum   + v_now
!         w_sum   = w_sum   + w_now
!         p_sum   = p_sum   + p_now
!       enddo
!       rho_target_x_1d(j) = rho_sum*over_nz_physical
!       u_target_x_1d(j)   = u_sum*over_nz_physical
!       v_target_x_1d(j)   = v_sum*over_nz_physical
!       w_target_x_1d(j)   = w_sum*over_nz_physical
!       p_target_x_1d(j)   = p_sum*over_nz_physical
!     enddo

!     !===================================================================
!     ! xスポンジだけを処理する。全領域を走査せず、対象となる末端nsp_x点に
!     ! カーネル範囲を限定してGPU上の不要なメモリアクセスを避ける。
!     !===================================================================
!     !$cuf kernel do(3)<<<*,*>>>
!     do k = 1, nz
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

!             rho_blend = (1.d0-sigma_x)*rho_now + sigma_x*rho_target_x_1d(j)
!             u_blend   = (1.d0-sigma_x)*u_now   + sigma_x*u_target_x_1d(j)
!             v_blend   = (1.d0-sigma_x)*v_now   + sigma_x*v_target_x_1d(j)
!             w_blend   = (1.d0-sigma_x)*w_now   + sigma_x*w_target_x_1d(j)
!             p_blend   = (1.d0-sigma_x)*p_now   + sigma_x*p_target_x_1d(j)
!             Q(i,1,j,k) = rho_blend * over_jacobian_tab(i,j)
!             Q(i,2,j,k) = rho_blend*u_blend * over_jacobian_tab(i,j)
!             Q(i,3,j,k) = rho_blend*v_blend * over_jacobian_tab(i,j)
!             Q(i,4,j,k) = rho_blend*w_blend * over_jacobian_tab(i,j)
!             Q(i,5,j,k) = (p_blend/(gamma-1.d0) &
!                          + 0.5d0*rho_blend*(u_blend**2+v_blend**2+w_blend**2)) * over_jacobian_tab(i,j)

!           endif
!         enddo;enddo;enddo

!     !===================================================================
!     ! 下側yスポンジ。xスポンジの後に適用して角部を外部一様流へ戻す。
!     !===================================================================
!     !$cuf kernel do(3)<<<*,*>>>
!     do k = 1, nz
!       do j = 2, 1+nsp_y
!         do i = 2, nx-1
!           sigma_y = sigma_y_1d(j)
!           if (sigma_y > 0.d0) then
!             rho_now = Q(i,1,j,k) * jacobian(i,j)
!             u_now   = Q(i,2,j,k) / Q(i,1,j,k)
!             v_now   = Q(i,3,j,k) / Q(i,1,j,k)
!             w_now   = Q(i,4,j,k) / Q(i,1,j,k)
!             p_now   = (gamma-1.d0) * (Q(i,5,j,k)*jacobian(i,j) &
!                       - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
!             rho_init = rho_target_1d(j)
!             u_init   = u_target_1d(j)
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
!     ! 上側yスポンジ。
!     !===================================================================
!     !$cuf kernel do(3)<<<*,*>>>
!     do k = 1, nz
!       do j = ny-nsp_y, ny-1
!         do i = 2, nx-1
!           sigma_y = sigma_y_1d(j)
!           if (sigma_y > 0.d0) then
!             rho_now = Q(i,1,j,k) * jacobian(i,j)
!             u_now   = Q(i,2,j,k) / Q(i,1,j,k)
!             v_now   = Q(i,3,j,k) / Q(i,1,j,k)
!             w_now   = Q(i,4,j,k) / Q(i,1,j,k)
!             p_now   = (gamma-1.d0) * (Q(i,5,j,k)*jacobian(i,j) &
!                       - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
!             rho_init = rho_target_1d(j)
!             u_init   = u_target_1d(j)
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
!     ! y方向境界(j=1, j=ny): 1点ゴースト、物理量で1次外挿してJ補正
!     !===================================================================
!     !$cuf kernel do(2)<<<*,*>>>
!     do k = 1, nz
!       do i = 1, nx
!         do l = 1, 5
!           Q(i,l,ny,k) = ( 2.0d0*(Q(i,l,ny-1,k)*jacobian(i,ny-1)) &
!                              - (Q(i,l,ny-2,k)*jacobian(i,ny-2)) ) * over_jacobian_tab(i,ny)

!           Q(i,l,1,k) = ( 2.0d0*(Q(i,l,2,k)*jacobian(i,2)) &
!                             - (Q(i,l,3,k)*jacobian(i,3)) ) * over_jacobian_tab(i,1)
!         enddo;enddo;enddo 
!     !===================================================================
!     ! z方向境界: 周期境界(3点ゴースト)
!     !===================================================================
!     !$cuf kernel do(2)<<<*,*>>>
!     do j = 1, ny
!       do i = 1, nx
!         do l = 1, 5
!           Q(i,l,j,1) = Q(i,l,j,nz-5)
!           Q(i,l,j,2) = Q(i,l,j,nz-4)
!           Q(i,l,j,3) = Q(i,l,j,nz-3)
!           Q(i,l,j,nz-2) = Q(i,l,j,4)
!           Q(i,l,j,nz-1) = Q(i,l,j,5)
!           Q(i,l,j,nz)   = Q(i,l,j,6)
!         enddo;enddo;enddo 

!   end subroutine set_bc

!   subroutine set_bc_mut(nx, ny, nz, mut, qc2)
!     integer, intent(in), value     :: nx, ny, nz
!     real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
!     call set_bc_mut_common(nx, ny, nz, mut, qc2)
!   end subroutine set_bc_mut

! end module set
