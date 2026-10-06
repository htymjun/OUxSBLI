module set
  use mod_globals , only : nx, ny, nz, Lx, Ly, Lz, gamma, R,dt,nt, &
                           Lx_main, Lx_buf,Ly_main, Ly_buf,Lz_main,Lz_buf,nx_main, nx_buf,ny_main, ny_buf, nz_main, nz_buf, &
                           npx, nx_global
  use set_bc_common
  use set_coordinate
  implicit none

  
  ! グリッド生成用の内部ルーチン(他モジュールから見せない)
  private :: partitioned_axis, stretched_widths, buffer_faces, calc_smooth_stretch_ratio

  real(8), device :: rho_target(ny,nz), u_target(ny,nz), inlet_envelope(ny,nz)
  real(8), device :: sigma_x_1d(nx), sigma_y_1d(ny), sigma_z_1d(nz)
  real(8), device :: rho_target_x(ny,nz), u_target_x(ny,nz)
  real(8), device :: v_target_x(ny,nz), w_target_x(ny,nz), p_target_x(ny,nz)
  real(8), device :: over_jacobian_tab(nx,ny)
  integer, save :: inlet_bc_calls = 0

  ! Non-uniform z without touching the shared solver.  main.f90 builds the 2-D
  ! Jacobian with dz(1) (set_Jacobian_xy3) and calc_R weights the x/y faces with
  ! the local width hz(k) = (dz(k-1)+dz(k))/2, so the update of plane k equals
  ! s(k) = hz(k)/dz(1) times the finite-volume update dt*(dE/hx + dF/hy + dG/hz).
  ! Every RK4 stage output is Q0 - (linear combination of residuals), with Q0 the
  ! state at the start of the step, so set_bc restores the exact update with
  ! Q = Q0 + (Q - Q0)/s(k) before any boundary/sponge treatment.
  real(8), device :: inv_z_update_scale(nz)
  real(8), allocatable, device :: q0_1(:,:,:), q0_2(:,:,:), q0_3(:,:,:), q0_4(:,:,:), q0_5(:,:,:)
  real(8), save :: dx_saved(nx-1), dy_saved(ny-1), dz_saved(nz-1)
  integer, parameter :: nsp_x = nx_buf, nsp_y = 7

  ! ---- x方向MPI分割(mod_globalsのnpx)。偶数(計算)ランクごとに1スラブ、奇数(I/O)ランクは相方と同じスラブ。
  ! global i = local i + i_offset。接続側のゴーストは3面で、set_bcの最後(exchange_x)で隣スラブの
  ! 内部3面を受け取る。共有カーネルは配列端から3面目以降の面/セルを高次で計算するので、
  ! 内部セルは分割なしと同じステンシルになる。
  integer, save :: my_slab = 0          ! このランクのスラブ番号(0..npx-1)
  integer, save :: i_offset = 0         ! ローカル→全体のx番号オフセット
  integer, save :: comm_compute         ! 計算ランクだけのコミュニケータ
  ! 配列端のヤコビアンは隣の列のコピー(set_Jacobian_xy3)なので、端のゴースト面だけ
  ! QJ = Q/J を真のJとの比で換算して受け取る(一様部では1)。
  real(8), save :: x_fac_lo = 1.d0, x_fac_hi = 1.d0
  real(8), allocatable, device :: xs_lo_d(:), xs_hi_d(:), xr_lo_d(:), xr_hi_d(:)
  ! ホスト側の送受信バッファはピン留めメモリ(ページ可能メモリより転送が速い)
  real(8), allocatable, pinned :: xs_lo_h(:), xs_hi_h(:), xr_lo_h(:), xr_hi_h(:)

  ! ---- パッシブスカラー(混合分率 xi)。mod_globalsのscalar_onで有効化 ----
  real(8), allocatable, device :: sc_phi(:,:,:)       ! xi (時刻n)
  real(8), allocatable, device :: sc_phis(:,:,:)      ! RK中間段
  real(8), allocatable, device :: sc_phib(:,:,:)      ! RK 段2の結果
  real(8), allocatable, device :: sc_gam(:,:,:)       ! rho*D = mu/Sc
  real(8), device  :: sc_mx(nx), sc_my(ny), sc_mz(nz)           ! セル中心のメトリック 1/h = 2/(d(i-1)+d(i))
  real(8), device  :: sc_idx(nx-1), sc_idy(ny-1), sc_idz(nz-1)  ! 1/(中心間隔)
  real(8), device  :: sc_in(ny,nz)                    ! 流入分布(初期の十字型分布と同じ)
  real(4), save    :: sc_xc(nx), sc_yc(ny), sc_zc(nz) ! VTK出力用のセル中心座標(流れ場のVTKと同じ)
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


!=====================================================================
! パッシブスカラー(混合分率 xi)。SL_spatial_inletと同じ方式を、このケースの十字型分布・
! 非周期で不等間隔のz・x分割に合わせて移植したもの。流れ場には一切影響しない。
!   xi = 1 : 流れ1(u1)側の流体、 xi = 0 : 流れ2(u2)側の流体
!   d(xi)/dt = -u.grad(xi) + (1/rho) div( (mu/Sc) grad(xi) )
! 流れ場の1ステップごと(最終RK段のset_bc)に、更新後の速度場でSSP-RK3により1ステップ進める。
! mod_globalsのscalar_on, Sc, scalar_output_everyで制御する。
!=====================================================================


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
! パッシブスカラー: 境界条件。流入は固定の十字型分布(先頭スラブ)、流出(最終スラブ)と
! y,zの両端は勾配ゼロ。zは周期ではない。xの接続面はsc_exchange_xで隣スラブから受け取る。
!=====================================================================
subroutine sc_bc(myrank, nx, ny, nz, f)
  integer, intent(in), value     :: myrank, nx, ny, nz
  real(8), intent(inout), device :: f(nx,ny,nz)
  integer :: i, j, k

  if (my_slab == 0) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        f(1,j,k)  = sc_in(j,k)
    enddo;enddo
  endif
  if (my_slab == npx-1) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        f(nx,j,k) = f(nx-1,j,k)
    enddo;enddo
  endif

  !$cuf kernel do(2)<<<*,*>>>
  do k = 2, nz-1
    do i = 1, nx
      f(i,1,k)  = f(i,2,k)
      f(i,ny,k) = f(i,ny-1,k)
  enddo;enddo

  !$cuf kernel do(2)<<<*,*>>>
  do j = 1, ny
    do i = 1, nx
      f(i,j,1)  = f(i,j,2)
      f(i,j,nz) = f(i,j,nz-1)
  enddo;enddo

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
  do k = 1, nz
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

end subroutine sc_calc_gam


!=====================================================================
! パッシブスカラー: 1セルの右辺  d(xi)/dt = -u.grad(xi) + (1/rho) div(rho*D grad(xi))
! 移流はKorenリミッタ付き3次風上(端は添字クランプで1次風上に落ちる)、拡散は2次中心。
! 近傍のxiは最初にまとめて読み込む(互いに独立な読み出しをまとめて発行してレイテンシを重ねる)。
! 読み込んだ値は界面値と拡散項の計算で使い回す。全セルを計算する(省略はしない)。
!=====================================================================
attributes(device) function sc_rhs_cell(nx, ny, nz, i, j, k, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam) result(r)
  integer, intent(in), value  :: nx, ny, nz, i, j, k
  real(8), intent(in), device :: jacobian(nx,ny)
  real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
  real(8), intent(in), device :: f(nx,ny,nz), gam(nx,ny,nz)
  real(8) :: r
  integer :: im2, ip2, jm2, jp2, km2, kp2
  real(8) :: fc, fim1, fip1, fim2, fip2, fjm1, fjp1, fjm2, fjp2, fkm1, fkp1, fkm2, fkp2
  real(8) :: rho, u, v, w, adv, dif, d, dz_adv

  im2 = max(i-2, 1);  ip2 = min(i+2, nx)
  jm2 = max(j-2, 1);  jp2 = min(j+2, ny)
  km2 = max(k-2, 1);  kp2 = min(k+2, nz)

  fc   = f(i,j,k)
  fim1 = f(i-1,j,k);  fip1 = f(i+1,j,k);  fim2 = f(im2,j,k);  fip2 = f(ip2,j,k)
  fjm1 = f(i,j-1,k);  fjp1 = f(i,j+1,k);  fjm2 = f(i,jm2,k);  fjp2 = f(i,jp2,k)
  fkm1 = f(i,j,k-1);  fkp1 = f(i,j,k+1);  fkm2 = f(i,j,km2);  fkp2 = f(i,j,kp2)

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
  ! --- z (不等間隔格子: 面値の再構成はx,yと同じ、メトリックは局所の格子幅) ---
  adv = adv + w*dz_adv*sc_mz(k)
  dif = dif + sc_mz(k)*( 0.5d0*(gam(i,j,k)+gam(i,j,k+1))*(fkp1-fc)*sc_idz(k) &
                       - 0.5d0*(gam(i,j,k)+gam(i,j,k-1))*(fc-fkm1)*sc_idz(k-1) )

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

  !$cuf kernel do(3)<<<*,(32,4,1)>>>
  do k = 2, nz-1
    do j = 2, ny-1
      do i = 2, nx-1
        out(i,j,k) = ca*phi0(i,j,k) + cb*( f(i,j,k) + dt*sc_rhs_cell(nx, ny, nz, i, j, k, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam) )
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

  !$cuf kernel do(3)<<<*,(32,4,1)>>>
  do k = 2, nz-1
    do j = 2, ny-1
      do i = 2, nx-1
        phi(i,j,k) = ca*phi(i,j,k) + cb*( f(i,j,k) + dt*sc_rhs_cell(nx, ny, nz, i, j, k, jacobian, Q_1, Q_2, Q_3, Q_4, f, gam) )
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
        f(i,j,k) = sc_in(j,k)
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
! ParaViewで直接開ける。ゴースト面も流れ場と同じく含む(xの接続部は隣スラブと重複)。
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

  allocate(h(nx,ny,nz), h4(nx,ny,nz))
  h  = sc_phi
  h4 = real(h, 4)
  fmin = minval(h(2:nx-1,2:ny-1,2:nz-1))
  fmax = maxval(h(2:nx-1,2:ny-1,2:nz-1))
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
  print '(1x,a,i0,a,a,a,es13.6,a,es13.6)', "myrank is ", myrank, "  scalar written: ", trim(filename), &
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
! なければ流入分布(十字型のtanh)を全域に与える。
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
    if (mod(myrank,2) == 1) then
      write(rank_dir, "(a, i0)") "mkdir -p data/", myrank
      call execute_command_line(trim(rank_dir), wait=.true., exitstat=ierr)
    endif
  endif

  call set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
                            nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)
  call set_profile_and_sponges(nx, ny, nz, yc, zc)
  call set_reference_volume(nx, ny, nz, dx, dy, dz)

  ! パッシブスカラー用のメトリック(セル中心間隔から)と出力用の座標
  block
    real(8) :: mx_h(nx), my_h(ny), mz_h(nz)
    integer :: ii
    do ii = 2, nx-1; mx_h(ii) = 2.d0/(dx(ii-1) + dx(ii)); enddo
    mx_h(1) = mx_h(2);  mx_h(nx) = mx_h(nx-1)
    do ii = 2, ny-1; my_h(ii) = 2.d0/(dy(ii-1) + dy(ii)); enddo
    my_h(1) = my_h(2);  my_h(ny) = my_h(ny-1)
    do ii = 2, nz-1; mz_h(ii) = 2.d0/(dz(ii-1) + dz(ii)); enddo
    mz_h(1) = mz_h(2);  mz_h(nz) = mz_h(nz-1)
    sc_mx = mx_h;  sc_my = my_h;  sc_mz = mz_h
    sc_idx = 1.d0/dx;  sc_idy = 1.d0/dy;  sc_idz = 1.d0/dz
    sc_xc = real(xc, 4);  sc_yc = real(yc, 4);  sc_zc = real(zc, 4)
  end block

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
! Reference volume shared with the solver.  The Jacobian is rebuilt with
! the same routine main.f90 uses (set_Jacobian_xy3, z width = dz(1)), so
! over_jacobian_tab matches the Jacobian passed to set_bc bit for bit.
! inv_z_update_scale(k) = dz(1)/hz(k) undoes the z-dependent factor that
! the 2-D Jacobian leaves in the shared update (see the module header).
!=====================================================================
subroutine set_reference_volume(nx, ny, nz, dx, dy, dz)
  use mod_globals, only : stretched_z_correction
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: dx(nx-1), dy(ny-1), dz(nz-1)
  real(8) :: jac(nx,ny), inv_scale_h(nz)
  integer :: k

  dx_saved = dx
  dy_saved = dy
  dz_saved = dz
  call set_Jacobian_xy3(nx, ny, nz, dx, dy, dz, jac)
  over_jacobian_tab = 1.0d0 / jac

  inv_scale_h = 1.d0
  if (stretched_z_correction) then
    do k = 2, nz-1
      inv_scale_h(k) = dz(1) / (0.5d0*(dz(k-1) + dz(k)))
    enddo
  endif
  inv_z_update_scale = inv_scale_h

  print *, "stretched_z_correction =", stretched_z_correction
  print *, "z update factor hz(k)/dz(1): min =", 1.d0/maxval(inv_scale_h(2:nz-1)), &
           " max =", 1.d0/minval(inv_scale_h(2:nz-1))
end subroutine set_reference_volume

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
! Main-domain faces of a transverse axis on [0, length], symmetric about
! the centre: uniform section of n_uniform points in the middle, and on
! each side a stretched section whose widths start at the uniform spacing
! and grow with a quintic smoothstep (stretched_widths).
!=====================================================================
subroutine partitioned_axis(n, length, uniform_length, n_uniform, coord)
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  integer, intent(in) :: n, n_uniform
  real(8), intent(in) :: length, uniform_length
  real(8), intent(out) :: coord(n)
  real(8) :: h, width(n-1), side_length
  integer :: n_side, first, j

  if (n < 2 .or. n_uniform < 2 .or. n_uniform > n) error stop 'Invalid uniform point count'
  if (.not.ieee_is_finite(length) .or. .not.ieee_is_finite(uniform_length)) &
    error stop 'Nonfinite mesh length'
  if (uniform_length <= 0.d0 .or. uniform_length > length) error stop 'Invalid uniform length'
  if (mod(n-n_uniform,2) /= 0) error stop 'Symmetric mesh requires even n-n_uniform'
  h = uniform_length / dble(n_uniform-1)

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


!=====================================================================
! Buffer faces: continue an axis from the junction face faces(i0) over
! nbuf cells in direction dir (+1 / -1).  The first width is dr and the
! growth rate ramps (cosine) from r0 to the r1 that makes the buffer
! length lbuf (calc_smooth_stretch_ratio).  On return dr is the width
! following the last buffer cell, used for the ghost faces.
!=====================================================================
subroutine buffer_faces(faces, i0, dir, nbuf, r0, lbuf, dr, r1)
  real(8), intent(inout) :: faces(:)
  integer, intent(in)    :: i0, dir, nbuf
  real(8), intent(in)    :: r0, lbuf
  real(8), intent(inout) :: dr
  real(8), intent(out)   :: r1
  real(8) :: pi, eta, r_local
  integer :: t

  pi = acos(-1.d0)
  r1 = calc_smooth_stretch_ratio(dr, r0, nbuf, lbuf)
  do t = 0, nbuf-1
    eta = dble(t) / dble(max(nbuf-1,1))
    r_local = r0 + (r1 - r0) * (1.0d0 - cos(pi*eta)) / 2.0d0
    if (dir > 0) then
      faces(i0+t+1) = faces(i0+t) + dr
    else
      faces(i0-t-1) = faces(i0-t) - dr
    endif
    dr = dr * r_local
  enddo
end subroutine buffer_faces


!=====================================================================
! Faces, cell centres and spacings.
!   x  : uniform main domain + outflow buffer
!   y/z: transverse_faces (uniform centre, smoothstep stretch, buffers)
! One ghost cell per side; the ghost width repeats the adjacent width.
!=====================================================================
subroutine set_grid_main_buffer(nx, ny, nz, Lx_main, Lx_buf,Ly_main, Ly_buf,Lz,&
  nx_main, nx_buf,ny_main, ny_buf,xc, yc, zc, dx, dy, dz)

  use mod_globals, only : Ly_uniform, ny_uniform, Lz_uniform, nz_uniform
  integer, intent(in) :: nx, ny, nz
  real(8), intent(in) :: Lx_main, Lx_buf, Ly_main, Ly_buf, Lz
  integer, intent(in) :: nx_main, nx_buf, ny_main, ny_buf
  real(8), intent(out) :: xc(nx), yc(ny), zc(nz)
  real(8), intent(out) :: dx(nx-1), dy(ny-1), dz(nz-1)

  ! x は全体格子(nxg点)で生成し、このスラブの範囲を i_offset で切り出す。
  integer, parameter :: nxg = nx_global
  real(8) :: x(nxg+1), xcg(nxg), dxg(nxg-1), y(ny+1), z(nz+1)
  real(8) :: dr, r1_x
  integer :: i, j, k, ig

  !=================================================================
  ! x: uniform main domain + outflow buffer (inflow end has no buffer)
  !=================================================================
  if (nx_main < 3 .or. Lx_main <= 0.d0) error stop 'Invalid uniform x mesh'
  if (nxg /= nx_main+nx_buf+2) error stop 'nx_global must be nx_main+nx_buf+2'
  if (i_offset + nx > nxg) error stop 'x slab exceeds the global grid'
  do i = 2, nx_main+1
    x(i) = Lx_main*dble(i-2)/dble(nx_main-1)
  enddo
  dr = x(nx_main+1) - x(nx_main)
  call buffer_faces(x, nx_main+1, +1, nx_buf, 1.d0, Lx_buf, dr, r1_x)
  x(1)     = x(2) - (x(3)-x(2))
  x(nxg)   = x(nxg-1) + dr
  x(nxg+1) = x(nxg)   + dr

  print *, "x(nx-1) - x(nx_main+1) =", x(nxg-1) - x(nx_main+1), " (目標:", Lx_buf, ")"
  print *, "xバッファ 開始伸び率=", 1.d0, " 終端伸び率=", r1_x

  !=================================================================
  ! y/z: same generator
  !=================================================================
  call transverse_faces(ny, ny_main, ny_buf, Ly_main, Ly_buf, Ly_uniform, ny_uniform, y)
  call transverse_faces(nz, nz_main, nz_buf, Lz_main, Lz_buf, Lz_uniform, nz_uniform, z)

  !=================================================================
  ! Cell centres and centre-to-centre spacings
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
  do j = 1, ny-1; dy(j) = yc(j+1)-yc(j); enddo
  do k = 1, nz-1; dz(k) = zc(k+1)-zc(k); enddo

  !=================================================================
  ! 確認出力
  !=================================================================
  print *, "=== 格子生成確認 ==="
  print *, "x: 主計算領域 ", x(2), "~", x(nx_main+1), " バッファ ~", x(nxg-1)
  print *, "dx_main(一様) =", x(3)-x(2), " dx_buf_max =", x(nxg-1)-x(nxg-2)
  if (npx > 1) print *, "x slab", my_slab, " of", npx, ": global i =", i_offset+1, "~", i_offset+nx, &
                        " x =", xc(1), "~", xc(nx)
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
  !   主領域 = 中心一様区間 + 外側伸長(partitioned_axis)
  !   バッファ = 接続点の格子幅と伸び率を引き継ぐコサインランプ伸長(buffer_faces)
  !=====================================================================
  subroutine transverse_faces(n, nmain, nbuf, lmain, lbuf, luni, nuni, faces)
    integer, intent(in) :: n, nmain, nbuf, nuni
    real(8), intent(in) :: lmain, lbuf, luni
    real(8), intent(out) :: faces(n+1)
    real(8) :: dr, r0, r1
    if (nmain < 4 .or. nbuf < 2 .or. n /= nmain+2*nbuf+2) error stop 'Invalid transverse mesh'

    ! 主領域: 0..lmain で生成し、下バッファ長さ分だけ平行移動
    call partitioned_axis(nmain, lmain, luni, nuni, faces(nbuf+2:nbuf+nmain+1))
    faces(nbuf+2:nbuf+nmain+1) = faces(nbuf+2:nbuf+nmain+1) + lbuf

    ! 下バッファ: 接続点の格子幅と、その外側への局所伸び率を引き継ぐ
    dr = faces(nbuf+3) - faces(nbuf+2)
    r0 = dr / (faces(nbuf+4) - faces(nbuf+3))
    call buffer_faces(faces, nbuf+2, -1, nbuf, r0, lbuf, dr, r1)

    ! 上バッファ
    dr = faces(nbuf+nmain+1) - faces(nbuf+nmain)
    r0 = dr / (faces(nbuf+nmain) - faces(nbuf+nmain-1))
    call buffer_faces(faces, nbuf+nmain+1, +1, nbuf, r0, lbuf, dr, r1)

    ! ゴースト
    faces(1)   = faces(2) - (faces(3)-faces(2))
    faces(n)   = faces(n-1) + dr
    faces(n+1) = faces(n)   + dr
  end subroutine transverse_faces


  !=====================================================================
  ! Cell-centre cross profile (same cross_profile as set_init) and the
  ! sponge weights, cached on the device for set_bc.
  !=====================================================================
  subroutine set_profile_and_sponges(nx, ny, nz, yc, zc)
    use mod_globals, only : u1, u2
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: yc(ny), zc(nz)
    real(8) :: rho_h(ny,nz), u_h(ny,nz), env_h(ny,nz)
    real(8) :: sx_h(nx), sy_h(ny), sz_h(nz), temp_h
    integer :: ii, jj, kk
    do kk = 1, nz
      do jj = 1, ny
        call cross_profile(yc(jj), zc(kk), u_h(jj,kk), temp_h, rho_h(jj,kk), env_h(jj,kk))
      enddo
    enddo
    do ii = 1, nx
      sx_h(ii) = sponge_weight(dble(ii+i_offset-(nx_global-1-nsp_x))/dble(nsp_x), sigma_max_x)
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
    ! パッシブスカラーの流入分布: 速度の十字型分布と同じ形。u1側で1、u2側で0。
    sc_in = (u_h - u2)/(u1 - u2)
    sigma_x_1d = sx_h
    sigma_y_1d = sy_h
    sigma_z_1d = sz_h
  end subroutine set_profile_and_sponges

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
    use mod_globals, only : p, stretched_z_correction
    integer, intent(in) :: myrank, nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    ! Same layout as main.f90 / pre_calc: Q(i,j,k,variable)
    real(8), intent(out) :: Q(nx,ny,nz,5)
    integer :: i, j, k, l
    real(8) :: vel, temp, rho, env
    real(8), allocatable :: jac(:,:), qj_h(:,:,:)
    do k = 1, nz
      do j = 1, ny
        call cross_profile(y(j), z(k), vel, temp, rho, env)
        do i = 1, nx
          Q(i,j,k,1) = rho
          Q(i,j,k,2) = rho*vel
          Q(i,j,k,3) = 0.d0
          Q(i,j,k,4) = 0.d0
          Q(i,j,k,5) = p/(gamma-1.d0)+0.5d0*rho*vel**2
        enddo
      enddo
    enddo

    ! Q0 (state at the start of the first step) for the stretched-z
    ! correction in set_bc: QJ exactly as pre_calc forms it, Q / Jacobian.
    if (stretched_z_correction) then
      allocate(jac(nx,ny), qj_h(nx,ny,nz))
      call set_Jacobian_xy3(nx, ny, nz, dx_saved, dy_saved, dz_saved, jac)
      allocate(q0_1(nx,ny,nz), q0_2(nx,ny,nz), q0_3(nx,ny,nz), q0_4(nx,ny,nz), q0_5(nx,ny,nz))
      do l = 1, 5
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              qj_h(i,j,k) = Q(i,j,k,l) / jac(i,j)
        enddo;enddo;enddo
        select case (l)
        case (1); q0_1 = qj_h
        case (2); q0_2 = qj_h
        case (3); q0_3 = qj_h
        case (4); q0_4 = qj_h
        case (5); q0_5 = qj_h
        end select
      enddo
      deallocate(jac, qj_h)
    endif
  end subroutine set_init

subroutine set_bc(myrank, nx, ny, nz, Jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)
  use mod_constant, only : id_accuracy
  use mod_globals, only : u1, rho1, u2, rho2, p, amp, dt, gamma, Ly, Lz, Lx,T1,T2,delta_bl, step_offset, &
                          inlet_fluctuation_rms, inlet_random_seed, inlet_rk_stages, stretched_z_correction, scalar_on
  use set_coordinate
  integer, intent(in), value            :: myrank, nx, ny, nz
  real(8), intent(in), device           :: jacobian(nx,ny)
  real(8), intent(inout), device        :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
  !real(8), intent(in), device, optional :: Qre(:)
  real(8) :: u, v, w, u_init, rho_init
  integer :: i, j, k, l, i_target, inlet_frame, rk_stage

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
    ! RK4 stage of this call (1..4): RungeKutta calls set_bc once per stage.
    rk_stage = mod(inlet_bc_calls, inlet_rk_stages) + 1

    !===================================================================
    ! Stretched z: Q = Q0 + (Q - Q0)*dz(1)/hz(k) on the points calc_step
    ! updates, before any boundary or sponge treatment (module header).
    !===================================================================
    if (stretched_z_correction) then
      if (inlet_rk_stages /= 4) error stop "stretched_z_correction requires RK=4 (inlet_rk_stages=4)"
      if (.not. allocated(q0_1)) error stop "stretched_z_correction requires set_init (RESTART is not supported)"
      ! Separate routine: with Q0 as dummy arguments the compiler may assume
      ! no aliasing with Q and launches a full grid (inline it ran as 1 block).
      call apply_stretched_z(nx, ny, nz, inv_z_update_scale, q0_1, q0_2, q0_3, q0_4, q0_5, &
                             Q_1, Q_2, Q_3, Q_4, Q_5)
    endif

    inlet_frame = inlet_bc_calls / inlet_rk_stages
    inlet_bc_calls = inlet_bc_calls + 1

    if (my_slab == 0) then
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

        Q_1(i,j,k) = rho_init*over_jacobian_tab(i,j)
        Q_2(i,j,k) = rho_init*u*over_jacobian_tab(i,j)
        Q_3(i,j,k) = rho_init*v*over_jacobian_tab(i,j)
        Q_4(i,j,k) = rho_init*w*over_jacobian_tab(i,j)
        Q_5(i,j,k) = (p/(gamma-1.d0) + 0.5d0*rho_init*(u**2+v**2+w**2))*over_jacobian_tab(i,j)
      enddo;enddo
    endif

    ! x流出スポンジ(目標面の取得を含む)は最終スラブだけが持つ。
    if (my_slab == npx-1) then
    !===================================================================
    ! Capture the full yz plane before damping.  A z average would erase
    ! the reversed high/low streams.  Keep each (j,k) target independently.
    !===================================================================
    i_target = nx - 1 - nsp_x
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        rho_now = Q_1(i_target,j,k)*jacobian(i_target,j)
        u_now = Q_2(i_target,j,k)/Q_1(i_target,j,k)
        v_now = Q_3(i_target,j,k)/Q_1(i_target,j,k)
        w_now = Q_4(i_target,j,k)/Q_1(i_target,j,k)
        p_now = (gamma-1.d0)*(Q_5(i_target,j,k)*jacobian(i_target,j) &
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
            rho_now = Q_1(i,j,k) * jacobian(i,j)
            u_now   = Q_2(i,j,k) / Q_1(i,j,k)
            v_now   = Q_3(i,j,k) / Q_1(i,j,k)
            w_now   = Q_4(i,j,k) / Q_1(i,j,k)
            p_now   = (gamma-1.d0) * (Q_5(i,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))

            rho_blend = (1.d0-sigma_x)*rho_now + sigma_x*rho_target_x(j,k)
            u_blend   = (1.d0-sigma_x)*u_now   + sigma_x*u_target_x(j,k)
            v_blend   = (1.d0-sigma_x)*v_now   + sigma_x*v_target_x(j,k)
            w_blend   = (1.d0-sigma_x)*w_now   + sigma_x*w_target_x(j,k)
            p_blend   = (1.d0-sigma_x)*p_now   + sigma_x*p_target_x(j,k)
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
    ! Apply y/z sponge weights together to preserve symmetry at corners.
    ! Both sides relax to the same yz-dependent cross profile as the inlet.
    !===================================================================
    !$cuf kernel do(3)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          sigma_y = 1.d0-(1.d0-sigma_y_1d(j))*(1.d0-sigma_z_1d(k))
          if (sigma_y > 0.d0) then
            rho_now = Q_1(i,j,k) * jacobian(i,j)
            u_now   = Q_2(i,j,k) / Q_1(i,j,k)
            v_now   = Q_3(i,j,k) / Q_1(i,j,k)
            w_now   = Q_4(i,j,k) / Q_1(i,j,k)
            p_now   = (gamma-1.d0) * (Q_5(i,j,k)*jacobian(i,j) &
                      - 0.5d0*rho_now*(u_now**2+v_now**2+w_now**2))
            rho_init = rho_target(j,k)
            u_init   = u_target(j,k)
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
    ! x = nx 流出境界(1点ゴースト): 物理量で1次外挿してJ補正(最終スラブのみ)
    !===================================================================
    if (my_slab == npx-1) then
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do l = 1, 5
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
        enddo;enddo;enddo 
    endif
    !===================================================================
    ! y方向境界(j=1, j=ny): 0次外挿。物理保存量をコピーしてJ補正。
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1, nz
      do i = 1, nx
          Q_1(i,ny,k) = Q_1(i,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)
          Q_2(i,ny,k) = Q_2(i,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)
          Q_3(i,ny,k) = Q_3(i,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)
          Q_4(i,ny,k) = Q_4(i,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)
          Q_5(i,ny,k) = Q_5(i,ny-1,k)*jacobian(i,ny-1)*over_jacobian_tab(i,ny)

          Q_1(i,1,k) = Q_1(i,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
          Q_2(i,1,k) = Q_2(i,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
          Q_3(i,1,k) = Q_3(i,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
          Q_4(i,1,k) = Q_4(i,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
          Q_5(i,1,k) = Q_5(i,2,k)*jacobian(i,2)*over_jacobian_tab(i,1)
        enddo;enddo
    !===================================================================
    ! z: zeroth-order extrapolation copies the adjacent interior state.
    ! Jacobian is independent of z, so its factors cancel here.
    !===================================================================
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
          Q_1(i,j,1) = Q_1(i,j,2)
          Q_2(i,j,1) = Q_2(i,j,2)
          Q_3(i,j,1) = Q_3(i,j,2)
          Q_4(i,j,1) = Q_4(i,j,2)
          Q_5(i,j,1) = Q_5(i,j,2)

          Q_1(i,j,nz) = Q_1(i,j,nz-1)
          Q_2(i,j,nz) = Q_2(i,j,nz-1)
          Q_3(i,j,nz) = Q_3(i,j,nz-1)
          Q_4(i,j,nz) = Q_4(i,j,nz-1)
          Q_5(i,j,nz) = Q_5(i,j,nz-1)
      enddo
    enddo

    ! x方向の接続面: 隣スラブの内部3面をゴーストに受け取る(npx=1なら何もしない)。
    call exchange_x(myrank, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)

    ! After the last stage this array is Q^(n+1), i.e. Q0 of the next step.
    ! (交換の後に取るので、ゴースト面のQ0も隣スラブの値になる。)
    if (stretched_z_correction .and. rk_stage == inlet_rk_stages) then
      q0_1 = Q_1; q0_2 = Q_2; q0_3 = Q_3; q0_4 = Q_4; q0_5 = Q_5
    endif

    ! パッシブスカラー(混合分率)。流れ場には一切影響しない。
    if (scalar_on) call sc_step(myrank, nx, ny, nz, jacobian, Q_1, Q_2, Q_3, Q_4, Q_5)

  end subroutine set_bc

  !=====================================================================
  ! x方向の袖交換(非周期: 両端のスラブは外側に相手を持たない)。
  ! 右隣へ内部の最後の3面(nx-5..nx-3)、左隣へ最初の3面(4..6)を送り、
  ! ゴースト(1..3, nx-2..nx)に受け取る。y/zのゴーストも含む全(j,k)。
  ! GPU上で1次元バッファに詰め、ホスト経由でMPI_SENDRECVする。
  !=====================================================================
  subroutine exchange_x(myrank, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
    use mpi
    integer, intent(in), value     :: myrank, nx, ny, nz
    real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer :: left, right, n, ierr
    integer :: istat(MPI_STATUS_SIZE)

    if (npx == 1) return
    n = 3*ny*nz*5
    call alloc_exchange_x(ny, nz)
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
  !> Q/J を真のJ/ローカルのJ の比 fac で換算する(module header参照)。
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

  !> Q = Q0 + (Q - Q0)*fz(k) on the points calc_step updates (see module header).
  subroutine apply_stretched_z(nx, ny, nz, fz, P_1, P_2, P_3, P_4, P_5, Q_1, Q_2, Q_3, Q_4, Q_5)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: fz(nz)
    real(8), intent(in), device    :: P_1(nx,ny,nz), P_2(nx,ny,nz), P_3(nx,ny,nz), P_4(nx,ny,nz), P_5(nx,ny,nz)
    real(8), intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer :: i, j, k
    real(8) :: f
    !$cuf kernel do(3)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          f = fz(k)
          Q_1(i,j,k) = P_1(i,j,k) + (Q_1(i,j,k) - P_1(i,j,k))*f
          Q_2(i,j,k) = P_2(i,j,k) + (Q_2(i,j,k) - P_2(i,j,k))*f
          Q_3(i,j,k) = P_3(i,j,k) + (Q_3(i,j,k) - P_3(i,j,k))*f
          Q_4(i,j,k) = P_4(i,j,k) + (Q_4(i,j,k) - P_4(i,j,k))*f
          Q_5(i,j,k) = P_5(i,j,k) + (Q_5(i,j,k) - P_5(i,j,k))*f
    enddo;enddo;enddo
  end subroutine apply_stretched_z

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
