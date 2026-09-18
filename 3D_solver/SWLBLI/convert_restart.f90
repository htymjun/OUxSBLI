program convert_restart
  implicit none
  integer, parameter :: nx = 257, ny = 257
  integer, parameter :: nz_old = 9,  nz_int_old = nz_old - 6   ! = 3
  integer, parameter :: nz_new = 140, nz_int_new = nz_new - 6   ! = 12
  real(8), allocatable :: Q_old(:,:,:,:), Q_new(:,:,:,:)
  real(8), allocatable :: Q_interior(:,:,:,:)
  integer i, j, k, m, k_old
  real(8) zf

  allocate(Q_old(nx,ny,nz_old,5))
  allocate(Q_interior(nx,ny,nz_int_old,5))
  allocate(Q_new(nx,ny,nz_new,5))

  ! --- 旧ファイルを読み込み ---
  open(10, file="recal/Q00001.dat", form="unformatted", access="stream", status="old")
  read(10) Q_old
  close(10)

  ! --- ゴーストを除いた内部平面だけ取り出す (k=4..nz_old-3) ---
  Q_interior(:,:,1:nz_int_old,:) = Q_old(:,:,4:nz_old-3,:)

  ! --- 最近傍複製でnz_int_newへアップサンプリング ---
  ! nz_int_new / nz_int_old が整数倍でなくても、各新平面が対応する
  ! 旧平面を「最も近い」ものとして選ぶので端数があっても動作する
  do k = 1, nz_int_new
    zf = dble(k-1) * dble(nz_int_old) / dble(nz_int_new)   ! 0 <= zf < nz_int_old
    k_old = mod(nint(zf), nz_int_old) + 1                   ! 最近傍のインデックス（周期的）
    do m = 1, 5
      do j = 1, ny
        do i = 1, nx
          Q_new(i,j,k+3,m) = Q_interior(i,j,k_old,m)
        enddo
      enddo
    enddo
  enddo

  ! --- ゴースト平面を埋める（set_bc_cyclic_zと同じロジック） ---
  Q_new(:,:,1,:) = Q_new(:,:,nz_new-5,:)
  Q_new(:,:,2,:) = Q_new(:,:,nz_new-4,:)
  Q_new(:,:,3,:) = Q_new(:,:,nz_new-3,:)
  Q_new(:,:,nz_new-2,:) = Q_new(:,:,4,:)
  Q_new(:,:,nz_new-1,:) = Q_new(:,:,5,:)
  Q_new(:,:,nz_new,:)   = Q_new(:,:,6,:)

  ! --- 新しいリスタートファイルとして保存 ---
  open(20, file="recal/Q00001_new.dat", form="unformatted", access="stream", status="replace")
  write(20) Q_new
  close(20)
end program convert_restart