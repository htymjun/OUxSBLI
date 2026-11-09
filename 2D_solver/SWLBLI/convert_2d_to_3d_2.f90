program convert_2d_to_3d_2
  implicit none
  ! DNS
  integer, parameter :: nx = 2700
  integer, parameter :: ny = 480
  integer, parameter :: nz_new = 35
  integer, parameter :: n_compute_ranks = 1   ! ★ 3D計算で使う計算ランク数（GPU数）に合わせる

  ! --- 擾乱パラメータ ---
  real(8), parameter :: eps = 1.d-14        ! 擾乱振幅

  real(8), allocatable :: Q2d(:,:,:), Q3d(:,:,:,:)
  real(8), allocatable :: noise(:,:,:,:)
  integer k, irank
  character(len=40) filename

  allocate(Q2d(nx,ny,4))
  allocate(Q3d(nx,ny,nz_new,5))
  allocate(noise(nx,ny,nz_new,5))

  ! --- 2次元のリスタートファイルを読み込み（全ランク共通の元データ） ---
  open(10, file="recal/Q00001.dat", form="unformatted", access="stream", status="old")
  read(10) Q2d
  close(10)
  !q = REAL(Q2d, kind=4)

  ! --- 計算ランクごとにファイルを生成 ---
  do irank = 0, n_compute_ranks - 1

    ! 全z平面（ゴースト含む）に同じ値を複製し、rho*w=0を挿入
    do k = 1, nz_new
      Q3d(:,:,k,1) = Q2d(:,:,1)   ! rho
      Q3d(:,:,k,2) = Q2d(:,:,2)   ! rho*u
      Q3d(:,:,k,3) = Q2d(:,:,3)   ! rho*v
      Q3d(:,:,k,4) = 0.d0         ! rho*w = 0
      Q3d(:,:,k,5) = Q2d(:,:,4)   ! E
    enddo

    ! --- 全保存量に丸め誤差程度のランダム擾乱を追加 ---
    ! スパン方向の一様性を崩すため、複製後に全領域で点ごとに独立な乱数を与える
    call random_number(noise)
    Q3d = Q3d + eps * noise

    ! --- このランク用のリスタートファイルとして保存 ---
    write(filename, "(a, i5.5, a)") "recal/Q_3d", irank+1, ".dat"
    open(20, file=filename, form="unformatted", access="stream", status="replace")
    write(20) Q3d
    close(20)
  enddo
end program convert_2d_to_3d_2
