module calc_para
  use mpi
  use cudafor
  implicit none
  interface exchange
    module procedure exchange_cyclic, exchange_rescale
  end interface exchange
contains
  subroutine flatten(nx, ny, nz, overlap, Q, Q1d_left, Q1d_right)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(out), device :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q1d_left(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)  = Q(l,overlap+i,j+1,k+3)
            Q1d_right(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l) = Q(l,nx-2*overlap+i,j+1,k+3)
    enddo;enddo;enddo;enddo
  end subroutine flatten

  subroutine flatten_left(nx, ny, nz, overlap, Q, Q1d_left)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q1d_left(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)  = Q(l,overlap+i,j+1,k+3)
    enddo;enddo;enddo;enddo
  end subroutine flatten_left

  subroutine flatten_right(nx, ny, nz, overlap, Q, Q1d_right)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q1d_right(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l) = Q(l,nx-2*overlap+i,j+1,k+3)
    enddo;enddo;enddo;enddo
  end subroutine flatten_right
  
  subroutine flatten_rescale(nx, ny, nz, nre, overlap, Q, Q1d_right)
    integer, intent(in), value   :: nx, ny, nz, nre, overlap
    real(8), intent(in), device  :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q1d_right(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l) = Q(l,nre-2*overlap+i,j+1,k+3)
    enddo;enddo;enddo;enddo
  end subroutine flatten_rescale
  
  subroutine reconstruct(nx, ny, nz, overlap, Q1d_left, Q1d_right, Q)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(in), device  :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(out), device :: Q(5,nx,ny,nz)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q(l,i,j+1,k+3)            =  Q1d_left(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)
            Q(l,nx-overlap+i,j+1,k+3) = Q1d_right(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct

  subroutine reconstruct_left(nx, ny, nz, overlap, Q1d_left, Q)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q1d_left(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(out), device :: Q(5,nx,ny,nz)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q(l,i,j+1,k+3) = Q1d_left(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_left
  
  subroutine reconstruct_right(nx, ny, nz, overlap, Q1d_right, Q)
    integer, intent(in), value   :: nx, ny, nz, overlap
    real(8), intent(in), device  :: Q1d_right(overlap*(ny-2)*(nz-6)*5)
    real(8), intent(out), device :: Q(5,nx,ny,nz)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          do l = 1, 5
            Q(l,nx-overlap+i,j+1,k+3) = Q1d_right(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_right

  subroutine reconstruct_sbli_inlet(nx, ny1, ny2, nz, overlap, Q1d, Q)
    integer, intent(in), value     :: nx, ny1, ny2, nz, overlap
    real(8), intent(in), device    :: Q1d(overlap*(ny1-2)*(nz-6)*5)
    real(8), intent(inout), device :: Q(5,nx,ny2,nz)
    integer i, j, k, l, ni, nj, nk
    ni = overlap
    nj = ny1-2
    nk = nz-6
    !$cuf kernel do(2)<<<*,(32,4)>>>
    do k = 1, nk
      do j = 1, nj
        do i = 1, overlap
          do l = 1, 5
            Q(l,i,j+1,k+3) = Q1d(ni*nj*5*(k-1)+ni*5*(j-1)+5*(i-1)+l)
    enddo;enddo;enddo;enddo
  end subroutine reconstruct_sbli_inlet

  subroutine exchange_cyclic(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
    integer(kind=2), intent(in), value :: id_rescale
    integer, intent(in), value         :: myrank, nranks, overlap, nx, ny, nz
    real(8), intent(inout), device     :: QJ(5,nx,ny,nz) ! Q / Jacobian
    integer rank1, rank2, ierr, ireq4(4), istat(MPI_STATUS_SIZE), istat4(MPI_STATUS_SIZE,4)
    real(8), dimension(overlap*(ny-2)*(nz-6)*5)         :: Qs_left,   Qs_right,   Qr_left,   Qr_right
    real(8), dimension(overlap*(ny-2)*(nz-6)*5), device :: Qs1d_left, Qs1d_right, Qr1d_left, Qr1d_right
    integer j, k, ni, nj, nk

    if (2 <= myrank .and. myrank <= nranks-4) then
      rank1 = myrank-2
      rank2 = myrank+2
    elseif (myrank == 0 .and. 4 <= nranks) then
      rank1 = nranks-2
      rank2 = myrank+2
    elseif (myrank == nranks-2 .and. 4 <= nranks) then
      rank1 = myrank-2
      rank2 = 0
    endif

    call flatten(nx, ny, nz, overlap, QJ, Qs1d_left, Qs1d_right)

    Qs_left    = Qs1d_left
    call MPI_SENDRECV(Qs_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, &
                      Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
    Qr1d_right = Qr_right

    Qs_right   = Qs1d_right
    call MPI_SENDRECV(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, &
                      Qr_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
    Qr1d_left  = Qr_left

    call reconstruct(nx, ny, nz, overlap, Qr1d_left, Qr1d_right, QJ)
  end subroutine exchange_cyclic

  subroutine exchange_rescale(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
    integer(kind=4), intent(in), value :: id_rescale
    integer, intent(in), value         :: myrank, nranks, overlap, nx, ny, nz
    real(8), intent(inout), device     :: QJ(5,nx,ny,nz) ! Q / Jacobian
    integer rank1, rank2, stat, ierr, ireq4(4), istat(MPI_STATUS_SIZE), istat4(MPI_STATUS_SIZE,4)
    real(8), dimension(overlap*(ny-2)*(nz-6)*5)         :: Qs_left,   Qs_right,   Qr_left,   Qr_right
    real(8), dimension(overlap*(ny-2)*(nz-6)*5), device :: Qs1d_left, Qs1d_right, Qr1d_left, Qr1d_right
    integer j, k, ni, nj, nk

    if (2 <= myrank .and. myrank <= nranks-4) then
      rank1 = myrank-2
      rank2 = myrank+2

      call flatten(nx, ny, nz, overlap, QJ, Qs1d_left, Qs1d_right)

      stat = cudaMemcpy(Qs_left, Qs1d_left, 5*overlap*(ny-2)*(nz-6), cudaMemcpyDeviceToHost)
      call MPI_SENDRECV(Qs_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, &
                        Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
      stat = cudaMemcpyAsync(Qr1d_right, Qr_right, 5*overlap*(ny-2)*(nz-6), cudaMemcpyHostToDevice, 1)

      stat = cudaMemcpy(Qs_right, Qs1d_right, 5*overlap*(ny-2)*(nz-6), cudaMemcpyDeviceToHost)
      call MPI_SENDRECV(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, &
                        Qr_left,  5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
      stat = cudaMemcpyAsync(Qr1d_left,  Qr_left,  5*overlap*(ny-2)*(nz-6), cudaMemcpyHostToDevice, 2)
      stat = cudaDeviceSynchronize()

      call reconstruct(nx, ny, nz, overlap, Qr1d_left, Qr1d_right, QJ)
    elseif (myrank == 0 .and. 4 <= nranks) then
      rank2 = myrank+2

      call flatten_right(nx, ny, nz, overlap, QJ, Qs1d_right)

      call MPI_RECV(Qr_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, istat, ierr)
      stat = cudaMemcpyAsync(Qr1d_right, Qr_right, 5*overlap*(ny-2)*(nz-6), cudaMemcpyHostToDevice, 1)

      stat = cudaMemcpy(Qs_right, Qs1d_right, 5*overlap*(ny-2)*(nz-6), cudaMemcpyDeviceToHost)
      call MPI_SEND(Qs_right, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank2, 0, MPI_COMM_WORLD, ierr)

      stat = cudaDeviceSynchronize()
      call reconstruct_right(nx, ny, nz, overlap, Qr1d_right, QJ)
    elseif (myrank == nranks-2 .and. 4 <= nranks) then
      rank1 = myrank-2

      call flatten_left(nx, ny, nz, overlap, QJ, Qs1d_left)

      stat = cudaMemcpy(Qs_left, Qs1d_left, 5*overlap*(ny-2)*(nz-6), cudaMemcpyDeviceToHost)
      call MPI_SEND(Qs_left, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, ierr)

      call MPI_RECV(Qr_left, 5*overlap*(ny-2)*(nz-6), MPI_REAL8, rank1, 0, MPI_COMM_WORLD, istat, ierr)
      stat = cudaMemcpy(Qr1d_left, Qr_left, 5*overlap*(ny-2)*(nz-6), cudaMemcpyHostToDevice)

      call reconstruct_left(nx, ny, nz, overlap, Qr1d_left, QJ)
    endif
  end subroutine exchange_rescale
end module calc_para

