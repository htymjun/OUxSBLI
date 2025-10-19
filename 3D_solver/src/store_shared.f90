module store_shared
  use mod_globals, only : threads, blocks
  implicit none
contains
  attributes(device) subroutine store_shared_visc_2nd(nx, ny, nz, it, jt, kt, Q, mud, u, v, w, mu)
    integer, intent(in), value :: nx, ny, nz, it, jt, kt
    real(8), intent(in)        :: Q(5,nx,ny,nz), mud(nx,ny,nz)
    real(8), intent(inout)     ::  u(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8), intent(inout)     ::  v(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8), intent(inout)     ::  w(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8), intent(inout)     :: mu(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    integer i, j, k, ii, jj, kk, i_base, j_base, k_base
    i_base = (blockIdx%x-1)*blockDim%x + 1
    j_base = (blockIdx%y-1)*blockDim%y + 1
    k_base = (blockIdx%z-1)*blockDim%z + 1
    do kk = kt-1, threads%z+1, blockDim%z
      do jj = jt-1, threads%y+1, blockDim%y
        do ii = it-1, threads%x+1, blockDim%x
          i = i_base + ii
          j = j_base + jj
          k = k_base + kk
          if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
             u(ii,jj,kk) = Q(2,i,j,k)
             v(ii,jj,kk) = Q(3,i,j,k)
             w(ii,jj,kk) = Q(4,i,j,k)
            mu(ii,jj,kk) = mud(i,j,k)
          endif
    enddo;enddo;enddo
    call syncthreads()
  end subroutine store_shared_visc_2nd


  attributes(device) subroutine store_shared_visc_4th(nx, ny, nz, it, jt, kt, Q, mud, u, v, w, mu)
    integer, intent(in), value :: nx, ny, nz, it, jt, kt
    real(8), intent(in)        :: Q(5,nx,ny,nz), mud(nx,ny,nz)
    real(8), intent(inout)     ::  u(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8), intent(inout)     ::  v(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8), intent(inout)     ::  w(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8), intent(inout)     :: mu(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    integer i, j, k, ii, jj, kk, i_base, j_base, k_base
    i_base = (blockIdx%x-1)*blockDim%x + 1
    j_base = (blockIdx%y-1)*blockDim%y + 1
    k_base = (blockIdx%z-1)*blockDim%z + 1
    do kk = kt-2, threads%z+2, blockDim%z
      do jj = jt-2, threads%y+2, blockDim%y
        do ii = it-2, threads%x+2, blockDim%x
          i = i_base + ii
          j = j_base + jj
          k = k_base + kk
          if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
             u(ii,jj,kk) = Q(2,i,j,k)
             v(ii,jj,kk) = Q(3,i,j,k)
             w(ii,jj,kk) = Q(4,i,j,k)
            mu(ii,jj,kk) = mud(i,j,k)
          endif
    enddo;enddo;enddo
    call syncthreads()
  end subroutine store_shared_visc_4th
end module store_shared

