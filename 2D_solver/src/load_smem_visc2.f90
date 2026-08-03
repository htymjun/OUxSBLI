module load_smem_visc2
  use mod_globals, only : threadsEv, threadsFv
  implicit none
  private
  public load_smem_visc2_x, load_smem_visc2_y
contains
  attributes(device) subroutine load_smem_visc2_x(it, jt, j, idx, nx, ny, Q_2, Q_3, u, v)
    integer, intent(in), value              :: it         !< local idx for x direction
    integer, intent(in), value              :: jt         !< local idx for y direction
    integer, intent(in), value              :: j          !< global idx for y direction
    integer, intent(in), value              :: idx        !< index for shared memory
    integer, intent(in), value              :: nx         !< number of grid points in x direction
    integer, intent(in), value              :: ny         !< number of grid points in y direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny) !< u
    real(8), intent(in), device, contiguous :: Q_3(nx,ny) !< v
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    real(8), intent(inout) :: u(0:sx*sy-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy-1) !< attribute(shared)
    integer i_base, ii, i, idx_l, offset_y
    i_base   = (blockIdx%x-1)*blockDim%x
    offset_y = (jt-1)*sx
    do ii = it, threadsEv%x+1, blockDim%x
      i = i_base + ii
      idx_l = (ii-1) + offset_y
      if (1 <= i .and. i <= nx .and. j <= ny) then
        u(idx_l) = Q_2(i,j)
        v(idx_l) = Q_3(i,j)
      else
        u(idx_l) = 0.d0
        v(idx_l) = 0.d0
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc2_x


  attributes(device) subroutine load_smem_visc2_y(it, jt, i, idx, nx, ny, Q_2, Q_3, u, v)
    integer, intent(in), value              :: it         !< local idx for x direction
    integer, intent(in), value              :: jt         !< local idx for y direction
    integer, intent(in), value              :: i          !< global idx for x direction
    integer, intent(in), value              :: idx        !< index for shared memory
    integer, intent(in), value              :: nx         !< number of grid points in x direction
    integer, intent(in), value              :: ny         !< number of grid points in y direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny) !< u
    real(8), intent(in), device, contiguous :: Q_3(nx,ny) !< v
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    real(8), intent(inout) :: u(0:sx*sy-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy-1) !< attribute(shared)
    integer j_base, jj, j, idx_l, offset_x
    j_base   = (blockIdx%y-1)*blockDim%y
    offset_x = (it-1)*sy
    do jj = jt, threadsFv%y+1, blockDim%y
      j = j_base + jj
      idx_l = (jj-1) + offset_x
      if (i <= nx .and. 1 <= j .and. j <= ny) then
        u(idx_l) = Q_2(i,j)
        v(idx_l) = Q_3(i,j)
      else
        u(idx_l) = 0.d0
        v(idx_l) = 0.d0
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc2_y
end module load_smem_visc2
