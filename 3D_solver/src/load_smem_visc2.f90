module load_smem_visc2
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  implicit none
  private
  public load_smem_visc2_x, load_smem_visc2_y, load_smem_visc2_z, load_smem_visc2_z_koff
contains
  attributes(device) subroutine load_smem_visc2_x(it, jt, kt, j, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: idx           !< index for shared memory
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), intent(inout) :: u(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: w(0:sx*sy*sz-1) !< attribute(shared)
    integer i_base, ii, i, idx_l, offset_yz
    logical :: jk_in_range
    i_base    = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    jk_in_range = (j <= ny .and. k <= nz)
    do ii = it, threadsEv%x+1, blockDim%x
      i = i_base + ii
      idx_l = (ii-1) + offset_yz
      if (1 <= i .and. i <= nx .and. jk_in_range) then
        u(idx_l) = Q_2(i,j,k)
        v(idx_l) = Q_3(i,j,k)
        w(idx_l) = Q_4(i,j,k)
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc2_x


  attributes(device) subroutine load_smem_visc2_y(it, jt, kt, i, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: idx           !< index for shared memory
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), intent(inout) :: u(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: w(0:sx*sy*sz-1) !< attribute(shared)
    integer j_base, jj, j, idx_l, offset_xz
    logical :: ik_in_range
    j_base    = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    ik_in_range = (i <= nx .and. k <= nz)
    do jj = jt, threadsFv%y+1, blockDim%y
      j = j_base + jj
      idx_l = (jj-1) + offset_xz
      if (ik_in_range .and. 1 <= j .and. j <= ny) then
        u(idx_l) = Q_2(i,j,k)
        v(idx_l) = Q_3(i,j,k)
        w(idx_l) = Q_4(i,j,k)
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc2_y


  attributes(device) subroutine load_smem_visc2_z(it, jt, kt, i, j, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: idx           !< index for shared memory
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), intent(inout) :: u(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: w(0:sx*sy*sz-1) !< attribute(shared)
    integer k_base, kk, k, idx_l, offset_xy
    logical :: ij_in_range
    k_base    = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    ij_in_range = (i <= nx .and. j <= ny)
    do kk = kt, threadsGv%z+1, blockDim%z
      k = k_base + kk
      idx_l = (kk-1) + offset_xy
      if (ij_in_range .and. 1 <= k .and. k <= nz) then
        u(idx_l) = Q_2(i,j,k)
        v(idx_l) = Q_3(i,j,k)
        w(idx_l) = Q_4(i,j,k)
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc2_z


  !> Like load_smem_visc2_z but k_base is shifted by k_lo-1 for koff kernel launches
  attributes(device) subroutine load_smem_visc2_z_koff(it, jt, kt, i, j, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w, k_lo)
    integer, intent(in), value              :: it, jt, kt
    integer, intent(in), value              :: i, j, idx
    integer, intent(in), value              :: nx, ny, nz
    integer, intent(in), value              :: k_lo
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), intent(inout) :: u(0:sx*sy*sz-1)
    real(8), intent(inout) :: v(0:sx*sy*sz-1)
    real(8), intent(inout) :: w(0:sx*sy*sz-1)
    integer k_base, kk, k, idx_l, offset_xy
    logical :: ij_in_range
    k_base    = (blockIdx%z-1)*blockDim%z + k_lo - 1
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    ij_in_range = (i <= nx .and. j <= ny)
    do kk = kt, threadsGv%z+1, blockDim%z
      k = k_base + kk
      idx_l = (kk-1) + offset_xy
      if (ij_in_range .and. 1 <= k .and. k <= nz) then
        u(idx_l) = Q_2(i,j,k)
        v(idx_l) = Q_3(i,j,k)
        w(idx_l) = Q_4(i,j,k)
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc2_z_koff
end module load_smem_visc2
