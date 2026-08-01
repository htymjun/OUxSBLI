module load_smem_visc2_curv
  use wmma
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  implicit none
  private
  public load_smem_visc2_curv_x, load_smem_visc2_curv_y, load_smem_visc2_curv_z
contains
  attributes(device) subroutine load_smem_visc2_curv_x(it, jt, kt, j, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: idx           !< index for shared memory
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables (rho*u)
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables (rho*v)
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables (rho*w)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), intent(inout) :: u(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: w(0:sx*sy*sz-1) !< attribute(shared)
    integer i_base, ii, i, idx_l, offset_yz
    i_base    = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    do ii = it, threadsEv%x+1, blockDim%x
      i = i_base + ii
      idx_l = (ii-1) + offset_yz
      if (1 <= i .and. i <= nx .and. j <= ny .and. k <= nz) then
        call pipelineMemcpyAsync(u(idx_l), Q_2(i,j,k))
        call pipelineMemcpyAsync(v(idx_l), Q_3(i,j,k))
        call pipelineMemcpyAsync(w(idx_l), Q_4(i,j,k))
      else
        u(idx_l) = 0.d0
        v(idx_l) = 0.d0
        w(idx_l) = 0.d0
      endif
    enddo
    call pipelineCommit()
    call pipelineWaitPrior(0)
    call syncthreads()
  end subroutine load_smem_visc2_curv_x


  attributes(device) subroutine load_smem_visc2_curv_y(it, jt, kt, i, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: idx           !< index for shared memory
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables (rho*u)
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables (rho*v)
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables (rho*w)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), intent(inout) :: u(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: w(0:sx*sy*sz-1) !< attribute(shared)
    integer j_base, jj, j, idx_l, offset_xz
    j_base    = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    do jj = jt, threadsFv%y+1, blockDim%y
      j = j_base + jj
      idx_l = (jj-1) + offset_xz
      if (i <= nx .and. 1 <= j .and. j <= ny .and. k <= nz) then
        call pipelineMemcpyAsync(u(idx_l), Q_2(i,j,k))
        call pipelineMemcpyAsync(v(idx_l), Q_3(i,j,k))
        call pipelineMemcpyAsync(w(idx_l), Q_4(i,j,k))
      else
        u(idx_l) = 0.d0
        v(idx_l) = 0.d0
        w(idx_l) = 0.d0
      endif
    enddo
    call pipelineCommit()
    call pipelineWaitPrior(0)
    call syncthreads()
  end subroutine load_smem_visc2_curv_y


  attributes(device) subroutine load_smem_visc2_curv_z(it, jt, kt, i, j, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: idx           !< index for shared memory
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables (rho*u)
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables (rho*v)
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables (rho*w)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), intent(inout) :: u(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: v(0:sx*sy*sz-1) !< attribute(shared)
    real(8), intent(inout) :: w(0:sx*sy*sz-1) !< attribute(shared)
    integer k_base, kk, k, idx_l, offset_xy
    k_base    = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    do kk = kt, threadsGv%z+1, blockDim%z
      k = k_base + kk
      idx_l = (kk-1) + offset_xy
      if (i <= nx .and. j <= ny .and. 1 <= k .and. k <= nz) then
        call pipelineMemcpyAsync(u(idx_l), Q_2(i,j,k))
        call pipelineMemcpyAsync(v(idx_l), Q_3(i,j,k))
        call pipelineMemcpyAsync(w(idx_l), Q_4(i,j,k))
      else
        u(idx_l) = 0.d0
        v(idx_l) = 0.d0
        w(idx_l) = 0.d0
      endif
    enddo
    call pipelineCommit()
    call pipelineWaitPrior(0)
    call syncthreads()
  end subroutine load_smem_visc2_curv_z
end module load_smem_visc2_curv
