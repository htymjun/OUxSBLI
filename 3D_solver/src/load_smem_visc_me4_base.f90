module load_smem_visc_me4_base
  use mod_globals, only : threadsEv, threadsFv, threadsGv
  use mod_constant, only : two_third, one_twelfth
  implicit none
  private
  public load_smem_visc_me4_base_x, load_smem_visc_me4_base_y, load_smem_visc_me4_base_z, load_smem_visc_me4_base_z_koff
contains
  attributes(device) subroutine load_smem_visc_me4_base_x(it, jt, kt, j, k, &
                                                          nx, ny, nz, inv_dy_j, inv_dz_k, Q_2, Q_3, Q_4, u, v, w, uy, vy, uz, wz)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), value              :: inv_dy_j      !< inverse grid spacing in y (1/dy) at j
    real(8), intent(in), value              :: inv_dz_k      !< inverse grid spacing in z (1/dz) at k
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsEv%x + 2*io_v + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: uy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: uz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer i_base, ii, i, idx, offset_yz
    logical :: jk_in_range, jk_in_range_grad_y, jk_in_range_grad_z
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1)*sx + (kt-1)*sx*sy
    jk_in_range        = (j <= ny .and. k <= nz)
    jk_in_range_grad_y = (3 <= j .and. j <= ny-2 .and. k <= nz)
    jk_in_range_grad_z = (j <= ny .and. 3 <= k .and. k <= nz-2)
    do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
      i = i_base + ii
      idx = ii + offset_yz
      if (1 <= i .and. i <= nx .and. jk_in_range) then
        u(idx) = Q_2(i,j,k)
        v(idx) = Q_3(i,j,k)
        w(idx) = Q_4(i,j,k)
      endif
      if (1 <= i .and. i <= nx .and. jk_in_range_grad_y) then
        vy(idx) = (two_third * (-Q_3(i,j-1,k) + Q_3(i,j+1,k)) - one_twelfth * (-Q_3(i,j-2,k) + Q_3(i,j+2,k))) * inv_dy_j
        uy(idx) = (two_third * (-Q_2(i,j-1,k) + Q_2(i,j+1,k)) - one_twelfth * (-Q_2(i,j-2,k) + Q_2(i,j+2,k))) * inv_dy_j
      endif
      if (1 <= i .and. i <= nx .and. jk_in_range_grad_z) then
        wz(idx) = (two_third * (-Q_4(i,j,k-1) + Q_4(i,j,k+1)) - one_twelfth * (-Q_4(i,j,k-2) + Q_4(i,j,k+2))) * inv_dz_k
        uz(idx) = (two_third * (-Q_2(i,j,k-1) + Q_2(i,j,k+1)) - one_twelfth * (-Q_2(i,j,k-2) + Q_2(i,j,k+2))) * inv_dz_k
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc_me4_base_x


  attributes(device) subroutine load_smem_visc_me4_base_y(it, jt, kt, i, k, &
                                                          nx, ny, nz, inv_dx_i, inv_dz_k, Q_2, Q_3, Q_4, u, v, w, ux, vx, vz, wz)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: k             !< global idx for z direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), value              :: inv_dx_i      !< inverse grid spacing in x (1/dx) at i
    real(8), intent(in), value              :: inv_dz_k      !< inverse grid spacing in z (1/dz) at k
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 2*io_v + 1
    integer, parameter :: sz = threadsFv%z
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: ux(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vx(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wz(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer j_base, jj, j, idx, offset_xz
    logical :: ik_in_range, ik_in_range_grad_x, ik_in_range_grad_z
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1)*sy + (kt-1)*sy*sx
    ik_in_range        = (i <= nx .and. k <= nz)
    ik_in_range_grad_x = (3 <= i .and. i <= nx-2 .and. k <= nz)
    ik_in_range_grad_z = (i <= nx .and. 3 <= k .and. k <= nz-2)
    do jj = jt-io_v, threadsFv%y+io_v+1, blockDim%y
      j = j_base + jj
      idx = jj + offset_xz
      if (ik_in_range .and. 1 <= j .and. j <= ny) then
        u(idx) = Q_2(i,j,k)
        v(idx) = Q_3(i,j,k)
        w(idx) = Q_4(i,j,k)
      endif
      if (ik_in_range_grad_x .and. 1 <= j .and. j <= ny) then
        ux(idx) = (two_third * (-Q_2(i-1,j,k) + Q_2(i+1,j,k)) - one_twelfth * (-Q_2(i-2,j,k) + Q_2(i+2,j,k))) * inv_dx_i
        vx(idx) = (two_third * (-Q_3(i-1,j,k) + Q_3(i+1,j,k)) - one_twelfth * (-Q_3(i-2,j,k) + Q_3(i+2,j,k))) * inv_dx_i
      endif
      if (ik_in_range_grad_z .and. 1 <= j .and. j <= ny) then
        vz(idx) = (two_third * (-Q_3(i,j,k-1) + Q_3(i,j,k+1)) - one_twelfth * (-Q_3(i,j,k-2) + Q_3(i,j,k+2))) * inv_dz_k
        wz(idx) = (two_third * (-Q_4(i,j,k-1) + Q_4(i,j,k+1)) - one_twelfth * (-Q_4(i,j,k-2) + Q_4(i,j,k+2))) * inv_dz_k
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc_me4_base_y


  attributes(device) subroutine load_smem_visc_me4_base_z(it, jt, kt, i, j, &
                                                          nx, ny, nz, inv_dx_i, inv_dy_j, Q_2, Q_3, Q_4, u, v, w, ux, wx, vy, wy)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), value              :: inv_dx_i      !< inverse grid spacing in x (1/dx) at i
    real(8), intent(in), value              :: inv_dy_j      !< inverse grid spacing in y (1/dy) at j
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 2*io_v + 1
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: ux(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wx(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer k_base, kk, k, idx, offset_xy
    logical :: ij_in_range, ij_in_range_grad_x, ij_in_range_grad_y
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    ij_in_range        = (i <= nx .and. j <= ny)
    ij_in_range_grad_x = (3 <= i .and. i <= nx-2 .and. j <= ny)
    ij_in_range_grad_y = (i <= nx .and. 3 <= j .and. j <= ny-2)
    do kk = kt-io_v, threadsGv%z+io_v+1, blockDim%z
      k = k_base + kk
      idx = kk + offset_xy
      if (ij_in_range .and. 1 <= k .and. k <= nz) then
        u(idx) = Q_2(i,j,k)
        v(idx) = Q_3(i,j,k)
        w(idx) = Q_4(i,j,k)
      endif
      if (ij_in_range_grad_x .and. 1 <= k .and. k <= nz) then
        ux(idx) = (two_third * (-Q_2(i-1,j,k) + Q_2(i+1,j,k)) - one_twelfth * (-Q_2(i-2,j,k) + Q_2(i+2,j,k))) * inv_dx_i
        wx(idx) = (two_third * (-Q_4(i-1,j,k) + Q_4(i+1,j,k)) - one_twelfth * (-Q_4(i-2,j,k) + Q_4(i+2,j,k))) * inv_dx_i
      endif
      if (ij_in_range_grad_y .and. 1 <= k .and. k <= nz) then
        vy(idx) = (two_third * (-Q_3(i,j-1,k) + Q_3(i,j+1,k)) - one_twelfth * (-Q_3(i,j-2,k) + Q_3(i,j+2,k))) * inv_dy_j
        wy(idx) = (two_third * (-Q_4(i,j-1,k) + Q_4(i,j+1,k)) - one_twelfth * (-Q_4(i,j-2,k) + Q_4(i,j+2,k))) * inv_dy_j
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc_me4_base_z


  !> Like load_smem_visc4_z but k_base is shifted by k_lo-1 for koff kernel launches
  attributes(device) subroutine load_smem_visc_me4_base_z_koff(it, jt, kt, i, j, &
                                                               nx, ny, nz, inv_dx_i, inv_dy_j, Q_2, Q_3, Q_4, u, v, w, ux, wx, vy, wy, k_lo)
    integer, intent(in), value              :: it            !< local idx for x direction
    integer, intent(in), value              :: jt            !< local idx for y direction
    integer, intent(in), value              :: kt            !< local idx for z direction
    integer, intent(in), value              :: i             !< global idx for x direction
    integer, intent(in), value              :: j             !< global idx for y direction
    integer, intent(in), value              :: nx            !< number of grid points in x direction
    integer, intent(in), value              :: ny            !< number of grid points in y direction
    integer, intent(in), value              :: nz            !< number of grid points in z direction
    real(8), intent(in), value              :: inv_dx_i      !< inverse grid spacing in x (1/dx) at i
    real(8), intent(in), value              :: inv_dy_j      !< inverse grid spacing in y (1/dy) at j
    real(8), intent(in), device, contiguous :: Q_2(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_3(nx,ny,nz) !< conservative variables
    real(8), intent(in), device, contiguous :: Q_4(nx,ny,nz) !< conservative variables
    integer, intent(in), value              :: k_lo          !< z-index lower bound (global, 1-based)
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 2*io_v + 1
    real(8), intent(inout) ::  u(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  v(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) ::  w(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared)
    real(8), intent(inout) :: ux(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wx(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: vy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    real(8), intent(inout) :: wy(-(io_v-1):sx*sy*sz-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer k_base, kk, k, idx, offset_xy
    logical :: ij_in_range, ij_in_range_grad_x, ij_in_range_grad_y
    k_base = (blockIdx%z-1)*blockDim%z + k_lo - 1
    offset_xy = (jt-1)*sz + (it-1)*sz*sy
    ij_in_range        = (i <= nx .and. j <= ny)
    ij_in_range_grad_x = (3 <= i .and. i <= nx-2 .and. j <= ny)
    ij_in_range_grad_y = (i <= nx .and. 3 <= j .and. j <= ny-2)
    do kk = kt-io_v, threadsGv%z+io_v+1, blockDim%z
      k = k_base + kk
      idx = kk + offset_xy
      if (ij_in_range .and. 1 <= k .and. k <= nz) then
        u(idx) = Q_2(i,j,k)
        v(idx) = Q_3(i,j,k)
        w(idx) = Q_4(i,j,k)
      endif
      if (ij_in_range_grad_x .and. 1 <= k .and. k <= nz) then
        ux(idx) = (two_third * (-Q_2(i-1,j,k) + Q_2(i+1,j,k)) - one_twelfth * (-Q_2(i-2,j,k) + Q_2(i+2,j,k))) * inv_dx_i
        wx(idx) = (two_third * (-Q_4(i-1,j,k) + Q_4(i+1,j,k)) - one_twelfth * (-Q_4(i-2,j,k) + Q_4(i+2,j,k))) * inv_dx_i
      endif
      if (ij_in_range_grad_y .and. 1 <= k .and. k <= nz) then
        vy(idx) = (two_third * (-Q_3(i,j-1,k) + Q_3(i,j+1,k)) - one_twelfth * (-Q_3(i,j-2,k) + Q_3(i,j+2,k))) * inv_dy_j
        wy(idx) = (two_third * (-Q_4(i,j-1,k) + Q_4(i,j+1,k)) - one_twelfth * (-Q_4(i,j-2,k) + Q_4(i,j+2,k))) * inv_dy_j
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc_me4_base_z_koff
end module load_smem_visc_me4_base
