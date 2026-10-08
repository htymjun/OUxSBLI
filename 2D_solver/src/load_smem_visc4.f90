module load_smem_visc4
  use mod_precision
  use mod_globals, only : threadsEv, threadsFv
  use mod_constant, only : two_third => two_third_visc, one_twelfth => one_twelfth_visc
  implicit none
  private
  public load_smem_visc4_x, load_smem_visc4_y
contains
  !> Loads u, v tiles and the E-direction's own cross-derivative uy=du/dy.
  !> vy=dv/dy (the gradient this direction doesn't own) is precomputed once
  !> per RK stage by calc_div and read directly from DRAM by calc_Ev4 instead.
  attributes(device) subroutine load_smem_visc4_x(it, jt, j, &
                                                  nx, ny, inv_dy, Q_2, Q_3, u, v, uy)
    integer, intent(in), value              :: it           !< local idx for x direction
    integer, intent(in), value              :: jt           !< local idx for y direction
    integer, intent(in), value              :: j            !< global idx for y direction
    integer, intent(in), value              :: nx           !< number of grid points in x direction
    integer, intent(in), value              :: ny           !< number of grid points in y direction
    real(kd_arr), intent(in), device, contiguous :: inv_dy(ny-1) !< inverse grid spacing in y (1/dy)
    real(kd_arr), intent(in), device, contiguous :: Q_2(nx,ny)   !< rho*u
    real(kd_arr), intent(in), device, contiguous :: Q_3(nx,ny)   !< rho*v
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsEv%x + 2*io_v + 1
    integer, parameter :: sy = threadsEv%y
    real(kd_visc), intent(inout) ::  u(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(kd_visc), intent(inout) ::  v(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(kd_visc), intent(inout) :: uy(-(io_v-1):sx*sy-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer i_base, ii, i, idx, offset_y
    i_base = (blockIdx%x-1)*blockDim%x
    offset_y = (jt-1)*sx
    do ii = it-io_v, threadsEv%x+io_v+1, blockDim%x
      i = i_base + ii
      idx = ii + offset_y
      if (1 <= i .and. i <= nx .and. j <= ny) then
        u(idx) = Q_2(i,j)
        v(idx) = Q_3(i,j)
      endif
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2) then
        uy(idx) = (two_third * (-real(Q_2(i,j-1), kd_visc) + real(Q_2(i,j+1), kd_visc)) - one_twelfth * (-real(Q_2(i,j-2), kd_visc) + real(Q_2(i,j+2), kd_visc))) * real(inv_dy(j), kd_visc)
      else
        uy(idx) = 0._kd_visc
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc4_x


  !> Loads u, v tiles and the F-direction's own cross-derivative vx=dv/dx.
  !> ux=du/dx (the gradient this direction doesn't own) is precomputed once
  !> per RK stage by calc_div and read directly from DRAM by calc_Fv4 instead.
  attributes(device) subroutine load_smem_visc4_y(it, jt, i, &
                                                  nx, ny, inv_dx, Q_2, Q_3, u, v, vx)
    integer, intent(in), value              :: it           !< local idx for x direction
    integer, intent(in), value              :: jt           !< local idx for y direction
    integer, intent(in), value              :: i            !< global idx for x direction
    integer, intent(in), value              :: nx           !< number of grid points in x direction
    integer, intent(in), value              :: ny           !< number of grid points in y direction
    real(kd_arr), intent(in), device, contiguous :: inv_dx(nx-1) !< inverse grid spacing in x (1/dx)
    real(kd_arr), intent(in), device, contiguous :: Q_2(nx,ny)   !< rho*u
    real(kd_arr), intent(in), device, contiguous :: Q_3(nx,ny)   !< rho*v
    integer, parameter :: io_v = 2
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 2*io_v + 1
    real(kd_visc), intent(inout) ::  u(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(kd_visc), intent(inout) ::  v(-(io_v-1):sx*sy-io_v) !< attribute(shared)
    real(kd_visc), intent(inout) :: vx(-(io_v-1):sx*sy-io_v) !< attribute(shared) gradient 4th-order accuracy
    integer j_base, jj, j, idx, offset_x
    j_base = (blockIdx%y-1)*blockDim%y
    offset_x = (it-1)*sy
    do jj = jt-io_v, threadsFv%y+io_v+1, blockDim%y
      j = j_base + jj
      idx = jj + offset_x
      if (i <= nx .and. 1 <= j .and. j <= ny) then
        u(idx) = Q_2(i,j)
        v(idx) = Q_3(i,j)
      endif
      if (3 <= i .and. i <= nx-2 .and. 1 <= j .and. j <= ny) then
        vx(idx) = (two_third * (-real(Q_3(i-1,j), kd_visc) + real(Q_3(i+1,j), kd_visc)) - one_twelfth * (-real(Q_3(i-2,j), kd_visc) + real(Q_3(i+2,j), kd_visc))) * real(inv_dx(i), kd_visc)
      else
        vx(idx) = 0._kd_visc
      endif
    enddo
    call syncthreads()
  end subroutine load_smem_visc4_y
end module load_smem_visc4
