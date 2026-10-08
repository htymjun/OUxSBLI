!> Module for du/dx, dv/dy at cell-center
!> Precomputed once per RK stage (VISC_ORDER>2) and reused by calc_Ev4/calc_Fv4
!> instead of each direction independently recomputing the other's gradient in
!> shared memory -- mirrors 3D_solver/src/calc_div.f90's ux/vy/wz split.
module calc_div
  use mod_precision
  use mod_constant, only : one_third => one_third_visc
  implicit none
  real(kd_visc), parameter :: one_24 = 1._kd_visc / 24._kd_visc
contains
  !> du/dx at cell-center (4th-order, interior stencil)
  subroutine calc_div_ux_4_in(nx, ny, xix, Q_2, ux)
    integer, intent(in), value                     :: nx, ny
    real(kd_arr), intent(in), dimension(nx-1), device   :: xix ! 1 / dx
    real(kd_arr), intent(in), dimension(nx,ny), device  :: Q_2
    real(kd_arr), intent(out), dimension(nx,ny), device :: ux
    integer, parameter :: io_v = 1
    integer i, j
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny
      do i = io_v+2, nx-io_v-1
        ux(i,j) = (one_third * (-real(Q_2(i-1,j), kd_visc) + real(Q_2(i+1,j), kd_visc)) &
                    - one_24 * (-real(Q_2(i-2,j), kd_visc) + real(Q_2(i+2,j), kd_visc))) * (real(xix(i-1), kd_visc) + real(xix(i), kd_visc))
      enddo
    enddo
  end subroutine calc_div_ux_4_in

  !> dv/dy at cell-center (4th-order, interior stencil)
  subroutine calc_div_vy_4_in(nx, ny, etay, Q_3, vy)
    integer, intent(in), value                     :: nx, ny
    real(kd_arr), intent(in), dimension(ny-1), device   :: etay ! 1 / dy
    real(kd_arr), intent(in), dimension(nx,ny), device  :: Q_3
    real(kd_arr), intent(out), dimension(nx,ny), device :: vy
    integer, parameter :: io_v = 1
    integer i, j
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = io_v+2, ny-io_v-1
      do i = 1, nx
        vy(i,j) = (one_third * (-real(Q_3(i,j-1), kd_visc) + real(Q_3(i,j+1), kd_visc)) &
                    - one_24 * (-real(Q_3(i,j-2), kd_visc) + real(Q_3(i,j+2), kd_visc))) * (real(etay(j-1), kd_visc) + real(etay(j), kd_visc))
      enddo
    enddo
  end subroutine calc_div_vy_4_in
end module calc_div
