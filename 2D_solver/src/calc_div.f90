!> Module for du/dx, dv/dy at cell-center
!> Precomputed once per RK stage (VISC_ORDER>2) and reused by calc_Ev4/calc_Fv4
!> instead of each direction independently recomputing the other's gradient in
!> shared memory -- mirrors 3D_solver/src/calc_div.f90's ux/vy/wz split.
module calc_div
  use mod_constant, only : one_third
  implicit none
  real(8), parameter :: one_24 = 1.d0 / 24.d0
contains
  !> du/dx at cell-center (4th-order, interior stencil)
  subroutine calc_div_ux_4_in(nx, ny, xix, Q_2, ux)
    integer, intent(in), value                     :: nx, ny
    real(8), intent(in), dimension(nx-1), device   :: xix ! 1 / dx
    real(8), intent(in), dimension(nx,ny), device  :: Q_2
    real(8), intent(out), dimension(nx,ny), device :: ux
    integer, parameter :: io_v = 1
    integer i, j
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny
      do i = io_v+2, nx-io_v-1
        ux(i,j) = (one_third * (-Q_2(i-1,j) + Q_2(i+1,j)) &
                    - one_24 * (-Q_2(i-2,j) + Q_2(i+2,j))) * (xix(i-1) + xix(i))
      enddo
    enddo
  end subroutine calc_div_ux_4_in

  !> dv/dy at cell-center (4th-order, interior stencil)
  subroutine calc_div_vy_4_in(nx, ny, etay, Q_3, vy)
    integer, intent(in), value                     :: nx, ny
    real(8), intent(in), dimension(ny-1), device   :: etay ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device  :: Q_3
    real(8), intent(out), dimension(nx,ny), device :: vy
    integer, parameter :: io_v = 1
    integer i, j
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = io_v+2, ny-io_v-1
      do i = 1, nx
        vy(i,j) = (one_third * (-Q_3(i,j-1) + Q_3(i,j+1)) &
                    - one_24 * (-Q_3(i,j-2) + Q_3(i,j+2))) * (etay(j-1) + etay(j))
      enddo
    enddo
  end subroutine calc_div_vy_4_in
end module calc_div
