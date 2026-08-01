!> Module for div(u) at cell-center
module calc_div
  use cudafor
  use mod_constant, only : one_third
  implicit none
  real(8), parameter :: one_24  = 1.d0 / 24.d0
  real(8), parameter :: one_120 = 1.d0 / 120.d0
contains
  !> Compute div(u) at cell-center (2nd-order), over the full [1,nx]x[1,ny]x[1,nz]
  !! domain -- matches calc_quantities_T_3D's coverage, since consumers read
  !! ux/vy/wz near the array edges (e.g. Ev's interior branch reads vy at
  !! i-io_v, which can reach index 1).
  subroutine calc_div_2(nx, ny, nz, xix, etay, zetaz, Q_2, Q_3, Q_4, ux, vy, wz)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: xix   ! 1 / dx
    real(8), intent(in), dimension(ny-1), device      :: etay  ! 1 / dy
    real(8), intent(in), dimension(nz-1), device      :: zetaz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_2, Q_3, Q_4
    real(8), intent(out), dimension(nx,ny,nz), device :: ux, vy, wz
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (2 <= i .and. i <= nx-1) then
            ux(i,j,k) = 0.25d0 * (-Q_2(i-1,j,k) + Q_2(i+1,j,k)) * (xix(i-1) + xix(i))
          endif
          if (2 <= j .and. j <= ny-1) then
            vy(i,j,k) = 0.25d0 * (-Q_3(i,j-1,k) + Q_3(i,j+1,k)) * (etay(j-1) + etay(j))
          endif
          if (2 <= k .and. k <= nz-1) then
            wz(i,j,k) = 0.25d0 * (-Q_4(i,j,k-1) + Q_4(i,j,k+1)) * (zetaz(k-1) + zetaz(k))
          endif
        enddo
      enddo
    enddo
  end subroutine calc_div_2

  ! ux/vy/wz (4th- and 6th-order) are each computed by an independent,
  ! single-component kernel with a directly-bounded loop in its own margin
  ! direction (no if-branching, no other component's arithmetic live at the
  ! same time) -- a prior combined-kernel version carried all three
  ! components' index arithmetic and Q reads in one register footprint,
  ! gated by 3 separate if-statements per thread; splitting them dropped
  ! register pressure substantially (see calc_flux_base.f90.fypp's 3 calls).
  ! Only the swept direction's own margin matters (see calc_div_2's docstring
  ! for why: the other two array dimensions are read at whatever index the
  ! thread already owns, never validity-checked).

  !> du/dx at cell-center (4th-order, interior stencil)
  subroutine calc_div_ux_4_in(nx, ny, nz, xix, Q_2, ux)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: xix   ! 1 / dx
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_2
    real(8), intent(out), dimension(nx,ny,nz), device :: ux
    integer, parameter :: io_v = 1
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = 1, ny
        do i = io_v+2, nx-io_v-1
          ux(i,j,k) = (one_third * (-Q_2(i-1,j,k) + Q_2(i+1,j,k)) &
                        - one_24 * (-Q_2(i-2,j,k) + Q_2(i+2,j,k))) * (xix(i-1) + xix(i))
        enddo
      enddo
    enddo
  end subroutine calc_div_ux_4_in

  !> dv/dy at cell-center (4th-order, interior stencil)
  subroutine calc_div_vy_4_in(nx, ny, nz, etay, Q_3, vy)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(ny-1), device      :: etay  ! 1 / dy
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_3
    real(8), intent(out), dimension(nx,ny,nz), device :: vy
    integer, parameter :: io_v = 1
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = io_v+2, ny-io_v-1
        do i = 1, nx
          vy(i,j,k) = (one_third * (-Q_3(i,j-1,k) + Q_3(i,j+1,k)) &
                        - one_24 * (-Q_3(i,j-2,k) + Q_3(i,j+2,k))) * (etay(j-1) + etay(j))
        enddo
      enddo
    enddo
  end subroutine calc_div_vy_4_in

  !> dw/dz at cell-center (4th-order, interior stencil)
  subroutine calc_div_wz_4_in(nx, ny, nz, zetaz, Q_4, wz)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nz-1), device      :: zetaz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_4
    real(8), intent(out), dimension(nx,ny,nz), device :: wz
    integer, parameter :: io_v = 1
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = io_v+2, nz-io_v-1
      do j = 1, ny
        do i = 1, nx
          wz(i,j,k) = (one_third * (-Q_4(i,j,k-1) + Q_4(i,j,k+1)) &
                        - one_24 * (-Q_4(i,j,k-2) + Q_4(i,j,k+2))) * (zetaz(k-1) + zetaz(k))
        enddo
      enddo
    enddo
  end subroutine calc_div_wz_4_in

  !> du/dx at cell-center (6th-order, interior stencil)
  subroutine calc_div_ux_6_in(nx, ny, nz, xix, Q_2, ux)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device      :: xix   ! 1 / dx
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_2
    real(8), intent(out), dimension(nx,ny,nz), device :: ux
    integer, parameter :: io_v = 2
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = 1, ny
        do i = io_v+2, nx-io_v-1
          ux(i,j,k) = (one_120 * (-Q_2(i-3,j,k) + Q_2(i+3,j,k)) &
                      + 0.075d0 * (Q_2(i-2,j,k) - Q_2(i+2,j,k)) &
                     + 0.375d0 * (-Q_2(i-1,j,k) + Q_2(i+1,j,k))) * (xix(i-1) + xix(i))
        enddo
      enddo
    enddo
  end subroutine calc_div_ux_6_in

  !> dv/dy at cell-center (6th-order, interior stencil)
  subroutine calc_div_vy_6_in(nx, ny, nz, etay, Q_3, vy)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(ny-1), device      :: etay  ! 1 / dy
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_3
    real(8), intent(out), dimension(nx,ny,nz), device :: vy
    integer, parameter :: io_v = 2
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = io_v+2, ny-io_v-1
        do i = 1, nx
          vy(i,j,k) = (one_120 * (-Q_3(i,j-3,k) + Q_3(i,j+3,k)) &
                      + 0.075d0 * (Q_3(i,j-2,k) - Q_3(i,j+2,k)) &
                     + 0.375d0 * (-Q_3(i,j-1,k) + Q_3(i,j+1,k))) * (etay(j-1) + etay(j))
        enddo
      enddo
    enddo
  end subroutine calc_div_vy_6_in

  !> dw/dz at cell-center (6th-order, interior stencil)
  subroutine calc_div_wz_6_in(nx, ny, nz, zetaz, Q_4, wz)
    integer, intent(in), value                        :: nx, ny, nz
    real(8), intent(in), dimension(nz-1), device      :: zetaz ! 1 / dz
    real(8), intent(in), dimension(nx,ny,nz), device  :: Q_4
    real(8), intent(out), dimension(nx,ny,nz), device :: wz
    integer, parameter :: io_v = 2
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = io_v+2, nz-io_v-1
      do j = 1, ny
        do i = 1, nx
          wz(i,j,k) = (one_120 * (-Q_4(i,j,k-3) + Q_4(i,j,k+3)) &
                      + 0.075d0 * (Q_4(i,j,k-2) - Q_4(i,j,k+2)) &
                     + 0.375d0 * (-Q_4(i,j,k-1) + Q_4(i,j,k+1))) * (zetaz(k-1) + zetaz(k))
        enddo
      enddo
    enddo
  end subroutine calc_div_wz_6_in
end module calc_div
