!> Module for div(u) at cell-center
module calc_div
  use mod_precision
  use cudafor
  use mod_constant, only : one_third => one_third_visc
  implicit none
  real(kd_visc), parameter :: one_24  = 1._kd_visc / 24._kd_visc
  real(kd_visc), parameter :: one_120 = 1._kd_visc / 120._kd_visc
contains
  !> Compute div(u) at cell-center (2nd-order), over the full [1,nx]x[1,ny]x[1,nz]
  !! domain -- matches calc_quantities_T_3D's coverage, since consumers read
  !! ux/vy/wz near the array edges (e.g. Ev's interior branch reads vy at
  !! i-io_v, which can reach index 1).
  subroutine calc_div_2(nx, ny, nz, xix, etay, zetaz, Q_2, Q_3, Q_4, ux, vy, wz)
    integer, intent(in), value                        :: nx, ny, nz
    real(kd_arr), intent(in), dimension(nx-1), device      :: xix   ! 1 / dx
    real(kd_arr), intent(in), dimension(ny-1), device      :: etay  ! 1 / dy
    real(kd_arr), intent(in), dimension(nz-1), device      :: zetaz ! 1 / dz
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_2, Q_3, Q_4
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: ux, vy, wz
    integer i, j, k
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (2 <= i .and. i <= nx-1) then
            ux(i,j,k) = 0.25_kd_visc * (-real(Q_2(i-1,j,k), kd_visc) + real(Q_2(i+1,j,k), kd_visc)) * (real(xix(i-1), kd_visc) + real(xix(i), kd_visc))
          endif
          if (2 <= j .and. j <= ny-1) then
            vy(i,j,k) = 0.25_kd_visc * (-real(Q_3(i,j-1,k), kd_visc) + real(Q_3(i,j+1,k), kd_visc)) * (real(etay(j-1), kd_visc) + real(etay(j), kd_visc))
          endif
          if (2 <= k .and. k <= nz-1) then
            wz(i,j,k) = 0.25_kd_visc * (-real(Q_4(i,j,k-1), kd_visc) + real(Q_4(i,j,k+1), kd_visc)) * (real(zetaz(k-1), kd_visc) + real(zetaz(k), kd_visc))
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
  !
  ! Every kernel sweeps only the z planes [k_lo, k_hi]. The non-COMMZ path
  ! passes the whole range (1..nz for ux/vy, io_v+2..nz-io_v-1 for wz); the
  ! COMMZ path computes the planes whose stencil is complete before the
  ! z-halo exchange and the remaining planes afterwards (calc_flux_base).
  ! The caller keeps wz's range inside the stencil-safe interior; an empty
  ! range is a no-op.

  !> du/dx at cell-center (4th-order, interior stencil)
  subroutine calc_div_ux_4_in(nx, ny, nz, xix, Q_2, ux, k_lo, k_hi)
    integer, intent(in), value                        :: nx, ny, nz, k_lo, k_hi
    real(kd_arr), intent(in), dimension(nx-1), device      :: xix   ! 1 / dx
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_2
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: ux
    integer, parameter :: io_v = 1
    integer i, j, k
    if (k_hi < k_lo) return
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = io_v+2, nx-io_v-1
          ux(i,j,k) = (one_third * (-real(Q_2(i-1,j,k), kd_visc) + real(Q_2(i+1,j,k), kd_visc)) &
                        - one_24 * (-real(Q_2(i-2,j,k), kd_visc) + real(Q_2(i+2,j,k), kd_visc))) * (real(xix(i-1), kd_visc) + real(xix(i), kd_visc))
        enddo
      enddo
    enddo
  end subroutine calc_div_ux_4_in

  !> dv/dy at cell-center (4th-order, interior stencil)
  subroutine calc_div_vy_4_in(nx, ny, nz, etay, Q_3, vy, k_lo, k_hi)
    integer, intent(in), value                        :: nx, ny, nz, k_lo, k_hi
    real(kd_arr), intent(in), dimension(ny-1), device      :: etay  ! 1 / dy
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_3
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: vy
    integer, parameter :: io_v = 1
    integer i, j, k
    if (k_hi < k_lo) return
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = io_v+2, ny-io_v-1
        do i = 1, nx
          vy(i,j,k) = (one_third * (-real(Q_3(i,j-1,k), kd_visc) + real(Q_3(i,j+1,k), kd_visc)) &
                        - one_24 * (-real(Q_3(i,j-2,k), kd_visc) + real(Q_3(i,j+2,k), kd_visc))) * (real(etay(j-1), kd_visc) + real(etay(j), kd_visc))
        enddo
      enddo
    enddo
  end subroutine calc_div_vy_4_in

  !> dw/dz at cell-center (4th-order, interior stencil); caller keeps
  !> [k_lo, k_hi] inside [io_v+2, nz-io_v-1] = [3, nz-2]
  subroutine calc_div_wz_4_in(nx, ny, nz, zetaz, Q_4, wz, k_lo, k_hi)
    integer, intent(in), value                        :: nx, ny, nz, k_lo, k_hi
    real(kd_arr), intent(in), dimension(nz-1), device      :: zetaz ! 1 / dz
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_4
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: wz
    integer i, j, k
    if (k_hi < k_lo) return
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          wz(i,j,k) = (one_third * (-real(Q_4(i,j,k-1), kd_visc) + real(Q_4(i,j,k+1), kd_visc)) &
                        - one_24 * (-real(Q_4(i,j,k-2), kd_visc) + real(Q_4(i,j,k+2), kd_visc))) * (real(zetaz(k-1), kd_visc) + real(zetaz(k), kd_visc))
        enddo
      enddo
    enddo
  end subroutine calc_div_wz_4_in

  !> du/dx at cell-center (6th-order, interior stencil)
  subroutine calc_div_ux_6_in(nx, ny, nz, xix, Q_2, ux, k_lo, k_hi)
    integer, intent(in), value                        :: nx, ny, nz, k_lo, k_hi
    real(kd_arr), intent(in), dimension(nx-1), device      :: xix   ! 1 / dx
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_2
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: ux
    integer, parameter :: io_v = 2
    integer i, j, k
    if (k_hi < k_lo) return
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = io_v+2, nx-io_v-1
          ux(i,j,k) = (one_120 * (-real(Q_2(i-3,j,k), kd_visc) + real(Q_2(i+3,j,k), kd_visc)) &
                      + 0.075_kd_visc * (real(Q_2(i-2,j,k), kd_visc) - real(Q_2(i+2,j,k), kd_visc)) &
                     + 0.375_kd_visc * (-real(Q_2(i-1,j,k), kd_visc) + real(Q_2(i+1,j,k), kd_visc))) * (real(xix(i-1), kd_visc) + real(xix(i), kd_visc))
        enddo
      enddo
    enddo
  end subroutine calc_div_ux_6_in

  !> dv/dy at cell-center (6th-order, interior stencil)
  subroutine calc_div_vy_6_in(nx, ny, nz, etay, Q_3, vy, k_lo, k_hi)
    integer, intent(in), value                        :: nx, ny, nz, k_lo, k_hi
    real(kd_arr), intent(in), dimension(ny-1), device      :: etay  ! 1 / dy
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_3
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: vy
    integer, parameter :: io_v = 2
    integer i, j, k
    if (k_hi < k_lo) return
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = io_v+2, ny-io_v-1
        do i = 1, nx
          vy(i,j,k) = (one_120 * (-real(Q_3(i,j-3,k), kd_visc) + real(Q_3(i,j+3,k), kd_visc)) &
                      + 0.075_kd_visc * (real(Q_3(i,j-2,k), kd_visc) - real(Q_3(i,j+2,k), kd_visc)) &
                     + 0.375_kd_visc * (-real(Q_3(i,j-1,k), kd_visc) + real(Q_3(i,j+1,k), kd_visc))) * (real(etay(j-1), kd_visc) + real(etay(j), kd_visc))
        enddo
      enddo
    enddo
  end subroutine calc_div_vy_6_in

  !> dw/dz at cell-center (6th-order, interior stencil); caller keeps
  !> [k_lo, k_hi] inside [io_v+2, nz-io_v-1] = [4, nz-3]
  subroutine calc_div_wz_6_in(nx, ny, nz, zetaz, Q_4, wz, k_lo, k_hi)
    integer, intent(in), value                        :: nx, ny, nz, k_lo, k_hi
    real(kd_arr), intent(in), dimension(nz-1), device      :: zetaz ! 1 / dz
    real(kd_arr), intent(in), dimension(nx,ny,nz), device  :: Q_4
    real(kd_arr), intent(out), dimension(nx,ny,nz), device :: wz
    integer i, j, k
    if (k_hi < k_lo) return
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          wz(i,j,k) = (one_120 * (-real(Q_4(i,j,k-3), kd_visc) + real(Q_4(i,j,k+3), kd_visc)) &
                      + 0.075_kd_visc * (real(Q_4(i,j,k-2), kd_visc) - real(Q_4(i,j,k+2), kd_visc)) &
                     + 0.375_kd_visc * (-real(Q_4(i,j,k-1), kd_visc) + real(Q_4(i,j,k+1), kd_visc))) * (real(zetaz(k-1), kd_visc) + real(zetaz(k), kd_visc))
        enddo
      enddo
    enddo
  end subroutine calc_div_wz_6_in
end module calc_div
