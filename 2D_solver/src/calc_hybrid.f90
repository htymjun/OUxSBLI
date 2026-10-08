!> Module for shock detection and hybrid scheme support
!> Computes Ducros sensor for automatic scheme switching between KEEP and SLAU
module calc_hybrid
  use mod_precision
  use cudafor
  use mod_globals, only : sp
  use mod_constant, only : gamma => gamma_conv
  implicit none
contains
  !> Compute Ducros shock sensor for hybrid scheme
  !> Uses ratio of dilatation (divergence) to vorticity to detect shocks
  !> Values closer to 1 indicate shock regions, close to 0 indicates smooth flow
  attributes(global) subroutine calc_Ducros(nx, ny, dx, dy, Q_2, Q_3, fd)
    integer, intent(in), value                      :: nx, ny
    real(kd_arr), intent(in), dimension(nx-1), device    :: dx ! 1 / dx
    real(kd_arr), intent(in), dimension(ny-1), device    :: dy ! 1 / dy
    real(kd_arr), intent(in), dimension(nx,ny), device   :: Q_2 ! u
    real(kd_arr), intent(in), dimension(nx,ny), device   :: Q_3 ! v
    real(sp), intent(out), device                   :: fd(nx,ny)
    integer i, j
    real(kd_arr) dudx, dudy, dvdx, dvdy
    real(kd_arr) dx_tmp, dy_tmp
    real(sp) div, rot
    real(sp), parameter :: eps = 1.0e-12_sp
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1

    if (nx-1 < i .or. ny-1 < j) return
    dx_tmp = 0.25_kd_arr * (dx(i-1) + dx(i))
    dy_tmp = 0.25_kd_arr * (dy(j-1) + dy(j))
    dudx = (-Q_2(i-1,j) + Q_2(i+1,j)) * dx_tmp
    dvdx = (-Q_3(i-1,j) + Q_3(i+1,j)) * dx_tmp
    dudy = (-Q_2(i,j-1) + Q_2(i,j+1)) * dy_tmp
    dvdy = (-Q_3(i,j-1) + Q_3(i,j+1)) * dy_tmp
    ! Ducros shock sensor: detector based on dilatation vs. vorticity
    div = real(dudx + dvdy, kind=sp) ! Divergence: ∇·u  !div = real(dudx + dvdy, kind=sp)
    
    ! Vorticity vector: ω = ∇ × u
    rot = real(dvdx - dudy, kind=sp)
    
    ! Sensor: f_d = (∇·u)² / [(∇·u)² + (∇×u)²]
    ! Returns ~1 in shocks (high compression), ~0 in smooth vortical flows
    fd(i,j) = (div**2) / (div**2 + rot**2 + eps)

    fd(i,j) = min(1.0_sp, fd(i,j))

    ! boundary
    ! x direction
    if (i == 2) then
      fd(1,j) = fd(2,j)
    elseif (i == nx-1) then
      fd(nx,j) = fd(nx-1,j)
    endif
    ! y direction
    if (j == 2) then
      fd(i,1) = fd(i,2)
    elseif (j == ny-1) then
      fd(i,ny) = fd(i,ny-1)
    endif
  end subroutine calc_Ducros


  pure attributes(device) function Albada(e, rho) result(phi)
    real(kd_conv), intent(in), dimension(4), device :: e, rho
    real(kd_conv) :: d1, d2, d3, phim, phip, phi
    real(kd_conv), parameter :: eps = 1.e-16_kd_conv
    d1   = -e(1) / rho(1) + e(2) / rho(2)
    d2   = -e(2) / rho(2) + e(3) / rho(3)
    d3   = -e(3) / rho(3) + e(4) / rho(4)
    phip = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phim = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi  = max(min(1._kd_conv - min(phim, phip), 1._kd_conv), 0._kd_conv)
  end function Albada


  pure attributes(device) function wiggle_detector(phi) result(ans)
    real(kd_conv), intent(in) :: phi(4)
    real(kd_conv) ans, phi1, phi2
    phi1 = (-phi(1) + phi(2)) * (-phi(2) + phi(3))
    phi2 = (-phi(3) + phi(4)) * (-phi(2) + phi(3))
    ans  = 0.5_kd_conv * (1._kd_conv - sign(1._kd_conv, min(phi1, phi2)))
  end function wiggle_detector
end module calc_hybrid
