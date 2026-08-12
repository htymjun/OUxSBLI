!> Module for shock detection and hybrid scheme support
!> Computes Ducros sensor for automatic scheme switching between KEEP and SLAU
module calc_hybrid
  use cudafor
  use mod_globals, only : gamma, sp
  implicit none
contains
  !> Compute Ducros shock sensor for hybrid scheme
  !> Uses ratio of dilatation (divergence) to vorticity to detect shocks
  !> Values closer to 1 indicate shock regions, close to 0 indicates smooth flow
  attributes(global) subroutine calc_Ducros(nx, ny, dx, dy, Q_2, Q_3, fd)
    integer, intent(in), value                      :: nx, ny
    real(8), intent(in), dimension(nx-1), device    :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device    :: dy ! 1 / dy
    real(8), intent(in), dimension(nx,ny), device   :: Q_2 ! u
    real(8), intent(in), dimension(nx,ny), device   :: Q_3 ! v
    real(sp), intent(out), device                   :: fd(nx,ny)
    integer i, j
    real(8) dudx, dudy, dvdx, dvdy
    real(8) dx_tmp, dy_tmp
    real(sp) div, rot
    real(sp), parameter :: eps = 1.0e-12_sp
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1

    if (nx-1 < i .or. ny-1 < j) return
    dx_tmp = 0.25d0 * (dx(i-1) + dx(i))
    dy_tmp = 0.25d0 * (dy(j-1) + dy(j))
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
    real(8), intent(in), dimension(4), device :: e, rho
    real(8) :: d1, d2, d3, phim, phip, phi
    real(8), parameter :: eps = 1.d-16
    d1   = -e(1) / rho(1) + e(2) / rho(2)
    d2   = -e(2) / rho(2) + e(3) / rho(3)
    d3   = -e(3) / rho(3) + e(4) / rho(4)
    phip = (d2 * d1 + d1**2) / (d2**2 + d1**2 + eps)
    phim = (d2 * d3 + d3**2) / (d2**2 + d3**2 + eps)
    phi  = max(min(1.d0 - min(phim, phip), 1.d0), 0.d0)
  end function Albada


  pure attributes(device) function wiggle_detector(phi) result(ans)
    real(8), intent(in) :: phi(4)
    real(8) ans, phi1, phi2
    phi1 = (-phi(1) + phi(2)) * (-phi(2) + phi(3))
    phi2 = (-phi(3) + phi(4)) * (-phi(2) + phi(3))
    ans  = 0.5d0 * (1.d0 - sign(1.d0, min(phi1, phi2)))
  end function wiggle_detector
end module calc_hybrid
