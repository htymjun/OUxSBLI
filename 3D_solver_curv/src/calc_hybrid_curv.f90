!> Module for shock detection and hybrid scheme support in curvilinear grids
!> Computes Ducros sensor using physical velocity gradients via chain rule with inverse metrics
module calc_hybrid_curv
  use cudafor
  use mod_globals, only : gamma, sp
  implicit none
contains
  !> Compute Ducros shock sensor for curvilinear grids
  !> Uses ratio of dilatation (divergence) to vorticity to detect shocks
  !> Physical velocity gradients computed from computational derivatives via chain rule
  attributes(global) subroutine calc_Ducros_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, Q_2, Q_3, Q_4, fd)
    integer, intent(in), value                :: nx, ny, nz
    real(8), intent(in), value                :: dz
    real(8), intent(in), device, contiguous   :: xi_x(nx,ny), xi_y(nx,ny)
    real(8), intent(in), device, contiguous   :: eta_x(nx,ny), eta_y(nx,ny)
    real(8), intent(in), device, contiguous   :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous   :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous   :: Q_4(nx,ny,nz)
    real(sp), intent(out), device, contiguous :: fd(nx,ny,nz)
    integer i, j, k
    real(8) :: dudxi, dvdxi, dwdxi, dudeta, dvdeta, dwdeta, dudz, dvdz, dwdz
    real(8) :: dudx, dudy, dvdx, dvdy, dwdx, dwdy
    real(sp) :: div, rot(3)
    real(sp), parameter :: eps = 1.0e-12_sp
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    ! Computational space derivatives (Δξ=Δη=1)
    dudxi  = 0.5d0 * (Q_2(i+1,j,k) - Q_2(i-1,j,k))
    dvdxi  = 0.5d0 * (Q_3(i+1,j,k) - Q_3(i-1,j,k))
    dwdxi  = 0.5d0 * (Q_4(i+1,j,k) - Q_4(i-1,j,k))
    dudeta = 0.5d0 * (Q_2(i,j+1,k) - Q_2(i,j-1,k))
    dvdeta = 0.5d0 * (Q_3(i,j+1,k) - Q_3(i,j-1,k))
    dwdeta = 0.5d0 * (Q_4(i,j+1,k) - Q_4(i,j-1,k))
    dudz   = 0.5d0 * (Q_2(i,j,k+1) - Q_2(i,j,k-1)) / dz
    dvdz   = 0.5d0 * (Q_3(i,j,k+1) - Q_3(i,j,k-1)) / dz
    dwdz   = 0.5d0 * (Q_4(i,j,k+1) - Q_4(i,j,k-1)) / dz
    ! Physical space derivatives via chain rule
    ! du/dx = ∂u/∂ξ * ∂ξ/∂x + ∂u/∂η * ∂η/∂x = xi_x(i,j)*dudxi + eta_x(i,j)*dudeta
    ! du/dy = ∂u/∂ξ * ∂ξ/∂y + ∂u/∂η * ∂η/∂y = xi_y(i,j)*dudxi + eta_y(i,j)*dudeta
    dudx = xi_x(i,j)*dudxi + eta_x(i,j)*dudeta
    dudy = xi_y(i,j)*dudxi + eta_y(i,j)*dudeta
    dvdx = xi_x(i,j)*dvdxi + eta_x(i,j)*dvdeta
    dvdy = xi_y(i,j)*dvdxi + eta_y(i,j)*dvdeta
    dwdx = xi_x(i,j)*dwdxi + eta_x(i,j)*dwdeta
    dwdy = xi_y(i,j)*dwdxi + eta_y(i,j)*dwdeta
    ! Ducros shock sensor: detector based on dilatation vs. vorticity
    div    = real(dudx + dvdy + dwdz, sp) ! Divergence: ∇·u
    rot(1) = real(dwdy - dvdz, sp)        ! ω_x = dw/dy - dv/dz
    rot(2) = real(dudz - dwdx, sp)        ! ω_y = du/dz - dw/dx
    rot(3) = real(dvdx - dudy, sp)        ! ω_z = dv/dx - du/dy
    ! Sensor: f_d = (∇·u)² / [(∇·u)² + (∇×u)²]
    ! Returns ~1 in shocks (high compression), ~0 in smooth vortical flows
    fd(i,j,k) = min(1.0_sp, (div**2) / (div**2 + rot(1)**2 + rot(2)**2 + rot(3)**2 + eps))
    ! xi is periodic on the O-grid: ghost cells wrap to opposite interior cells
    if (i == 2)    fd(1,j,k)  = fd(nx-1,j,k)
    if (i == nx-1) fd(nx,j,k) = fd(2,j,k)
    if (j == 2)    fd(i,1,k)  = fd(i,2,k)
    if (j == ny-1) fd(i,ny,k) = fd(i,ny-1,k)
    if (k == 2)    fd(i,j,1)  = fd(i,j,2)
    if (k == nz-1) fd(i,j,nz) = fd(i,j,nz-1)
  end subroutine calc_Ducros_curv
end module calc_hybrid_curv
