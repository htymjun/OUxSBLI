module calc_hybrid
  use cudafor
  implicit none
contains
  attributes(global) subroutine calc_Ducros(nx, ny, nz, dx, dy, dz, Q, fd)
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device       :: dx ! 1 / dx
    real(8), intent(in), dimension(ny-1), device       :: dy ! 1 / dy
    real(8), intent(in), dimension(nz-1), device       :: dz ! 1 / dz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(out), dimension(nx,ny,nz), device  :: fd
    integer i, j, k
    real(8) dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(8) div, rot(3)
    real(8) dx_tmp, dy_tmp, dz_tmp
    real(8), parameter :: eps = 1.d-16
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1 
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    dx_tmp = 0.25d0 * (dx(i-1) + dx(i))
    dy_tmp = 0.25d0 * (dy(j-1) + dy(j))
    dz_tmp = 0.25d0 * (dz(k-1) + dz(k))
    dudx = (-Q(2,i-1,j,k) + Q(2,i+1,j,k)) * dx_tmp
    dvdx = (-Q(3,i-1,j,k) + Q(3,i+1,j,k)) * dx_tmp
    dwdx = (-Q(4,i-1,j,k) + Q(4,i+1,j,k)) * dx_tmp
    dudy = (-Q(2,i,j-1,k) + Q(2,i,j+1,k)) * dy_tmp
    dvdy = (-Q(3,i,j-1,k) + Q(3,i,j+1,k)) * dy_tmp
    dwdy = (-Q(4,i,j-1,k) + Q(4,i,j+1,k)) * dy_tmp
    dudz = (-Q(2,i,j,k-1) + Q(2,i,j,k+1)) * dz_tmp
    dvdz = (-Q(3,i,j,k-1) + Q(3,i,j,k+1)) * dz_tmp
    dwdz = (-Q(4,i,j,k-1) + Q(4,i,j,k+1)) * dz_tmp
    div = dudx + dvdy + dwdz
    rot(1) = dwdy - dvdz
    rot(2) = dudz - dwdx
    rot(3) = dvdx - dudy
    fd(i,j,k) = (div**2) / (div**2 + (rot(1)**2 + rot(2)**2 + rot(3)**2) + eps)

    fd(i,j,k) = min(1.d0, fd(i,j,k))

    ! boundary
    ! x direction
    if (i == 2) then
      fd(1,j,k) = fd(2,j,k)
    elseif (i == nx-1) then
      fd(nx,j,k) = fd(nx-1,j,k)
    endif
    ! y direction
    if (j == 2) then
      fd(i,1,k) = fd(i,2,k)
    elseif (j == ny-1) then
      fd(i,ny,k) = fd(i,ny-1,k)
    endif
    ! z direction
    if (k == 2) then
      fd(i,j,1) = fd(i,j,2)
    elseif (k == nz-1) then
      fd(i,j,nz) = fd(i,j,nz-1)
    endif
  end subroutine calc_Ducros


  attributes(device) function wiggle_detector(phi) result(ans)
    real(8), intent(in), device :: phi(4)
    real(8) ans, phi1, phi2
    phi1 = (-phi(1) + phi(2)) * (-phi(2) + phi(3))
    phi2 = (-phi(3) + phi(4)) * (-phi(2) + phi(3))
    ans  = 0.5d0 * (1.d0 - sign(1.d0, min(phi1, phi2)))
  end function wiggle_detector
end module calc_hybrid

