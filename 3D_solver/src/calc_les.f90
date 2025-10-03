module calc_les
  use mod_constant, only : one_third
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_differential(u,v,w,dx,dy,dz,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w
    real(8), intent(in), value                    :: dx, dy, dz ! 1 /dx, 1 / dy, 1 / dz
    real(8), intent(out) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    dudx = 0.5d0 * (-u(1,2,2) + u(3,2,2)) * dx
    dudy = 0.5d0 * (-u(2,1,2) + u(2,3,2)) * dy
    dudz = 0.5d0 * (-u(2,2,1) + u(2,2,3)) * dz
    dvdx = 0.5d0 * (-v(1,2,2) + v(3,2,2)) * dx
    dvdy = 0.5d0 * (-v(2,1,2) + v(2,3,2)) * dy
    dvdz = 0.5d0 * (-v(2,2,1) + v(2,2,3)) * dz
    dwdx = 0.5d0 * (-w(1,2,2) + w(3,2,2)) * dx
    dwdy = 0.5d0 * (-w(2,1,2) + w(2,3,2)) * dy
    dwdz = 0.5d0 * (-w(2,2,1) + w(2,2,3)) * dz
  end subroutine calc_differential


  !$dir inline
  attributes(device) function SMS(u, v, w, uh, vh, wh, dx, dy, dz, qc2) result(nut)
    real(8), intent(in), dimension(3,3,3), device :: u, v, w, uh, vh, wh
    real(8), intent(in), value                    :: dx, dy, dz, qc2 ! 1 / dx, 1 / dy, 1 / dz
    real(8) nut, ftheta, vor2, vorh2, vord2, S2
    real(8), parameter :: Cm = 0.06d0, alpha = 0.5d0
    real(8), parameter :: over_tan = 1.d0 / (tan(10.d0 * acos(-1.d0) / 180.d0))**2
    block
      real(8) vor(3), vorh(3), vord(3)
      block
        real(8) dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        call calc_differential(u, v, w, dx, dy, dz, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
        vor(1) = dwdy - dvdz
        vor(2) = dudz - dwdx
        vor(3) = dvdx - dudy
        S2     = 2.d0 * (dudx**2 + dvdy**2 + dwdz**2) + ((dvdx + dudy)**2 + (dwdy + dvdz)**2 + (dudz + dwdx)**2)
      end block
      vor2    = vor(1)**2  + vor(2)**2  + vor(3)**2
      block
        real(8) dudxh, dudyh, dudzh, dvdxh, dvdyh, dvdzh, dwdxh, dwdyh, dwdzh
        call calc_differential(uh, vh, wh, dx, dy, dz, dudxh, dudyh, dudzh, dvdxh, dvdyh, dvdzh, dwdxh, dwdyh, dwdzh)
        vorh(1) = dwdyh - dvdzh
        vorh(2) = dudzh - dwdxh
        vorh(3) = dvdxh - dudyh
      end block
      vorh2   = vorh(1)**2 + vorh(2)**2 + vorh(3)**2
      vord(:) = vor(:) - vorh(:)
      vord2   = vord(1)**2 + vord(2)**2 + vord(3)**2
    end block
    block
      real(8) svor2, tan2
      svor2 = 2.d0 * sqrt(vor2) * sqrt(vorh2) + vorh2 + vor2 - vord2
      if (abs(svor2) < 1.d-30) then
        tan2 = 0.d0
      else
        tan2 = (2.d0 * sqrt(vorh2) * sqrt(vor2) - vorh2 - vor2 + vord2) / svor2
      endif
      block
        real(8) rtheta
        rtheta = tan2 * over_tan
        ftheta = min(1.d0, rtheta*rtheta)
      end block
    end block
    ! aspect ratio
    block
      real(8) a1, a2, delta
      a1    = min(dx, dy, dz) / dx
      a2    = min(dx, dy, dz) / dy
      delta = cosh(sqrt(4.d0 * ((log(a1))**2 - log(a1) * log(a2) + (log(a2))**2) / 27.d0)) * (dx * dy * dz)**(-one_third) 
      nut   = ftheta * Cm * (S2**(0.5d0 * alpha)) * (qc2**(0.5d0 * (1.d0 - alpha))) * (delta**(1.d0 + alpha))
    end block
  end function SMS


  attributes(device) function stride_filter(x) result(xh)
    real(8), intent(in), device :: x(5,5,5)
    integer i, j, k
    real(8) xh(3,3,3), x0
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          x0 = x(i+1,j+1,k+1)
          xh(i,j,k) = x0 + 0.25d0 * (x(i,j+1,k+1) - 2.d0 * x0 + x(i+2,j+1,k+1) &
                                   + x(i+1,j,k+1) - 2.d0 * x0 + x(i+1,j+2,k+1) &
                                   + x(i+1,j+1,k) - 2.d0 * x0 + x(i+1,j+1,k+2))
    enddo;enddo;enddo
  end function stride_filter


  attributes(global) subroutine calc_mut(nx, ny, nz, dx, dy, dz, Q, mut, qc2)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(in), device  :: Q(5,nx,ny,nz)
    real(8), intent(out), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    real(8), dimension(3,3,3), device :: u3, v3, w3, uh, vh, wh
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (i <= nx-1 .and. j <= ny-1 .and. k <= nz-1) return
    u3 = Q(2,i-1:i+1,j-1:j+1,k-1:k+1)
    v3 = Q(3,i-1:i+1,j-1:j+1,k-1:k+1)
    w3 = Q(4,i-1:i+1,j-1:j+1,k-1:k+1)
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8) u5(5,5,5)
        u5 = Q(2,i-2:i+2,j-2:j+2,k-2:k+2)
        uh = stride_filter(u5)
      end block
      block
        real(8) v5(5,5,5)
        v5 = Q(3,i-2:i+2,j-2:j+2,k-2:k+2)
        vh = stride_filter(v5)
      end block
      block
        real(8) w5(5,5,5)
        w5 = Q(4,i-2:i+2,j-2:j+2,k-2:k+2)
        wh = stride_filter(w5)
      end block
    else
      uh = u3
      vh = v3
      wh = w3
    endif
    qc2(i,j,k) = 0.5d0 * ((u3(2,2,2) - uh(2,2,2))**2 + (v3(2,2,2) - vh(2,2,2))**2 + (w3(2,2,2) - wh(2,2,2))**2)
    mut(i,j,k) = Q(1,i,j,k) * SMS(u3, v3, w3, uh, vh, wh, dx(i), dy(j), dz(k), qc2(i,j,k))
  end subroutine calc_mut
end module calc_les

