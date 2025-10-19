module calc_visc
  use cudafor
  use mod_globals, only : threads, blocks
  use mod_constant, only : Cp_over_Pr, one_third, two_third, four_third, one_twelfth
  use store_shared
  implicit none
  real(8), parameter :: over_144 = 1.d0 / 144.d0
contains
  attributes(global) subroutine calc_visc_Laplacian2nd(nx, ny, nz, over_dxs, over_dys, over_dzs, over_Jacobian, Q, T, mud, Rv)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: over_dxs(nx-1)
    real(8), intent(in), device  :: over_dys(ny-1)
    real(8), intent(in), device  :: over_dzs(nz-1)
    real(8), intent(in), device  :: over_Jacobian(nx,ny)
    real(8), intent(in), device  :: Q(5,nx,ny,nz), T(nx,ny,nz), mud(nx,ny,nz)
    real(8), intent(out), device :: Rv(4,nx-2,ny-2,nz-2)
    real(8), shared ::  u(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8), shared ::  v(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8), shared ::  w(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8), shared :: mu(0:threads%x+1,0:threads%y+1,0:threads%z+1)
    real(8) mx, my, mz, ux, vy, wz, ux_vy_wz, vy_wz_ux, wz_ux_vy, vx_uy, wy_vz, uz_wx
    real(8) over_dx, over_dy, over_dz, over_dx2, over_dy2, over_dz2
    integer i, j, k, it, jt, kt
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    call store_shared_visc_2nd(nx, ny, nz, it, jt, kt, Q, mud, u, v, w, mu)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    over_dx = over_dxs(i)
    over_dy = over_dys(j)
    over_dz = over_dzs(k)
    over_dx2 = over_dx*over_dx
    over_dy2 = over_dy*over_dy
    over_dz2 = over_dz*over_dz
    mx = 0.5d0 * (-mu(it-1,jt,kt) + mu(it+1,jt,kt)) * over_dx
    my = 0.5d0 * (-mu(it,jt-1,kt) + mu(it,jt+1,kt)) * over_dy
    mz = 0.5d0 * (-mu(it,jt,kt-1) + mu(it,jt,kt+1)) * over_dz
    ux = 0.5d0 * (-u(it-1,jt,kt) + u(it+1,jt,kt)) * over_dx
    vy = 0.5d0 * (-v(it,jt-1,kt) + v(it,jt+1,kt)) * over_dy
    wz = 0.5d0 * (-w(it,jt,kt-1) + w(it,jt,kt+1)) * over_dz
    ux_vy_wz = 2.d0 * ux - vy - wz
    vy_wz_ux = 2.d0 * vy - wz - ux
    wz_ux_vy = 2.d0 * wz - ux - vy
    vx_uy = 0.5d0 * (-v(it-1,jt,kt) + v(it+1,jt,kt)) * over_dx &
          + 0.5d0 * (-u(it,jt-1,kt) + u(it,jt+1,kt)) * over_dy
    wy_vz = 0.5d0 * (-w(it,jt-1,kt) + w(it,jt+1,kt)) * over_dy &
          + 0.5d0 * (-v(it,jt,kt-1) + v(it,jt,kt+1)) * over_dz
    uz_wx = 0.5d0 * (-u(it,jt,kt-1) + u(it,jt,kt+1)) * over_dz &
          + 0.5d0 * (-w(it-1,jt,kt) + w(it+1,jt,kt)) * over_dx
    block
      real(8) d2vover_dxover_dy, d2wover_dxover_dz, d2uover_dx2, d2uover_dy2, d2uover_dz2
      d2vover_dxover_dy = 0.25d0 * (v(it-1,jt-1,kt) - v(it-1,jt+1,kt) - v(it+1,jt-1,kt) + v(it+1,jt+1,kt)) * over_dx * over_dy
      d2wover_dxover_dz = 0.25d0 * (w(it-1,jt,kt-1) - w(it-1,jt,kt+1) - w(it+1,jt,kt-1) + w(it+1,jt,kt+1)) * over_dx * over_dz
      d2uover_dx2  = (u(it-1,jt,kt) - 2.d0 * u(it,jt,kt) + u(it+1,jt,kt)) * over_dx2
      d2uover_dy2  = (u(it,jt-1,kt) - 2.d0 * u(it,jt,kt) + u(it,jt+1,kt)) * over_dy2
      d2uover_dz2  = (u(it,jt,kt-1) - 2.d0 * u(it,jt,kt) + u(it,jt,kt+1)) * over_dz2
      Rv(1,i-1,j-1,k-1) = two_third * mx * ux_vy_wz + my * vx_uy + mz * uz_wx &
                          + mu(it,jt,kt) * (four_third * d2uover_dx2 + d2uover_dy2 + d2uover_dz2 + one_third * (d2vover_dxover_dy + d2wover_dxover_dz))
    end block
    block
      real(8) d2wover_dyover_dz, d2uover_dyover_dx, d2vover_dy2, d2vover_dz2, d2vover_dx2
      d2wover_dyover_dz = 0.25d0 * (w(it,jt-1,kt-1) - w(it,jt-1,kt+1) - w(it,jt+1,kt-1) + w(it,jt+1,kt+1)) * over_dy * over_dz
      d2uover_dyover_dx = 0.25d0 * (u(it-1,jt-1,kt) - u(it+1,jt-1,kt) - u(it-1,jt+1,kt) + u(it+1,jt+1,kt)) * over_dy * over_dx
      d2vover_dy2  = (v(it,jt-1,kt) - 2.d0 * v(it,jt,kt) + v(it,jt+1,kt)) * over_dy2
      d2vover_dz2  = (v(it,jt,kt-1) - 2.d0 * v(it,jt,kt) + v(it,jt,kt+1)) * over_dz2
      d2vover_dx2  = (v(it-1,jt,kt) - 2.d0 * v(it,jt,kt) + v(it+1,jt,kt)) * over_dx2
      Rv(2,i-1,j-1,k-1) = two_third * my * vy_wz_ux + mz * wy_vz + mx * vx_uy &
                          + mu(it,jt,kt) * (four_third * d2vover_dy2 + d2vover_dz2 + d2vover_dx2 + one_third * (d2wover_dyover_dz + d2uover_dyover_dx))
    end block
    block
      real(8) d2uover_dzover_dx, d2vover_dzover_dy, d2wover_dz2, d2wover_dx2, d2wover_dy2
      d2uover_dzover_dx = 0.25d0 * (u(it-1,jt,kt-1) - u(it+1,jt,kt-1) - u(it-1,jt,kt+1) + u(it+1,jt,kt+1)) * over_dz * over_dx
      d2vover_dzover_dy = 0.25d0 * (v(it,jt-1,kt-1) - v(it,jt+1,kt-1) - v(it,jt-1,kt+1) + v(it,jt+1,kt+1)) * over_dz * over_dy
      d2wover_dz2  = (w(it,jt,kt-1) - 2.d0 * w(it,jt,kt) + w(it,jt,kt+1)) * over_dz2
      d2wover_dx2  = (w(it-1,jt,kt) - 2.d0 * w(it,jt,kt) + w(it+1,jt,kt)) * over_dx2
      d2wover_dy2  = (w(it,jt-1,kt) - 2.d0 * w(it,jt,kt) + w(it,jt+1,kt)) * over_dy2
      Rv(3,i-1,j-1,k-1) = two_third * mz * wz_ux_vy + mx * uz_wx + my * wy_vz &
                          + mu(it,jt,kt) * (four_third * d2wover_dz2 + d2wover_dx2 + d2wover_dy2 + one_third * (d2uover_dzover_dx + d2vover_dzover_dy))
    end block
    block
      real(8) d2Tdx2, d2Tdy2, d2Tdz2, dTdx, dTdy, dTdz, dmdx, dmdy, dmdz
      d2Tdx2 = (T(i-1,j,k) - 2.d0 * T(i,j,k) + T(i+1,j,k)) * over_dx2
      d2Tdy2 = (T(i,j-1,k) - 2.d0 * T(i,j,k) + T(i,j+1,k)) * over_dy2
      d2Tdz2 = (T(i,j,k-1) - 2.d0 * T(i,j,k) + T(i,j,k+1)) * over_dz2
      dTdx   = 0.5d0 * (-T(i-1,j,k) + T(i+1,j,k)) * over_dx
      dTdy   = 0.5d0 * (-T(i,j-1,k) + T(i,j+1,k)) * over_dy
      dTdz   = 0.5d0 * (-T(i,j,k-1) + T(i,j,k+1)) * over_dz
      dmdx   = 0.5d0 * (-mu(it-1,jt,kt) + mu(it+1,jt,kt)) * over_dx
      dmdy   = 0.5d0 * (-mu(it,jt-1,kt) + mu(it,jt+1,kt)) * over_dy
      dmdz   = 0.5d0 * (-mu(it,jt,kt-1) + mu(it,jt,kt+1)) * over_dz
      Rv(4,i-1,j-1,k-1) = u(it,jt,kt)*Rv(1,i-1,j-1,k-1) + v(it,jt,kt)*Rv(2,i-1,j-1,k-1) + w(it,jt,kt)*Rv(3,i-1,j-1,k-1) &
                        + mu(it,jt,kt) * (two_third * (ux * ux_vy_wz + vy * vy_wz_ux + wz * wz_ux_vy) &
                        + vx_uy*vx_uy + wy_vz*wy_vz + uz_wx*uz_wx) &
                        + Cp_over_Pr * (mu(it,jt,kt) * (d2Tdx2 + d2Tdy2 + d2Tdz2) + dmdx * dTdx + dmdy * dTdy + dmdz * dTdz)
    end block
    Rv(:,i-1,j-1,k-1) = Rv(:,i-1,j-1,k-1) * over_Jacobian(i,j)
  end subroutine calc_visc_Laplacian2nd

  
  attributes(global) subroutine calc_visc_Laplacian4th(nx, ny, nz, over_dxs, over_dys, over_dzs, over_Jacobian, Q, T, mud, Rv)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: over_dxs(nx-1)
    real(8), intent(in), device  :: over_dys(ny-1)
    real(8), intent(in), device  :: over_dzs(nz-1)
    real(8), intent(in), device  :: over_Jacobian(nx,ny)
    real(8), intent(in), device  :: Q(5,nx,ny,nz), T(nx,ny,nz), mud(nx,ny,nz)
    real(8), intent(out), device :: Rv(4,nx-2,ny-2,nz-2)
    real(8), shared ::  u(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8), shared ::  v(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8), shared ::  w(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8), shared :: mu(-1:threads%x+2,-1:threads%y+2,-1:threads%z+2)
    real(8) mx, my, mz, ux, vy, wz, ux_vy_wz, vy_wz_ux, wz_ux_vy, vx_uy, wy_vz, uz_wx
    real(8) over_dx, over_dy, over_dz, over_dx2, over_dy2, over_dz2
    integer i, j, k, it, jt, kt
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    call store_shared_visc_4th(nx, ny, nz, it, jt, kt, Q, mud, u, v, w, mu)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    over_dx = over_dxs(i)
    over_dy = over_dys(j)
    over_dz = over_dzs(k)
    over_dx2 = over_dx*over_dx
    over_dy2 = over_dy*over_dy
    over_dz2 = over_dz*over_dz
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      mx = one_twelfth * (mu(it-2,jt,kt) - 8.d0 * (mu(it-1,jt,kt) - mu(it+1,jt,kt)) - mu(it+2,jt,kt)) * over_dx
      my = one_twelfth * (mu(it,jt-2,kt) - 8.d0 * (mu(it,jt-1,kt) - mu(it,jt+1,kt)) - mu(it,jt+2,kt)) * over_dy
      mz = one_twelfth * (mu(it,jt,kt-2) - 8.d0 * (mu(it,jt,kt-1) - mu(it,jt,kt+1)) - mu(it,jt,kt+2)) * over_dz
      ux = one_twelfth * ( u(it-2,jt,kt) - 8.d0 * ( u(it-1,jt,kt) -  u(it+1,jt,kt)) -  u(it+2,jt,kt)) * over_dx
      vy = one_twelfth * ( v(it,jt-2,kt) - 8.d0 * ( v(it,jt-1,kt) -  v(it,jt+1,kt)) -  v(it,jt+2,kt)) * over_dy
      wz = one_twelfth * ( w(it,jt,kt-2) - 8.d0 * ( w(it,jt,kt-1) -  w(it,jt,kt+1)) -  w(it,jt,kt+2)) * over_dz
      ux_vy_wz = 2.d0 * ux - vy - wz
      vy_wz_ux = 2.d0 * vy - wz - ux
      wz_ux_vy = 2.d0 * wz - ux - vy
      vx_uy = one_twelfth * (v(it-2,jt,kt) - 8.d0 * (v(it-1,jt,kt) - v(it+1,jt,kt)) - v(it+2,jt,kt)) * over_dx &
            + one_twelfth * (u(it,jt-2,kt) - 8.d0 * (u(it,jt-1,kt) - u(it,jt+1,kt)) - u(it,jt+2,kt)) * over_dy
      wy_vz = one_twelfth * (w(it,jt-2,kt) - 8.d0 * (w(it,jt-1,kt) - w(it,jt+1,kt)) - w(it,jt+2,kt)) * over_dy &
            + one_twelfth * (v(it,jt,kt-2) - 8.d0 * (v(it,jt,kt-1) - v(it,jt,kt+1)) - v(it,jt,kt+2)) * over_dz
      uz_wx = one_twelfth * (u(it,jt,kt-2) - 8.d0 * (u(it,jt,kt-1) - u(it,jt,kt+1)) - u(it,jt,kt+2)) * over_dz &
            + one_twelfth * (w(it-2,jt,kt) - 8.d0 * (w(it-1,jt,kt) - w(it+1,jt,kt)) - w(it+2,jt,kt)) * over_dx
      block
        real(8) d2vdxdy, d2wdxdz, d2udx2, d2udy2, d2udz2
        d2vdxdy = (         (v(it-2,jt-2,kt) - 8.d0 * (v(it-1,jt-2,kt) - v(it+1,jt-2,kt)) - v(it+2,jt-2,kt)) &
                 & - 8.d0 * (v(it-2,jt-1,kt) - 8.d0 * (v(it-1,jt-1,kt) - v(it+1,jt-1,kt)) - v(it+2,jt-1,kt)) &
                 & + 8.d0 * (v(it-2,jt+1,kt) - 8.d0 * (v(it-1,jt+1,kt) - v(it+1,jt+1,kt)) - v(it+2,jt+1,kt)) &
                 &        - (v(it-2,jt+2,kt) - 8.d0 * (v(it-1,jt+2,kt) - v(it+1,jt+2,kt)) - v(it+2,jt+2,kt))) * over_144 * over_dx * over_dy
        d2wdxdz = (         (w(it-2,jt,kt-2) - 8.d0 * (w(it-1,jt,kt-2) - w(it+1,jt,kt-2)) - w(it+2,jt,kt-2)) &
                 & - 8.d0 * (w(it-2,jt,kt-1) - 8.d0 * (w(it-1,jt,kt-1) - w(it+1,jt,kt-1)) - w(it+2,jt,kt-1)) &
                 & + 8.d0 * (w(it-2,jt,kt+1) - 8.d0 * (w(it-1,jt,kt+1) - w(it+1,jt,kt+1)) - w(it+2,jt,kt+1)) &
                 &        - (w(it-2,jt,kt+2) - 8.d0 * (w(it-1,jt,kt+2) - w(it+1,jt,kt+2)) - w(it+2,jt,kt+2))) * over_144 * over_dx * over_dz
        d2udx2  = one_twelfth * (-u(it-2,jt,kt) + 16.d0 * u(it-1,jt,kt) - 30.d0 * u(it,jt,kt) + 16.d0 * u(it+1,jt,kt) - u(it+2,jt,kt)) * over_dx2
        d2udy2  = one_twelfth * (-u(it,jt-2,kt) + 16.d0 * u(it,jt-1,kt) - 30.d0 * u(it,jt,kt) + 16.d0 * u(it,jt+1,kt) - u(it,jt+2,kt)) * over_dy2
        d2udz2  = one_twelfth * (-u(it,jt,kt-2) + 16.d0 * u(it,jt,kt-1) - 30.d0 * u(it,jt,kt) + 16.d0 * u(it,jt,kt+1) - u(it,jt,kt+2)) * over_dz2
        Rv(1,i-1,j-1,k-1) = two_third * mx * ux_vy_wz + my * vx_uy + mz * uz_wx &
                            + mu(it,jt,kt) * (four_third * d2udx2 + d2udy2 + d2udz2 + one_third * (d2vdxdy + d2wdxdz))
      end block
      block
        real(8) d2wdydz, d2udydx, d2vdy2, d2vdz2, d2vdx2
        d2wdydz = (         (w(it,jt-2,kt-2) - 8.d0 * (w(it,jt-1,kt-2) - w(it,jt+1,kt-2)) - w(it,jt+2,kt-2)) &
                 & - 8.d0 * (w(it,jt-2,kt-1) - 8.d0 * (w(it,jt-1,kt-1) - w(it,jt+1,kt-1)) - w(it,jt+2,kt-1)) &
                 & + 8.d0 * (w(it,jt-2,kt+1) - 8.d0 * (w(it,jt-1,kt+1) - w(it,jt+1,kt+1)) - w(it,jt+2,kt+1)) &
                 &        - (w(it,jt-2,kt+2) - 8.d0 * (w(it,jt-1,kt+2) - w(it,jt+1,kt+2)) - w(it,jt+2,kt+2))) * over_144 * over_dy * over_dz
        d2udydx = (         (u(it-2,jt-2,kt) - 8.d0 * (u(it-2,jt-1,kt) - u(it-2,jt+1,kt)) - u(it-2,jt+2,kt)) &
                 & - 8.d0 * (u(it-1,jt-2,kt) - 8.d0 * (u(it-1,jt-1,kt) - u(it-1,jt+1,kt)) - u(it-1,jt+2,kt)) &
                 & + 8.d0 * (u(it+1,jt-2,kt) - 8.d0 * (u(it+1,jt-1,kt) - u(it+1,jt+1,kt)) - u(it+1,jt+2,kt)) &
                 &        - (u(it+2,jt-2,kt) - 8.d0 * (u(it+2,jt-1,kt) - u(it+2,jt+1,kt)) - u(it+2,jt+2,kt))) * over_144 * over_dy * over_dx
        d2vdy2  = one_twelfth * (-v(it,jt-2,kt) + 16.d0 * v(it,jt-1,kt) - 30.d0 * v(it,jt,kt) + 16.d0 * v(it,jt+1,kt) - v(it,jt+2,kt)) * over_dy2
        d2vdz2  = one_twelfth * (-v(it,jt,kt-2) + 16.d0 * v(it,jt,kt-1) - 30.d0 * v(it,jt,kt) + 16.d0 * v(it,jt,kt+1) - v(it,jt,kt+2)) * over_dz2
        d2vdx2  = one_twelfth * (-v(it-2,jt,kt) + 16.d0 * v(it-1,jt,kt) - 30.d0 * v(it,jt,kt) + 16.d0 * v(it+1,jt,kt) - v(it+2,jt,kt)) * over_dx2
        Rv(2,i-1,j-1,k-1) = two_third * my * vy_wz_ux + mz * wy_vz + mx * vx_uy &
                            + mu(it,jt,kt) * (four_third * d2vdy2 + d2vdz2 + d2vdx2 + one_third * (d2wdydz + d2udydx))
      end block
      block
        real(8) d2udzdx, d2vdzdy, d2wdz2, d2wdx2, d2wdy2
        d2udzdx = (         (u(it-2,jt,kt-2) - 8.d0 * (u(it-1,jt,kt-2) - u(it+1,jt,kt-2)) - u(it+2,jt,kt-2)) &
                 & - 8.d0 * (u(it-2,jt,kt-1) - 8.d0 * (u(it-1,jt,kt-1) - u(it+1,jt,kt-1)) - u(it+2,jt,kt-1)) &
                 & + 8.d0 * (u(it-2,jt,kt+1) - 8.d0 * (u(it-1,jt,kt+1) - u(it+1,jt,kt+1)) - u(it+2,jt,kt+1)) &
                 &        - (u(it-2,jt,kt+2) - 8.d0 * (u(it-1,jt,kt+2) - u(it+1,jt,kt+2)) - u(it+2,jt,kt+2))) * over_144 * over_dz * over_dx
        d2vdzdy = (         (v(it,jt-2,kt-2) - 8.d0 * (v(it,jt-1,kt-2) - v(it,jt+1,kt-2)) - v(it,jt+2,kt-2)) &
                 & - 8.d0 * (v(it,jt-2,kt-1) - 8.d0 * (v(it,jt-1,kt-1) - v(it,jt+1,kt-1)) - v(it,jt+2,kt-1)) &
                 & + 8.d0 * (v(it,jt-2,kt+1) - 8.d0 * (v(it,jt-1,kt+1) - v(it,jt+1,kt+1)) - v(it,jt+2,kt+1)) &
                 &        - (v(it,jt-2,kt+2) - 8.d0 * (v(it,jt-1,kt+2) - v(it,jt+1,kt+2)) - v(it,jt+2,kt+2))) * over_144 * over_dz * over_dy
        d2wdz2  = one_twelfth * (-w(it,jt,kt-2) + 16.d0 * w(it,jt,kt-1) - 30.d0 * w(it,jt,kt) + 16.d0 * w(it,jt,kt+1) - w(it,jt,kt+2)) * over_dz2
        d2wdx2  = one_twelfth * (-w(it-2,jt,kt) + 16.d0 * w(it-1,jt,kt) - 30.d0 * w(it,jt,kt) + 16.d0 * w(it+1,jt,kt) - w(it+2,jt,kt)) * over_dx2
        d2wdy2  = one_twelfth * (-w(it,jt-2,kt) + 16.d0 * w(it,jt-1,kt) - 30.d0 * w(it,jt,kt) + 16.d0 * w(it,jt+1,kt) - w(it,jt+2,kt)) * over_dy2
        Rv(3,i-1,j-1,k-1) = two_third * mz * wz_ux_vy + mx * uz_wx + my * wy_vz &
                            + mu(it,jt,kt) * (four_third * d2wdz2 + d2wdx2 + d2wdy2 + one_third * (d2udzdx + d2vdzdy))
      end block
      block
        real(8) d2Tdx2, d2Tdy2, d2Tdz2, dTdx, dTdy, dTdz, dmdx, dmdy, dmdz
        d2Tdx2 = one_twelfth * (-T(i-2,j,k) + 16.d0 * T(i-1,j,k) - 30.d0 * T(i,j,k) + 16.d0 * T(i+1,j,k) - T(i+2,j,k)) * over_dx2
        d2Tdy2 = one_twelfth * (-T(i,j-2,k) + 16.d0 * T(i,j-1,k) - 30.d0 * T(i,j,k) + 16.d0 * T(i,j+1,k) - T(i,j+2,k)) * over_dy2
        d2Tdz2 = one_twelfth * (-T(i,j,k-2) + 16.d0 * T(i,j,k-1) - 30.d0 * T(i,j,k) + 16.d0 * T(i,j,k+1) - T(i,j,k+2)) * over_dz2
        dTdx   = one_twelfth * ( T(i-2,j,k) - 8.d0 * ( T(i-1,j,k) -  T(i+1,j,k)) -  T(i+2,j,k)) * over_dx
        dTdy   = one_twelfth * ( T(i,j-2,k) - 8.d0 * ( T(i,j-1,k) -  T(i,j+1,k)) -  T(i,j+2,k)) * over_dy
        dTdz   = one_twelfth * ( T(i,j,k-2) - 8.d0 * ( T(i,j,k-1) -  T(i,j,k+1)) -  T(i,j,k+2)) * over_dz
        dmdx   = one_twelfth * (mu(it-2,jt,kt) - 8.d0 * (mu(it-1,jt,kt) - mu(it+1,jt,kt)) - mu(it+2,jt,kt)) * over_dx
        dmdy   = one_twelfth * (mu(it,jt-2,kt) - 8.d0 * (mu(it,jt-1,kt) - mu(it,jt+1,kt)) - mu(it,jt+2,kt)) * over_dy
        dmdz   = one_twelfth * (mu(it,jt,kt-2) - 8.d0 * (mu(it,jt,kt-1) - mu(it,jt,kt+1)) - mu(it,jt,kt+2)) * over_dz
        Rv(4,i-1,j-1,k-1) = u(it,jt,kt)*Rv(1,i-1,j-1,k-1) + v(it,jt,kt)*Rv(2,i-1,j-1,k-1) + w(it,jt,kt)*Rv(3,i-1,j-1,k-1) &
                          + mu(it,jt,kt) * (two_third * (ux * ux_vy_wz + vy * vy_wz_ux + wz * wz_ux_vy) &
                          + vx_uy*vx_uy + wy_vz*wy_vz + uz_wx*uz_wx) &
                          + Cp_over_Pr * (mu(it,jt,kt) * (d2Tdx2 + d2Tdy2 + d2Tdz2) + dmdx * dTdx + dmdy * dTdy + dmdz * dTdz)
      end block
    else
      mx = 0.5d0 * (-mu(it-1,jt,kt) + mu(it+1,jt,kt)) * over_dx
      my = 0.5d0 * (-mu(it,jt-1,kt) + mu(it,jt+1,kt)) * over_dy
      mz = 0.5d0 * (-mu(it,jt,kt-1) + mu(it,jt,kt+1)) * over_dz
      ux = 0.5d0 * ( -u(it-1,jt,kt) +  u(it+1,jt,kt)) * over_dx
      vy = 0.5d0 * ( -v(it,jt-1,kt) +  v(it,jt+1,kt)) * over_dy
      wz = 0.5d0 * ( -w(it,jt,kt-1) +  w(it,jt,kt+1)) * over_dz
      ux_vy_wz = 2.d0 * ux - vy - wz
      vy_wz_ux = 2.d0 * vy - wz - ux
      wz_ux_vy = 2.d0 * wz - ux - vy
      vx_uy = 0.5d0 * (-v(it-1,jt,kt) + v(it+1,jt,kt)) * over_dx &
            + 0.5d0 * (-u(it,jt-1,kt) + u(it,jt+1,kt)) * over_dy
      wy_vz = 0.5d0 * (-w(it,jt-1,kt) + w(it,jt+1,kt)) * over_dy &
            + 0.5d0 * (-v(it,jt,kt-1) + v(it,jt,kt+1)) * over_dz
      uz_wx = 0.5d0 * (-u(it,jt,kt-1) + u(it,jt,kt+1)) * over_dz &
            + 0.5d0 * (-w(it-1,jt,kt) + w(it+1,jt,kt)) * over_dx
      block
        real(8) d2vdxdy, d2wdxdz, d2udx2, d2udy2, d2udz2
        d2vdxdy = 0.25d0 * (v(it-1,jt-1,kt) - v(it-1,jt+1,kt) - v(it+1,jt-1,kt) + v(it+1,jt+1,kt)) * over_dx * over_dy
        d2wdxdz = 0.25d0 * (w(it-1,jt,kt-1) - w(it-1,jt,kt+1) - w(it+1,jt,kt-1) + w(it+1,jt,kt+1)) * over_dx * over_dz
        d2udx2  = (u(it-1,jt,kt) - 2.d0 * u(it,jt,kt) + u(it+1,jt,kt)) * over_dx2
        d2udy2  = (u(it,jt-1,kt) - 2.d0 * u(it,jt,kt) + u(it,jt+1,kt)) * over_dy2
        d2udz2  = (u(it,jt,kt-1) - 2.d0 * u(it,jt,kt) + u(it,jt,kt+1)) * over_dz2
        Rv(1,i-1,j-1,k-1) = two_third * mx * ux_vy_wz + my * vx_uy + mz * uz_wx &
                            + mu(it,jt,kt) * (four_third * d2udx2 + d2udy2 + d2udz2 + one_third * (d2vdxdy + d2wdxdz))
      end block
      block
        real(8) d2wdydz, d2udydx, d2vdy2, d2vdz2, d2vdx2
        d2wdydz = 0.25d0 * (w(it,jt-1,kt-1) - w(it,jt-1,kt+1) - w(it,jt+1,kt-1) + w(it,jt+1,kt+1)) * over_dy * over_dz
        d2udydx = 0.25d0 * (u(it-1,jt-1,kt) - u(it+1,jt-1,kt) - u(it-1,jt+1,kt) + u(it+1,jt+1,kt)) * over_dy * over_dx
        d2vdy2  = (v(it,jt-1,kt) - 2.d0 * v(it,jt,kt) + v(it,jt+1,kt)) * over_dy2
        d2vdz2  = (v(it,jt,kt-1) - 2.d0 * v(it,jt,kt) + v(it,jt,kt+1)) * over_dz2
        d2vdx2  = (v(it-1,jt,kt) - 2.d0 * v(it,jt,kt) + v(it+1,jt,kt)) * over_dx2
        Rv(2,i-1,j-1,k-1) = two_third * my * vy_wz_ux + mz * wy_vz + mx * vx_uy &
                            + mu(it,jt,kt) * (four_third * d2vdy2 + d2vdz2 + d2vdx2 + one_third * (d2wdydz + d2udydx))
      end block
      block
        real(8) d2udzdx, d2vdzdy, d2wdz2, d2wdx2, d2wdy2
        d2udzdx = 0.25d0 * (u(it-1,jt,kt-1) - u(it+1,jt,kt-1) - u(it-1,jt,kt+1) + u(it+1,jt,kt+1)) * over_dz * over_dx
        d2vdzdy = 0.25d0 * (v(it,jt-1,kt-1) - v(it,jt+1,kt-1) - v(it,jt-1,kt+1) + v(it,jt+1,kt+1)) * over_dz * over_dy
        d2wdz2  = (w(it,jt,kt-1) - 2.d0 * w(it,jt,kt) + w(it,jt,kt+1)) * over_dz2
        d2wdx2  = (w(it-1,jt,kt) - 2.d0 * w(it,jt,kt) + w(it+1,jt,kt)) * over_dx2
        d2wdy2  = (w(it,jt-1,kt) - 2.d0 * w(it,jt,kt) + w(it,jt+1,kt)) * over_dy2
        Rv(3,i-1,j-1,k-1) = two_third * mz * wz_ux_vy + mx * uz_wx + my * wy_vz &
                            + mu(it,jt,kt) * (four_third * d2wdz2 + d2wdx2 + d2wdy2 + one_third * (d2udzdx + d2vdzdy))
      end block
      block
        real(8) d2Tdx2, d2Tdy2, d2Tdz2, dTdx, dTdy, dTdz, dmdx, dmdy, dmdz
        d2Tdx2 = (T(i-1,j,k) - 2.d0 * T(i,j,k) + T(i+1,j,k)) * over_dx2
        d2Tdy2 = (T(i,j-1,k) - 2.d0 * T(i,j,k) + T(i,j+1,k)) * over_dy2
        d2Tdz2 = (T(i,j,k-1) - 2.d0 * T(i,j,k) + T(i,j,k+1)) * over_dz2
        dTdx   = 0.5d0 * (-T(i-1,j,k) + T(i+1,j,k)) * over_dx
        dTdy   = 0.5d0 * (-T(i,j-1,k) + T(i,j+1,k)) * over_dy
        dTdz   = 0.5d0 * (-T(i,j,k-1) + T(i,j,k+1)) * over_dz
        dmdx   = 0.5d0 * (-mu(it-1,jt,kt) + mu(it+1,jt,kt)) * over_dx
        dmdy   = 0.5d0 * (-mu(it,jt-1,kt) + mu(it,jt+1,kt)) * over_dy
        dmdz   = 0.5d0 * (-mu(it,jt,kt-1) + mu(it,jt,kt+1)) * over_dz
        Rv(4,i-1,j-1,k-1) = u(it,jt,kt)*Rv(1,i-1,j-1,k-1) + v(it,jt,kt)*Rv(2,i-1,j-1,k-1) + w(it,jt,kt)*Rv(3,i-1,j-1,k-1) &
                          + mu(it,jt,kt) * (two_third * (ux * ux_vy_wz + vy * vy_wz_ux + wz * wz_ux_vy) &
                          + vx_uy*vx_uy + wy_vz*wy_vz + uz_wx*uz_wx) &
                          + Cp_over_Pr * (mu(it,jt,kt) * (d2Tdx2 + d2Tdy2 + d2Tdz2) + dmdx * dTdx + dmdy * dTdy + dmdz * dTdz)
      end block
    endif
    Rv(:,i-1,j-1,k-1) = Rv(:,i-1,j-1,k-1) * over_Jacobian(i,j)
  end subroutine calc_visc_Laplacian4th
end module calc_visc

