module calc_steps
  use cudafor
  use mod_globals, only : nx, ny, nz, dt, threads
  use mod_constant, only : one_sixth, one_third
  implicit none
  interface calc_R
    module procedure calc_R_Euler, calc_R_NS
  end interface
  
  interface calc_step1
    module procedure calc_step1_Euler, calc_step1_NS
  end interface
  
  interface calc_step
    module procedure calc_step_Euler, calc_step_NS
  end interface
  
  interface calc_step2
    module procedure calc_step2_Euler, calc_step2_NS
  end interface

  interface calc_step3
    module procedure calc_step3_Euler, calc_step3_NS
  end interface
  
  interface calc_step4
    module procedure calc_step4_Euler, calc_step4_NS
  end interface
contains
  !$dir inline
  attributes(device) function calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv) result(R)
    integer, intent(in), value  :: i, j, k, nx, ny, nz
    real(8), intent(in), device :: dx(nx-1)
    real(8), intent(in), device :: dy(ny-1)
    real(8), intent(in), device :: dz(nz-1)
    real(8), intent(in), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device :: Rv(1,1,1,1)
    real(8) R(5), dxdy, dydz, dzdx
    integer l
    dydz = 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dzdx = 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
    dxdy = 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
    do l = 1, 5
      R(l) = dt * &
      &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
      & + dzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dxdy * (-G(l,i,j,k) + G(l,i,j,k+1)))
    enddo
  end function calc_R_base_Euler


  !$dir inline
  attributes(device) function calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv) result(R)
    integer, intent(in), value  :: i, j, k, nx, ny, nz
    real(8), intent(in), device :: dx(nx-1)
    real(8), intent(in), device :: dy(ny-1)
    real(8), intent(in), device :: dz(nz-1)
    real(8), intent(in), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device :: Rv(4,nx-2,ny-2,nz-2)
    real(8) R(5), dxdy, dydz, dzdx
    integer l
    dydz = 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dzdx = 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
    dxdy = 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
    R(1) = dt * &
    &  (dydz * (-E(1,i,j,k) + E(1,i+1,j,k)) &
    & + dzdx * (-F(1,i,j,k) + F(1,i,j+1,k)) &
    & + dxdy * (-G(1,i,j,k) + G(1,i,j,k+1)))
    do l = 2, 5
      R(l) = dt * &
      &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
      & + dzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dxdy * (-G(l,i,j,k) + G(l,i,j,k+1)) - Rv(l-1,i,j,k))
    enddo
  end function calc_R_base_NS


  attributes(global) subroutine calc_R_Euler(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, R)
    integer(2), intent(in), value :: id_visc
    integer, intent(in), value    :: nx, ny, nz
    real(8), intent(in), device   :: dx(nx-1)
    real(8), intent(in), device   :: dy(ny-1)
    real(8), intent(in), device   :: dz(nz-1)
    real(8), intent(in), device   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device   :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device   :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device   :: Rv(1,1,1,1)
    real(8), intent(out), device  :: R(5,nx-2,ny-2,nz-2)
    real(8) dydz, dzdx, dxdy
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    R(:,i,j,k) = calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
  end subroutine calc_R_Euler

  
  attributes(global) subroutine calc_R_NS(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, R)
    integer(4), intent(in), value :: id_visc
    integer, intent(in), value    :: nx, ny, nz
    real(8), intent(in), device   :: dx(nx-1)
    real(8), intent(in), device   :: dy(ny-1)
    real(8), intent(in), device   :: dz(nz-1)
    real(8), intent(in), device   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device   :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device   :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device   :: Rv(4,nx-2,ny-2,nz-2)
    real(8), intent(out), device  :: R(5,nx-2,ny-2,nz-2)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    R(:,i,j,k) = calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
  end subroutine calc_R_NS


  attributes(global) subroutine calc_step1_Euler(id_visc, nx, ny, nz, coef, dx, dy, dz, E, F, G, Rv, Q, Q2)
    integer(2), intent(in), value  :: id_visc
    integer, intent(in), value    :: nx, ny, nz
    real(8), intent(in), value    :: coef
    real(8), intent(in), device   :: dx(nx-1)
    real(8), intent(in), device   :: dy(ny-1)
    real(8), intent(in), device   :: dz(nz-1)
    real(8), intent(in), device   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device   :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device   :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device   :: Rv(1,1,1,1)
    real(8), intent(in), device   :: Q(5,nx,ny,nz)
    real(8), intent(out), device  :: Q2(5,nx,ny,nz)
    real(8) R
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Q2(:,i+1,j+1,k+1) = Q(:,i+1,j+1,k+1) - coef * calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
  end subroutine calc_step1_Euler


  attributes(global) subroutine calc_step1_NS(id_visc, nx, ny, nz, coef, dx, dy, dz, E, F, G, Rv, Q, Q2)
    integer(4), intent(in), value  :: id_visc
    integer, intent(in), value    :: nx, ny, nz
    real(8), intent(in), value    :: coef
    real(8), intent(in), device   :: dx(nx-1)
    real(8), intent(in), device   :: dy(ny-1)
    real(8), intent(in), device   :: dz(nz-1)
    real(8), intent(in), device   :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device   :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device   :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device   :: Rv(4,nx-2,ny-2,nz-2)
    real(8), intent(in), device   :: Q(5,nx,ny,nz)
    real(8), intent(out), device  :: Q2(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Q2(:,i+1,j+1,k+1) = Q(:,i+1,j+1,k+1) - coef * calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
  end subroutine calc_step1_NS

 
  attributes(global) subroutine calc_step_Euler(id_visc, nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Rv, Q, Q2, Rs)
    integer(2), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), value     :: coef1, coef2
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(1,1,1,1)
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(out), device   :: Q2(5,nx,ny,nz)
    real(8), intent(inout), device :: Rs(5,nx-2,ny-2,nz-2)
    real(8) R(5)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    R(:) = calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
    Q2(:,i+1,j+1,k+1) = Q(:,i+1,j+1,k+1) - coef1 * R
    Rs(:,i,j,k) = Rs(:,i,j,k) + coef2 * R
  end subroutine calc_step_Euler


  attributes(global) subroutine calc_step_NS(id_visc, nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Rv, Q, Q2, Rs)
    integer(4), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), value     :: coef1, coef2
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(4,nx-2,ny-2,nz-2)
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(out), device   :: Q2(5,nx,ny,nz)
    real(8), intent(inout), device :: Rs(5,nx-2,ny-2,nz-2)
    real(8) R(5)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    R(:) = calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
    Q2(:,i+1,j+1,k+1) = Q(:,i+1,j+1,k+1) - coef1 * R
    Rs(:,i,j,k) = Rs(:,i,j,k) + coef2 * R
  end subroutine calc_step_NS
 
 
  attributes(global) subroutine calc_step2_Euler(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, Qin, Qout)
    integer(2), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(1,1,1,1)
    real(8), intent(in), device    :: Qin(5,nx,ny,nz)
    real(8), intent(inout), device :: Qout(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Qout(:,i+1,j+1,k+1) = 0.25d0 * (3.d0 * Qin(:,i+1,j+1,k+1) + Qout(:,i+1,j+1,k+1) &
                          - calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv))
  end subroutine calc_step2_Euler


  attributes(global) subroutine calc_step2_NS(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, Qin, Qout)
    integer(4), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(4,nx-2,ny-2,nz-2)
    real(8), intent(in), device    :: Qin(5,nx,ny,nz)
    real(8), intent(inout), device :: Qout(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Qout(:,i+1,j+1,k+1) = 0.25d0 * (3.d0 * Qin(:,i+1,j+1,k+1) + Qout(:,i+1,j+1,k+1) &
                          - calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv))
  end subroutine calc_step2_NS
 

  attributes(global) subroutine calc_step3_Euler(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, Qin, Qout)
    integer(2), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(1,1,1,1)
    real(8), intent(in), device    :: Qin(5,nx,ny,nz)
    real(8), intent(inout), device :: Qout(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Qout(:,i+1,j+1,k+1) = (2.d0 * Qin(:,i+1,j+1,k+1) + Qout(:,i+1,j+1,k+1) &
                          - calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)) * one_third
  end subroutine calc_step3_Euler


  attributes(global) subroutine calc_step3_NS(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, Qin, Qout)
    integer(4), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(4,nx-2,ny-2,nz-2)
    real(8), intent(in), device    :: Qin(5,nx,ny,nz)
    real(8), intent(inout), device :: Qout(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Qout(:,i+1,j+1,k+1) = (2.d0 * Qin(:,i+1,j+1,k+1) + Qout(:,i+1,j+1,k+1) &
                          - calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)) * one_third
  end subroutine calc_step3_NS
 
 
  attributes(global) subroutine calc_step4_Euler(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, Rs, Q)
    integer(2), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(1,1,1,1)
    real(8), intent(inout), device :: Rs(5,nx-2,ny-2,nz-2)
    real(8), intent(inout), device :: Q(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Rs(:,i,j,k) = Rs(:,i,j,k) + calc_R_base_Euler(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
    Q(:,i+1,j+1,k+1) = Q(:,i+1,j+1,k+1) - Rs(:,i,j,k) * one_sixth
    Rs(:,i,j,k) = 0.d0
  end subroutine calc_step4_Euler


  attributes(global) subroutine calc_step4_NS(id_visc, nx, ny, nz, dx, dy, dz, E, F, G, Rv, Rs, Q)
    integer(4), intent(in), value  :: id_visc
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1)
    real(8), intent(in), device    :: dy(ny-1)
    real(8), intent(in), device    :: dz(nz-1)
    real(8), intent(in), device    :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(in), device    :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(in), device    :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(in), device    :: Rv(4,nx-2,ny-2,nz-2)
    real(8), intent(inout), device :: Rs(5,nx-2,ny-2,nz-2)
    real(8), intent(inout), device :: Q(5,nx,ny,nz)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    Rs(:,i,j,k) = Rs(:,i,j,k) + calc_R_base_NS(i, j, k, nx, ny, nz, dx, dy, dz, E, F, G, Rv)
    Q(:,i+1,j+1,k+1) = Q(:,i+1,j+1,k+1) - Rs(:,i,j,k) * one_sixth
    Rs(:,i,j,k) = 0.d0
  end subroutine calc_step4_NS

 
  subroutine calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                     :: err
    integer i, j, k, l
    err = 0.d0
    !$cuf kernel do(3)<<<*,(32,2,2)>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          do l = 1, 5
            err = err + sqrt((R1(l,i,j,k) - R1_new(l,i,j,k)**2)) + sqrt((R2(l,i,j,k) - R2_new(l,i,j,k))**2)
    enddo;enddo;enddo;enddo
    err = err / (dble(nx - 2) * dble(ny - 2) * dble(nz - 2) * 5.d0)
  end subroutine calc_error


  attributes(global) subroutine calc_Gauss_step(nx, ny, nz, a1, a2, R1, R2, Q, Q2)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: a1, a2
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R1
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R2
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device      :: Q2
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    do l = 1, 5
      Q2(l,i,j,k) = Q(l,i,j,k) - (a1 * R1(l,i-1,j-1,k-1) + a2 * R2(l,i-1,j-1,k-1))
    enddo
  end subroutine calc_Gauss_step


  attributes(global) subroutine calc_Gauss_step_Q(nx, ny, nz, a1, a2, R1, R2, Q)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: a1, a2
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R1
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R2
    real(8), intent(inout), dimension(5,nx,ny,nz), device    :: Q
    integer i, j, k, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    do l = 1, 5
      Q(l,i,j,k) = Q(l,i,j,k) - (a1 * R1(l,i-1,j-1,k-1) + a2 * R2(l,i-1,j-1,k-1))
    enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

