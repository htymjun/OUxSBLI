module calc_steps
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth, one_third
  implicit none
contains
  subroutine calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,ny-1,nx-2,nz-2), device :: F
    real(8), intent(in), dimension(5,nz-1,ny-2,nx-2), device :: G
    real(8), intent(out), dimension(5,nx-2,ny-2,nz-2),device :: R
    real(8) dydz, dzdx, dxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dydz = dy(j) * dz(k)
          dzdx = dz(k) * dx(i)
          dxdy = dx(i) * dy(j)
          do l = 1, 5
            R(l,i,j,k) = dt * &
            &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dzdx * (-F(l,j,i,k) + F(l,j+1,i,k)) &
            & + dxdy * (-G(l,k,j,i) + G(l,k+1,j,i)))
    enddo;enddo;enddo;enddo
  end subroutine calc_R


  subroutine calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: coef
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,ny-1,nx-2,nz-2), device :: F
    real(8), intent(in), dimension(5,nz-1,ny-2,nx-2), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device      :: Q2
    real(8) R, dydz, dzdx, dxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dydz = dy(j) * dz(k)
          dzdx = dz(k) * dx(i)
          dxdy = dx(i) * dy(j)
          do l = 1, 5
            R = dt * &
            &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dzdx * (-F(l,j,i,k) + F(l,j+1,i,k)) &
            & + dxdy * (-G(l,k,j,i) + G(l,k+1,j,i)))
            Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step1
  
    
  subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), value                                  :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                :: dx
    real(8), intent(in), dimension(ny-1), device                :: dy
    real(8), intent(in), dimension(nz-1), device                :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device    :: E
    real(8), intent(in), dimension(5,ny-1,nx-2,nz-2), device    :: F
    real(8), intent(in), dimension(5,nz-1,ny-2,nx-2), device    :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device          :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device         :: Q2
    real(8), intent(inout), dimension(5,nx-2,ny-2,nz-2), device :: Rs
    real(8) R, dydz, dzdx, dxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dydz = dy(j) * dz(k)
          dzdx = dz(k) * dx(i)
          dxdy = dx(i) * dy(j)
          do l = 1, 5
            R = dt * &
            &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dzdx * (-F(l,j,i,k) + F(l,j+1,i,k)) &
            & + dxdy * (-G(l,k,j,i) + G(l,k+1,j,i)))
            Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef1 * R
            Rs(l,i,j,k) = Rs(l,i,j,k) + coef2 * R
    enddo;enddo;enddo;enddo
  end subroutine calc_step
  
    
  subroutine calc_step2(nx, ny, nz, dx, dy, dz, E, F, G, Qin, Qout)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,ny-1,nx-2,nz-2), device :: F
    real(8), intent(in), dimension(5,nz-1,ny-2,nx-2), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Qin
    real(8), intent(inout), dimension(5,nx,ny,nz), device    :: Qout
    real(8) R, dydz, dzdx, dxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dydz = dy(j) * dz(k)
          dzdx = dz(k) * dx(i)
          dxdy = dx(i) * dy(j)
          do l = 1, 5
            R = dt * &
            &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dzdx * (-F(l,j,i,k) + F(l,j+1,i,k)) &
            & + dxdy * (-G(l,k,j,i) + G(l,k+1,j,i)))
            Qout(l,i+1,j+1,k+1) = 0.25d0 * (3.d0 * Qin(l,i+1,j+1,k+1) + Qout(l,i+1,j+1,k+1) - R)
    enddo;enddo;enddo;enddo
  end subroutine calc_step2
  

  subroutine calc_step3(nx, ny, nz, dx, dy, dz, E, F, G, Qin, Qout)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,ny-1,nx-2,nz-2), device :: F
    real(8), intent(in), dimension(5,nz-1,ny-2,nx-2), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Qin
    real(8), intent(inout), dimension(5,nx,ny,nz), device    :: Qout
    real(8) R, dydz, dzdx, dxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dydz = dy(j) * dz(k)
          dzdx = dz(k) * dx(i)
          dxdy = dx(i) * dy(j)
          do l = 1, 5
            R = dt * &
            &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dzdx * (-F(l,j,i,k) + F(l,j+1,i,k)) &
            & + dxdy * (-G(l,k,j,i) + G(l,k+1,j,i)))
            Qout(l,i+1,j+1,k+1) = (2.d0 * Qin(l,i+1,j+1,k+1) + Qout(l,i+1,j+1,k+1) - 2.d0 * R) * one_third
    enddo;enddo;enddo;enddo
  end subroutine calc_step3
  
    
  subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                :: dx
    real(8), intent(in), dimension(ny-1), device                :: dy
    real(8), intent(in), dimension(nz-1), device                :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device    :: E
    real(8), intent(in), dimension(5,ny-1,nx-2,nz-2), device    :: F
    real(8), intent(in), dimension(5,nz-1,ny-2,nx-2), device    :: G
    real(8), intent(inout), dimension(5,nx-2,ny-2,nz-2), device :: Rs
    real(8), intent(inout), dimension(5,nx,ny,nz), device       :: Q
    real(8) R, dydz, dzdx, dxdy
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dydz = dy(j) * dz(k)
          dzdx = dz(k) * dx(i)
          dxdy = dx(i) * dy(j)
          do l = 1, 5
            R = dt * &
            &  (dydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            & + dzdx * (-F(l,j,i,k) + F(l,j+1,i,k)) &
            & + dxdy * (-G(l,k,j,i) + G(l,k+1,j,i)))
            Rs(l,i,j,k) = Rs(l,i,j,k) + R
            Q(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - Rs(l,i,j,k) * one_sixth
            Rs(l,i,j,k) = 0.d0
    enddo;enddo;enddo;enddo
  end subroutine calc_step4

  
  subroutine calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R1, R2, R1_new, R2_new
    real(8), intent(out)                                     :: err
    integer i, j, k, l
    err = 0.d0
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          do l = 1, 5
            err = err + sqrt((R1(l,i,j,k) - R1_new(l,i,j,k)**2)) + sqrt((R2(l,i,j,k) - R2_new(l,i,j,k))**2)
    enddo;enddo;enddo;enddo
    err = err / (dble(nx - 2) * dble(ny - 2) * dble(nz - 2) * 5.d0)
  end subroutine calc_error


  subroutine calc_Gauss_step(nx, ny, nz, a1, a2, R1, R2, Q, Q2)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: a1, a2
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R1
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R2
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device      :: Q2
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          do l = 1, 5
            Q2(l,i,j,k) = Q(l,i,j,k) - (a1 * R1(l,i-1,j-1,k-1) + a2 * R2(l,i-1,j-1,k-1))
    enddo;enddo;enddo;enddo
  end subroutine calc_Gauss_step


  subroutine calc_Gauss_step_Q(nx, ny, nz, a1, a2, R1, R2, Q)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: a1, a2
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R1
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-2), device :: R2
    real(8), intent(inout), dimension(5,nx,ny,nz), device    :: Q
    integer i, j, k, l
    !$cuf kernel do(3) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          do l = 1, 5
            Q(l,i,j,k) = Q(l,i,j,k) - (a1 * R1(l,i-1,j-1,k-1) + a2 * R2(l,i-1,j-1,k-1))
    enddo;enddo;enddo;enddo
  end subroutine calc_Gauss_step_Q
end module calc_steps

