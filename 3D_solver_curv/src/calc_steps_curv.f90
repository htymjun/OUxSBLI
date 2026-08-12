!> Step kernels for curvilinear grid.
!> Differences from calc_steps.f90:
!>   - E and F fluxes are area-scaled by |S_xi| / |S_eta|; dt_xi and dt_eta are scalars (dt*dz).
!>   - G flux is NOT area-scaled; dt_Szeta(nx-2,ny-2) = dt * J_2D(i,j) varies per cell.
module calc_steps_curv
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_R_curv(nx, ny, nz, i, j, k, dt_xi, dt_eta, dt_Szeta, E, F, G, R)
    integer, intent(in), value              :: nx, ny, nz, i, j, k
    real(8), intent(in), value              :: dt_xi    !< dt * dz  (uniform; E is area-scaled)
    real(8), intent(in), value              :: dt_eta   !< dt * dz  (uniform; F is area-scaled)
    real(8), intent(in), value              :: dt_Szeta !< dt * J_2D(i,j)
    real(8), intent(in), device, contiguous :: E(nx-1,ny-2,nz-2,5)
    real(8), intent(in), device, contiguous :: F(nx-2,ny-1,nz-2,5)
    real(8), intent(in), device, contiguous :: G(nx-2,ny-2,nz-1,5)
    real(8), intent(out), contiguous        :: R(5)
    R(1) = dt_xi * (-E(i,j,k,1) + E(i+1,j,k,1)) + dt_eta * (-F(i,j,k,1) + F(i,j+1,k,1)) + dt_Szeta * (-G(i,j,k,1) + G(i,j,k+1,1))
    R(2) = dt_xi * (-E(i,j,k,2) + E(i+1,j,k,2)) + dt_eta * (-F(i,j,k,2) + F(i,j+1,k,2)) + dt_Szeta * (-G(i,j,k,2) + G(i,j,k+1,2))
    R(3) = dt_xi * (-E(i,j,k,3) + E(i+1,j,k,3)) + dt_eta * (-F(i,j,k,3) + F(i,j+1,k,3)) + dt_Szeta * (-G(i,j,k,3) + G(i,j,k+1,3))
    R(4) = dt_xi * (-E(i,j,k,4) + E(i+1,j,k,4)) + dt_eta * (-F(i,j,k,4) + F(i,j+1,k,4)) + dt_Szeta * (-G(i,j,k,4) + G(i,j,k+1,4))
    R(5) = dt_xi * (-E(i,j,k,5) + E(i+1,j,k,5)) + dt_eta * (-F(i,j,k,5) + F(i,j+1,k,5)) + dt_Szeta * (-G(i,j,k,5) + G(i,j,k+1,5))
  end subroutine calc_R_curv


  attributes(global) subroutine calc_step1_curv(nx, ny, nz, coef, dt_xi, dt_eta, dt_Szeta, E, F, G, &
                                                 Q_1, Q_2, Q_3, Q_4, Q_5, Q2_1, Q2_2, Q2_3, Q2_4, Q2_5)
    integer, intent(in), value               :: nx, ny, nz
    real(8), intent(in), value               :: coef, dt_xi, dt_eta
    real(8), intent(in), device, contiguous  :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous  :: E(nx-1,ny-2,nz-2,5)
    real(8), intent(in), device, contiguous  :: F(nx-2,ny-1,nz-2,5)
    real(8), intent(in), device, contiguous  :: G(nx-2,ny-2,nz-1,5)
    real(8), intent(in), device, contiguous  :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(out), device, contiguous :: Q2_1(nx,ny,nz), Q2_2(nx,ny,nz), Q2_3(nx,ny,nz), Q2_4(nx,ny,nz), Q2_5(nx,ny,nz)
    real(8) R(5)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, coef*dt_xi, coef*dt_eta, coef*dt_Szeta(i,j), E, F, G, R)
    Q2_1(i+1,j+1,k+1) = Q_1(i+1,j+1,k+1) - R(1)
    Q2_2(i+1,j+1,k+1) = Q_2(i+1,j+1,k+1) - R(2)
    Q2_3(i+1,j+1,k+1) = Q_3(i+1,j+1,k+1) - R(3)
    Q2_4(i+1,j+1,k+1) = Q_4(i+1,j+1,k+1) - R(4)
    Q2_5(i+1,j+1,k+1) = Q_5(i+1,j+1,k+1) - R(5)
  end subroutine calc_step1_curv


  attributes(global) subroutine calc_step2_3_curv(nx, ny, nz, coef1, coef2, coef3, coef4_inv, dt_xi, dt_eta, dt_Szeta, E, F, G, &
                                                    Qin_1, Qin_2, Qin_3, Qin_4, Qin_5, Qout_1, Qout_2, Qout_3, Qout_4, Qout_5)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), value                 :: coef1, coef2, coef3, coef4_inv, dt_xi, dt_eta
    real(8), intent(in), device, contiguous    :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous    :: E(nx-1,ny-2,nz-2,5)
    real(8), intent(in), device, contiguous    :: F(nx-2,ny-1,nz-2,5)
    real(8), intent(in), device, contiguous    :: G(nx-2,ny-2,nz-1,5)
    real(8), intent(in), device, contiguous    :: Qin_1(nx,ny,nz), Qin_2(nx,ny,nz), Qin_3(nx,ny,nz), Qin_4(nx,ny,nz), Qin_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Qout_1(nx,ny,nz), Qout_2(nx,ny,nz), Qout_3(nx,ny,nz), Qout_4(nx,ny,nz), Qout_5(nx,ny,nz)
    real(8) R(5)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, coef3*dt_xi, coef3*dt_eta, coef3*dt_Szeta(i,j), E, F, G, R)
    Qout_1(i+1,j+1,k+1) = (coef1*Qin_1(i+1,j+1,k+1) + coef2*Qout_1(i+1,j+1,k+1) - R(1)) * coef4_inv
    Qout_2(i+1,j+1,k+1) = (coef1*Qin_2(i+1,j+1,k+1) + coef2*Qout_2(i+1,j+1,k+1) - R(2)) * coef4_inv
    Qout_3(i+1,j+1,k+1) = (coef1*Qin_3(i+1,j+1,k+1) + coef2*Qout_3(i+1,j+1,k+1) - R(3)) * coef4_inv
    Qout_4(i+1,j+1,k+1) = (coef1*Qin_4(i+1,j+1,k+1) + coef2*Qout_4(i+1,j+1,k+1) - R(4)) * coef4_inv
    Qout_5(i+1,j+1,k+1) = (coef1*Qin_5(i+1,j+1,k+1) + coef2*Qout_5(i+1,j+1,k+1) - R(5)) * coef4_inv
  end subroutine calc_step2_3_curv


  attributes(global) subroutine calc_step_curv(nx, ny, nz, coef1, coef2, dt_xi, dt_eta, dt_Szeta, E, F, G, &
                                                Q_1, Q_2, Q_3, Q_4, Q_5, Q2_1, Q2_2, Q2_3, Q2_4, Q2_5, Rs)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), value                 :: coef1, coef2, dt_xi, dt_eta
    real(8), intent(in), device, contiguous    :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous    :: E(nx-1,ny-2,nz-2,5)
    real(8), intent(in), device, contiguous    :: F(nx-2,ny-1,nz-2,5)
    real(8), intent(in), device, contiguous    :: G(nx-2,ny-2,nz-1,5)
    real(8), intent(in), device, contiguous    :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(out), device, contiguous   :: Q2_1(nx,ny,nz), Q2_2(nx,ny,nz), Q2_3(nx,ny,nz), Q2_4(nx,ny,nz), Q2_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Rs(nx-2,ny-2,nz-2,5)
    real(8) R(5)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, dt_xi, dt_eta, dt_Szeta(i,j), E, F, G, R)
    Q2_1(i+1,j+1,k+1) = Q_1(i+1,j+1,k+1) - coef1 * R(1)
    Rs(i,j,k,1) = Rs(i,j,k,1) + coef2 * R(1)
    Q2_2(i+1,j+1,k+1) = Q_2(i+1,j+1,k+1) - coef1 * R(2)
    Rs(i,j,k,2) = Rs(i,j,k,2) + coef2 * R(2)
    Q2_3(i+1,j+1,k+1) = Q_3(i+1,j+1,k+1) - coef1 * R(3)
    Rs(i,j,k,3) = Rs(i,j,k,3) + coef2 * R(3)
    Q2_4(i+1,j+1,k+1) = Q_4(i+1,j+1,k+1) - coef1 * R(4)
    Rs(i,j,k,4) = Rs(i,j,k,4) + coef2 * R(4)
    Q2_5(i+1,j+1,k+1) = Q_5(i+1,j+1,k+1) - coef1 * R(5)
    Rs(i,j,k,5) = Rs(i,j,k,5) + coef2 * R(5)
  end subroutine calc_step_curv


  attributes(global) subroutine calc_step4_curv(nx, ny, nz, dt_xi, dt_eta, dt_Szeta, E, F, G, Rs, &
                                                 Q_1, Q_2, Q_3, Q_4, Q_5)
    integer, intent(in), value                 :: nx, ny, nz
    real(8), intent(in), value                 :: dt_xi, dt_eta
    real(8), intent(in), device, contiguous    :: dt_Szeta(nx-2,ny-2)
    real(8), intent(in), device, contiguous    :: E(nx-1,ny-2,nz-2,5)
    real(8), intent(in), device, contiguous    :: F(nx-2,ny-1,nz-2,5)
    real(8), intent(in), device, contiguous    :: G(nx-2,ny-2,nz-1,5)
    real(8), intent(inout), device, contiguous :: Rs(nx-2,ny-2,nz-2,5)
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8) R(5)
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    call calc_R_curv(nx, ny, nz, i, j, k, dt_xi, dt_eta, dt_Szeta(i,j), E, F, G, R)
    R(1) = Rs(i,j,k,1) + R(1)
    Q_1(i+1,j+1,k+1) = Q_1(i+1,j+1,k+1) - R(1) * one_sixth
    Rs(i,j,k,1) = 0.d0
    R(2) = Rs(i,j,k,2) + R(2)
    Q_2(i+1,j+1,k+1) = Q_2(i+1,j+1,k+1) - R(2) * one_sixth
    Rs(i,j,k,2) = 0.d0
    R(3) = Rs(i,j,k,3) + R(3)
    Q_3(i+1,j+1,k+1) = Q_3(i+1,j+1,k+1) - R(3) * one_sixth
    Rs(i,j,k,3) = 0.d0
    R(4) = Rs(i,j,k,4) + R(4)
    Q_4(i+1,j+1,k+1) = Q_4(i+1,j+1,k+1) - R(4) * one_sixth
    Rs(i,j,k,4) = 0.d0
    R(5) = Rs(i,j,k,5) + R(5)
    Q_5(i+1,j+1,k+1) = Q_5(i+1,j+1,k+1) - R(5) * one_sixth
    Rs(i,j,k,5) = 0.d0
  end subroutine calc_step4_curv
end module calc_steps_curv
