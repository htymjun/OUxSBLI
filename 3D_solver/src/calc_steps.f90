module calc_steps
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  use libm
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_R(nx, ny, nz, i, j, k, dtdxdy, dtdydz, dtdzdx, E, F, G, R)
    integer, intent(in), value              :: nx                  !< number of grid points in x direction
    integer, intent(in), value              :: ny                  !< number of grid points in y direction
    integer, intent(in), value              :: nz                  !< number of grid points in z direction
    integer, intent(in), value              :: i, j, k             !< index
    real(8), intent(in), value              :: dtdxdy              !< dt * Sxy
    real(8), intent(in), value              :: dtdydz              !< dt * Syz
    real(8), intent(in), value              :: dtdzdx              !< dt * Szx
    real(8), intent(in), device, contiguous :: E(nx-1,ny-2,nz-2,5) !< Flux in x direction
    real(8), intent(in), device, contiguous :: F(nx-2,ny-1,nz-2,5) !< Flux in y direction
    real(8), intent(in), device, contiguous :: G(nx-2,ny-2,nz-1,5) !< Flux in z direction
    real(8), intent(out), contiguous        :: R(5)
    real(8) v0, v1
    ! x direction
    !R(1) = dtdydz * (-E(i,j,k,1) + E(i+1,j,k,1))
    v0 = E(i,j,k,1); v1 = E(i+1,j,k,1); R(1) = dtdydz * (v1 - v0)
    v0 = E(i,j,k,2); v1 = E(i+1,j,k,2); R(2) = dtdydz * (v1 - v0)
    v0 = E(i,j,k,3); v1 = E(i+1,j,k,3); R(3) = dtdydz * (v1 - v0)
    v0 = E(i,j,k,4); v1 = E(i+1,j,k,4); R(4) = dtdydz * (v1 - v0)
    v0 = E(i,j,k,5); v1 = E(i+1,j,k,5); R(5) = dtdydz * (v1 - v0)
    ! y direction
    !R(1) = R(1) + dtdzdx * (-F(i,j,k,1) + F(i,j+1,k,1))
    v0 = F(i,j,k,1); v1 = F(i,j+1,k,1); R(1) = fma(dtdzdx, v1 - v0, R(1))
    v0 = F(i,j,k,2); v1 = F(i,j+1,k,2); R(2) = fma(dtdzdx, v1 - v0, R(2))
    v0 = F(i,j,k,3); v1 = F(i,j+1,k,3); R(3) = fma(dtdzdx, v1 - v0, R(3))
    v0 = F(i,j,k,4); v1 = F(i,j+1,k,4); R(4) = fma(dtdzdx, v1 - v0, R(4))
    v0 = F(i,j,k,5); v1 = F(i,j+1,k,5); R(5) = fma(dtdzdx, v1 - v0, R(5))
    ! z direction
    !R(1) = R(1) + dtdxdy * (-G(i,j,k,1) + G(i,j,k+1,1))
    v0 = G(i,j,k,1); v1 = G(i,j,k+1,1); R(1) = fma(dtdxdy, v1 - v0, R(1))
    v0 = G(i,j,k,2); v1 = G(i,j,k+1,2); R(2) = fma(dtdxdy, v1 - v0, R(2))
    v0 = G(i,j,k,3); v1 = G(i,j,k+1,3); R(3) = fma(dtdxdy, v1 - v0, R(3))
    v0 = G(i,j,k,4); v1 = G(i,j,k+1,4); R(4) = fma(dtdxdy, v1 - v0, R(4))
    v0 = G(i,j,k,5); v1 = G(i,j,k+1,5); R(5) = fma(dtdxdy, v1 - v0, R(5))
  end subroutine calc_R


  !> CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step1(nx, ny, nz, coef, dtdxdy, dtdydz, dtdzdx, E, F, G, &
                                            Q_1, Q_2, Q_3, Q_4, Q_5, Q2_1, Q2_2, Q2_3, Q2_4, Q2_5)
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), value               :: coef                !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous  :: dtdxdy(nx-2,ny-2)   !< dt * Sxy
    real(8), intent(in), device, contiguous  :: dtdydz(ny-2,nz-2)   !< dt * Syz
    real(8), intent(in), device, contiguous  :: dtdzdx(nx-2,nz-2)   !< dt * Szx
    real(8), intent(in), device, contiguous  :: E(nx-1,ny-2,nz-2,5) !< Flux in x direction
    real(8), intent(in), device, contiguous  :: F(nx-2,ny-1,nz-2,5) !< Flux in y direction
    real(8), intent(in), device, contiguous  :: G(nx-2,ny-2,nz-1,5) !< Flux in z direction
    real(8), intent(in), device, contiguous  :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)  !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(in), device, contiguous  :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(out), device, contiguous :: Q2_1(nx,ny,nz), Q2_2(nx,ny,nz), Q2_3(nx,ny,nz) !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous :: Q2_4(nx,ny,nz), Q2_5(nx,ny,nz)
    real(8) R(5), coef_dtdxdy, coef_dtdydz, coef_dtdzdx
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== Conservative Update via TVD RK3: Stage 1 ==========
    ! Q^(1) = Q^n - (coef) * dt/vol * (Flux_divergence)
    coef_dtdxdy = coef * dtdxdy(i,j)
    coef_dtdydz = coef * dtdydz(j,k)
    coef_dtdzdx = coef * dtdzdx(i,k)
    call calc_R(nx, ny, nz, i, j, k, coef_dtdxdy, coef_dtdydz, coef_dtdzdx, E, F, G, R)
    ! Q2 is write-once here and only consumed by the next kernel launch: __stcs
    call __stcs(Q2_1(i+1,j+1,k+1), Q_1(i+1,j+1,k+1) - R(1))
    call __stcs(Q2_2(i+1,j+1,k+1), Q_2(i+1,j+1,k+1) - R(2))
    call __stcs(Q2_3(i+1,j+1,k+1), Q_3(i+1,j+1,k+1) - R(3))
    call __stcs(Q2_4(i+1,j+1,k+1), Q_4(i+1,j+1,k+1) - R(4))
    call __stcs(Q2_5(i+1,j+1,k+1), Q_5(i+1,j+1,k+1) - R(5))
  end subroutine calc_step1


  !> CUDA Fortran kernel for 1st~3rd step of 4-4 Runge-Kutta
  attributes(global) subroutine calc_step(nx, ny, nz, coef1, coef2, dtdxdy, dtdydz, dtdzdx, E, F, G, &
                                           Q_1, Q_2, Q_3, Q_4, Q_5, Q2_1, Q2_2, Q2_3, Q2_4, Q2_5, Rs)
    integer, intent(in), value                 :: nx                   !< number of grid points in x direction
    integer, intent(in), value                 :: ny                   !< number of grid points in y direction
    integer, intent(in), value                 :: nz                   !< number of grid points in z direction
    real(8), intent(in), value                 :: coef1                !< coefficient for Runge-Kutta
    real(8), intent(in), value                 :: coef2                !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous    :: dtdxdy(nx-2,ny-2)    !< dt * Sxy
    real(8), intent(in), device, contiguous    :: dtdydz(ny-2,nz-2)    !< dt * Syz
    real(8), intent(in), device, contiguous    :: dtdzdx(nx-2,nz-2)    !< dt * Szx
    real(8), intent(in), device, contiguous    :: E(nx-1,ny-2,nz-2,5)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(nx-2,ny-1,nz-2,5)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(nx-2,ny-2,nz-1,5)  !< Flux in z direction
    real(8), intent(in), device, contiguous    :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)   !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(out), device, contiguous   :: Q2_1(nx,ny,nz), Q2_2(nx,ny,nz), Q2_3(nx,ny,nz) !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous   :: Q2_4(nx,ny,nz), Q2_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Rs(nx-2,ny-2,nz-2,5) !< accumulation for 4-4 Runge-Kutta
    real(8) R(5), dtdxdy_tmp, dtdydz_tmp, dtdzdx_tmp
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== Conservative Update via 4-4 RK: Stage 1-3 ==========
    ! For stage 1-3: Q^(s) = Q^(s-1) - coef1 * dt/vol * Flux_div + accumulate in Rs
    ! coef2 applies weighting to residual for final 4th stage assembly
    dtdxdy_tmp = dtdxdy(i,j)
    dtdydz_tmp = dtdydz(j,k)
    dtdzdx_tmp = dtdzdx(i,k)
    call calc_R(nx, ny, nz, i, j, k, dtdxdy_tmp, dtdydz_tmp, dtdzdx_tmp, E, F, G, R)
    ! Q2 is write-once here, not re-read until the next kernel launch: __stcs
    call __stcs(Q2_1(i+1,j+1,k+1), fma(-coef1, R(1), Q_1(i+1,j+1,k+1)))
    Rs(i,j,k,1) = fma(coef2, R(1), Rs(i,j,k,1))         ! Accumulate weighted residual
    call __stcs(Q2_2(i+1,j+1,k+1), fma(-coef1, R(2), Q_2(i+1,j+1,k+1)))
    Rs(i,j,k,2) = fma(coef2, R(2), Rs(i,j,k,2))
    call __stcs(Q2_3(i+1,j+1,k+1), fma(-coef1, R(3), Q_3(i+1,j+1,k+1)))
    Rs(i,j,k,3) = fma(coef2, R(3), Rs(i,j,k,3))
    call __stcs(Q2_4(i+1,j+1,k+1), fma(-coef1, R(4), Q_4(i+1,j+1,k+1)))
    Rs(i,j,k,4) = fma(coef2, R(4), Rs(i,j,k,4))
    call __stcs(Q2_5(i+1,j+1,k+1), fma(-coef1, R(5), Q_5(i+1,j+1,k+1)))
    Rs(i,j,k,5) = fma(coef2, R(5), Rs(i,j,k,5))
  end subroutine calc_step


  !> CUDA Fortran kernel for 2nd & 3rd step of 3-3 TVD Runge-Kutta
  !> TVD RK3 Stage 2 & 3: Q^(n+1) = (α*Q^n + β*Q^(*) - γ*R)/(α+β)
  attributes(global) subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4_inv, dtdxdy, dtdydz, dtdzdx, E, F, G, &
                                              Qin_1, Qin_2, Qin_3, Qin_4, Qin_5, Qout_1, Qout_2, Qout_3, Qout_4, Qout_5)
    integer, intent(in), value                 :: nx                  !< number of grid points in x direction
    integer, intent(in), value                 :: ny                  !< number of grid points in y direction
    integer, intent(in), value                 :: nz                  !< number of grid points in z direction
    real(8), intent(in), value                 :: coef1               !< α coefficient (weight of original Q^n)
    real(8), intent(in), value                 :: coef2               !< β coefficient (weight of Q^(*))
    real(8), intent(in), value                 :: coef3               !< γ coefficient (weight of flux residual)
    real(8), intent(in), value                 :: coef4_inv           !< 1/(α+β) normalization factor
    real(8), intent(in), device, contiguous    :: dtdxdy(nx-2,ny-2)   !< dt * Sxy
    real(8), intent(in), device, contiguous    :: dtdydz(ny-2,nz-2)   !< dt * Syz
    real(8), intent(in), device, contiguous    :: dtdzdx(nx-2,nz-2)   !< dt * Szx
    real(8), intent(in), device, contiguous    :: E(nx-1,ny-2,nz-2,5) !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(nx-2,ny-1,nz-2,5) !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(nx-2,ny-2,nz-1,5) !< Flux in z direction
    real(8), intent(in), device, contiguous    :: Qin_1(nx,ny,nz), Qin_2(nx,ny,nz), Qin_3(nx,ny,nz)    !< Q^n (original from previous step)
    real(8), intent(in), device, contiguous    :: Qin_4(nx,ny,nz), Qin_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Qout_1(nx,ny,nz), Qout_2(nx,ny,nz), Qout_3(nx,ny,nz)  !< Q^(*) on input, Q^(n+1) on output
    real(8), intent(inout), device, contiguous :: Qout_4(nx,ny,nz), Qout_5(nx,ny,nz)
    real(8) R(5), coef3_dtdxdy, coef3_dtdydz, coef3_dtdzdx
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== TVD RK3 Stage 2 & 3 Update ==========
    ! Q^(n+1) = (α·Q^n + β·Q^(*) - γ·dt/vol·∇·F) / (α+β)
    ! Stage 2: α=3/4, β=1/4 (from Q^n and Q^(1)), coef4 = 1.d0 (compiler eliminates this division)
    ! Stage 3: α=1/3, β=2/3 (from Q^n and Q^(2)), coef4 = 3.d0 (requires division or inversion)
    coef3_dtdxdy = coef3 * dtdxdy(i,j)
    coef3_dtdydz = coef3 * dtdydz(j,k)
    coef3_dtdzdx = coef3 * dtdzdx(i,k)
    call calc_R(nx, ny, nz, i, j, k, coef3_dtdxdy, coef3_dtdydz, coef3_dtdzdx, E, F, G, R)
    ! Convex combination: weighted average of Qin and Qout minus scaled residual
    ! Qout is write-once here, not re-read until the next kernel launch: __stcs
    block
      real(8) qin_val, qout_val, val
      qin_val  = Qin_1(i+1,j+1,k+1)
      qout_val = Qout_1(i+1,j+1,k+1)
      val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(1))) * coef4_inv
      call __stcs(Qout_1(i+1,j+1,k+1), val)
      qin_val  = Qin_2(i+1,j+1,k+1)
      qout_val = Qout_2(i+1,j+1,k+1)
      val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(2))) * coef4_inv
      call __stcs(Qout_2(i+1,j+1,k+1), val)
      qin_val  = Qin_3(i+1,j+1,k+1)
      qout_val = Qout_3(i+1,j+1,k+1)
      val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(3))) * coef4_inv
      call __stcs(Qout_3(i+1,j+1,k+1), val)
      qin_val  = Qin_4(i+1,j+1,k+1)
      qout_val = Qout_4(i+1,j+1,k+1)
      val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(4))) * coef4_inv
      call __stcs(Qout_4(i+1,j+1,k+1), val)
      qin_val  = Qin_5(i+1,j+1,k+1)
      qout_val = Qout_5(i+1,j+1,k+1)
      val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(5))) * coef4_inv
      call __stcs(Qout_5(i+1,j+1,k+1), val)
    end block
  end subroutine calc_step2_3


  !> CUDA Fortran kernel for 4th step of 4-4 Runge-Kutta
  !> Final RK4 Stage: Q^n+1 = Q^n - (1/6)·∑(R_ᵢ) where R_ᵢ indexed over 4 stages
  attributes(global) subroutine calc_step4(nx, ny, nz, dtdxdy, dtdydz, dtdzdx, E, F, G, Rs, &
                                            Q_1, Q_2, Q_3, Q_4, Q_5)
    integer, intent(in), value                 :: nx                   !< number of grid points in x direction
    integer, intent(in), value                 :: ny                   !< number of grid points in y direction
    integer, intent(in), value                 :: nz                   !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dtdxdy(nx-2,ny-2)    !< dt * Sxy
    real(8), intent(in), device, contiguous    :: dtdydz(ny-2,nz-2)    !< dt * Syz
    real(8), intent(in), device, contiguous    :: dtdzdx(nx-2,nz-2)    !< dt * Szx
    real(8), intent(in), device, contiguous    :: E(nx-1,ny-2,nz-2,5)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(nx-2,ny-1,nz-2,5)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(nx-2,ny-2,nz-1,5)  !< Flux in z direction
    real(8), intent(inout), device, contiguous :: Rs(nx-2,ny-2,nz-2,5) !< accumulated residuals from stages 1-3
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)  !< Q^n on input, Q^n+1 on output
    real(8), intent(inout), device, contiguous :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8) R(5), dtdxdy_tmp, dtdydz_tmp, dtdzdx_tmp
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    ! ========== 4-4 RK Final Assembly ==========
    ! Compute 4th stage residual and accumulate with previous stages
    ! Final update: Q^(n+1) = Q^n - (one_sixth) * (R1 + 2*R2 + 2*R3 + R4)
    ! one_sixth ≈ 1/6 is the standard RK4 weight
    dtdxdy_tmp = dtdxdy(i,j)
    dtdydz_tmp = dtdydz(j,k)
    dtdzdx_tmp = dtdzdx(i,k)
    call calc_R(nx, ny, nz, i, j, k, dtdxdy_tmp, dtdydz_tmp, dtdzdx_tmp, E, F, G, R)
    ! Accumulate 4th stage residual (not multiplied by coefficient yet)
    R(1) = Rs(i,j,k,1) + R(1)
    ! Q is write-once here, not re-read until the next kernel launch: __stcs
    call __stcs(Q_1(i+1,j+1,k+1), fma(-one_sixth, R(1), Q_1(i+1,j+1,k+1)))
    Rs(i,j,k,1) = 0.d0
    R(2) = Rs(i,j,k,2) + R(2)
    call __stcs(Q_2(i+1,j+1,k+1), fma(-one_sixth, R(2), Q_2(i+1,j+1,k+1)))
    Rs(i,j,k,2) = 0.d0
    R(3) = Rs(i,j,k,3) + R(3)
    call __stcs(Q_3(i+1,j+1,k+1), fma(-one_sixth, R(3), Q_3(i+1,j+1,k+1)))
    Rs(i,j,k,3) = 0.d0
    R(4) = Rs(i,j,k,4) + R(4)
    call __stcs(Q_4(i+1,j+1,k+1), fma(-one_sixth, R(4), Q_4(i+1,j+1,k+1)))
    Rs(i,j,k,4) = 0.d0
    R(5) = Rs(i,j,k,5) + R(5)
    call __stcs(Q_5(i+1,j+1,k+1), fma(-one_sixth, R(5), Q_5(i+1,j+1,k+1)))
    Rs(i,j,k,5) = 0.d0
  end subroutine calc_step4
end module calc_steps
