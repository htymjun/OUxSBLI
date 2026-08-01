module calc_steps
  use cudafor
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  use libm
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_R(nx, ny, i, j, dtdx, dtdy, E, F, R)
    integer, intent(in), value              :: nx             !< number of grid points in x direction
    integer, intent(in), value              :: ny             !< number of grid points in y direction
    integer, intent(in), value              :: i, j           !< index
    real(8), intent(in), value              :: dtdx           !< dt * Sx, dz = 1
    real(8), intent(in), value              :: dtdy           !< dt * Sy, dz = 1
    real(8), intent(in), device, contiguous :: E(4,nx-1,ny-2) !< Flux in x direction
    real(8), intent(in), device, contiguous :: F(4,nx-2,ny-1) !< Flux in y direction
    real(8), intent(out), contiguous        :: R(4)
    real(8) v0, v1
    ! x direction
    !R(1) = dtdy * (-E(1,i,j) + E(1,i+1,j))
    v0 = E(1,i,j); v1 = E(1,i+1,j); R(1) = dtdy * (v1 - v0)
    v0 = E(2,i,j); v1 = E(2,i+1,j); R(2) = dtdy * (v1 - v0)
    v0 = E(3,i,j); v1 = E(3,i+1,j); R(3) = dtdy * (v1 - v0)
    v0 = E(4,i,j); v1 = E(4,i+1,j); R(4) = dtdy * (v1 - v0)
    ! y direction
    !R(1) = R(1) + dtdx * (-F(1,i,j) + F(1,i,j+1))
    v0 = F(1,i,j); v1 = F(1,i,j+1); R(1) = fma(dtdx, v1 - v0, R(1))
    v0 = F(2,i,j); v1 = F(2,i,j+1); R(2) = fma(dtdx, v1 - v0, R(2))
    v0 = F(3,i,j); v1 = F(3,i,j+1); R(3) = fma(dtdx, v1 - v0, R(3))
    v0 = F(4,i,j); v1 = F(4,i,j+1); R(4) = fma(dtdx, v1 - v0, R(4))
  end subroutine calc_R


  !> CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step1(nx, ny, coef, dtdx, dtdy, E, F, Q, Q2)
    integer, intent(in), value               :: nx              !< number of grid points in x direction
    integer, intent(in), value               :: ny              !< number of grid points in y direction
    real(8), intent(in), value               :: coef            !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous  :: dtdx(nx-2)      !< dt * Sx, dz = 1
    real(8), intent(in), device, contiguous  :: dtdy(ny-2)      !< dt * Sy, dz = 1
    real(8), intent(in), device, contiguous  :: E(4,nx-1,ny-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous  :: F(4,nx-2,ny-1)  !< Flux in y direction
    real(8), intent(in), device, contiguous  :: Q(nx,4,ny)      !< present Q(rho, rhou, rhov, E) / Jacobian
    real(8), intent(out), device, contiguous :: Q2(nx,4,ny)     !< next    Q(rho, rhou, rhov, E) / Jacobian
    real(8) R(4), coef_dtdx, coef_dtdy
    integer i, j, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    if (nx-2 < i .or. ny-2 < j) return
    ! ========== Conservative Update via TVD RK3: Stage 1 ==========
    ! Q^(1) = Q^n - (coef) * dt/vol * (Flux_divergence)
    coef_dtdx = coef * dtdx(i)
    coef_dtdy = coef * dtdy(j)
    call calc_R(nx, ny, i, j, coef_dtdx, coef_dtdy, E, F, R)
    do l = 1, 4  ! Loop over all conserved variables (rho, rhou, rhov, E)
      ! Q2 is write-once here and only consumed by the next kernel launch: __stcs
      call __stcs(Q2(i+1,l,j+1), Q(i+1,l,j+1) - R(l))
    enddo
  end subroutine calc_step1


  !> CUDA Fortran kernel for 1st~3rd step of 4-4 Runge-Kutta
  attributes(global) subroutine calc_step(nx, ny, coef1, coef2, dtdx, dtdy, E, F, Q, Q2, Rs)
    integer, intent(in), value                 :: nx              !< number of grid points in x direction
    integer, intent(in), value                 :: ny              !< number of grid points in y direction
    real(8), intent(in), value                 :: coef1           !< coefficient for Runge-Kutta
    real(8), intent(in), value                 :: coef2           !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous    :: dtdx(nx-2)      !< dt * Sx, dz = 1
    real(8), intent(in), device, contiguous    :: dtdy(ny-2)      !< dt * Sy, dz = 1
    real(8), intent(in), device, contiguous    :: E(4,nx-1,ny-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(4,nx-2,ny-1)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: Q(nx,4,ny)      !< present Q(rho, rhou, rhov, E) / Jacobian
    real(8), intent(out), device, contiguous   :: Q2(nx,4,ny)     !< next    Q(rho, rhou, rhov, E) / Jacobian
    real(8), intent(inout), device, contiguous :: Rs(nx-2,4,ny-2) !< accumulation for 4-4 Runge-Kutta
    real(8) R(4), dtdx_tmp, dtdy_tmp
    integer i, j, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    if (nx-2 < i .or. ny-2 < j) return
    ! ========== Conservative Update via 4-4 RK: Stage 1-3 ==========
    ! For stage 1-3: Q^(s) = Q^(s-1) - coef1 * dt/vol * Flux_div + accumulate in Rs
    ! coef2 applies weighting to residual for final 4th stage assembly
    dtdx_tmp = dtdx(i)
    dtdy_tmp = dtdy(j)
    call calc_R(nx, ny, i, j, dtdx_tmp, dtdy_tmp, E, F, R)
    do l = 1, 4
      !Q2(i+1,l,j+1) = Q(i+1,l,j+1) - coef1 * R(l) ! Intermediate Q for next stage
      ! Q2 is write-once here, not re-read until the next kernel launch: __stcs
      call __stcs(Q2(i+1,l,j+1), fma(-coef1, R(l), Q(i+1,l,j+1)))
      Rs(i,l,j) = fma(coef2, R(l), Rs(i,l,j))     ! Accumulate weighted residual
    enddo
  end subroutine calc_step
 

  !> CUDA Fortran kernel for 2nd & 3rd step of 3-3 TVD Runge-Kutta
  !> TVD RK3 Stage 2 & 3: Q^(n+1) = (α*Q^n + β*Q^(*) - γ*R)/(α+β)
  attributes(global) subroutine calc_step2_3(nx, ny, coef1, coef2, coef3, coef4_inv, dtdx, dtdy, E, F, Qin, Qout)
    integer, intent(in), value                 :: nx              !< number of grid points in x direction
    integer, intent(in), value                 :: ny              !< number of grid points in y direction
    real(8), intent(in), value                 :: coef1           !< α coefficient (weight of original Q^n)
    real(8), intent(in), value                 :: coef2           !< β coefficient (weight of Q^(*))
    real(8), intent(in), value                 :: coef3           !< γ coefficient (weight of flux residual)
    real(8), intent(in), value                 :: coef4_inv       !< 1/(α+β) normalization factor
    real(8), intent(in), device, contiguous    :: dtdx(nx-2)      !< dt * Sx, dz = 1
    real(8), intent(in), device, contiguous    :: dtdy(ny-2)      !< dt * Sy, dz = 1
    real(8), intent(in), device, contiguous    :: E(4,nx-1,ny-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(4,nx-2,ny-1)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: Qin(nx,4,ny)    !< Q^n (original from previous step)
    real(8), intent(inout), device, contiguous :: Qout(nx,4,ny)   !< Q^(*) on input, Q^(n+1) on output
    real(8) R(4), coef3_dtdx, coef3_dtdy
    integer i, j, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    if (nx-2 < i .or. ny-2 < j) return
    ! ========== TVD RK3 Stage 2 & 3 Update ==========
    ! Q^(n+1) = (α·Q^n + β·Q^(*) - γ·dt/vol·∇·F) / (α+β)
    ! Stage 2: α=3/4, β=1/4 (from Q^n and Q^(1)), coef4 = 1.d0 (compiler eliminates this division)
    ! Stage 3: α=1/3, β=2/3 (from Q^n and Q^(2)), coef4 = 3.d0 (requires division or inversion)
    coef3_dtdx = coef3 * dtdx(i)
    coef3_dtdy = coef3 * dtdy(j)
    call calc_R(nx, ny, i, j, coef3_dtdx, coef3_dtdy, E, F, R)
    do l = 1, 4  ! All conserved variables
      ! Convex combination: weighted average of Qin and Qout minus scaled residual
      !Qout(i+1,l,j+1) = (coef1 * Qin(i+1,l,j+1) + coef2 * Qout(i+1,l,j+1) - R(l)) * coef4_inv
      ! Qout is write-once here, not re-read until the next kernel launch: __stcs
      block
        real(8) qin_val, qout_val, val
        qin_val  = Qin(i+1,l,j+1)
        qout_val = Qout(i+1,l,j+1)
        val      = fma(coef1, qin_val, fma(coef2, qout_val, -R(l))) * coef4_inv
        call __stcs(Qout(i+1,l,j+1), val)
      end block
    enddo
  end subroutine calc_step2_3
  
 
  !> CUDA Fortran kernel for 4th step of 4-4 Runge-Kutta
  !> Final RK4 Stage: Q^n+1 = Q^n - (1/6)·∑(R_ᵢ) where R_ᵢ indexed over 4 stages
  attributes(global) subroutine calc_step4(nx, ny, dtdx, dtdy, E, F, Rs, Q)
    integer, intent(in), value                 :: nx              !< number of grid points in x direction
    integer, intent(in), value                 :: ny              !< number of grid points in y direction
    real(8), intent(in), device, contiguous    :: dtdx(nx-2)      !< dt * Sxy
    real(8), intent(in), device, contiguous    :: dtdy(ny-2)      !< dt * Syz
    real(8), intent(in), device, contiguous    :: E(4,nx-1,ny-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(4,nx-2,ny-1)  !< Flux in y direction
    real(8), intent(inout), device, contiguous :: Rs(nx-2,4,ny-2) !< accumulated residuals from stages 1-3
    real(8), intent(inout), device, contiguous :: Q(nx,4,ny)      !< Q^n on input, Q^n+1 on output
    real(8) R(4), dtdx_tmp, dtdy_tmp
    integer i, j, l
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    if (nx-2 < i .or. ny-2 < j) return
    ! ========== 4-4 RK Final Assembly ==========
    ! Compute 4th stage residual and accumulate with previous stages
    ! Final update: Q^(n+1) = Q^n - (one_sixth) * (R1 + 2*R2 + 2*R3 + R4)
    ! one_sixth ≈ 1/6 is the standard RK4 weight
    dtdx_tmp = dtdx(i)
    dtdy_tmp = dtdy(j)
    call calc_R(nx, ny, i, j, dtdx_tmp, dtdy_tmp, E, F, R)
    do l = 1, 4 ! All conserved variables
      ! Accumulate 4th stage residual (not multiplied by coefficient yet)
      R(l) = Rs(i,l,j) + R(l)
      ! Apply full RK4 update with (1/6) weighting to final solution
      !Q(i+1,l,j+1) = Q(i+1,l,j+1) - R(l) * one_sixth
      ! Q is write-once here, not re-read until the next kernel launch: __stcs
      call __stcs(Q(i+1,l,j+1), fma(-one_sixth, R(l), Q(i+1,l,j+1)))
      ! Clear residual accumulator for next time step
      Rs(i,l,j) = 0.d0
    enddo
  end subroutine calc_step4
end module calc_steps
