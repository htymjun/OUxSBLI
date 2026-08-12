module calc_steps
  use cudafor
  use wmma
  use mod_globals, only : dt, threads
  use mod_constant, only : one_sixth
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_R_EF(nx, ny, nz, i, j, k, dtdydz, dtdzdx, E, F, R)
    integer, intent(in), value              :: nx                  !< number of grid points in x direction
    integer, intent(in), value              :: ny                  !< number of grid points in y direction
    integer, intent(in), value              :: nz                  !< number of grid points in z direction
    integer, intent(in), value              :: i, j, k             !< index
    real(8), intent(in), value              :: dtdydz              !< dt * Syz
    real(8), intent(in), value              :: dtdzdx              !< dt * Szx
    real(8), intent(in), device, contiguous :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device, contiguous :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    real(8), intent(out), contiguous        :: R(5)
    ! x direction
    R(1) = dtdydz * (-E(1,i,j,k) + E(1,i+1,j,k)) 
    R(2) = dtdydz * (-E(2,i,j,k) + E(2,i+1,j,k)) 
    R(3) = dtdydz * (-E(3,i,j,k) + E(3,i+1,j,k)) 
    R(4) = dtdydz * (-E(4,i,j,k) + E(4,i+1,j,k)) 
    R(5) = dtdydz * (-E(5,i,j,k) + E(5,i+1,j,k)) 
    ! y direction
    R(1) = R(1) + dtdzdx * (-F(1,i,j,k) + F(1,i,j+1,k))
    R(2) = R(2) + dtdzdx * (-F(2,i,j,k) + F(2,i,j+1,k))
    R(3) = R(3) + dtdzdx * (-F(3,i,j,k) + F(3,i,j+1,k))
    R(4) = R(4) + dtdzdx * (-F(4,i,j,k) + F(4,i,j+1,k))
    R(5) = R(5) + dtdzdx * (-F(5,i,j,k) + F(5,i,j+1,k))
  end subroutine calc_R_EF


  !$dir inline
  attributes(device) subroutine calc_R_G(it, jt, kt, dtdxdy, G, R)
    integer, intent(in), value       :: it, jt, kt                           !< index
    real(8), intent(in), value       :: dtdxdy                               !< dt * Sxy
    real(8), intent(in), contiguous  :: G(5,threads%x,threads%z+1,threads%y) !< Flux in z direction (shared)
    real(8), intent(out), contiguous :: R(5)
    ! z direction
    R(1) = R(1) + dtdxdy * (-G(1,it,kt,jt) + G(1,it,kt+1,jt))
    R(2) = R(2) + dtdxdy * (-G(2,it,kt,jt) + G(2,it,kt+1,jt))
    R(3) = R(3) + dtdxdy * (-G(3,it,kt,jt) + G(3,it,kt+1,jt))
    R(4) = R(4) + dtdxdy * (-G(4,it,kt,jt) + G(4,it,kt+1,jt))
    R(5) = R(5) + dtdxdy * (-G(5,it,kt,jt) + G(5,it,kt+1,jt))
  end subroutine calc_R_G


  !$dir inline
  attributes(device) subroutine load_smem_G(it, jt, kt, i, j, nx, ny, nz, G, Gs)
    integer, intent(in), value              :: it, jt, kt                            !< local index
    integer, intent(in), value              :: i, j                                  !< global index (x,y)
    integer, intent(in), value              :: nx, ny, nz
    real(8), intent(in), device, contiguous :: G(5,nx-2,ny-2,nz-1)                   !< Flux in z direction (global)
    real(8), intent(inout), contiguous      :: Gs(5,threads%x,threads%z+1,threads%y) !< Flux in z direction (shared)
    integer k_base, kk, k
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt, threads%z+1, blockDim%z
      k = k_base + kk
      if (i <= nx-2 .and. j <= ny-2 .and. k <= nz-1) then
        call pipelineMemcpyAsync(Gs(1,it,kk,jt), G(1,i,j,k))
        call pipelineMemcpyAsync(Gs(2,it,kk,jt), G(2,i,j,k))
        call pipelineMemcpyAsync(Gs(3,it,kk,jt), G(3,i,j,k))
        call pipelineMemcpyAsync(Gs(4,it,kk,jt), G(4,i,j,k))
        call pipelineMemcpyAsync(Gs(5,it,kk,jt), G(5,i,j,k))
      endif
    enddo
    call pipelineCommit()
  end subroutine load_smem_G


  !$dir inline
  attributes(device) subroutine load_smem_G_cooperative(it, jt, kt, nx, ny, nz, G, Gs)
    integer, intent(in), value               :: it, jt, kt
    integer, intent(in), value               :: nx, ny, nz
    real(8), intent(in), device, contiguous  :: G(5,nx-2,ny-2,nz-1)
    real(8), intent(inout), contiguous       :: Gs(5,threads%x,threads%z+1,threads%y)
    integer :: tid, block_size, total_elements
    integer :: e_idx, l, i_loc, k_loc, j_loc
    integer :: i_glob, j_glob, k_glob
    integer :: k_base
    tid = (it - 1) + (kt - 1)*blockDim%x + (jt - 1)*blockDim%x*blockDim%z
    block_size = blockDim%x * blockDim%y * blockDim%z
    ! (5 × blockDim%x × (blockDim%z + 1) × blockDim%y)
    total_elements = 5 * blockDim%x * (blockDim%z + 1) * blockDim%y
    k_base = (blockIdx%z - 1) * blockDim%z
    do e_idx = tid, total_elements - 1, block_size
       ! decode 1D idx, (e_idx) to 4D idx, Gs(l, i_loc, k_loc, j_loc)
       l     = mod(e_idx, 5) + 1
       i_loc = mod(e_idx / 5, blockDim%x) + 1
       k_loc = mod(e_idx / (5 * blockDim%x), blockDim%z + 1) + 1
       j_loc = (e_idx / (5 * blockDim%x * (blockDim%z + 1))) + 1
       ! interpret to 4D idx
       i_glob = (blockIdx%x - 1)*blockDim%x + i_loc
       j_glob = (blockIdx%y - 1)*blockDim%y + j_loc
       k_glob = k_base + k_loc
       if (i_glob <= nx-2 .and. j_glob <= ny-2 .and. k_glob <= nz-1) then
          ! coalesced cooperative copy
          call pipelineMemcpyAsync(Gs(l, i_loc, k_loc, j_loc), G(l, i_glob, j_glob, k_glob))
       endif
    enddo
    call pipelineCommit()
  end subroutine load_smem_G_cooperative


  !> CUDA Fortran kernel for 1st step of 3-3 TVD Runge-Kutta
  attributes(global) subroutine calc_step1(nx, ny, nz, coef, dtdxdy, dtdydz, dtdzdx, E, F, G, Q, Q2)
    integer, intent(in), value               :: nx                  !< number of grid points in x direction
    integer, intent(in), value               :: ny                  !< number of grid points in y direction
    integer, intent(in), value               :: nz                  !< number of grid points in z direction
    real(8), intent(in), value               :: coef                !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous  :: dtdxdy(nx-2,ny-2)   !< dt * Sxy
    real(8), intent(in), device, contiguous  :: dtdydz(ny-2,nz-2)   !< dt * Syz
    real(8), intent(in), device, contiguous  :: dtdzdx(nx-2,nz-2)   !< dt * Szx
    real(8), intent(in), device, contiguous  :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device, contiguous  :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    real(8), intent(in), device, contiguous  :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    real(8), intent(in), device, contiguous  :: Q(nx,5,ny,nz)       !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous :: Q2(nx,5,ny,nz)      !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), shared :: Gs(5,threads%x,threads%z+1,threads%y)
    real(8) R(5), coef_dtdxdy, coef_dtdydz, coef_dtdzdx
    integer i, j, k, l, it, jt, kt
    logical valid
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt
    valid = (i <= nx-2 .and. j <= ny-2 .and. k <= nz-2)
    ! ========== load to shared memory ==========
    !call load_smem_G(it, jt, kt, i, j, nx, ny, nz, G, Gs)
    call load_smem_G_cooperative(it, jt, kt, nx, ny, nz, G, Gs)
    if (valid) then
      ! ========== Conservative Update via TVD RK3: Stage 1 ==========
      ! Q^(1) = Q^n - (coef) * dt/vol * (Flux_divergence)
      coef_dtdxdy = coef * dtdxdy(i,j)
      coef_dtdydz = coef * dtdydz(j,k)
      coef_dtdzdx = coef * dtdzdx(i,k)
      call calc_R_EF(nx, ny, nz, i, j, k, coef_dtdydz, coef_dtdzdx, E, F, R)
    endif
    call pipelineWaitPrior(0)
    call syncthreads() ! Gs(:,it,jt,kt+1) may have been loaded by another thread; all threads must reach this barrier
    if (.not. valid) return
    call calc_R_G(it, jt, kt, coef_dtdxdy, Gs, R)
    do l = 1, 5  ! Loop over all conserved variables (rho, rhou, rhov, rhow, E)
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - R(l)
    enddo
  end subroutine calc_step1


  !> CUDA Fortran kernel for 1st~3rd step of 4-4 Runge-Kutta
  attributes(global) subroutine calc_step(nx, ny, nz, coef1, coef2, dtdxdy, dtdydz, dtdzdx, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                 :: nx                   !< number of grid points in x direction
    integer, intent(in), value                 :: ny                   !< number of grid points in y direction
    integer, intent(in), value                 :: nz                   !< number of grid points in z direction
    real(8), intent(in), value                 :: coef1                !< coefficient for Runge-Kutta
    real(8), intent(in), value                 :: coef2                !< coefficient for Runge-Kutta
    real(8), intent(in), device, contiguous    :: dtdxdy(nx-2,ny-2)    !< dt * Sxy
    real(8), intent(in), device, contiguous    :: dtdydz(ny-2,nz-2)    !< dt * Syz
    real(8), intent(in), device, contiguous    :: dtdzdx(nx-2,nz-2)    !< dt * Szx
    real(8), intent(in), device, contiguous    :: E(5,nx-1,ny-2,nz-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(5,nx-2,ny-1,nz-2)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(5,nx-2,ny-2,nz-1)  !< Flux in z direction
    real(8), intent(in), device, contiguous    :: Q(nx,5,ny,nz)        !< present Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(out), device, contiguous   :: Q2(nx,5,ny,nz)       !< next    Q(rho, rhou, rhov, rhow, E) / Jacobian
    real(8), intent(inout), device, contiguous :: Rs(nx-2,5,ny-2,nz-2) !< accumulation for 4-4 Runge-Kutta
    real(8), shared :: Gs(5,threads%x,threads%z+1,threads%y)
    real(8) R(5), dtdxdy_tmp, dtdydz_tmp, dtdzdx_tmp
    integer i, j, k, l, it, jt, kt
    logical valid
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt
    valid = (i <= nx-2 .and. j <= ny-2 .and. k <= nz-2)
    ! ========== load to shared memory ==========
    !call load_smem_G(it, jt, kt, i, j, nx, ny, nz, G, Gs)
    call load_smem_G_cooperative(it, jt, kt, nx, ny, nz, G, Gs)
    if (valid) then
      ! ========== Conservative Update via 4-4 RK: Stage 1-3 ==========
      ! For stage 1-3: Q^(s) = Q^(s-1) - coef1 * dt/vol * Flux_div + accumulate in Rs
      ! coef2 applies weighting to residual for final 4th stage assembly
      dtdxdy_tmp = dtdxdy(i,j)
      dtdydz_tmp = dtdydz(j,k)
      dtdzdx_tmp = dtdzdx(i,k)
      call calc_R_EF(nx, ny, nz, i, j, k, dtdydz_tmp, dtdzdx_tmp, E, F, R)
    endif
    call pipelineWaitPrior(0)
    call syncthreads() ! Gs(:,it,jt,kt+1) may have been loaded by another thread; all threads must reach this barrier
    if (.not. valid) return
    call calc_R_G(it, jt, kt, dtdxdy_tmp, Gs, R)
    do l = 1, 5
      Q2(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - coef1 * R(l) ! Intermediate Q for next stage
      Rs(i,l,j,k) = Rs(i,l,j,k) + coef2 * R(l)            ! Accumulate weighted residual
    enddo
  end subroutine calc_step
 

  !> CUDA Fortran kernel for 2nd & 3rd step of 3-3 TVD Runge-Kutta
  !> TVD RK3 Stage 2 & 3: Q^(n+1) = (α*Q^n + β*Q^(*) - γ*R)/(α+β)
  attributes(global) subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4_inv, dtdxdy, dtdydz, dtdzdx, E, F, G, Qin, Qout)
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
    real(8), intent(in), device, contiguous    :: E(5,nx-1,ny-2,nz-2) !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(5,nx-2,ny-1,nz-2) !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(5,nx-2,ny-2,nz-1) !< Flux in z direction
    real(8), intent(in), device, contiguous    :: Qin(nx,5,ny,nz)     !< Q^n (original from previous step)
    real(8), intent(inout), device, contiguous :: Qout(nx,5,ny,nz)    !< Q^(*) on input, Q^(n+1) on output
    real(8), shared :: Gs(5,threads%x,threads%z+1,threads%y)
    real(8) R(5), coef3_dtdxdy, coef3_dtdydz, coef3_dtdzdx
    integer i, j, k, l, it, jt, kt
    logical valid
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt
    valid = (i <= nx-2 .and. j <= ny-2 .and. k <= nz-2)
    ! ========== load to shared memory ==========
    !call load_smem_G(it, jt, kt, i, j, nx, ny, nz, G, Gs)
    call load_smem_G_cooperative(it, jt, kt, nx, ny, nz, G, Gs)
    if (valid) then
      ! ========== TVD RK3 Stage 2 & 3 Update ==========
      ! Q^(n+1) = (α·Q^n + β·Q^(*) - γ·dt/vol·∇·F) / (α+β)
      ! Stage 2: α=3/4, β=1/4 (from Q^n and Q^(1)), coef4 = 1.d0 (compiler eliminates this division)
      ! Stage 3: α=1/3, β=2/3 (from Q^n and Q^(2)), coef4 = 3.d0 (requires division or inversion)
      coef3_dtdxdy = coef3 * dtdxdy(i,j)
      coef3_dtdydz = coef3 * dtdydz(j,k)
      coef3_dtdzdx = coef3 * dtdzdx(i,k)
      call calc_R_EF(nx, ny, nz, i, j, k, coef3_dtdydz, coef3_dtdzdx, E, F, R)
    endif
    call pipelineWaitPrior(0)
    call syncthreads() ! Gs(:,it,jt,kt+1) may have been loaded by another thread; all threads must reach this barrier
    if (.not. valid) return
    call calc_R_G(it, jt, kt, coef3_dtdxdy, Gs, R)
    do l = 1, 5  ! All conserved variables
      ! Convex combination: weighted average of Qin and Qout minus scaled residual
      Qout(i+1,l,j+1,k+1) = (coef1 * Qin(i+1,l,j+1,k+1) + coef2 * Qout(i+1,l,j+1,k+1) - R(l)) * coef4_inv
    enddo
  end subroutine calc_step2_3
 
 
  !> CUDA Fortran kernel for 4th step of 4-4 Runge-Kutta
  !> Final RK4 Stage: Q^n+1 = Q^n - (1/6)·∑(R_ᵢ) where R_ᵢ indexed over 4 stages
  attributes(global) subroutine calc_step4(nx, ny, nz, dtdxdy, dtdydz, dtdzdx, E, F, G, Rs, Q)
    integer, intent(in), value                 :: nx                   !< number of grid points in x direction
    integer, intent(in), value                 :: ny                   !< number of grid points in y direction
    integer, intent(in), value                 :: nz                   !< number of grid points in z direction
    real(8), intent(in), device, contiguous    :: dtdxdy(nx-2,ny-2)    !< dt * Sxy
    real(8), intent(in), device, contiguous    :: dtdydz(ny-2,nz-2)    !< dt * Syz
    real(8), intent(in), device, contiguous    :: dtdzdx(nx-2,nz-2)    !< dt * Szx
    real(8), intent(in), device, contiguous    :: E(5,nx-1,ny-2,nz-2)  !< Flux in x direction
    real(8), intent(in), device, contiguous    :: F(5,nx-2,ny-1,nz-2)  !< Flux in y direction
    real(8), intent(in), device, contiguous    :: G(5,nx-2,ny-2,nz-1)  !< Flux in z direction
    real(8), intent(inout), device, contiguous :: Rs(nx-2,5,ny-2,nz-2) !< accumulated residuals from stages 1-3
    real(8), intent(inout), device, contiguous :: Q(nx,5,ny,nz)        !< Q^n on input, Q^n+1 on output
    real(8), shared :: Gs(5,threads%x,threads%z+1,threads%y)
    real(8) R(5), dtdxdy_tmp, dtdydz_tmp, dtdzdx_tmp
    integer i, j, k, l, it, jt, kt
    logical valid
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt
    valid = (i <= nx-2 .and. j <= ny-2 .and. k <= nz-2)
    ! ========== load to shared memory ==========
    !call load_smem_G(it, jt, kt, i, j, nx, ny, nz, G, Gs)
    call load_smem_G_cooperative(it, jt, kt, nx, ny, nz, G, Gs)
    if (valid) then
      ! ========== 4-4 RK Final Assembly ==========
      ! Compute 4th stage residual and accumulate with previous stages
      ! Final update: Q^(n+1) = Q^n - (one_sixth) * (R1 + 2*R2 + 2*R3 + R4)
      ! one_sixth ≈ 1/6 is the standard RK4 weight
      dtdxdy_tmp = dtdxdy(i,j)
      dtdydz_tmp = dtdydz(j,k)
      dtdzdx_tmp = dtdzdx(i,k)
      call calc_R_EF(nx, ny, nz, i, j, k, dtdydz_tmp, dtdzdx_tmp, E, F, R)
    endif
    call pipelineWaitPrior(0)
    call syncthreads() ! Gs(:,it,jt,kt+1) may have been loaded by another thread; all threads must reach this barrier
    if (.not. valid) return
    call calc_R_G(it, jt, kt, dtdxdy_tmp, Gs, R)
    do l = 1, 5 ! All conserved variables
      ! Accumulate 4th stage residual (not multiplied by coefficient yet)
      R(l) = Rs(i,l,j,k) + R(l)
      ! Apply full RK4 update with (1/6) weighting to final solution
      Q(i+1,l,j+1,k+1) = Q(i+1,l,j+1,k+1) - R(l) * one_sixth
      ! Clear residual accumulator for next time step
      Rs(i,l,j,k) = 0.d0
    enddo
  end subroutine calc_step4

end module calc_steps
