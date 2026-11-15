#ifndef REAL_T
#define REAL_T real(8)
#endif
module calc_physical_quantities
  use cudafor
#ifdef PREC_DF
  use mod_prec
#endif
  use mod_globals, only : gamma, R
  use mod_constant, only : gamma_1, mu0_T0_S_over_T0_2_3
  implicit none
contains
  subroutine calc_quantities_2D(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, Q_1, Q_2, Q_3, Q_4, T)
    integer, intent(in), value                 :: nx, ny
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny) ! Q / Jacobian
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny), Q_2(nx,ny), Q_3(nx,ny), Q_4(nx,ny)
    real(8), intent(inout), device, contiguous :: T(nx,ny)
    integer i, j
    real(8) :: over_Q1, rho, u, v, p
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny
      do i = 1, nx
        over_Q1  = 1.d0 / QJ_1(i,j)
        rho      = Jacobian(i,j) * QJ_1(i,j)
        u        = QJ_2(i,j) * over_Q1
        v        = QJ_3(i,j) * over_Q1
        p        = gamma_1 * (Jacobian(i,j) * QJ_4(i,j) - 0.5d0 * rho * (u*u + v*v))
        Q_1(i,j) = rho
        Q_2(i,j) = u
        Q_3(i,j) = v
        Q_4(i,j) = p
        T(i,j)   = p / (R * rho)
    enddo;enddo
  end subroutine calc_quantities_2D


  subroutine calc_quantities_T_2D(nx, ny, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, Q_1, Q_2, Q_3, Q_4, T, mu)
    integer, intent(in), value                 :: nx, ny
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ_1(nx,ny), QJ_2(nx,ny), QJ_3(nx,ny), QJ_4(nx,ny) ! Q / Jacobian
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny), Q_2(nx,ny), Q_3(nx,ny), Q_4(nx,ny)
    real(8), intent(inout), device, contiguous :: T(nx,ny)
    real(8), intent(inout), device, contiguous :: mu(nx,ny)
    integer i, j
    real(8) :: over_Q1, rho, u, v, p, temp
    !$cuf kernel do(2) <<<*,(32,4)>>>
    do j = 1, ny
      do i = 1, nx
        over_Q1  = 1.d0 / QJ_1(i,j)
        rho      = Jacobian(i,j) * QJ_1(i,j)
        u        = QJ_2(i,j) * over_Q1
        v        = QJ_3(i,j) * over_Q1
        p        = gamma_1 * (Jacobian(i,j) * QJ_4(i,j) - 0.5d0 * rho * (u*u + v*v))
        Q_1(i,j) = rho
        Q_2(i,j) = u
        Q_3(i,j) = v
        Q_4(i,j) = p
        temp     = p / (R * rho)
        T(i,j)   = temp
        mu(i,j)  = mu0_T0_S_over_T0_2_3 / (temp + 110.4d0) * (temp * sqrt(temp)) !mu(i,j)  = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo;enddo
  end subroutine calc_quantities_T_2D


  subroutine calc_quantities_3D(nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, &
                                 Q_1, Q_2, Q_3, Q_4, Q_5, T, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    REAL_T, intent(in), device, contiguous     :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz) ! Q / Jacobian
    REAL_T, intent(in), device, contiguous    :: QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: T(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p
#ifdef PREC_DF
    call prim_df<<<dim3((nx+31)/32, (ny+3)/4, (k_hi-k_lo+2)/2), dim3(32,4,2)>>>(nx, ny, nz, k_lo, k_hi, Jacobian, &
         QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Q_1, Q_2, Q_3, Q_4, Q_5, T, T, .false.)
#else
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          over_Q1  = 1.d0 / QJ_1(i,j,k)
          rho      = Jacobian(i,j) * QJ_1(i,j,k)
          u        = QJ_2(i,j,k) * over_Q1
          v        = QJ_3(i,j,k) * over_Q1
          w        = QJ_4(i,j,k) * over_Q1
          p        = gamma_1 * (Jacobian(i,j) * QJ_5(i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q_1(i,j,k) = rho
          Q_2(i,j,k) = u
          Q_3(i,j,k) = v
          Q_4(i,j,k) = w
          Q_5(i,j,k) = p
          T(i,j,k)   = p / (R * rho)
    enddo;enddo;enddo
#endif
  end subroutine calc_quantities_3D


  subroutine calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, &
                                   Q_1, Q_2, Q_3, Q_4, Q_5, T, mu, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    REAL_T, intent(in), device, contiguous     :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz) ! Q / Jacobian
    REAL_T, intent(in), device, contiguous    :: QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: T(nx,ny,nz)
    REAL_T, intent(inout), device, contiguous :: mu(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p, temp
#ifdef PREC_DF
    call prim_df<<<dim3((nx+31)/32, (ny+3)/4, (k_hi-k_lo+2)/2), dim3(32,4,2)>>>(nx, ny, nz, k_lo, k_hi, Jacobian, &
         QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Q_1, Q_2, Q_3, Q_4, Q_5, T, mu, .true.)
#else
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = k_lo, k_hi
      do j = 1, ny
        do i = 1, nx
          over_Q1  = 1.d0 / QJ_1(i,j,k)
          rho      = Jacobian(i,j) * QJ_1(i,j,k)
          u        = QJ_2(i,j,k) * over_Q1
          v        = QJ_3(i,j,k) * over_Q1
          w        = QJ_4(i,j,k) * over_Q1
          p        = gamma_1 * (Jacobian(i,j) * QJ_5(i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q_1(i,j,k) = rho
          Q_2(i,j,k) = u
          Q_3(i,j,k) = v
          Q_4(i,j,k) = w
          Q_5(i,j,k) = p
          temp       = p / (R * rho)
          T(i,j,k)   = temp
          mu(i,j,k)  = mu0_T0_S_over_T0_2_3 / (temp + 110.4d0) * (temp * sqrt(temp)) !mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo;enddo;enddo
#endif
  end subroutine calc_quantities_T_3D

#ifdef PREC_DF


  !> PREC = 'DF': CUF kernels cannot call the fltflt operators (nvfortran 26.5), so
  !> calc_quantities_3D/_T_3D run as this kernel (same formulas, mu only if visc).
  attributes(global) subroutine prim_df(nx, ny, nz, k_lo, k_hi, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, &
                                        Q_1, Q_2, Q_3, Q_4, Q_5, T, mu, visc)
    integer, intent(in), value    :: nx, ny, nz, k_lo, k_hi
    logical, intent(in), value    :: visc
    real(8), intent(in), device   :: Jacobian(nx,ny)
    REAL_T, intent(in), device    :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    REAL_T, intent(inout), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    REAL_T, intent(inout), device :: T(nx,ny,nz), mu(nx,ny,nz)
    REAL_T over_Q1, rho, u, v, w, p, temp
    integer i, j, k
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + k_lo - 1
    if (i > nx .or. j > ny .or. k > k_hi) return
    over_Q1 = 1.d0 / QJ_1(i,j,k)
    rho     = Jacobian(i,j) * QJ_1(i,j,k)
    u       = QJ_2(i,j,k) * over_Q1
    v       = QJ_3(i,j,k) * over_Q1
    w       = QJ_4(i,j,k) * over_Q1
    p       = gamma_1 * (Jacobian(i,j) * QJ_5(i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
    temp    = p / (R * rho)
    Q_1(i,j,k) = rho; Q_2(i,j,k) = u; Q_3(i,j,k) = v; Q_4(i,j,k) = w; Q_5(i,j,k) = p; T(i,j,k) = temp
    if (visc) mu(i,j,k) = mu0_T0_S_over_T0_2_3 / (temp + 110.4d0) * (temp * sqrt(temp))
  end subroutine prim_df
#endif
end module calc_physical_quantities

