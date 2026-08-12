module calc_physical_quantities
  use cudafor
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
        mu(i,j)  = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo;enddo
  end subroutine calc_quantities_T_2D


  !> 1D: no Jacobian (uniform grid) -- decode primitive rho/u/p and Sutherland's
  !> mu(T) directly from the conservative Q(rho, rho*u, E)
  subroutine calc_quantities_T_1D(nx, Q, ruvwp, T, mu)
    integer, intent(in), value               :: nx
    real(8), intent(in), device, contiguous  :: Q(nx,3)
    real(8), intent(out), device, contiguous :: ruvwp(nx,3)
    real(8), intent(out), device, contiguous :: T(nx)
    real(8), intent(out), device, contiguous :: mu(nx)
    integer i
    real(8) :: rho, u, p, temp
    !$cuf kernel do(1) <<<*,128>>>
    do i = 1, nx
      rho        = Q(i,1)
      u          = Q(i,2) / rho
      p          = gamma_1 * (Q(i,3) - 0.5d0 * rho * u*u)
      ruvwp(i,1) = rho
      ruvwp(i,2) = u
      ruvwp(i,3) = p
      temp       = p / (R * rho)
      T(i)       = temp
      mu(i)      = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo
  end subroutine calc_quantities_T_1D


  subroutine calc_quantities_3D(nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, &
                                 Q_1, Q_2, Q_3, Q_4, Q_5, T, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, contiguous    :: QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: T(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p
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
  end subroutine calc_quantities_3D


  subroutine calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, &
                                   Q_1, Q_2, Q_3, Q_4, Q_5, T, mu, k_lo, k_hi)
    integer, intent(in), value                 :: nx, ny, nz, k_lo, k_hi
    real(8), intent(in), device, contiguous    :: Jacobian(nx,ny)
    real(8), intent(in), device, contiguous    :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz) ! Q / Jacobian
    real(8), intent(in), device, contiguous    :: QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: T(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: mu(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p, temp
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
          mu(i,j,k)  = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * (temp * sqrt(temp))
    enddo;enddo;enddo
  end subroutine calc_quantities_T_3D
end module calc_physical_quantities

