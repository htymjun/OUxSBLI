module calc_physical_quantities
  use cudafor
  use mod_globals, only : gamma, R, dim => dimension
  use mod_constant, only : gamma_1, mu0_T0_S_over_T0_2_3
  implicit none
contains
  subroutine calc_quantities_2D(nx,ny,Jacobian,QJ,rho,u,v,p)
    integer, intent(in), value   :: nx, ny
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(4,nx,ny) ! Q / Jacobian
    real(8), intent(out), device :: rho(nx,ny), u(nx,ny), v(nx,ny), p(nx,ny)
    integer i, j
    real(8) :: over_Q1
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        over_Q1  = 1.d0 / QJ(1,i,j)
        rho(i,j) = Jacobian(i,j) * QJ(1,i,j)
        u(i,j)   = QJ(2,i,j) * over_Q1
        v(i,j)   = QJ(3,i,j) * over_Q1
        p(i,j)   = gamma_1 * (Jacobian(i,j) * QJ(4,i,j)- 0.5d0 * rho(i,j) * (u(i,j)*u(i,j) + v(i,j)*v(i,j)))
    enddo;enddo
  end subroutine calc_quantities_2D


  subroutine calc_quantities_3D(nx, ny, nz, Jacobian, QJ, Q, T)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: Q(5,nx,ny,nz), T(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(1,i,j,k)
          rho        = Jacobian(i,j) * QJ(1,i,j,k)
          u          = QJ(2,i,j,k) * over_Q1
          v          = QJ(3,i,j,k) * over_Q1
          w          = QJ(4,i,j,k) * over_Q1
          p          = gamma_1 * (Jacobian(i,j) * QJ(5,i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q(1,i,j,k) = rho
          Q(2,i,j,k) = u
          Q(3,i,j,k) = v
          Q(4,i,j,k) = w
          Q(5,i,j,k) = p
          T(i,j,k)   = p / (R * rho)
    enddo;enddo;enddo
  end subroutine calc_quantities_3D
  

  subroutine calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, Q, T, mu)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p, temp
    !$cuf kernel do(3) <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(1,i,j,k)
          rho        = Jacobian(i,j) * QJ(1,i,j,k)
          u          = QJ(2,i,j,k) * over_Q1
          v          = QJ(3,i,j,k) * over_Q1
          w          = QJ(4,i,j,k) * over_Q1
          p          = gamma_1 * (Jacobian(i,j) * QJ(5,i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q(1,i,j,k) = rho
          Q(2,i,j,k) = u
          Q(3,i,j,k) = v
          Q(4,i,j,k) = w
          Q(5,i,j,k) = p
          temp       = p / (R * rho)
          T(i,j,k)   = temp
          mu(i,j,k)  = mu0_T0_S_over_T0_2_3 / (temp + 111.d0) * temp**1.5d0
    enddo;enddo;enddo
  end subroutine calc_quantities_T_3D
end module calc_physical_quantities

