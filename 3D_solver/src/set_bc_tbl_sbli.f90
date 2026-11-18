module set_bc_tbl_sbli
  use mod_globals, only : R, gamma
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1, Cp
  use mod_prec
  use cudafor
  implicit none
contains
  subroutine set_bc_Neumann_tbl_top_down(nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: nx, ny, nz, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    REAL_T, intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    ! hand-written kernel instead of a CUF loop: with PREC = 'DF' a CUF loop cannot hold
    ! fltflt assignments or dble() of the state
    call bc_Neumann_tbl_kernel<<<dim3((iend-istart+1+31)/32, (nz-2*offset+3)/4, 1), dim3(32, 4, 1)>>>( &
        nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
  end subroutine set_bc_Neumann_tbl_top_down


  attributes(global) subroutine bc_Neumann_tbl_kernel(nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    use mod_globals, only : Taw, rf, u0, p0
    integer, intent(in), value     :: nx, ny, nz, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    REAL_T, intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8) :: Jacobian_tmp, over_QJ1, p_wall, rhob, ub, vb, wb, pb
    integer i, k
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + istart - 1
    k = (blockIdx%y - 1) * blockDim%y + threadIdx%y + offset
    if (i > iend .or. k > nz-offset) return
    Jacobian_tmp = 1.d0 / Jacobian(i,ny)
    over_QJ1 = 1.d0 / dble(QJ_1(i,ny-1,k))
    rhob = dble(QJ_1(i,ny-1,k)) * Jacobian(i,ny-1)
    ub   = dble(QJ_2(i,ny-1,k)) * over_QJ1
    vb   = dble(QJ_3(i,ny-1,k)) * over_QJ1
    wb   = dble(QJ_4(i,ny-1,k)) * over_QJ1
    pb   = gamma_1 * (dble(QJ_5(i,ny-1,k)) * Jacobian(i,ny-1) - 0.5d0 * rhob * (ub**2 + vb**2 + wb**2))
    QJ_1(i,ny,k) = rhob * Jacobian_tmp
    QJ_2(i,ny,k) = rhob * ub * Jacobian_tmp
    QJ_3(i,ny,k) = rhob * vb * Jacobian_tmp
    QJ_4(i,ny,k) = rhob * wb * Jacobian_tmp
    QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (ub**2 + vb**2 + wb**2)) * Jacobian_tmp
    ! NoSlip
    QJ_1(i,1,k) = QJ_1(i,2,k)
    QJ_2(i,1,k) = 0.d0
    QJ_3(i,1,k) = 0.d0
    QJ_4(i,1,k) = 0.d0
    p_wall = gamma_1 * (dble(QJ_5(i,2,k)) - 0.5d0 * (dble(QJ_2(i,2,k))**2 + dble(QJ_3(i,2,k))**2 + dble(QJ_4(i,2,k))**2) / dble(QJ_1(i,2,k)))
    QJ_5(i,1,k) = p_wall * over_gamma_1
  end subroutine bc_Neumann_tbl_kernel


  subroutine set_bc_Riemann_tbl_top_down(nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: nx, ny, nz, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    REAL_T, intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    ! hand-written kernel instead of a CUF loop: with PREC = 'DF' a CUF loop cannot hold
    ! fltflt assignments or dble() of the state
    call bc_Riemann_tbl_kernel<<<dim3((iend-istart+1+31)/32, (nz-2*offset+3)/4, 1), dim3(32, 4, 1)>>>( &
        nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
  end subroutine set_bc_Riemann_tbl_top_down


  attributes(global) subroutine bc_Riemann_tbl_kernel(nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    use mod_globals, only : Taw, rf, u0, p0
    integer, intent(in), value     :: nx, ny, nz, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    REAL_T, intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), parameter :: T       = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0    = p0 / (R * T)
    real(8), parameter :: c0      = sqrt(gamma * p0 / rho0)
    real(8), parameter :: over_c0 = 1.d0 / c0
    real(8), parameter :: v0      = 0.d0
    real(8) :: Jacobian_tmp, p_wall, rhoin, pin, cin, vin, Rp, Rm, rhob, vb, cb, pb
    integer i, k
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + istart - 1
    k = (blockIdx%y - 1) * blockDim%y + threadIdx%y + offset
    if (i > iend .or. k > nz-offset) return
    Jacobian_tmp = 1.d0 / Jacobian(i,ny)
    pin   = gamma_1 * (dble(QJ_5(i,ny-1,k)) - 0.5d0 * (dble(QJ_2(i,ny-1,k))**2 + dble(QJ_3(i,ny-1,k))**2 + dble(QJ_4(i,ny-1,k))**2) &
            / dble(QJ_1(i,ny-1,k))) * Jacobian(i,ny-1)
    rhoin = dble(QJ_1(i,ny-1,k)) * Jacobian(i,ny-1)
    cin   = sqrt(gamma * pin / rhoin)
    vin   = dble(QJ_3(i,ny-1,k)) / dble(QJ_1(i,ny-1,k))
    Rp   = vin + 2.d0 * cin * over_gamma_1
    Rm   = v0  - 2.d0 * c0  * over_gamma_1
    vb   = 0.5d0 * (Rp + Rm)
    cb   = 0.25d0 * gamma_1 * (Rp - Rm)
    rhob = (cb * over_c0)**(2.d0 * over_gamma_1) * rho0
    pb   = (rhob * cb**2) * over_gamma
    QJ_1(i,ny,k) = rhob * Jacobian_tmp
    QJ_2(i,ny,k) = rhob * u0 * Jacobian_tmp
    QJ_3(i,ny,k) = rhob * vb * Jacobian_tmp
    QJ_4(i,ny,k) = 0.d0
    QJ_5(i,ny,k) = (pb * over_gamma_1 + 0.5d0 * rhob * (u0**2 + vb**2)) * Jacobian_tmp
    ! NoSlip
    QJ_1(i,1,k) = QJ_1(i,2,k)
    QJ_2(i,1,k) = 0.d0
    QJ_3(i,1,k) = 0.d0
    QJ_4(i,1,k) = 0.d0
    p_wall = gamma_1 * (dble(QJ_5(i,2,k)) - 0.5d0 * (dble(QJ_2(i,2,k))**2 + dble(QJ_3(i,2,k))**2 + dble(QJ_4(i,2,k))**2) / dble(QJ_1(i,2,k)))
    QJ_5(i,1,k) = p_wall * over_gamma_1
  end subroutine bc_Riemann_tbl_kernel
end module set_bc_tbl_sbli
