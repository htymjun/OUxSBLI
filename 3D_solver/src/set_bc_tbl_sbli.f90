module set_bc_tbl_sbli
  use mod_globals, only : R, gamma
  use mod_constant, only : gamma_1, over_gamma, over_gamma_1, Cp
  implicit none
contains
  subroutine set_bc_Neumann_tbl_top_down(nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    use mod_globals, only : Taw, rf, u0, p0
    integer, intent(in), value     :: nx, ny, nz, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8) :: Jacobian_tmp, over_QJ1, p_wall, rhob, ub, vb, wb, pb
    integer i, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1+offset, nz-offset
      do i = istart, iend
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        over_QJ1 = 1.d0 / QJ_1(i,ny-1,k)
        rhob = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
        ub   = QJ_2(i,ny-1,k) * over_QJ1
        vb   = QJ_3(i,ny-1,k) * over_QJ1
        wb   = QJ_4(i,ny-1,k) * over_QJ1
        pb   = gamma_1 * (QJ_5(i,ny-1,k) * Jacobian(i,ny-1) - 0.5d0 * rhob * (ub**2 + vb**2 + wb**2))
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
        p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
        QJ_5(i,1,k) = p_wall * over_gamma_1
    enddo;enddo
  end subroutine set_bc_Neumann_tbl_top_down


  subroutine set_bc_Riemann_tbl_top_down(nx, ny, nz, offset, istart, iend, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    use mod_globals, only : Taw, rf, u0, p0
    integer, intent(in), value     :: nx, ny, nz, offset, istart, iend
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), parameter :: T       = Taw - rf * u0**2 / (2.d0 * Cp)
    real(8), parameter :: rho0    = p0 / (R * T)
    real(8), parameter :: c0      = sqrt(gamma * p0 / rho0)
    real(8), parameter :: over_c0 = 1.d0 / c0
    real(8) :: Jacobian_tmp, p_wall, rhoin, pin, cin, vin, Rp, Rm, rhob, vb, cb, pb, v0 = 0.d0
    integer i, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 1+offset, nz-offset
      do i = istart, iend
        Jacobian_tmp = 1.d0 / Jacobian(i,ny)
        pin   = gamma_1 * (QJ_5(i,ny-1,k) - 0.5d0 * (QJ_2(i,ny-1,k)**2 + QJ_3(i,ny-1,k)**2 + QJ_4(i,ny-1,k)**2) &
                / QJ_1(i,ny-1,k)) * Jacobian(i,ny-1)
        rhoin = QJ_1(i,ny-1,k) * Jacobian(i,ny-1)
        cin   = sqrt(gamma * pin / rhoin)
        vin   = QJ_3(i,ny-1,k) / QJ_1(i,ny-1,k)
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
        p_wall = gamma_1 * (QJ_5(i,2,k) - 0.5d0 * (QJ_2(i,2,k)**2 + QJ_3(i,2,k)**2 + QJ_4(i,2,k)**2) / QJ_1(i,2,k))
        QJ_5(i,1,k) = p_wall * over_gamma_1
    enddo;enddo
  end subroutine set_bc_Riemann_tbl_top_down
end module set_bc_tbl_sbli
