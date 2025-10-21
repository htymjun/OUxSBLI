module calc_visc_common
  use mod_constant, only : Cp, Cp_over_Pr, one_third
  implicit none
contains
  !dir$ inline
  attributes(device) function interpolation6_scalar(a1, a2, a3, a4, a5, a6) result(ans)
    real(8), intent(in), value :: a1, a2, a3, a4, a5, a6
    real(8) ans(3)
    ans(1) = 0.0625d0 * (-a1 + 9.d0 * (a2 + a3) -a4)
    ans(2) = 0.0625d0 * (-a2 + 9.d0 * (a3 + a4) -a5)
    ans(3) = 0.0625d0 * (-a3 + 9.d0 * (a4 + a5) -a6)
  end function interpolation6_scalar

  !$dir inline
  attributes(device) function dy5(a1, a2, a3, a4, d) result(ans)
    real(8), intent(in), value :: a1, a2, a3, a4, d
    real(8) ans
    ans = (2.d0 * (-a2 + a3) - 0.25d0 * (-a1 + a4)) * d * one_third
  end function dy5

  !dir$ inline
  attributes(device) function dx6(a, dx) result(ans)
    real(8), intent(in), device :: a(6)
    real(8), intent(in), value  :: dx
    real(8) ans(3)
    ans(:) = 0.125d0 * (9.d0 * (-a(2:4) + a(3:5)) - (-a(1:3) + a(4:6)) * one_third) * dx
  end function dx6

  !dir$ inline
  attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8) ans
    ans = 0.125d0 * (26.d0 * a(2) - (a(1) + a(3))) * one_third
  end function flux4

  attributes(device) subroutine tauxx_4(mu, ux, vy, wz, u1, u2, u3, u4, u5, u6, txx, utxx)
    real(8), intent(in), dimension(3), device :: mu, ux, vy, wz
    real(8), intent(in), value :: u1, u2, u3, u4, u5, u6
    real(8), intent(out)  :: txx, utxx
    real(8), dimension(3) :: tmp
    tmp(:) = 2.d0 * mu(:) * (2.d0 * ux(:) - vy(:) - wz(:)) * one_third ! tau
    txx    = flux4(tmp(:))
    tmp(:) = interpolation6_scalar(u1, u2, u3, u4, u5, u6) * tmp(:) ! u tau
    utxx   = flux4(tmp(:))
  end subroutine tauxx_4

  attributes(device) subroutine tauxy_4(mu, uy, vx, v1, v2, v3, v4, v5, v6, txy, vtxy)
    real(8), intent(in), dimension(3), device :: mu, uy, vx
    real(8), intent(in), value :: v1, v2, v3, v4, v5, v6
    real(8), intent(out)  :: txy, vtxy
    real(8), dimension(3) :: tmp
    tmp(:) = mu(:) * (uy(:) + vx(:)) ! tau
    txy    = flux4(tmp(:))
    tmp(:) = interpolation6_scalar(v1, v2, v3, v4, v5, v6) * tmp(:) ! vtau
    vtxy   = flux4(tmp(:))
  end subroutine tauxy_4

  !dir$ inline
  attributes(device) function heat_conduction6(mu, T, dx) result(ans)
    real(8), intent(in), device :: mu(3), T(6)
    real(8), intent(in), value  :: dx
    real(8), dimension(3) ::  kTx
    real(8) :: ans
    kTx(:) = Cp_over_Pr * mu(:) * dx6(T(:), dx)
    ans = flux4(kTx(:))
  end function heat_conduction6
end module calc_visc_common

