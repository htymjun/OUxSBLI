module calc_term
  use mod_globals, only : dimension
  use mod_constant, only : one_third, one_twelfth, one_sixty
  use calc_common_dim
  implicit none
contains
  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function Phi4(a) result(ans)
    real(8), intent(in), device :: a(4)
    real(8) ans(3)
    ans(1) = 0.5d0 * (a(2) + a(3))
    ans(2) = 0.5d0 * (a(2) + a(4))
    ans(3) = 0.5d0 * (a(1) + a(3))
  end function Phi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function Phi6(a) result(ans)
    real(8), intent(in), device :: a(6)
    real(8) ans(6)
    ans(1) = 0.5d0 * (a(3) + a(4))
    ans(2) = 0.5d0 * (a(3) + a(5))
    ans(3) = 0.5d0 * (a(2) + a(4))
    ans(4) = 0.5d0 * (a(3) + a(6))
    ans(5) = 0.5d0 * (a(2) + a(5))
    ans(6) = 0.5d0 * (a(1) + a(4))
  end function Phi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoPhi4(rho, u) result(ans)
    real(8), intent(in), dimension(4), device :: rho, u
    real(8) ans(3)
    ans(1) = 0.25d0 * (rho(2) + rho(3)) * (u(2) + u(3))
    ans(2) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(3) = 0.25d0 * (rho(1) + rho(3)) * (u(1) + u(3))
  end function RhoPhi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !$dir inline
  attributes(device) function RhoPhi6(rho, u) result(ans)
    real(8), intent(in), dimension(6), device :: rho, u
    real(8) ans(6)
    ans(1) = 0.25d0 * (rho(3) + rho(4)) * (u(3) + u(4))
    ans(2) = 0.25d0 * (rho(3) + rho(5)) * (u(3) + u(5))
    ans(3) = 0.25d0 * (rho(2) + rho(4)) * (u(2) + u(4))
    ans(4) = 0.25d0 * (rho(3) + rho(6)) * (u(3) + u(6))
    ans(5) = 0.25d0 * (rho(2) + rho(5)) * (u(2) + u(5))
    ans(6) = 0.25d0 * (rho(1) + rho(4)) * (u(1) + u(4))
  end function RhoPhi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoPhiU4(rhou, ph) result(ans)
    real(8), intent(in), dimension(3), device :: rhou
    real(8), intent(in), dimension(4), device :: ph
    real(8) ans(3)
    ans(1) = rhou(1) * 0.5d0 * (ph(2) + ph(3))
    ans(2) = rhou(2) * 0.5d0 * (ph(2) + ph(4))
    ans(3) = rhou(3) * 0.5d0 * (ph(1) + ph(3))
  end function RhoPhiU4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoPhiU6(rhou, ph) result(ans)
    real(8), intent(in), dimension(6), device :: rhou, ph
    real(8) ans(6)
    ans(1) = rhou(1) * 0.5d0 * (ph(3) + ph(4))
    ans(2) = rhou(2) * 0.5d0 * (ph(3) + ph(5))
    ans(3) = rhou(3) * 0.5d0 * (ph(2) + ph(4))
    ans(4) = rhou(4) * 0.5d0 * (ph(3) + ph(6))
    ans(5) = rhou(5) * 0.5d0 * (ph(2) + ph(5))
    ans(6) = rhou(6) * 0.5d0 * (ph(1) + ph(4))
  end function RhoPhiU6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoPhiP_Rho4(rhou, p, rho) result(ans)
    real(8), intent(in), dimension(3), device :: rhou
    real(8), intent(in), dimension(4), device :: p, rho
    real(8) ans(3)
    ans(1) = rhou(1) * 0.5d0 * (p(2) / rho(2) + p(3) / rho(3))
    ans(2) = rhou(2) * 0.5d0 * (p(2) / rho(2) + p(4) / rho(4))
    ans(3) = rhou(3) * 0.5d0 * (p(1) / rho(1) + p(3) / rho(3))
  end function RhoPhiP_Rho4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoPhiP_Rho6(rhou, p, rho) result(ans)
    real(8), intent(in), dimension(6), device :: rhou, p, rho
    real(8) ans(6)
    ans(1) = rhou(1) * 0.5d0 * (p(3) / rho(3) + p(4) / rho(4))
    ans(2) = rhou(2) * 0.5d0 * (p(3) / rho(3) + p(5) / rho(5))
    ans(3) = rhou(3) * 0.5d0 * (p(2) / rho(2) + p(4) / rho(4))
    ans(4) = rhou(4) * 0.5d0 * (p(3) / rho(3) + p(6) / rho(6))
    ans(5) = rhou(5) * 0.5d0 * (p(2) / rho(2) + p(5) / rho(5))
    ans(6) = rhou(6) * 0.5d0 * (p(1) / rho(1) + p(4) / rho(4))
  end function RhoPhiP_Rho6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoUPhiPhi4(rhou, u, v, w) result(ans)
    real(8), intent(in), dimension(3), device :: rhou
    real(8), intent(in), dimension(4), device :: u, v, w
    real(8) ans(3)
    ans(1) = rhou(1) * 0.5d0 * (u(2) * u(3) + v(2) * v(3) + w(2) * w(3))
    ans(2) = rhou(2) * 0.5d0 * (u(2) * u(4) + v(2) * v(4) + w(2) * w(4))
    ans(3) = rhou(3) * 0.5d0 * (u(1) * u(3) + v(1) * v(3) + w(1) * w(3))
  end function RhoUPhiPhi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function RhoUPhiPhi6(rhou, u, v, w) result(ans)
    real(8), intent(in), dimension(6), device :: rhou, u, v, w
    real(8) ans(6)
    ans(1) = rhou(1) * 0.5d0 * (u(3) * u(4) + v(3) * v(4) + w(3) * w(4))
    ans(2) = rhou(2) * 0.5d0 * (u(3) * u(5) + v(3) * v(5) + w(3) * w(5))
    ans(3) = rhou(3) * 0.5d0 * (u(2) * u(4) + v(2) * v(4) + w(2) * w(4))
    ans(4) = rhou(4) * 0.5d0 * (u(3) * u(6) + v(3) * v(6) + w(3) * w(6))
    ans(5) = rhou(5) * 0.5d0 * (u(2) * u(5) + v(2) * v(5) + w(2) * w(5))
    ans(6) = rhou(6) * 0.5d0 * (u(1) * u(4) + v(1) * v(4) + w(1) * w(4))
  end function RhoUPhiPhi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function PhiPsi4(ph, psi) result(ans)
    real(8), intent(in), dimension(4), device :: ph, psi
    real(8) ans(3)
    ans(1) = 0.5d0 * (ph(2) * psi(3) + ph(3) * psi(2))
    ans(2) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(3) = 0.5d0 * (ph(1) * psi(3) + ph(3) * psi(1))
  end function PhiPsi4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function PhiPsi6(ph, psi) result(ans)
    real(8), intent(in), dimension(6), device :: ph, psi
    real(8) ans(6)
    ans(1) = 0.5d0 * (ph(3) * psi(4) + ph(4) * psi(3))
    ans(2) = 0.5d0 * (ph(3) * psi(5) + ph(5) * psi(3))
    ans(3) = 0.5d0 * (ph(2) * psi(4) + ph(4) * psi(2))
    ans(4) = 0.5d0 * (ph(3) * psi(6) + ph(6) * psi(3))
    ans(5) = 0.5d0 * (ph(2) * psi(5) + ph(5) * psi(2))
    ans(6) = 0.5d0 * (ph(1) * psi(4) + ph(4) * psi(1))
  end function PhiPsi6

  !KEEP 4th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$dir inline
  attributes(device) function Flux4(ph) result(ans)
    real(8), intent(in), dimension(3), device :: ph
    real(8) :: ans
    ans = 2.d0 * ((2.d0 * one_third) * ph(1) - (ph(2) + ph(3)) * one_twelfth)
  end function Flux4

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !$dir inline
  attributes(device) function Flux6(ph) result(ans)
    real(8), intent(in), dimension(6), device :: ph
    real(8) :: ans
    ans = 2.d0 * (0.75d0 * ph(1) - 3.d0 * (ph(2) + ph(3)) * 0.05d0 &
          + (ph(4) + ph(5) + ph(6)) * one_sixty)
  end function Flux6
end module

