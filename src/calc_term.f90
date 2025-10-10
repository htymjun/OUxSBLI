module calc_term
  use mod_constant, only : four_third, one_thirty, one_sixth
  implicit none
contains
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

  !KEEP 6th!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !$dir inline
  attributes(device) function Flux6(ph) result(ans)
    real(8), intent(in), dimension(6), device :: ph
    real(8) :: ans
    ans = 1.5d0 * ph(1) - 0.3d0 * (ph(2) + ph(3)) + (ph(4) + ph(5) + ph(6)) * one_thirty
  end function Flux6
end module

