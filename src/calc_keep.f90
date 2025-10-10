module calc_keep
  use mod_globals, only : dim => dimension, gamma
  use mod_constant, only : R_over_gamma_1, one_twelfth, two_third, seven_twelfth
  use calc_term
  implicit none
contains
  !KEEP energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 6th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  attributes(device) function  Et6(rho, p, u, v, w, uu, RhoV) result(Et)
    real(8), intent(in), dimension(6) :: rho, p, u, v, w, uu, RhoV
    real(8) Et(6)
    Et(:) = RhoPhiP_Rho6(RhoV(:), p(:), rho(:)) * R_over_gamma_1 ! internal energy
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi6(uu(:), p(:)) ! pressure diffusion
  end function Et6

  !KEEP main!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  attributes(device) function KEEP2(rho, u, v, w, uu, p, T, Normal) result(F)
    real(8), intent(in), dimension(2) :: rho, u, v, w, uu, p, T
    real(8), intent(in), dimension(5) :: Normal
    real(8) F(5)
    F(1) = 0.25d0 * (rho(1) + rho(2)) * (uu(1) + uu(2))
    F(2) = 0.5d0 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)) * Normal(2))
    F(3) = 0.5d0 * (F(1) * (v(1) + v(2)) + (p(1) + p(2)) * Normal(3))
    F(4) = 0.5d0 * (F(1) * (w(1) + w(2)) + (p(1) + p(2)) * Normal(4))
    F(5) = F(1) * 0.5d0 * (T(1) + T(2)) * R_over_gamma_1 ! internal energy
    F(5) = F(5) + 0.5d0 * (uu(1) * p(2) + uu(2) * p(1)) ! pressure diffusion
    F(5) = F(5) + 0.5d0 * F(1) * (u(1) * u(2) + v(1) * v(2) + w(1) * w(2)) ! kinetic energy
  end function KEEP2


  !$dir inline
  attributes(device) function KEEP4(rho, u, v, w, uu, p, T, Normal) result(F)
    real(8), intent(in), dimension(4)  :: rho, u, v, w, uu, p, T
    real(8), intent(in), dimension(5)  :: Normal
    real(8) F(5), RhoV1, RhoV2, RhoV3, RhoV1_RhoV2, RhoV1_RhoV3 
    real(8) u1, u2, u3, u4, v1, v2, v3, v4, w1, w2, w3, w4
    RhoV1 = 0.25d0 * (rho(2) + rho(3)) * (uu(2) + uu(3))
    RhoV2 = 0.25d0 * (rho(2) + rho(4)) * (uu(2) + uu(4))
    RhoV3 = 0.25d0 * (rho(1) + rho(3)) * (uu(1) + uu(3))
    F(1)  = four_third * RhoV1 - (RhoV2 + RhoV3) * one_sixth
    RhoV1 = two_third * RhoV1
    RhoV2 = one_twelfth * RhoV2
    RhoV3 = one_twelfth * RhoV3
    RhoV1_RhoV2 = RhoV1 - RhoV2
    RhoV1_RhoV3 = RhoV1 - RhoV3
    u1 = u(1); u2 = u(2); u3 = u(3); u4 = u(4)
    v1 = v(1); v2 = v(2); v3 = v(3); v4 = v(4)
    w1 = w(1); w2 = w(2); w3 = w(3); w4 = w(4)
    block
      real(8) pres
      pres = -one_twelfth * p(1) + seven_twelfth * (p(2) + p(3)) - one_twelfth * p(4)
      F(2) = -RhoV3 * u1 + RhoV1_RhoV2 * u2 + RhoV1_RhoV3 * u3 - RhoV2 * u4 + pres * Normal(2)
      F(3) = -RhoV3 * v1 + RhoV1_RhoV2 * v2 + RhoV1_RhoV3 * v3 - RhoV2 * v4 + pres * Normal(3)
      F(4) = -RhoV3 * w1 + RhoV1_RhoV2 * w2 + RhoV1_RhoV3 * w3 - RhoV2 * w4 + pres * Normal(4)
    end block
    block
      real(8) ene
      ! internal energy
      ene  = (-RhoV3 * T(1) + RhoV1_RhoV2 * T(2) + RhoV1_RhoV3 * T(3) - RhoV2 * T(4)) * R_over_gamma_1
      ! kinetic energy
      ene  = ene + ((RhoV1 * (u2*u3 + v2*v3 + w2*w3)) &
                  - (RhoV2 * (u2*u4 + v2*v4 + w2*w4) &
                   + RhoV3 * (u1*u3 + v1*v3 + w1*w3)))
      ! pressure diffusion
      F(5) = ene + (two_third * (uu(2)*p(3)+uu(3)*p(2)) &
                - one_twelfth * (uu(2)*p(4)+uu(4)*p(2) &
                                +uu(1)*p(3)+uu(3)*p(1)))
    end block
  end function KEEP4


  attributes(device) function KEEP6(rho, u, v, w, uu, p, T, Normal) result(F)
    real(8), intent(in), dimension(6)     :: rho, u, v, w, uu, p, T
    real(8), intent(in), dimension(dim+2) :: Normal
    real(8) F(dim+2), RhoV(6)
    RhoV(:) = RhoPhi6(rho(:), uu(:))
    F(1)    = Flux6(RhoV(:))
    block
      real(8) RhoVV_P(6)
      RhoVV_P(:) = RhoPhiU6(RhoV(:), u(:)) + Phi6(p(:)) * Normal(2)
      F(2)       = Flux6(RhoVV_P(:))
      RhoVV_P(:) = RhoPhiU6(RhoV(:), v(:)) + Phi6(p(:)) * Normal(3)
      F(3)       = Flux6(RhoVV_P(:))
      RhoVV_P(:) = RhoPhiU6(RhoV(:), w(:)) + Phi6(p(:)) * Normal(4)
      F(4)       = Flux6(RhoVV_P(:))
    end block
    block
      real(8) Energy(6)
      Energy(:) = Et6(rho, p, u, v, w, uu, RhoV)
      F(dim+2)  = Flux6(Energy(:))
    end block
  end function KEEP6
end module calc_keep

