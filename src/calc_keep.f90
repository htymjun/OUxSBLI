module calc_keep
  use mod_globals, only : dim => dimension, gamma
  use mod_constant, only : over_gamma_1
  use calc_common_dim
  use calc_term
  use calc_mat
  implicit none
  interface Et2
    module procedure EtKEEP2, EtKEEPPE2, EtKEP2
  end interface

  interface Et4
    module procedure EtKEEP4, EtKEEPPE4, EtKEP4
  end interface

  interface Et6
    module procedure EtKEEP6, EtKEEPPE6, EtKEP6
  end interface
contains
  !KEEP energy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2nd-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function EtKEEP2(id_keep, C, rho, p, u, v, w, uu) result(Et)
    integer(kind=2), intent(in), value :: id_keep
    real(8), intent(in), value         :: C
    real(8), intent(in), dimension(2)  :: rho, p, u, v, w, uu
    real(8) Et
    Et = C * 0.5d0 * (p(1) / rho(1) + p(2) / rho(2)) * over_gamma_1 ! internal energy
    Et = Et + 0.5d0 * (uu(1) * p(2) + uu(2) * p(1)) ! pressure diffusion
    Et = Et + 0.5d0 * C * (u(1) * u(2) + v(1) * v(2) + w(1) * w(2)) ! kinetic energy
  end function EtKEEP2

  attributes(device) function EtKEEPPE2(id_keep, C, rho, p, u, v, w, uu) result(Et)
    integer(kind=4), intent(in), value :: id_keep
    real(8), intent(in), value         :: C
    real(8), intent(in), dimension(2)  :: rho, p, u, v, w, uu
    real(8) :: Et
    Et = 0.25d0 * (uu(1) + uu(2)) * (p(1) + p(2)) * over_gamma_1 ! internal energy
    Et = Et + 0.5d0 * (uu(1) * p(2) + uu(2) * p(1)) ! pressure diffusion
    Et = Et + 0.5d0 * C * (u(1) * u(2) + v(1) * v(2) + w(1) * w(2)) ! kinetic energy
  end function EtKEEPPE2

  attributes(device) function EtKEP2(id_keep, C, rho, p, u, v, w, uu) result(Et)
    integer(kind=8), intent(in), value :: id_keep
    real(8), intent(in), value         :: C
    real(8), intent(in), dimension(2)  :: rho, p, u, v, w, uu
    real(8) :: Et
    Et = C * 0.5d0 * gamma * over_gamma_1 * (p(1) / rho(1) + p(2) / rho(2)) ! enthalpy
    Et = Et + 0.5d0 * C * (u(1) * u(2) + v(1) * v(2) + w(1) * w(2)) ! kinetic energy
  end function EtKEP2

  ! 4th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function  EtKEEP4(id_keep, rho, p, u, v, w, uu, RhoV) result(Et)
    integer(kind=2), intent(in), value :: id_keep
    real(8), intent(in), dimension(4)  :: rho, p, u, v, w, uu
    real(8), intent(in), dimension(3)  :: RhoV
    real(8) Et(3)
    Et(:) = RhoPhiP_Rho4(RhoV(:), p(:), rho(:)) * over_gamma_1 ! internal energy
    Et(:) = Et(:) + RhoUPhiPhi4(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi4(uu(:), p(:)) ! pressure diffusion
  end function EtKEEP4

  attributes(device) function EtKEEPPE4(id_keep, rho, p, u, v, w, uu, RhoV) result(Et)
    integer(kind=4), intent(in), value :: id_keep
    real(8), intent(in), dimension(4)  :: rho, p, u, v, w, uu
    real(8), intent(in), dimension(3)  :: RhoV
    real(8) Et(3)
    block
      real(8) Vm(3)
      Vm(:) = Phi4(uu(:))
      Et(:) = RhoPhiU4(Vm(:), p(:)) * over_gamma_1 ! internal energy    
    end block
    Et(:) = Et(:) + RhoUPhiPhi4(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi4(uu(:), p(:)) ! pressure diffusion
  end function EtKEEPPE4

  attributes(device) function EtKEP4(id_keep, rho, p, u, v, w, uu, RhoV) result(Et)
    integer(kind=8), intent(in), value :: id_keep
    real(8), intent(in), dimension(4)  :: rho, p, u, v, w, uu
    real(8), intent(in), dimension(3)  :: RhoV
    real(8) Et(3)
    Et(:) = RhoPhiP_Rho4(RhoV(:), p(:), rho(:)) * gamma * over_gamma_1 ! enthalpy
    Et(:) = Et(:) + RhoUPhiPhi4(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
  end function EtKEP4

  ! 6th-order accuracy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function  EtKEEP6(id_keep, rho, p, u, v, w, uu, RhoV) result(Et)
    integer(kind=2), intent(in), value :: id_keep
    real(8), intent(in), dimension(6)  :: rho, p, u, v, w, uu, RhoV
    real(8) Et(6)
    Et(:) = RhoPhiP_Rho6(RhoV(:), p(:), rho(:)) * over_gamma_1 ! internal energy
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi6(uu(:), p(:)) ! pressure diffusion
  end function EtKEEP6

  attributes(device) function EtKEEPPE6(id_keep, rho, p, u, v, w, uu, RhoV) result(Et)
    integer(kind=4), intent(in), value :: id_keep
    real(8), intent(in), dimension(6)  :: rho, p, u, v, w, uu, RhoV
    real(8) Et(6)
    block
      real(8) Vm(6)
      Vm(:) = Phi6(uu(:))
      Et(:) = RhoPhiU6(Vm(:), p(:)) * over_gamma_1 ! internal energy    
    end block
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
    Et(:) = Et(:) + PhiPsi6(uu(:), p(:)) ! pressure diffusion
  end function EtKEEPPE6

  attributes(device) function EtKEP6(id_keep, rho, p, u, v, w, uu, RhoV) result(Et)
    integer(kind=8), intent(in), value :: id_keep
    real(8), intent(in), dimension(6)  :: rho, p, u, v, w, uu, RhoV
    real(8) Et(6)
    Et(:) = RhoPhiP_Rho6(RhoV(:), p(:), rho(:)) * gamma * over_gamma_1 ! enthalpy
    Et(:) = Et(:) + RhoUPhiPhi6(RhoV(:), u(:), v(:), w(:)) ! kinetic energy
  end function EtKEP6

  !KEEP main!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  attributes(device) function KEEP2(rho, u, v, w, uu, p, Normal) result(F)
    use mod_globals, only : id_keep
    real(8), intent(in), dimension(2)     :: rho, u, v, w, uu, p
    real(8), intent(in), dimension(dim+2) :: Normal
    real(8) F(dim+2)
    F(1) = 0.25d0 * (rho(1) + rho(2)) * (uu(1) + uu(2))
    F(2) = 0.5d0 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)) * Normal(2))
    F(3) = 0.5d0 * (F(1) * (v(1) + v(2)) + (p(1) + p(2)) * Normal(3))
    F(4) = 0.5d0 * (F(1) * (w(1) + w(2)) + (p(1) + p(2)) * Normal(4))
    F(5) = Et2(id_keep,F(1),rho,p,u,v,w,uu)
  end function KEEP2

  attributes(device) function KEEP4(rho, u, v, w, uu, p, Normal) result(F)
    use mod_globals, only : id_keep
    real(8), intent(in), dimension(4)     :: rho, u, v, w, uu, p
    real(8), intent(in), dimension(dim+2) :: Normal
    real(8) F(dim+2), RhoV(3)
    RhoV(:) = RhoPhi4(rho(:), uu(:))
    F(1)    = Flux4(RhoV(:))
    block
      real(8) RhoVV_P(3)
      RhoVV_P(:) = RhoPhiU4(RhoV(:), u(:)) + Phi4(p(:)) * Normal(2)
      F(2)       = Flux4(RhoVV_P(:))
      RhoVV_P(:) = RhoPhiU4(RhoV(:), v(:)) + Phi4(p(:)) * Normal(3)
      F(3)       = Flux4(RhoVV_P(:))
      RhoVV_P(:) = RhoPhiU4(RhoV(:), w(:)) + Phi4(p(:)) * Normal(4)
      F(4)       = Flux4(RhoVV_P(:))
    end block
    block
      real(8) Energy(3)
      Energy(:) =  Et4(id_keep, rho, p, u, v, w, uu, RhoV)
      F(dim+2)  = Flux4(Energy(:))
    end block
  end function KEEP4

  attributes(device) function KEEP6(rho, u, v, w, uu, p, Normal) result(F)
    use mod_globals, only : id_keep
    real(8), intent(in), dimension(6)     :: rho, u, v, w, uu, p
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
      Energy(:) = Et6(id_keep, rho, p, u, v, w, uu, RhoV)
      F(dim+2)  = Flux6(Energy(:))
    end block
  end function KEEP6
end module calc_keep

