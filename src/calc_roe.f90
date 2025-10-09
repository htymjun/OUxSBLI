module calc_roe
  use mod_constant, only : gamma_1, over_gamma_1, gamma_over_gamma_1
  implicit none
contains
  !$dir inline
  attributes(device) subroutine calc_RARinv_fast(u, v, w, H, c, dQ, F)
    real(8), intent(in), value :: u, v, w, H, c
    real(8), intent(in)        :: dQ(5)
    real(8), intent(inout)     :: F(5)
    real(8) q2, R1, R2, R3, R4, R5
    q2 = 0.5d0 * (u*u + v*v + w*w)
    block
      real(8) oc, b1, b2, dQ1, dQ2, dQ3, dQ4, dQ5
      dQ1 = dQ(1); dQ2 = dQ(2); dQ3 = dQ(3); dQ4 = dQ(4); dQ5 = dQ(5)
      oc = 1.d0 / c
      b2 = gamma_1 * oc * oc
      b1 = q2 * b2
      R1 = 0.5d0*((b1+u*oc)*dQ1-(oc+b2*u)*dQ2-b2*v*dQ3-b2*w*dQ4+b2*dQ5)
      R2 = oc*( w*dQ1-dQ4)
      R3 = oc*(-v*dQ1+dQ3)
      R4 = (1.d0-b1)*dQ1+b2*(u*dQ2+v*dQ3+w*dQ4-dQ5)
      R5 = 0.5d0*((b1-u*oc)*dQ1+(oc-b2*u)*dQ2-b2*(v*dQ3+w*dQ4-dQ5))
    end block
    block
      real(8) a1, a2, a3, a2c, cu
      a1   = abs(u-c)
      a2   = abs(u)
      a3   = abs(u+c)
      a2c  = a2 * c
      cu   = c * u
      F(1) = F(1) - 0.5d0 * (a1*R1                            +a2*R4      +a3*R5)
      F(2) = F(2) - 0.5d0 * (a1*(u-c)*R1                    +a2*u*R4+a3*(u+c)*R5)
      F(3) = F(3) - 0.5d0 * (a1*v*R1                 +a2c*R3+a2*v*R4    +a3*v*R5)
      F(4) = F(4) - 0.5d0 * (a1*w*R1          -a2c*R2     +a2*w*R4      +a3*w*R5)
      F(5) = F(5) - 0.5d0 * (a1*(H-cu)*R1+a2*(-w*c*R2+v*c*R3+q2*R4)+a3*(H+cu)*R5)
    end block
  end subroutine calc_RARinv_fast


  !$dir inline
  attributes(device) subroutine calc_Roe_ave(rho, V, p, uroe, vroe, wroe, Hroe, croe)
    real(8), intent(in)  :: rho(2), V(2,3), p(2)
    real(8), intent(out) :: uroe, vroe, wroe, Hroe, croe
    real(8) rhol, rhor, rhol_rhor, H(2)
    H(:) = gamma_over_gamma_1 * p(:) / rho(:) + 0.5d0 * (V(:,1)*V(:,1) + V(:,2)*V(:,2) + V(:,3)*V(:,3))
    rhol = sqrt(rho(1))
    rhor = sqrt(rho(2))
    rhol_rhor = 1.d0 / (rhol + rhor)
    uroe = (rhol * V(1,1) + rhor * V(2,1)) * rhol_rhor
    vroe = (rhol * V(1,2) + rhor * V(2,2)) * rhol_rhor
    wroe = (rhol * V(1,3) + rhor * V(2,3)) * rhol_rhor
    Hroe = (rhol * H(1) + rhor * H(2)) * rhol_rhor
    croe = sqrt(gamma_1 * (Hroe - 0.5d0 * (uroe*uroe + vroe*vroe + wroe*wroe)))
  end subroutine calc_Roe_ave

 
  !$dir inline
  attributes(device) subroutine calc_Roe_common(rho, V, p, dQ, F)
    real(8), intent(in)  :: rho(2), V(2,3), p(2)
    real(8), intent(out) :: dQ(5), F(5)
    real(8) e(2)
    e(:)  = p(:) * over_gamma_1 + 0.5d0 * rho(:) * (V(:,1)*V(:,1) + V(:,2)*V(:,2) + V(:,3)*V(:,3))
    ! central difference term
    block
      real(8) rhou1, rhou2
      rhou1 = rho(1) * V(1,1)
      rhou2 = rho(2) * V(2,1)
      F(1)  = 0.5d0 * (rhou1 + rhou2)
      F(2)  = 0.5d0 * (rhou1 * V(1,1) + rhou2 * V(2,1) + (p(1) + p(2)))
      F(3)  = 0.5d0 * (rhou1 * V(1,2) + rhou2 * V(2,2))
      F(4)  = 0.5d0 * (rhou1 * V(1,3) + rhou2 * V(2,3))
      F(5)  = 0.5d0 * ((e(1) + p(1)) * V(1,1) + (e(2) + p(2)) * V(2,1))
    end block
    dQ(1) = -rho(1) + rho(2)
    dQ(2) = -rho(1) * V(1,1) + rho(2) * V(2,1)
    dQ(3) = -rho(1) * V(1,2) + rho(2) * V(2,2)
    dQ(4) = -rho(1) * V(1,3) + rho(2) * V(2,3)
    dQ(5) = -e(1) + e(2)
  end subroutine calc_Roe_common

  
  attributes(device) function Roe(rho, V, p, Normal) result(E)
    real(8), intent(in) :: rho(2), V(2,3), p(2)
    real(8), intent(in) :: Normal(5)
    real(8) E(5), dQ(5), uroe, vroe, wroe, Hroe, croe
    call calc_Roe_common(rho, V, p, dQ, E)
    call calc_Roe_ave(rho, V, p, uroe, vroe, wroe, Hroe, croe)
    call calc_RARinv_fast(uroe, vroe, wroe, Hroe, croe, dQ, E)
  end function Roe
end module calc_roe

