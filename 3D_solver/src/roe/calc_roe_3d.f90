  pure attributes(device) subroutine calc_RARinv(u, v, w, H, c, dQ, F1, F2, F3, F4, F5)
    real(8), intent(in), value :: u, v, w, H, c
    real(8), intent(in)        :: dQ(5)
    real(8), intent(inout)     :: F1, F2, F3, F4, F5
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
      F1 = F1 - 0.5d0 * (a1*R1                            +a2*R4      +a3*R5)
      F2 = F2 - 0.5d0 * (a1*(u-c)*R1                    +a2*u*R4+a3*(u+c)*R5)
      F3 = F3 - 0.5d0 * (a1*v*R1                 +a2c*R3+a2*v*R4    +a3*v*R5)
      F4 = F4 - 0.5d0 * (a1*w*R1          -a2c*R2     +a2*w*R4      +a3*w*R5)
      F5 = F5 - 0.5d0 * (a1*(H-cu)*R1+a2*(-w*c*R2+v*c*R3+q2*R4)+a3*(H+cu)*R5)
    end block
  end subroutine calc_RARinv


  pure attributes(device) subroutine calc_Roe_ave(rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, k1, k2, uroe, vroe, wroe, Hroe, croe)
    real(8), intent(in)  :: rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, k1, k2
    real(8), intent(out) :: uroe, vroe, wroe, Hroe, croe
    real(8) :: gamma_over_gamma_1 = gamma * over_gamma_1
    real(8) rhol, rhor, rhol_rhor, H1, H2
    H1 = gamma_over_gamma_1 * p1 / rho1 + k1
    H2 = gamma_over_gamma_1 * p2 / rho2 + k2
    rhol = sqrt(rho1)
    rhor = sqrt(rho2)
    rhol_rhor = 1.d0 / (rhol + rhor)
    uroe = (rhol * u1 + rhor * u2) * rhol_rhor
    vroe = (rhol * v1 + rhor * v2) * rhol_rhor
    wroe = (rhol * w1 + rhor * w2) * rhol_rhor
    Hroe = (rhol * H1 + rhor * H2) * rhol_rhor
    croe = sqrt(gamma_1 * (Hroe - 0.5d0 * (uroe*uroe + vroe*vroe + wroe*wroe)))
  end subroutine calc_Roe_ave

 
  pure attributes(device) subroutine calc_Roe_common(rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, k1, k2, dQ, F1, F2, F3, F4, F5)
    use mod_constant, only : over_gamma_1
    real(8), intent(in)  :: rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, k1, k2
    real(8), intent(out) :: dQ(5), F1, F2, F3, F4, F5
    real(8) e1, e2
    e1 = p1 * over_gamma_1 + rho1 * k1
    e2 = p2 * over_gamma_1 + rho2 * k2
    ! central difference term
    block
      real(8) rhou1, rhou2
      rhou1 = rho1 * u1
      rhou2 = rho2 * u2
      F1  = 0.5d0 * (rhou1 + rhou2)
      F2  = 0.5d0 * (rhou1 * u1 + rhou2 * u2 + (p1 + p2))
      F3  = 0.5d0 * (rhou1 * v1 + rhou2 * v2)
      F4  = 0.5d0 * (rhou1 * w1 + rhou2 * w2)
      F5  = 0.5d0 * ((e1 + p1) * u1 + (e2 + p2) * u2)
    end block
    dQ(1) = -rho1 + rho2
    dQ(2) = -rho1 * u1 + rho2 * u2
    dQ(3) = -rho1 * v1 + rho2 * v2
    dQ(4) = -rho1 * w1 + rho2 * w2
    dQ(5) = -e1 + e2
  end subroutine calc_Roe_common

  
  pure attributes(device) subroutine Roe(rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, F1, F2, F3, F4, F5)
    real(8), intent(in)  :: rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2
    real(8), intent(out) :: F1, F2, F3, F4, F5
    real(8) dQ(5), uroe, vroe, wroe, Hroe, croe, k1, k2
    k1 = 0.5d0 * (u1*u1 + v1*v1 + w1*w1)
    k2 = 0.5d0 * (u2*u2 + v2*v2 + w2*w2)
    call calc_Roe_common(rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, k1, k2, dQ, F1, F2, F3, F4, F5)
    call calc_Roe_ave(rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, k1, k2, uroe, vroe, wroe, Hroe, croe)
    call calc_RARinv(uroe, vroe, wroe, Hroe, croe, dQ, F1, F2, F3, F4, F5)
  end subroutine Roe

