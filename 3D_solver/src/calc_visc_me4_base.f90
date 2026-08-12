  !> Pure device function: 4th-order accurate flux reconstruction from 3-point stencil
  !> Uses compact central difference: F(i+1/2) = (-F_i + 26*F_{i+1/2} - F_{i+1})/24
  !> Achieves O(dx^4) accuracy with implicit stencil via dispersion relation optimization
  attributes(device) function flux4(a) result(ans)
    real(8), intent(in) :: a(3) !< 3-point array of flux values
    real(8) ans                 !< 4th-order flux result (-a1 + 26*a2 - a3) / 24
    !ans = (-a(1) + 26.d0 * a(2) - a(3)) * one_24
    ans = one_24 * fma(26.d0, a(2), -(a(1) + a(3)))
  end function flux4

  !> Pure device subroutine: Compute diagonal stress tensor components via 4th-order stencils
  !> Diagonal: t_ii = (2/3)*mu*(2*u_i,i - u_j,j - u_k,k) [with bulk viscosity correction]
  !> Uses 6-point stencil for strain rates and 3-point for viscosity averaging
  !> Computes work term ut_ii = u_i * t_ii needed for energy equation viscous contribution
  attributes(device) subroutine calc_tau_straight(mu, u, vy, wz, d, t11, ut11)
    real(8), intent(in), contiguous :: mu(3) !< viscosity at 3 stencil points
    real(8), intent(in), contiguous :: u(6)  !< velocity u at 6-point stencil
    real(8), intent(in), contiguous :: vy(6) !< dv/dy at 6-point stencil
    real(8), intent(in), contiguous :: wz(6) !< dw/dz at 6-point stencil
    real(8), intent(in)             :: d     !< inverse grid spacing (1/dx or 1/dy or 1/dz)
    real(8), intent(out)            :: t11   !< stress tensor component t_11
    real(8), intent(out)            :: ut11  !< work term u * t_11
    real(8) tmp1, tmp2, tmp3
    real(8) ax, svy, swz
    !tmp1 = two_third * mu(1) * ((2.25d0 * (-u(2) + u(3)) - (-u(1) + u(4)) * one_twelfth) * d &
    !       - 0.0625d0 * (-vy(1) + 9.d0 * (vy(2) + vy(3)) - vy(4)) &
    !       - 0.0625d0 * (-wz(1) + 9.d0 * (wz(2) + wz(3)) - wz(4)))
    ax   = fma(2.25d0, u(3) - u(2), -one_twelfth * (u(4) - u(1))) * d
    svy  = fma(9.d0, vy(2) + vy(3), -(vy(1) + vy(4)))
    swz  = fma(9.d0, wz(2) + wz(3), -(wz(1) + wz(4)))
    tmp1 = two_third * mu(1) * fma(-0.0625d0, svy + swz, ax)
    !tmp2 = two_third * mu(2) * ((2.25d0 * (-u(3) + u(4)) - (-u(2) + u(5)) * one_twelfth) * d &
    !       - 0.0625d0 * (-vy(2) + 9.d0 * (vy(3) + vy(4)) - vy(5)) &
    !       - 0.0625d0 * (-wz(2) + 9.d0 * (wz(3) + wz(4)) - wz(5)))
    ax   = fma(2.25d0, u(4) - u(3), -one_twelfth * (u(5) - u(2))) * d
    svy  = fma(9.d0, vy(3) + vy(4), -(vy(2) + vy(5)))
    swz  = fma(9.d0, wz(3) + wz(4), -(wz(2) + wz(5)))
    tmp2 = two_third * mu(2) * fma(-0.0625d0, svy + swz, ax)
    !tmp3 = two_third * mu(3) * ((2.25d0 * (-u(4) + u(5)) - (-u(3) + u(6)) * one_twelfth) * d &
    !       - 0.0625d0 * (-vy(3) + 9.d0 * (vy(4) + vy(5)) - vy(6)) &
    !       - 0.0625d0 * (-wz(3) + 9.d0 * (wz(4) + wz(5)) - wz(6)))
    ax   = fma(2.25d0, u(5) - u(4), -one_twelfth * (u(6) - u(3))) * d
    svy  = fma(9.d0, vy(4) + vy(5), -(vy(3) + vy(6)))
    swz  = fma(9.d0, wz(4) + wz(5), -(wz(3) + wz(6)))
    tmp3 = two_third * mu(3) * fma(-0.0625d0, svy + swz, ax)
    !t11  = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
    t11  = one_24 * fma(26.d0, tmp2, -(tmp1 + tmp3))
    !tmp1 = 0.0625d0 * (-u(1) + 9.d0 * (u(2) + u(3)) - u(4)) * tmp1
    !tmp2 = 0.0625d0 * (-u(2) + 9.d0 * (u(3) + u(4)) - u(5)) * tmp2
    !tmp3 = 0.0625d0 * (-u(3) + 9.d0 * (u(4) + u(5)) - u(6)) * tmp3
    tmp1 = 0.0625d0 * fma(9.d0, u(2) + u(3), -(u(1) + u(4))) * tmp1
    tmp2 = 0.0625d0 * fma(9.d0, u(3) + u(4), -(u(2) + u(5))) * tmp2
    tmp3 = 0.0625d0 * fma(9.d0, u(4) + u(5), -(u(3) + u(6))) * tmp3
    !ut11 = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
    ut11 = one_24 * fma(26.d0, tmp2, -(tmp1 + tmp3))
  end subroutine calc_tau_straight


  !> Pure device subroutine: Compute diagonal stress with Smagorinsky LES turbulent viscosity
  !> Combines molecular + turbulent (SGS) viscosity: nu_total = nu + nu_t
  !> Turbulent part nu_t captures unresolved subgrid energy dissipation
  attributes(device) subroutine calc_tau_straight_LES(mu, mut, u, vy, wz, d, t11, ut11)
    real(8), intent(in), contiguous :: mu(3)  !< molecular viscosity at 3 stencil points
    real(8), intent(in), contiguous :: mut(3) !< turbulent viscosity at 3 stencil points
    real(8), intent(in), contiguous :: u(6)   !< velocity u at 6-point stencil
    real(8), intent(in), contiguous :: vy(6)  !< dv/dy at 6-point stencil
    real(8), intent(in), contiguous :: wz(6)  !< dw/dz at 6-point stencil
    real(8), intent(in)             :: d      !< inverse grid spacing
    real(8), intent(out)            :: t11    !< total stress (molecular + SGS)
    real(8), intent(out)            :: ut11   !< work term u * t_11
    real(8) tmp1(3), tmp2(3)
    real(8) ax, svy, swz
    !tmp1(:) = (2.25d0 * (-u(2:4) + u(3:5)) - (-u(1:3) + u(4:6)) * one_twelfth) * d
    !tmp1(:) = tmp1(:) - 0.0625d0 * (-vy(1:3) + 9.d0 * (vy(2:4) + vy(3:5)) - vy(4:6))
    !tmp1(:) = tmp1(:) - 0.0625d0 * (-wz(1:3) + 9.d0 * (wz(2:4) + wz(3:5)) - wz(4:6))
    ax      = fma(2.25d0, u(3) - u(2), -one_twelfth * (u(4) - u(1))) * d
    svy     = fma(9.d0, vy(2) + vy(3), -(vy(1) + vy(4)))
    swz     = fma(9.d0, wz(2) + wz(3), -(wz(1) + wz(4)))
    tmp1(1) = fma(-0.0625d0, svy + swz, ax)
    ax      = fma(2.25d0, u(4) - u(3), -one_twelfth * (u(5) - u(2))) * d
    svy     = fma(9.d0, vy(3) + vy(4), -(vy(2) + vy(5)))
    swz     = fma(9.d0, wz(3) + wz(4), -(wz(2) + wz(5)))
    tmp1(2) = fma(-0.0625d0, svy + swz, ax)
    ax      = fma(2.25d0, u(5) - u(4), -one_twelfth * (u(6) - u(3))) * d
    svy     = fma(9.d0, vy(4) + vy(5), -(vy(3) + vy(6)))
    swz     = fma(9.d0, wz(4) + wz(5), -(wz(3) + wz(6)))
    tmp1(3) = fma(-0.0625d0, svy + swz, ax)
    tmp2(:) = two_third * mu(:) * tmp1(:)
    !tmp2(1) = 0.0625d0 * (-u(1) + 9.d0 * (u(2) + u(3)) - u(4)) * tmp2(1)
    !tmp2(2) = 0.0625d0 * (-u(2) + 9.d0 * (u(3) + u(4)) - u(5)) * tmp2(2)
    !tmp2(3) = 0.0625d0 * (-u(3) + 9.d0 * (u(4) + u(5)) - u(6)) * tmp2(3)
    tmp2(1) = 0.0625d0 * fma(9.d0, u(2) + u(3), -(u(1) + u(4))) * tmp2(1)
    tmp2(2) = 0.0625d0 * fma(9.d0, u(3) + u(4), -(u(2) + u(5))) * tmp2(2)
    tmp2(3) = 0.0625d0 * fma(9.d0, u(4) + u(5), -(u(3) + u(6))) * tmp2(3)
    !ut11    = (-tmp2(1) + 26.d0 * tmp2(2) - tmp2(3)) * one_24
    ut11    = one_24 * fma(26.d0, tmp2(2), -(tmp2(1) + tmp2(3)))
    tmp2(:) = two_third * (mu(:) + mut(:)) * tmp1(:)
    !t11     = (-tmp2(1) + 26.d0 * tmp2(2) - tmp2(3)) * one_24
    t11     = one_24 * fma(26.d0, tmp2(2), -(tmp2(1) + tmp2(3)))
  end subroutine calc_tau_straight_LES


  !> Pure device subroutine: Compute shear (off-diagonal) stress tensor components
  !> Shear: t_ij = mu*(u_i,j + u_j,i) for i != j components
  !> 4th-order stencil preserves cross-derivatives symmetry (t_12 = t_21)
  attributes(device) subroutine calc_tau_cross(mu, v, uy, d, t12, vt12)
    real(8), intent(in), contiguous :: mu(3) !< viscosity at 3 stencil points
    real(8), intent(in), contiguous :: v(6)  !< velocity v at 6-point stencil
    real(8), intent(in), contiguous :: uy(6) !< du/dy at 6-point stencil
    real(8), intent(in)             :: d      !< inverse grid spacing
    real(8), intent(out)            :: t12    !< shear stress component t_12
    real(8), intent(out)            :: vt12   !< work term v * t_12
    real(8) tmp1, tmp2, tmp3
    real(8) av, suy
    !tmp1 = mu(1) * ((1.125d0 * (-v(2) + v(3)) - (-v(1) + v(4)) * one_24) * d &
    !                + 0.0625d0 * (-uy(1) + 9.d0 * (uy(2) + uy(3)) - uy(4)))
    !tmp2 = mu(2) * ((1.125d0 * (-v(3) + v(4)) - (-v(2) + v(5)) * one_24) * d &
    !                + 0.0625d0 * (-uy(2) + 9.d0 * (uy(3) + uy(4)) - uy(5)))
    !tmp3 = mu(3) * ((1.125d0 * (-v(4) + v(5)) - (-v(3) + v(6)) * one_24) * d &
    !                + 0.0625d0 * (-uy(3) + 9.d0 * (uy(4) + uy(5)) - uy(6)))
    av   = fma(1.125d0, v(3) - v(2), -one_24 * (v(4) - v(1))) * d
    suy  = fma(9.d0, uy(2) + uy(3), -(uy(1) + uy(4)))
    tmp1 = mu(1) * fma(0.0625d0, suy, av)
    av   = fma(1.125d0, v(4) - v(3), -one_24 * (v(5) - v(2))) * d
    suy  = fma(9.d0, uy(3) + uy(4), -(uy(2) + uy(5)))
    tmp2 = mu(2) * fma(0.0625d0, suy, av)
    av   = fma(1.125d0, v(5) - v(4), -one_24 * (v(6) - v(3))) * d
    suy  = fma(9.d0, uy(4) + uy(5), -(uy(3) + uy(6)))
    tmp3 = mu(3) * fma(0.0625d0, suy, av)
    !t12  = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
    t12  = one_24 * fma(26.d0, tmp2, -(tmp1 + tmp3))
    !tmp1 = 0.0625d0 * (-v(1) + 9.d0 * (v(2) + v(3)) - v(4)) * tmp1
    !tmp2 = 0.0625d0 * (-v(2) + 9.d0 * (v(3) + v(4)) - v(5)) * tmp2
    !tmp3 = 0.0625d0 * (-v(3) + 9.d0 * (v(4) + v(5)) - v(6)) * tmp3
    tmp1 = 0.0625d0 * fma(9.d0, v(2) + v(3), -(v(1) + v(4))) * tmp1
    tmp2 = 0.0625d0 * fma(9.d0, v(3) + v(4), -(v(2) + v(5))) * tmp2
    tmp3 = 0.0625d0 * fma(9.d0, v(4) + v(5), -(v(3) + v(6))) * tmp3
    !vt12 = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
    vt12 = one_24 * fma(26.d0, tmp2, -(tmp1 + tmp3))
  end subroutine calc_tau_cross


  !> Pure device subroutine: Compute shear stress with Smagorinsky LES turbulent model
  !> Off-diagonal components including both molecular and subgrid turbulent dissipation
  attributes(device) subroutine calc_tau_cross_LES(mu, mut, v, uy, d, t12, vt12)
    real(8), intent(in), contiguous :: mu(3)  !< molecular viscosity at 3 stencil points    
    real(8), intent(in), contiguous :: mut(3) !< turbulent viscosity at 3 stencil points
    real(8), intent(in), contiguous :: v(6)   !< velocity v at 6-point stencil
    real(8), intent(in), contiguous :: uy(6)  !< du/dy at 6-point stencil
    real(8), intent(in)             :: d      !< inverse grid spacing
    real(8), intent(out)            :: t12    !< total shear stress (molecular + SGS)
    real(8), intent(out)            :: vt12   !< work term v * t_12
    real(8) tmp1(3), tmp2(3)
    real(8) av, suy
    !tmp1(:) = (1.125d0 * (-v(2:4) + v(3:5)) - (-v(1:3) + v(4:6)) * one_24) * d
    !tmp1(:) = tmp1(:) + 0.0625d0 * (-uy(1:3) + 9.d0 * (uy(2:4) + uy(3:5)) - uy(4:6))
    av      = fma(1.125d0, v(3) - v(2), -one_24 * (v(4) - v(1))) * d
    suy     = fma(9.d0, uy(2) + uy(3), -(uy(1) + uy(4)))
    tmp1(1) = fma(0.0625d0, suy, av)
    av      = fma(1.125d0, v(4) - v(3), -one_24 * (v(5) - v(2))) * d
    suy     = fma(9.d0, uy(3) + uy(4), -(uy(2) + uy(5)))
    tmp1(2) = fma(0.0625d0, suy, av)
    av      = fma(1.125d0, v(5) - v(4), -one_24 * (v(6) - v(3))) * d
    suy     = fma(9.d0, uy(4) + uy(5), -(uy(3) + uy(6)))
    tmp1(3) = fma(0.0625d0, suy, av)
    tmp2(:) = mu(:) * tmp1(:)
    !tmp2(1) = 0.0625d0 * (-v(1) + 9.d0 * (v(2) + v(3)) - v(4)) * tmp2(1)
    !tmp2(2) = 0.0625d0 * (-v(2) + 9.d0 * (v(3) + v(4)) - v(5)) * tmp2(2)
    !tmp2(3) = 0.0625d0 * (-v(3) + 9.d0 * (v(4) + v(5)) - v(6)) * tmp2(3)
    tmp2(1) = 0.0625d0 * fma(9.d0, v(2) + v(3), -(v(1) + v(4))) * tmp2(1)
    tmp2(2) = 0.0625d0 * fma(9.d0, v(3) + v(4), -(v(2) + v(5))) * tmp2(2)
    tmp2(3) = 0.0625d0 * fma(9.d0, v(4) + v(5), -(v(3) + v(6))) * tmp2(3)
    !vt12    = (-tmp2(1) + 26.d0 * tmp2(2) - tmp2(3)) * one_24
    vt12    = one_24 * fma(26.d0, tmp2(2), -(tmp2(1) + tmp2(3)))
    tmp2(:) = (mu(:) + mut(:)) * tmp1(:)
    !t12     = (-tmp2(1) + 26.d0 * tmp2(2) - tmp2(3)) * one_24
    t12     = one_24 * fma(26.d0, tmp2(2), -(tmp2(1) + tmp2(3)))
  end subroutine calc_tau_cross_LES
