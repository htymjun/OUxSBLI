module calc_weno
  use mod_constant, only : one_sixth
  implicit none
contains
  !$dir inline
  attributes(device) function weno5z_left(v) result(vf)
    real(8), intent(in) :: v(5)
    real(8) :: vf, a0, a1, a2, s
    real(8), parameter :: eps = 1.0d-8, d0 = 0.1d0, d1 = 0.6d0, d2 = 0.3d0
    real(8), parameter :: thirteen_twelfth = 13.d0 / 12.d0
    block
      real(8) b0, b1, b2
      b0 = thirteen_twelfth * (v(1) - 2.d0 * v(2) + v(3))**2 + 0.25d0 * (v(1) - 4.d0 * v(2) + 3.d0 * v(3))**2
      b1 = thirteen_twelfth * (v(2) - 2.d0 * v(3) + v(4))**2 + 0.25d0 * (v(2) - v(4))**2
      b2 = thirteen_twelfth * (v(3) - 2.d0 * v(4) + v(5))**2 + 0.25d0 * (3.d0 * v(3) - 4.d0 * v(4) + v(5))**2
      block
        real(8) tau5
        tau5 = abs(b0 - b2)
        a0 = d0 * (1.d0 + (tau5 / (b0 + eps))**2)
        a1 = d1 * (1.d0 + (tau5 / (b1 + eps))**2)
        a2 = d2 * (1.d0 + (tau5 / (b2 + eps))**2)
        s  = a0 + a1 + a2
      end block
    end block
    block
      real(8) p0, p1, p2
      p0 = ( 2.d0 * v(1) - 7.d0 * v(2) + 11.d0 * v(3)) * one_sixth
      p1 = (-1.d0 * v(2) + 5.d0 * v(3) +  2.d0 * v(4)) * one_sixth
      p2 = ( 2.d0 * v(3) + 5.d0 * v(4) -  1.d0 * v(5)) * one_sixth
      vf = (a0 * p0 + a1 * p1 + a2 * p2) / s
    end block
  end function weno5z_left


  !$dir inline
  attributes(device) function weno5z_right(v) result(vf)
    real(8), intent(in) :: v(5)
    real(8) :: vf, a0, a1, a2, s
    real(8), parameter :: eps = 1.0d-8, d0 = 0.1d0, d1 = 0.6d0, d2 = 0.3d0
    real(8), parameter :: thirteen_twelfth = 13.d0 / 12.d0
    block
      real(8) b0, b1, b2
      b0 = thirteen_twelfth * (v(1) - 2.d0 * v(2) + v(3))**2 + 0.25d0 * (v(1) - 4.d0 * v(2) + 3.d0 * v(3))**2
      b1 = thirteen_twelfth * (v(2) - 2.d0 * v(3) + v(4))**2 + 0.25d0 * (v(2) - v(4))**2
      b2 = thirteen_twelfth * (v(3) - 2.d0 * v(4) + v(5))**2 + 0.25d0 * (3.d0 * v(3) - 4.d0 * v(4) + v(5))**2
      block
        real(8) tau5
        tau5 = abs(b0 - b2)
        a0 = d0 * (1.d0 + (tau5 / (b0 + eps))**2)
        a1 = d1 * (1.d0 + (tau5 / (b1 + eps))**2)
        a2 = d2 * (1.d0 + (tau5 / (b2 + eps))**2)
        s  = a0 + a1 + a2
      end block
    end block
    block
      real(8) p0, p1, p2
      p0 = (-1.d0 * v(1) + 5.d0 * v(2) +  2.d0 * v(3)) * one_sixth
      p1 = ( 2.d0 * v(2) + 5.d0 * v(3) -  1.d0 * v(4)) * one_sixth
      p2 = (11.d0 * v(3) - 7.d0 * v(4) +  2.d0 * v(5)) * one_sixth
      vf = (a0 * p0 + a1 * p1 + a2 * p2) / s
    end block
  end function weno5z_right
  

  !$dir inline
  attributes(device) subroutine calc_Reo_ave(rho, u, v, w, p, uroe, vroe, wroe, Hroe, croe)
    real(8), intent(in), dimension(2) :: rho, u, v, w, p
    real(8), intent(out)              :: uroe, vroe, wroe, Hroe, croe
    real(8) rhol, rhor, rhol_rhor, H(2)
    H(:) =
    rhol = sqrt(rho(1))
    rhor = sqrt(rho(2))
    rhol_rhor = 1.d0 / (rhol + rhor)
    uroe = (rhol * u(1) + rhor * u(2)) * rhol_rhor
    vroe = (rhol * v(1) + rhor * v(2)) * rhol_rhor
    wroe = (rhol * w(1) + rhor * w(2)) * rhol_rhor
    Hroe = (rhol * H(1) + rhor * H(2)) * rhol_rhor
    croe = sqrt(gamma_1 * (Hroe - 0.5d0 * (uroe**2 + vroe**2 + wroe**2)))
  end subroutine calc_Roe_ave

  
  !$dir inline
  attributes(device) subroutine calc_eigen_matrices(uroe, vroe, wroe, Hroe, croe, Lm, Rm, lam)
  end subroutine calc_eigen_matrices


  !$dir inline
  attributes(device) function weno5_z(f) result(f2)
    real(8), intent(in) :: f(6)
    real(8) f2(2)
    f2(1) = weno5z_left(f(1:5))
    f2(2) = weno5z_right(f(2:6))
  end function weno5_z

  
  !dir inline
  attributes(device) function weno3(f) result(f2)
    use mod_constant, only : one_third
    real(8), intent(in) :: f(4)
    real(8) f2(2), d(3)
    d(:)  = -f(1:3) + f(2:4)
    f2(1) = f(2) + 0.5d0 * (d(1) + 2.d0 * d(2)) * one_third ! 3rd non-TVD MUSCL
    f2(2) = f(3) - 0.5d0 * (d(3) + 2.d0 * d(2)) * one_third ! 3rd non-TVD MUSCL
  end function weno3


  attributes(device) subroutine weno_4points(rho, u, v, w, p, rho2, p2, V2)
    real(8), intent(in)          :: rho(4), u(4), v(4), w(4), p(4)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    rho2(:) = weno3(rho)
    V2(:,1) = weno3(u)
    V2(:,2) = weno3(v)
    V2(:,3) = weno3(w)
    p2(:)   = weno3(p)
  end subroutine weno_4points


  attributes(device) subroutine weno_6points(rho, u, v, w, p, rho2, p2, V2)
    real(8), intent(in)          :: rho(6), u(6), v(6), w(6), p(6)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    rho2(:) = weno5_z(rho)
    V2(:,1) = weno5_z(u)
    V2(:,2) = weno5_z(v)
    V2(:,3) = weno5_z(w)
    p2(:)   = weno5_z(p)
  end subroutine weno_6points
end module calc_weno

