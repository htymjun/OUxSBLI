  attributes(device) subroutine SLAU_common(rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, &
                                                 un1, un2, p1, p2, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    real(8), intent(in)  :: rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, un1, un2, p1, p2
    real(8), intent(out) :: c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm
    block
      real(8) cl, cr
      cl = sqrt(gamma * p1 * over_rho1)
      cr = sqrt(gamma * p2 * over_rho2)
      c  = 0.5d0 * (cl + cr)
    end block
    over_c = 1.d0 / c
    Mp  = un1 * over_c
    Mm  = un2 * over_c
    block
      real(8) g, one_g_Vt
      g   = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
      one_g_Vt = (1.d0 - g) * (rho1 * abs(un1) + rho2 * abs(un2)) / (rho1 + rho2)
      Vtp = one_g_Vt + g * abs(un1)
      Vtm = one_g_Vt + g * abs(un2)
    end block
    bp = merge(0.25d0 * (2.d0 - Mp) * (Mp + 1.d0)**2, &
               0.5d0 * (1.d0 + sign(1.d0, Mp)), &
               abs(Mp) < 1.d0)
    bm = merge(0.25d0 * (2.d0 + Mm) * (Mm - 1.d0)**2, &
               0.5d0 * (1.d0 + sign(1.d0, -Mm)), &
               abs(Mm) < 1.d0)
    dp = -p1 + p2
  end subroutine SLAU_common


  pure attributes(device) function phi(rho, k, p, over_rho) result(ans)
    use mod_globals, only : gamma
    use mod_constant, only : over_gamma_1
    real(8), intent(in) :: rho, k, p, over_rho
    real(8) :: gamma_over_gamma_1 = gamma * over_gamma_1
    real(8) ans
    ans = (p * gamma_over_gamma_1 + rho * k) * over_rho
  end function phi


  attributes(device) subroutine SLAU1(id_slau, rho1, rho2, u1, u2, v1, v2, &
                                           un1, un2, p1, p2, Norm, HR, F1, F2, F3, F4)
    integer(2), intent(in), value :: id_slau
    real(8), intent(in)           :: rho1, rho2, u1, u2, v1, v2, un1, un2, p1, p2, Norm(5)
    real(sp), intent(in), value   :: HR
    real(8), intent(out)          :: F1, F2, F3, F4
    real(8) c, over_c, Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass, mass1, mass2, over_rho1, over_rho2, k1, k2
    over_rho1 = 1.d0 / rho1
    over_rho2 = 1.d0 / rho2
    call SLAU_common(rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, un1, un2, &
                     p1, p2, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    k1 = 0.5d0 * (u1*u1 + v1*v1)
    k2 = 0.5d0 * (u2*u2 + v2*v2)
    M  = fmin(1.d0, sqrt(k1 + k2) * over_c)
    x  = (1.d0 - M) ** 2
    !mass = 0.25d0 * (rho1 * (un1 + Vtp) + rho2 * (un2 - Vtm) - x * dp * over_c)
    mass  = 0.25d0 * fma(rho1, un1 + Vtp, fma(rho2, un2 - Vtm, -x * dp * over_c))
    mass1 = mass + abs(mass)
    mass2 = mass - abs(mass)
    block
      ! pressure term
      real(8) pres
      !pres = 0.5d0 * ((p1 + p2) + (bp - bm) * (-dp) + (1.d0 - x) * (bp + bm - 1.d0) * (p1 + p2))
      !pres = 0.5d0 * ((p1 + p2) * ((1.d0 - x) * (bp + bm - 1.d0) + 1.d0) + (bm - bp) * dp)
      block
        real(8) p_sum, b_sum_m1, b_diff, term, term_p
        ! independent
        p_sum    = p1 + p2
        b_sum_m1 = bp + bm - 1.d0
        b_diff   = bm - bp
        ! FMA
        term     = fma(1.d0 - x, b_sum_m1, 1.d0)   ! (1.d0 - x) * b_sum_m1 + 1.d0
        term_p   = term * p_sum
        pres     = 0.5d0 * fma(dp, b_diff, term_p) ! dp * b_diff + term_p
      end block
      F1 = mass1       + mass2
      !F2 = mass1 * u1 + mass2 * u2 + pres * Norm(2)
      !F3 = mass1 * v1 + mass2 * v2 + pres * Norm(3)
      F2 = fma(mass1, u1, fma(mass2, u2, pres * Norm(2)))
      F3 = fma(mass1, v1, fma(mass2, v2, pres * Norm(3)))
      F4 = mass1 * phi(rho1, k1, p1, over_rho1) + mass2 * phi(rho2, k2, p2, over_rho2)
    end block
  end subroutine SLAU1


  attributes(device) subroutine HRSLAU2(id_slau, rho1, rho2, u1, u2, v1, v2, &
                                             un1, un2, p1, p2, Norm, HR, F1, F2, F3, F4)
    integer(4), intent(in), value :: id_slau
    real(8), intent(in)           :: rho1, rho2, u1, u2, v1, v2, un1, un2, p1, p2, Norm(5)
    real(sp), intent(in), value   :: HR
    real(8), intent(out)          :: F1, F2, F3, F4
    real(8) c, over_c, Mp, Mm
    real(8) Vtp, Vtm, dp, bp, bm, mass, mass1, mass2, Vec2, over_rho1, over_rho2, k1, k2
    over_rho1 = 1.d0 / rho1
    over_rho2 = 1.d0 / rho2
    call SLAU_common(rho1, rho2, over_rho1, over_rho2, u1, u2, v1, v2, un1, un2, &
                     p1, p2, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    k1   = 0.5d0 * (u1*u1 + v1*v1)
    k2   = 0.5d0 * (u2*u2 + v2*v2)
    Vec2 = sqrt(k1 + k2)
    block
      real(8) M, x
      M  = fmin(1.d0, Vec2 * over_c)
      x  = (1.d0 - M) ** 2
      !mass = 0.25d0 * (rho1 * (un1 + Vtp) + rho2 * (un2 - Vtm) - x * dp * over_c)
      mass  = 0.25d0 * fma(rho1, (un1 + Vtp), fma(rho2, (un2 - Vtm), -x * dp * over_c))
      mass1 = mass + abs(mass)
      mass2 = mass - abs(mass)
    end block
    block
      real(8) pres
      !pres = 0.5d0 * (p1 + p2 + (bp - bm) * (-dp) + HR * Vec2 * (bp + bm - 1.d0) * 0.5d0 * (rho1 + rho2) * c)
      !pres = 0.5d0 * (p1 + p2 + (bm - bp) * dp + HR * Vec2 * 0.5d0 * c * (bp + bm - 1.d0) * (rho1 + rho2))
      block
        real(8) p_sum, b_sum_m1, b_diff, rho_sum, coeff, term
        ! independent
        p_sum    = p1 + p2
        b_sum_m1 = bp + bm - 1.d0
        b_diff   = bm - bp
        rho_sum  = rho1 + rho2
        coeff    = HR * Vec2 * 0.5d0 * c * b_sum_m1
        ! FMA
        term     = fma(dp, b_diff, p_sum)            ! dp * b_diff + p_sum
        pres     = 0.5d0 * fma(coeff, rho_sum, term) ! 0.5d0 * (coeff * rho_sum + term)
      end block
      F1 = mass1       + mass2
      !F2 = mass1 * u1 + mass2 * u2 + pres * Norm(2)
      !F3 = mass1 * v1 + mass2 * v2 + pres * Norm(3)
      F2 = fma(mass1, u1, fma(mass2, u2, pres * Norm(2)))
      F3 = fma(mass1, v1, fma(mass2, v2, pres * Norm(3)))
      F4 = mass1 * phi(rho1, k1, p1, over_rho1) + mass2 * phi(rho2, k2, p2, over_rho2)
    end block
  end subroutine HRSLAU2
  