module calc_slau
  use mod_globals, only : dim => dimension, gamma
  use mod_constant, only : over_gamma_1
  implicit none
  interface SLAU
    module procedure SLAU1, HRSLAU2
  end interface SLAU
contains
  !$dir inline
  attributes(device) function vecsum(id, V) result(ans)
    integer, intent(in), value :: id
    real(8), intent(in)        :: V(2,dim)
    real(8) ans
    integer i
    ans = 0.d0
    do i = 1, dim
      ans = ans + V(id,i)**2
    enddo
  end function vecsum


  !$dir inline
  attributes(device) function q2(V) result(ans)
    real(8), intent(in) :: V(2,dim)
    real(8) ans
    integer i, j
    ans = 0.d0
    do j = 1, dim
      do i = 1, 2
        ans = ans + V(i,j)**2
    enddo;enddo
  end function q2


  !$dir inline
  attributes(device) subroutine SLAU_common(id, rho, p, V, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    integer, intent(in), value :: id
    real(8), intent(in)        :: rho(2), p(2)
    real(8), intent(in)        :: V(2,dim)
    real(8), intent(out)       :: c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm
    block
      real(8) cl, cr
      cl = sqrt(gamma * p(1) / rho(1))
      cr = sqrt(gamma * p(2) / rho(2))
      c  = 0.5d0 * (cl + cr)
    end block
    over_c = 1.d0 / c
    Mp  = V(1,id) * over_c
    Mm  = V(2,id) * over_c
    block
      real(8) g, one_g_Vt
      g   = -max(min(Mp, 0.d0), -1.d0) * min(max(Mm, 0.d0), 1.d0)
      one_g_Vt = (1.d0 - g) * (rho(1) * abs(V(1,id)) + rho(2) * abs(V(2,id))) / (rho(1) + rho(2))
      Vtp = one_g_Vt + g * abs(V(1,id))
      Vtm = one_g_Vt + g * abs(V(2,id))
    end block
    bp = merge(0.25d0 * (2.d0 - Mp) * (Mp + 1.d0)**2, &
               0.5d0 * (1.d0 + sign(1.d0, Mp)), &
               abs(Mp) < 1.d0)
    bm = merge(0.25d0 * (2.d0 + Mm) * (Mm - 1.d0)**2, &
               0.5d0 * (1.d0 + sign(1.d0, -Mm)), &
               abs(Mm) < 1.d0)
    dp = -p(1) + p(2)
  end subroutine SLAU_common


  !$dir inline
  attributes(device) function phil(rho, V, p) result(ans)
    real(8), intent(in) :: rho(2), V(2,dim), p(2)
    real(8) ans(dim+2)
    ans(1)       = 1.d0
    ans(2:dim+1) = V(1,:)
    ans(dim+2)   = (p(1) * over_gamma_1 + 0.5d0 * rho(1) * vecsum(1, V(:,:)) + p(1)) / rho(1)
  end function phil


  !$dir inline
  attributes(device) function phir(rho, V, p) result(ans)
    real(8), intent(in) :: rho(2), V(2,dim), p(2)
    real(8) ans(dim+2)
    ans(1)       = 1.d0
    ans(2:dim+1) = V(2,:)
    ans(dim+2)   = (p(2) * over_gamma_1 + 0.5d0 * rho(2) * vecsum(2, V(:,:)) + p(2)) / rho(2)
  end function phir


  attributes(device) function SLAU1(id_slau, id, rho, p, V, Norm, HR) result(F)
    integer(kind=2), intent(in), value :: id_slau
    integer, intent(in), value         :: id
    real(8), intent(in)                :: rho(2), p(2)
    real(8), intent(in)                :: V(2,dim)
    real(8), intent(in)                :: Norm(dim+2)
    real(8), intent(in), value         :: HR
    real(8) c, over_c, Mp, Mm, M, x
    real(8) Vtp, Vtm, dp, bp, bm, mass
    real(8), dimension(dim+2) :: F
    call SLAU_common(id, rho, p, V, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    M = min(1.d0, sqrt(0.5d0 * q2(V(:,:))) * over_c)
    x = (1.d0 - M) ** 2
    mass = 0.5d0 * (rho(1) * (V(1,id) + Vtp) + rho(2) * (V(2,id) - Vtm) - x * dp * over_c)
    block
      real(8) pres
      pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + (1.d0 - x) * (bp + bm - 1.d0) * (p(1) + p(2)))
      F(:) = 0.5d0 * ((mass + abs(mass)) * phil(rho, V, p) + (mass - abs(mass)) * phir(rho, V, p)) + pres * Norm(:)
    end block
  end function SLAU1


  attributes(device) function HRSLAU2(id_slau, id, rho, p, V, Norm, HR) result(F)
    integer(kind=4), intent(in), value :: id_slau
    integer, intent(in), value         :: id
    real(8), intent(in)                :: rho(2), p(2)
    real(8), intent(in)                :: V(2,dim)
    real(8), intent(in)                :: Norm(dim+2)
    real(8), intent(in), value         :: HR
    real(8) c, over_c, Mp, Mm
    real(8) Vtp, Vtm, dp, bp, bm, mass, V2
    real(8), dimension(dim+2) :: F
    call SLAU_common(id, rho, p, V, c, over_c, Mp, Mm, bp, bm, dp, Vtp, Vtm)
    V2 = sqrt(0.5d0 * q2(V(:,:)))
    block
      real(8) M, x
      M  = min(1.d0, V2 * over_c)
      x  = (1.d0 - M) ** 2
      mass = 0.5d0 * (rho(1) * (V(1,id) + Vtp) + rho(2) * (V(2,id) - Vtm) - x * dp * over_c)
    end block
    block
      real(8) pres
      pres = 0.5d0 * (p(1) + p(2) + (bp - bm) * (-dp) + HR * V2 * (bp + bm - 1.d0) * 0.5d0 * (rho(1) + rho(2)) * c)
      F(:) = 0.5d0 * ((mass + abs(mass)) * phil(rho, V, p) + (mass - abs(mass)) * phir(rho, V, p)) + pres * Norm(:)
    end block
  end function HRSLAU2
end module calc_slau

