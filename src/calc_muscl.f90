module calc_muscl
  use libm
  use mod_globals, only : sp
  use mod_constant, only : one_third, one_sixth, one_twelfth
  implicit none
  private
  public delta4, delta6
  interface minmod
    module procedure minmod2, minmod3
  end interface

  interface MUSCL3rd
    module procedure MUSCL3rdnonTVD, MUSCL3rdMinmod, MUSCL3rdThreshold
  end interface
  
  interface MUSCL4th
    module procedure MUSCL4thnonTVD, MUSCL4thTVD, MUSCL4thThreshold
  end interface
contains
  attributes(device) function minmod2(x, y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = copysign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod2


  attributes(device) function minmod3(x, y, z) result(ans)
    real(8), intent(in), value :: x, y, z
    real(8) :: ans, sgn
    sgn = copysign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y, sgn * z), 0.d0)
  end function minmod3


  attributes(device) function d33(d1, d2, d3) result(ans)
    real(8), intent(in), value :: d1, d2, d3 
    real(8) :: ans
    ans =          minmod(d1, 2.d0 * d2, 2.d0 * d3) &
          - 2.d0 * minmod(d2, 2.d0 * d1, 2.d0 * d3) &
                 + minmod(d3, 2.d0 * d1, 2.d0 * d2)
  end function d33


  attributes(device) subroutine MUSCL3rdnonTVD(id_tvd, sensor, a2, a3, d1, d2, d3, al, ar)
    integer(kind=2), intent(in) :: id_tvd
    real(sp), intent(in)        :: sensor
    real(8), intent(in)         :: a2, a3, d1, d2, d3
    real(8), intent(out)        :: al, ar
    al = a2 + fma(2.d0, d2, d1) * one_sixth
    ar = a3 - fma(2.d0, d2, d3) * one_sixth
  end subroutine MUSCL3rdnonTVD


  attributes(device) subroutine MUSCL3rdMinmod(id_tvd, sensor, a2, a3, d1, d2, d3, al, ar)
    integer(kind=4), intent(in) :: id_tvd
    real(sp), intent(in)        :: sensor
    real(8), intent(in)         :: a2, a3, d1, d2, d3
    real(8), intent(out)        :: al, ar
    real(8), constant :: b = (3.d0 - one_third) / (1.d0 - one_third)
    block
      real(8) dt1, dt2
      dt1 = minmod(d1, b * d2)
      dt2 = minmod(d2, b * d1)
      al  = a2 + fma(2.d0, dt2, dt1) * one_sixth
    end block
    block
      real(8) dt3, dt4
      dt3 = minmod(d3, b * d2)
      dt4 = minmod(d2, b * d3)
      ar  = a3 - fma(2.d0, dt4, dt3) * one_sixth
    end block
  end subroutine MUSCL3rdMinmod


  attributes(device) subroutine MUSCL3rdThreshold(id_tvd, sensor, a2, a3, d1, d2, d3, al, ar)
    use mod_globals, only : threshold
    integer(kind=8), intent(in) :: id_tvd
    real(sp), intent(in)        :: sensor
    real(8), intent(in)         :: a2, a3, d1, d2, d3
    real(8), intent(out)        :: al, ar
    integer(kind=2) :: id2 = 0_2
    integer(kind=4) :: id4 = 0_4
    if (sensor < threshold) then
      call MUSCL3rdnonTVD(id2, sensor, a2, a3, d1, d2, d3, al, ar)
    else
      call MUSCL3rdMinmod(id4, sensor, a2, a3, d1, d2, d3, al, ar)
    endif
  end subroutine MUSCL3rdThreshold


  attributes(device) subroutine MUSCL4thnonTVD(id_tvd, sensor, a2, a3, d1, d2, d3, d4, d5, al, ar)
    integer(kind=2), intent(in) :: id_tvd
    real(sp), intent(in)        :: sensor
    real(8), intent(in)         :: a2, a3, d1, d2, d3, d4, d5
    real(8), intent(out)        :: al, ar
    ! al = a2 + (-0.4d0 * d1 + 2.2d0 * d2 + 4.8d0 * d3 - 0.6d0 * d4) * one_twelfth
    ! ar = a3 - (-0.6d0 * d2 + 4.8d0 * d3 + 2.2d0 * d4 - 0.4d0 * d5) * one_twelfth
    real(8), parameter :: coef1 = - one_twelfth * 0.4d0
    real(8), parameter :: coef2 =   one_twelfth * 2.2d0
    real(8), parameter :: coef3 =   one_twelfth * 4.8d0
    real(8), parameter :: coef4 = - one_twelfth * 0.6d0
    !al = a2 + (coef1 * d1 + coef2 * d2 + coef3 * d3 + coef4 * d4)
    al = fma(coef1, d1, a2)
    al = fma(coef2, d2, al)
    al = fma(coef3, d3, al)
    al = fma(coef4, d4, al)
    !ar = a3 - (coef1 * d5 + coef2 * d4 + coef3 * d3 + coef4 * d2)
    ar = fma(-coef1, d5, a3)
    ar = fma(-coef2, d4, ar)
    ar = fma(-coef3, d3, ar)
    ar = fma(-coef4, d2, ar)
  end subroutine MUSCL4thnonTVD


  attributes(device) subroutine MUSCL4thTVD(id_tvd, sensor, a2, a3, d1, d2, d3, d4, d5, al, ar)
    integer(kind=4), intent(in) :: id_tvd
    real(sp), intent(in)        :: sensor
    real(8), intent(in)         :: a2, a3, d1, d2, d3, d4, d5
    real(8), intent(out)        :: al, ar
    real(8) delta2
    delta2 = d3 - d33(d2, d3, d4) * one_sixth
    block
      real(8) delta1, dl, dr
      delta1 = d2 - d33(d1, d2, d3) * one_sixth
      dl = minmod(delta1, 4.d0 * delta2)
      dr = minmod(delta2, 4.d0 * delta1)
      al = a2 + fma(2.d0, dr, dl) * one_sixth
    end block
    block
      real(8) delta3, dl, dr
      delta3 = d4 - d33(d3, d4, d5) * one_sixth
      dl = minmod(delta2, 4.d0 * delta3)
      dr = minmod(delta3, 4.d0 * delta2)
      ar = a3 - fma(2.d0, dl, dr) * one_sixth
    end block
  end subroutine MUSCL4thTVD


  attributes(device) subroutine MUSCL4thThreshold(id_tvd, sensor, a2, a3, d1, d2, d3, d4, d5, al, ar)
    use mod_globals, only : threshold
    integer(kind=8), intent(in) :: id_tvd
    real(sp), intent(in)        :: sensor
    real(8), intent(in)         :: a2, a3, d1, d2, d3, d4, d5
    real(8), intent(out)        :: al, ar
    real(8) alr(2)
    integer(kind=2) :: id2 = 0_2
    integer(kind=4) :: id4 = 0_4
    if (sensor < threshold) then
      call MUSCL4thnonTVD(id2, sensor, a2, a3, d1, d2, d3, d4, d5, al, ar)
    else
      call MUSCL4thTVD(id4, sensor, a2, a3, d1, d2, d3, d4, d5, al, ar)
    endif
  end subroutine MUSCL4thThreshold


  attributes(device) subroutine delta4(sensor, a, al, ar)
    use mod_globals, only : sp
    use mod_constant, only : id_tvd
    real(sp), intent(in), value     :: sensor
    real(8), intent(in), contiguous :: a(4)
    real(8), intent(out)            :: al, ar
    real(8) d1, d2, d3
    d1 = -a(1) + a(2)
    d2 = -a(2) + a(3)
    d3 = -a(3) + a(4)
    call MUSCL3rd(id_tvd, sensor, a(2), a(3), d1, d2, d3, al, ar)
  end subroutine delta4


  attributes(device) subroutine delta6(sensor, a, al, ar)
    use mod_globals, only : sp
    use mod_constant, only : id_tvd
    real(sp), intent(in), value     :: sensor
    real(8), intent(in), contiguous :: a(6)
    real(8), intent(out)            :: al, ar
    real(8) d1, d2, d3, d4, d5
    d1 = -a(1) + a(2)
    d2 = -a(2) + a(3)
    d3 = -a(3) + a(4)
    d4 = -a(4) + a(5)
    d5 = -a(5) + a(6)
    call MUSCL4th(id_tvd, sensor, a(3), a(4), d1, d2, d3, d4, d5, al, ar)
  end subroutine delta6
end module calc_muscl
