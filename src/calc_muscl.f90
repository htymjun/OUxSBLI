module calc_muscl
  use mod_constant, only : one_third, one_sixth, one_twelfth
  implicit none
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
  !dir$ inline
  attributes(device) function minmod2(x, y) result(ans)
    real(8), intent(in), value :: x, y
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y), 0.d0)
  end function minmod2

  !dir$ inline
  attributes(device) function minmod3(x, y, z) result(ans)
    real(8), intent(in), value :: x, y, z
    real(8) :: ans, sgn
    sgn = sign(1.d0, x)
    ans = sgn * max(min(abs(x), sgn * y, sgn * z), 0.d0)
  end function minmod3

  !dir$ inline
  attributes(device) function d33(d1, d2, d3) result(ans)
    real(8), intent(in), value :: d1, d2, d3 
    real(8) :: ans
    ans =          minmod(d1, 2.d0 * d2, 2.d0 * d3) &
          - 2.d0 * minmod(d2, 2.d0 * d1, 2.d0 * d3) &
                 + minmod(d3, 2.d0 * d1, 2.d0 * d2)
  end function d33

  !dir$ inline
  attributes(device) function MUSCL3rdnonTVD(id_tvd, sensor, a2, a3, d) result(alr)
    integer(kind=2), intent(in), value :: id_tvd
    real(8), intent(in), value         :: sensor, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8), constant :: b = (3.d0 - one_third) / (1.d0 - one_third)
    real(8) :: alr(2)
    alr(1) = a2 + 0.5d0 * (d(1) + 2.d0 * d(2)) * one_third
    alr(2) = a3 - 0.5d0 * (d(3) + 2.d0 * d(2)) * one_third
  end function MUSCL3rdnonTVD

  attributes(device) function MUSCL3rdMinmod(id_tvd, sensor, a2, a3, d) result(alr)
    integer(kind=4), intent(in), value :: id_tvd
    real(8), intent(in), value         :: sensor, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8), constant :: b = (3.d0 - one_third) / (1.d0 - one_third)
    real(8) :: alr(2)
    block
      real(8) dt1, dt2
      dt1 = minmod(d(1), b * d(2))
      dt2 = minmod(d(2), b * d(1))
      alr(1) = a2 + 0.5d0 * (dt1 + 2.d0 * dt2) * one_third
    end block
    block
      real(8) dt3, dt4
      dt3 = minmod(d(3), b * d(2))
      dt4 = minmod(d(2), b * d(3))
      alr(2) = a3 - 0.5d0 * (dt3 + 2.d0 * dt4) * one_third
    end block
  end function MUSCL3rdMinmod

  attributes(device) function MUSCL3rdThreshold(id_tvd, sensor, a2, a3, d) result(alr)
    use mod_globals, only : threshold
    integer(kind=8), intent(in), value :: id_tvd
    real(8), intent(in), value         :: sensor, a2, a3
    real(8), intent(in), device        :: d(3)
    real(8) alr(2)
    integer(kind=2) :: id2
    integer(kind=4) :: id4
    if (sensor < threshold) then
      alr = MUSCL3rdnonTVD(id2,sensor,a2,a3,d)
    else
      alr = MUSCL3rdMinmod(id4,sensor,a2,a3,d)
    endif
  end function MUSCL3rdThreshold

  !dir$ inline
  attributes(device) function MUSCL4thnonTVD(id_tvd, sensor, a2, a3, d) result(alr)
    integer(kind=2), intent(in), value :: id_tvd
    real(8), intent(in), value         :: sensor, a2, a3
    real(8), intent(in), device        :: d(5)
    real(8), constant :: phi = 1.d0 / 30.d0
    real(8) d3(3), alr(2)
    d3(:)  = d(1:3) - 2.d0 * d(2:4) + d(3:5)
    alr(1) = a2 + (2.d0 * d(2) - 12.d0 * phi * d3(1) &
                  + 4.d0 * d(3) - (1.d0 - 12.d0 * phi) * d3(2)) * one_twelfth
    alr(2) = a3 - (4.d0 * d(3) - (1.d0 - 12.d0 * phi) * d3(2) &
                  + 2.d0 * d(4) - 12.d0 * phi * d3(3)) * one_twelfth
  end function MUSCL4thnonTVD

  attributes(device) function MUSCL4thTVD(id_tvd, sensor, a2, a3, d) result(alr)
    integer(kind=4), intent(in), value :: id_tvd
    real(8), intent(in), value         :: sensor, a2, a3
    real(8), intent(in), device        :: d(5)
    real(8) delta2, alr(2)
    delta2 = d(3) - d33(d(2), d(3), d(4)) * one_sixth
    block
      real(8) delta1, dl, dr
      delta1 = d(2) - d33(d(1), d(2), d(3)) * one_sixth
      dl = minmod(delta1, 4.d0 * delta2)
      dr = minmod(delta2, 4.d0 * delta1)
      alr(1) = a2 + (dl + 2.d0 * dr) * one_sixth
    end block
    block
      real(8) delta3, dl, dr
      delta3 = d(4) - d33(d(3), d(4), d(5)) * one_sixth
      dl = minmod(delta2, 4.d0 * delta3)
      dr = minmod(delta3, 4.d0 * delta2)
      alr(2) = a3 - (dr + 2.d0 * dl) * one_sixth
    end block
  end function MUSCL4thTVD

  attributes(device) function MUSCL4thThreshold(id_tvd, sensor, a2, a3, d) result(alr)
    use mod_globals, only : threshold
    integer(kind=8), intent(in), value :: id_tvd
    real(8), intent(in), value         :: sensor, a2, a3
    real(8), intent(in), device        :: d(5)
    real(8) alr(2)
    integer(kind=2) :: id2
    integer(kind=4) :: id4
    if (sensor < threshold) then
      alr = MUSCL4thnonTVD(id2,sensor,a2,a3,d)
    else
      alr = MUSCL4thTVD(id4,sensor,a2,a3,d)
    endif
  end function MUSCL4thThreshold

  !dir$ inline
  attributes(device) function delta4(sensor, a) result(alr)
    use mod_globals, only : id_tvd
    real(8), intent(in), value :: sensor
    real(8), intent(in)        :: a(4)
    real(8) :: alr(2), d(3)
    d(:) = -a(1:3) + a(2:4)
    alr  = MUSCL3rd(id_tvd, sensor, a(2), a(3), d)
  end function delta4

  !dir$ inline
  attributes(device) function delta6(sensor, a) result(alr)
    use mod_globals, only : id_tvd
    real(8), intent(in), value :: sensor
    real(8), intent(in)        :: a(6)
    real(8) :: alr(2), d(5)
    d(:) = -a(1:5) + a(2:6)
    alr  = MUSCL4th(id_tvd, sensor, a(3), a(4), d)
  end function delta6

  attributes(device) subroutine calc_4points(sensor, rho, u, v, w, p, rho2, p2, V2)
    real(8), intent(in), value   :: sensor
    real(8), intent(in)          :: rho(4), u(4), v(4), w(4), p(4)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    rho2(:) = delta4(sensor, rho)
    V2(:,1) = delta4(sensor, u)
    V2(:,2) = delta4(sensor, v)
    V2(:,3) = delta4(sensor, w)
    p2(:)   = delta4(sensor, p)
  end subroutine calc_4points

  attributes(device) subroutine calc_6points(sensor, rho, u, v, w, p, rho2, p2, V2)
    real(8), intent(in), value   :: sensor
    real(8), intent(in)          :: rho(6), u(6), v(6), w(6), p(6)
    real(8), intent(out), device :: rho2(2), p2(2), V2(2,3)
    rho2(:) = delta6(sensor, rho)
    V2(:,1) = delta6(sensor, u)
    V2(:,2) = delta6(sensor, v)
    V2(:,3) = delta6(sensor, w)
    p2(:)   = delta6(sensor, p)
  end subroutine calc_6points
end module calc_muscl

