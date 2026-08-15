!> WENO-Z (Borges et al. 2008) 5th-order reconstruction -- dimension-agnostic
!> scalar-stencil math, the WENO sibling of src/calc_muscl.f90.fypp's
!> delta4/delta6. Only the 5-point (6-point-window) stencil is implemented
!> here; there is no WENO3, so 1D_solver only offers this at ORDER=6 (see
!> 1D_solver/src/calc_slau_kernel.f90.fypp).
module calc_weno
  implicit none
  private
  public delta6_weno
contains
  ! ------------------------- WENO-Z weights ---------------------------
  pure attributes(device) function weno5z_left(v1,v2,v3,v4,v5) result(vf)
    real(8), intent(in) :: v1,v2,v3,v4,v5
    real(8) :: vf, p0,p1,p2, b0,b1,b2, a0,a1,a2, s, tau5
    real(8), parameter :: eps = 1.0d-20, d0=1.0d0/10.0d0, d1=6.0d0/10.0d0, d2=3.0d0/10.0d0
    p0 = ( 2.0d0*v1 - 7.0d0*v2 + 11.0d0*v3)/6.0d0
    p1 = (-1.0d0*v2 + 5.0d0*v3 +  2.0d0*v4)/6.0d0
    p2 = ( 2.0d0*v3 + 5.0d0*v4 -  1.0d0*v5)/6.0d0
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    a0 = d0 * (1.0d0 + (tau5/(b0+eps))**2)
    a1 = d1 * (1.0d0 + (tau5/(b1+eps))**2)
    a2 = d2 * (1.0d0 + (tau5/(b2+eps))**2)
    s  = a0 + a1 + a2
    vf = (a0*p0 + a1*p1 + a2*p2) / s
  end function weno5z_left

  pure attributes(device) function weno5z_right(v1,v2,v3,v4,v5) result(vf)
    real(8), intent(in) :: v1,v2,v3,v4,v5
    real(8) :: vf, p0,p1,p2, b0,b1,b2, a0,a1,a2, s, tau5
    real(8), parameter :: eps = 1.0d-20, d0=1.0d0/10.0d0, d1=6.0d0/10.0d0, d2=3.0d0/10.0d0
    p0 = (-1.0d0*v1 + 5.0d0*v2 +  2.0d0*v3)/6.0d0
    p1 = ( 2.0d0*v2 + 5.0d0*v3 -  1.0d0*v4)/6.0d0
    p2 = (11.0d0*v3 - 7.0d0*v4 +  2.0d0*v5)/6.0d0
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    a0 = d0 * (1.0d0 + (tau5/(b0+eps))**2)
    a1 = d1 * (1.0d0 + (tau5/(b1+eps))**2)
    a2 = d2 * (1.0d0 + (tau5/(b2+eps))**2)
    s  = a0 + a1 + a2
    vf = (a0*p0 + a1*p1 + a2*p2) / s
  end function weno5z_right

  !> Matches calc_muscl's delta6(sensor, a, al, ar) call-site shape minus the
  !> unused sensor arg (WENO has no TVD-style limiter family to dispatch on):
  !> al/ar reconstructed at the face between a(3) and a(4). al's window is
  !> centred on a(3) (weno5z_left over a(1:5)); ar's window is centred on
  !> a(4) (weno5z_right over a(2:6)) -- the same 6-point span calc_muscl's
  !> delta6 already reads.
  pure attributes(device) subroutine delta6_weno(a, al, ar)
    real(8), intent(in), contiguous :: a(6)
    real(8), intent(out)            :: al, ar
    al = weno5z_left (a(1), a(2), a(3), a(4), a(5))
    ar = weno5z_right(a(2), a(3), a(4), a(5), a(6))
  end subroutine delta6_weno
end module calc_weno
