!> WENO-Z (Borges et al. 2008) 5th-order reconstruction -- dimension-agnostic
!> scalar-stencil math, the WENO sibling of src/calc_muscl.f90.fypp's
!> delta4/delta6. Only the 5-point (6-point-window) stencil is implemented
!> here; there is no WENO3, so 1D_solver only offers this at ORDER=6 (see
!> 1D_solver/src/calc_slau_kernel.f90.fypp).
!>
!> Division budget. `div.rn.f64` is not a hardware instruction: ptxas expands
!> each one into MUFU.RCP64H + ~9 DFMA plus a slow-path CALL, ~20 SASS
!> instructions. The textbook form of weno5z_* has seven divisions per call
!> (3 x /6, 3 x tau5/(b+eps), 1 x /s); cuobjdump on sm_89 measured 12 per
!> delta6_weno (the compiler CSEs two of the six /6 between the left and right
!> windows) out of 432 instructions, 109 of them DFMA -- i.e. over half the
!> function was division expansion. Two algebraic rewrites below cut that to
!> **two divisions per call, four per delta6_weno**:
!>
!>   1. The 1/6 on the candidate polynomials is a common factor of the WENO
!>      numerator, so it is folded into the single final normalisation.
!>   2. The three tau5/(b_i+eps) ratios share one reciprocal:
!>      tau5/c_i = tau5 * (prod_{j/=i} c_j) * invd,  invd = 1/(c0*c1*c2).
!>
!> Neither is bit-identical to the original -- (1) rounds the polynomials once
!> instead of twice (slightly *more* accurate) and (2) trades one correctly
!> rounded quotient for ~3 extra roundings on a quantity that is then squared
!> and used only as a smoothness weight. Both are far below the 1e-7 tolerance
!> report/check_vs_numpy.py applies.
module calc_weno
  implicit none
  private
  public delta6_weno, delta8_weno, delta10_weno
contains
  ! ------------------------- WENO-Z weights ---------------------------
  pure attributes(device) function weno5z_left(v1,v2,v3,v4,v5) result(vf)
    real(8), intent(in) :: v1,v2,v3,v4,v5
    real(8) :: vf, q0,q1,q2, b0,b1,b2, c0,c1,c2, a0,a1,a2, tau5
    real(8) :: c01,c12,c02, invd, tinv, r0,r1,r2
    real(8), parameter :: eps = 1.0d-20, d0=1.0d0/10.0d0, d1=6.0d0/10.0d0, d2=3.0d0/10.0d0
    ! Candidate polynomials scaled by 6; the factor is undone in the final
    ! division, where it costs one multiply instead of three divisions.
    q0 =  2.0d0*v1 - 7.0d0*v2 + 11.0d0*v3
    q1 = -1.0d0*v2 + 5.0d0*v3 +  2.0d0*v4
    q2 =  2.0d0*v3 + 5.0d0*v4 -  1.0d0*v5
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    ! Batched inversion. c_i >= eps = 1e-20 so the product is >= 1e-60, ~250
    ! decades above the FP64 min normal -- safe here, but NOT in the fltflt
    ! (double-float) twin, which carries FP32 exponent range; see
    ! 1D_solver/src/calc_slau_kernel.f90.fypp.
    c01 = c0*c1
    c12 = c1*c2
    c02 = c0*c2
    invd = 1.0d0 / (c01*c2)
    tinv = tau5 * invd
    r0 = tinv * c12
    r1 = tinv * c02
    r2 = tinv * c01
    a0 = d0 * (1.0d0 + r0*r0)
    a1 = d1 * (1.0d0 + r1*r1)
    a2 = d2 * (1.0d0 + r2*r2)
    vf = (a0*q0 + a1*q1 + a2*q2) / (6.0d0*(a0 + a1 + a2))
  end function weno5z_left

  !> Right-biased (v^+) reconstruction at the same face. The optimal linear
  !> weights are the MIRROR of the left-biased ones, not the same ones: the
  !> mirror about the face maps cell i-k to cell i+1+k, so this function's q0
  !> (the near-face candidate) pairs with weno5z_left's q2 and takes 3/10,
  !> while q2 (which reaches furthest from the face) takes 1/10.
  !>
  !> This was wrong until 2026-08-16 -- d was copied from weno5z_left in the
  !> same order, which silently made v^+ THIRD order instead of fifth while
  !> v^- stayed fifth. It survived because report/check_vs_numpy.py carried the
  !> identical mistake, so the two implementations agreed with each other; that
  !> check tests consistency between two codes, not correctness of the scheme.
  !> report/check_weno_order.py is the check that catches this, by measuring the
  !> convergence rate on a smooth solution.
  pure attributes(device) function weno5z_right(v1,v2,v3,v4,v5) result(vf)
    real(8), intent(in) :: v1,v2,v3,v4,v5
    real(8) :: vf, q0,q1,q2, b0,b1,b2, c0,c1,c2, a0,a1,a2, tau5
    real(8) :: c01,c12,c02, invd, tinv, r0,r1,r2
    real(8), parameter :: eps = 1.0d-20, d0=3.0d0/10.0d0, d1=6.0d0/10.0d0, d2=1.0d0/10.0d0
    q0 = -1.0d0*v1 + 5.0d0*v2 +  2.0d0*v3
    q1 =  2.0d0*v2 + 5.0d0*v3 -  1.0d0*v4
    q2 = 11.0d0*v3 - 7.0d0*v4 +  2.0d0*v5
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    c01 = c0*c1
    c12 = c1*c2
    c02 = c0*c2
    invd = 1.0d0 / (c01*c2)
    tinv = tau5 * invd
    r0 = tinv * c12
    r1 = tinv * c02
    r2 = tinv * c01
    a0 = d0 * (1.0d0 + r0*r0)
    a1 = d1 * (1.0d0 + r1*r1)
    a2 = d2 * (1.0d0 + r2*r2)
    vf = (a0*q0 + a1*q1 + a2*q2) / (6.0d0*(a0 + a1 + a2))
  end function weno5z_right


  ! ------------------------- WENO7-Z (r = 4) --------------------------
  ! GENERATED by 1D_solver/report/check_weno_order.py --gen 4 -- do not hand
  ! edit; regenerate instead. The candidate reconstructions, the optimal linear
  ! weights, the smoothness-indicator quadratic forms and the tau_7 combination
  ! are all solved from their defining conditions in exact rational arithmetic
  ! by that script, which then measures the convergence rate on a smooth
  ! solution (7.02 observed). Transcribing the published WENO7 beta tables by
  ! hand is the failure mode this avoids: a mistyped coefficient is invisible
  ! on a shock tube, which is first order at the discontinuity for every
  ! scheme. The generator is validated by regenerating r=3 and matching the
  ! hand-written WENO5-Z above term for term.
  !
  ! Same division-reduced shape as WENO5-Z: candidate polynomials carry a
  ! factor 12 that is undone in the single final division, and one batched
  ! reciprocal serves all four tau/(beta+eps) ratios -- 9 divisions -> 2.
  ! The batch product c0*c1*c2*c3 >= 1e-80, far above the FP64 min normal, but
  ! it would flush to zero in FP32/fltflt exponent range, so a double-float
  ! twin must keep the ratios separate (see calc_slau_kernel.f90.fypp).

  pure attributes(device) function weno7z_left(v1, v2, v3, v4, v5, v6, v7) result(vf)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(8) :: vf, q0, q1, q2, q3, b0, b1, b2, b3, c0, c1, c2, c3, a0, a1, a2, a3, tau
    real(8) :: invd, tinv, r0, r1, r2, r3, e01, e23
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=1.0d0/35.0d0, dd1=12.0d0/35.0d0, dd2=18.0d0/35.0d0, dd3=4.0d0/35.0d0
    ! candidate polynomials scaled by 12; undone in the final division
    q0 = - 3.0d0*v1 + 13.0d0*v2 - 23.0d0*v3 + 25.0d0*v4
    q1 = 1.0d0*v2 - 5.0d0*v3 + 13.0d0*v4 + 3.0d0*v5
    q2 = - 1.0d0*v3 + 7.0d0*v4 + 7.0d0*v5 - 1.0d0*v6
    q3 = 3.0d0*v4 + 13.0d0*v5 - 5.0d0*v6 + 1.0d0*v7
    b0 = v1*((547.0d0/240.0d0)*v1 - (647.0d0/40.0d0)*v2 + (2321.0d0/120.0d0)*v3 - (309.0d0/40.0d0)*v4) &
         + v2*((7043.0d0/240.0d0)*v2 - (8623.0d0/120.0d0)*v3 + (3521.0d0/120.0d0)*v4) &
         + v3*((11003.0d0/240.0d0)*v3 - (1567.0d0/40.0d0)*v4) &
         + v4*((2107.0d0/240.0d0)*v4)
    b1 = v2*((89.0d0/80.0d0)*v2 - (821.0d0/120.0d0)*v3 + (267.0d0/40.0d0)*v4 - (247.0d0/120.0d0)*v5) &
         + v3*((2843.0d0/240.0d0)*v3 - (2983.0d0/120.0d0)*v4 + (961.0d0/120.0d0)*v5) &
         + v4*((3443.0d0/240.0d0)*v4 - (1261.0d0/120.0d0)*v5) &
         + v5*((547.0d0/240.0d0)*v5)
    b2 = v3*((547.0d0/240.0d0)*v3 - (1261.0d0/120.0d0)*v4 + (961.0d0/120.0d0)*v5 - (247.0d0/120.0d0)*v6) &
         + v4*((3443.0d0/240.0d0)*v4 - (2983.0d0/120.0d0)*v5 + (267.0d0/40.0d0)*v6) &
         + v5*((2843.0d0/240.0d0)*v5 - (821.0d0/120.0d0)*v6) &
         + v6*((89.0d0/80.0d0)*v6)
    b3 = v4*((2107.0d0/240.0d0)*v4 - (1567.0d0/40.0d0)*v5 + (3521.0d0/120.0d0)*v6 - (309.0d0/40.0d0)*v7) &
         + v5*((11003.0d0/240.0d0)*v5 - (8623.0d0/120.0d0)*v6 + (2321.0d0/120.0d0)*v7) &
         + v6*((7043.0d0/240.0d0)*v6 - (647.0d0/40.0d0)*v7) &
         + v7*((547.0d0/240.0d0)*v7)
    tau = abs(1.0d0*b0 + 3.0d0*b1 - 3.0d0*b2 - 1.0d0*b3)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    c3 = b3 + eps
    e01 = c0*c1
    e23 = c2*c3
    invd = 1.0d0 / (e01*e23)
    tinv = tau * invd
    r0 = tinv * (c1*e23)
    r1 = tinv * (c0*e23)
    r2 = tinv * (e01*c3)
    r3 = tinv * (e01*c2)
    a0 = dd0 * (1.0d0 + r0*r0)
    a1 = dd1 * (1.0d0 + r1*r1)
    a2 = dd2 * (1.0d0 + r2*r2)
    a3 = dd3 * (1.0d0 + r3*r3)
    vf = (a0*q0 + a1*q1 + a2*q2 + a3*q3) / (12.0d0*(a0 + a1 + a2 + a3))
  end function weno7z_left

  pure attributes(device) function weno7z_right(v1, v2, v3, v4, v5, v6, v7) result(vf)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(8) :: vf, q0, q1, q2, q3, b0, b1, b2, b3, c0, c1, c2, c3, a0, a1, a2, a3, tau
    real(8) :: invd, tinv, r0, r1, r2, r3, e01, e23
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=4.0d0/35.0d0, dd1=18.0d0/35.0d0, dd2=12.0d0/35.0d0, dd3=1.0d0/35.0d0
    ! candidate polynomials scaled by 12; undone in the final division
    q3 = - 3.0d0*v7 + 13.0d0*v6 - 23.0d0*v5 + 25.0d0*v4
    q2 = 1.0d0*v6 - 5.0d0*v5 + 13.0d0*v4 + 3.0d0*v3
    q1 = - 1.0d0*v5 + 7.0d0*v4 + 7.0d0*v3 - 1.0d0*v2
    q0 = 3.0d0*v4 + 13.0d0*v3 - 5.0d0*v2 + 1.0d0*v1
    b3 = v7*((547.0d0/240.0d0)*v7 - (647.0d0/40.0d0)*v6 + (2321.0d0/120.0d0)*v5 - (309.0d0/40.0d0)*v4) &
         + v6*((7043.0d0/240.0d0)*v6 - (8623.0d0/120.0d0)*v5 + (3521.0d0/120.0d0)*v4) &
         + v5*((11003.0d0/240.0d0)*v5 - (1567.0d0/40.0d0)*v4) &
         + v4*((2107.0d0/240.0d0)*v4)
    b2 = v6*((89.0d0/80.0d0)*v6 - (821.0d0/120.0d0)*v5 + (267.0d0/40.0d0)*v4 - (247.0d0/120.0d0)*v3) &
         + v5*((2843.0d0/240.0d0)*v5 - (2983.0d0/120.0d0)*v4 + (961.0d0/120.0d0)*v3) &
         + v4*((3443.0d0/240.0d0)*v4 - (1261.0d0/120.0d0)*v3) &
         + v3*((547.0d0/240.0d0)*v3)
    b1 = v5*((547.0d0/240.0d0)*v5 - (1261.0d0/120.0d0)*v4 + (961.0d0/120.0d0)*v3 - (247.0d0/120.0d0)*v2) &
         + v4*((3443.0d0/240.0d0)*v4 - (2983.0d0/120.0d0)*v3 + (267.0d0/40.0d0)*v2) &
         + v3*((2843.0d0/240.0d0)*v3 - (821.0d0/120.0d0)*v2) &
         + v2*((89.0d0/80.0d0)*v2)
    b0 = v4*((2107.0d0/240.0d0)*v4 - (1567.0d0/40.0d0)*v3 + (3521.0d0/120.0d0)*v2 - (309.0d0/40.0d0)*v1) &
         + v3*((11003.0d0/240.0d0)*v3 - (8623.0d0/120.0d0)*v2 + (2321.0d0/120.0d0)*v1) &
         + v2*((7043.0d0/240.0d0)*v2 - (647.0d0/40.0d0)*v1) &
         + v1*((547.0d0/240.0d0)*v1)
    tau = abs(- 1.0d0*b0 - 3.0d0*b1 + 3.0d0*b2 + 1.0d0*b3)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    c3 = b3 + eps
    e01 = c0*c1
    e23 = c2*c3
    invd = 1.0d0 / (e01*e23)
    tinv = tau * invd
    r0 = tinv * (c1*e23)
    r1 = tinv * (c0*e23)
    r2 = tinv * (e01*c3)
    r3 = tinv * (e01*c2)
    a0 = dd0 * (1.0d0 + r0*r0)
    a1 = dd1 * (1.0d0 + r1*r1)
    a2 = dd2 * (1.0d0 + r2*r2)
    a3 = dd3 * (1.0d0 + r3*r3)
    vf = (a0*q0 + a1*q1 + a2*q2 + a3*q3) / (12.0d0*(a0 + a1 + a2 + a3))
  end function weno7z_right


  ! ------------------------- WENO9-Z (r = 5) --------------------------
  ! GENERATED by 1D_solver/report/check_weno_order.py --gen 5. Same provenance
  ! and same division-reduced shape as WENO7-Z above; verified at order 8.95
  ! (the convergence test hits the FP64 round-off floor by N=64 at this order,
  ! so 8.95 is read off the coarse end of the refinement -- see the
  ! ROUNDOFF_FLOOR note in that script). Candidate polynomials carry a factor
  ! 60; one batched reciprocal serves all five ratios, 11 divisions -> 2.

  pure attributes(device) function weno9z_left(v1, v2, v3, v4, v5, v6, v7, v8, v9) result(vf)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8) :: vf, q0, q1, q2, q3, q4, b0, b1, b2, b3, b4, c0, c1, c2, c3, c4, a0, a1, a2, a3, a4, tau
    real(8) :: invd, tinv, r0, r1, r2, r3, r4, e01, e012, e34, e234
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=1.0d0/126.0d0, dd1=10.0d0/63.0d0, dd2=10.0d0/21.0d0, dd3=20.0d0/63.0d0, dd4=5.0d0/126.0d0
    ! candidate polynomials scaled by 60; undone in the final division
    q0 = 12.0d0*v1 - 63.0d0*v2 + 137.0d0*v3 - 163.0d0*v4 + 137.0d0*v5
    q1 = - 3.0d0*v2 + 17.0d0*v3 - 43.0d0*v4 + 77.0d0*v5 + 12.0d0*v6
    q2 = 2.0d0*v3 - 13.0d0*v4 + 47.0d0*v5 + 27.0d0*v6 - 3.0d0*v7
    q3 = - 3.0d0*v4 + 27.0d0*v5 + 47.0d0*v6 - 13.0d0*v7 + 2.0d0*v8
    q4 = 12.0d0*v5 + 77.0d0*v6 - 43.0d0*v7 + 17.0d0*v8 - 3.0d0*v9
    b0 = v1*((11329.0d0/2520.0d0)*v1 - (208501.0d0/5040.0d0)*v2 + (121621.0d0/1680.0d0)*v3 - (288007.0d0/5040.0d0)*v4 + (86329.0d0/5040.0d0)*v5) &
         + v2*((482963.0d0/5040.0d0)*v2 - (142033.0d0/420.0d0)*v3 + (679229.0d0/2520.0d0)*v4 - (411487.0d0/5040.0d0)*v5) &
         + v3*((507131.0d0/1680.0d0)*v3 - (68391.0d0/140.0d0)*v4 + (252941.0d0/1680.0d0)*v5) &
         + v4*((1020563.0d0/5040.0d0)*v4 - (649501.0d0/5040.0d0)*v5) &
         + v5*((53959.0d0/2520.0d0)*v5)
    b1 = v2*((1727.0d0/1260.0d0)*v2 - (60871.0d0/5040.0d0)*v3 + (33071.0d0/1680.0d0)*v4 - (70237.0d0/5040.0d0)*v5 + (18079.0d0/5040.0d0)*v6) &
         + v3*((138563.0d0/5040.0d0)*v3 - (3229.0d0/35.0d0)*v4 + (168509.0d0/2520.0d0)*v5 - (88297.0d0/5040.0d0)*v6) &
         + v4*((135431.0d0/1680.0d0)*v4 - (25499.0d0/210.0d0)*v5 + (55051.0d0/1680.0d0)*v6) &
         + v5*((242723.0d0/5040.0d0)*v5 - (140251.0d0/5040.0d0)*v6) &
         + v6*((11329.0d0/2520.0d0)*v6)
    b2 = v3*((1727.0d0/1260.0d0)*v3 - (51001.0d0/5040.0d0)*v4 + (7547.0d0/560.0d0)*v5 - (38947.0d0/5040.0d0)*v6 + (8209.0d0/5040.0d0)*v7) &
         + v4*((104963.0d0/5040.0d0)*v4 - (24923.0d0/420.0d0)*v5 + (89549.0d0/2520.0d0)*v6 - (38947.0d0/5040.0d0)*v7) &
         + v5*((77051.0d0/1680.0d0)*v5 - (24923.0d0/420.0d0)*v6 + (7547.0d0/560.0d0)*v7) &
         + v6*((104963.0d0/5040.0d0)*v6 - (51001.0d0/5040.0d0)*v7) &
         + v7*((1727.0d0/1260.0d0)*v7)
    b3 = v4*((11329.0d0/2520.0d0)*v4 - (140251.0d0/5040.0d0)*v5 + (55051.0d0/1680.0d0)*v6 - (88297.0d0/5040.0d0)*v7 + (18079.0d0/5040.0d0)*v8) &
         + v5*((242723.0d0/5040.0d0)*v5 - (25499.0d0/210.0d0)*v6 + (168509.0d0/2520.0d0)*v7 - (70237.0d0/5040.0d0)*v8) &
         + v6*((135431.0d0/1680.0d0)*v6 - (3229.0d0/35.0d0)*v7 + (33071.0d0/1680.0d0)*v8) &
         + v7*((138563.0d0/5040.0d0)*v7 - (60871.0d0/5040.0d0)*v8) &
         + v8*((1727.0d0/1260.0d0)*v8)
    b4 = v5*((53959.0d0/2520.0d0)*v5 - (649501.0d0/5040.0d0)*v6 + (252941.0d0/1680.0d0)*v7 - (411487.0d0/5040.0d0)*v8 + (86329.0d0/5040.0d0)*v9) &
         + v6*((1020563.0d0/5040.0d0)*v6 - (68391.0d0/140.0d0)*v7 + (679229.0d0/2520.0d0)*v8 - (288007.0d0/5040.0d0)*v9) &
         + v7*((507131.0d0/1680.0d0)*v7 - (142033.0d0/420.0d0)*v8 + (121621.0d0/1680.0d0)*v9) &
         + v8*((482963.0d0/5040.0d0)*v8 - (208501.0d0/5040.0d0)*v9) &
         + v9*((11329.0d0/2520.0d0)*v9)
    tau = abs(1.0d0*b0 + 2.0d0*b1 - 6.0d0*b2 + 2.0d0*b3 + 1.0d0*b4)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    c3 = b3 + eps
    c4 = b4 + eps
    e01 = c0*c1
    e012 = e01*c2
    e34 = c3*c4
    e234 = c2*e34
    invd = 1.0d0 / (e01*e234)
    tinv = tau * invd
    r0 = tinv * (c1*e234)
    r1 = tinv * (c0*e234)
    r2 = tinv * (e01*e34)
    r3 = tinv * (e012*c4)
    r4 = tinv * (e012*c3)
    a0 = dd0 * (1.0d0 + r0*r0)
    a1 = dd1 * (1.0d0 + r1*r1)
    a2 = dd2 * (1.0d0 + r2*r2)
    a3 = dd3 * (1.0d0 + r3*r3)
    a4 = dd4 * (1.0d0 + r4*r4)
    vf = (a0*q0 + a1*q1 + a2*q2 + a3*q3 + a4*q4) / (60.0d0*(a0 + a1 + a2 + a3 + a4))
  end function weno9z_left

  pure attributes(device) function weno9z_right(v1, v2, v3, v4, v5, v6, v7, v8, v9) result(vf)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8) :: vf, q0, q1, q2, q3, q4, b0, b1, b2, b3, b4, c0, c1, c2, c3, c4, a0, a1, a2, a3, a4, tau
    real(8) :: invd, tinv, r0, r1, r2, r3, r4, e01, e012, e34, e234
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=5.0d0/126.0d0, dd1=20.0d0/63.0d0, dd2=10.0d0/21.0d0, dd3=10.0d0/63.0d0, dd4=1.0d0/126.0d0
    ! candidate polynomials scaled by 60; undone in the final division
    q4 = 12.0d0*v9 - 63.0d0*v8 + 137.0d0*v7 - 163.0d0*v6 + 137.0d0*v5
    q3 = - 3.0d0*v8 + 17.0d0*v7 - 43.0d0*v6 + 77.0d0*v5 + 12.0d0*v4
    q2 = 2.0d0*v7 - 13.0d0*v6 + 47.0d0*v5 + 27.0d0*v4 - 3.0d0*v3
    q1 = - 3.0d0*v6 + 27.0d0*v5 + 47.0d0*v4 - 13.0d0*v3 + 2.0d0*v2
    q0 = 12.0d0*v5 + 77.0d0*v4 - 43.0d0*v3 + 17.0d0*v2 - 3.0d0*v1
    b4 = v9*((11329.0d0/2520.0d0)*v9 - (208501.0d0/5040.0d0)*v8 + (121621.0d0/1680.0d0)*v7 - (288007.0d0/5040.0d0)*v6 + (86329.0d0/5040.0d0)*v5) &
         + v8*((482963.0d0/5040.0d0)*v8 - (142033.0d0/420.0d0)*v7 + (679229.0d0/2520.0d0)*v6 - (411487.0d0/5040.0d0)*v5) &
         + v7*((507131.0d0/1680.0d0)*v7 - (68391.0d0/140.0d0)*v6 + (252941.0d0/1680.0d0)*v5) &
         + v6*((1020563.0d0/5040.0d0)*v6 - (649501.0d0/5040.0d0)*v5) &
         + v5*((53959.0d0/2520.0d0)*v5)
    b3 = v8*((1727.0d0/1260.0d0)*v8 - (60871.0d0/5040.0d0)*v7 + (33071.0d0/1680.0d0)*v6 - (70237.0d0/5040.0d0)*v5 + (18079.0d0/5040.0d0)*v4) &
         + v7*((138563.0d0/5040.0d0)*v7 - (3229.0d0/35.0d0)*v6 + (168509.0d0/2520.0d0)*v5 - (88297.0d0/5040.0d0)*v4) &
         + v6*((135431.0d0/1680.0d0)*v6 - (25499.0d0/210.0d0)*v5 + (55051.0d0/1680.0d0)*v4) &
         + v5*((242723.0d0/5040.0d0)*v5 - (140251.0d0/5040.0d0)*v4) &
         + v4*((11329.0d0/2520.0d0)*v4)
    b2 = v7*((1727.0d0/1260.0d0)*v7 - (51001.0d0/5040.0d0)*v6 + (7547.0d0/560.0d0)*v5 - (38947.0d0/5040.0d0)*v4 + (8209.0d0/5040.0d0)*v3) &
         + v6*((104963.0d0/5040.0d0)*v6 - (24923.0d0/420.0d0)*v5 + (89549.0d0/2520.0d0)*v4 - (38947.0d0/5040.0d0)*v3) &
         + v5*((77051.0d0/1680.0d0)*v5 - (24923.0d0/420.0d0)*v4 + (7547.0d0/560.0d0)*v3) &
         + v4*((104963.0d0/5040.0d0)*v4 - (51001.0d0/5040.0d0)*v3) &
         + v3*((1727.0d0/1260.0d0)*v3)
    b1 = v6*((11329.0d0/2520.0d0)*v6 - (140251.0d0/5040.0d0)*v5 + (55051.0d0/1680.0d0)*v4 - (88297.0d0/5040.0d0)*v3 + (18079.0d0/5040.0d0)*v2) &
         + v5*((242723.0d0/5040.0d0)*v5 - (25499.0d0/210.0d0)*v4 + (168509.0d0/2520.0d0)*v3 - (70237.0d0/5040.0d0)*v2) &
         + v4*((135431.0d0/1680.0d0)*v4 - (3229.0d0/35.0d0)*v3 + (33071.0d0/1680.0d0)*v2) &
         + v3*((138563.0d0/5040.0d0)*v3 - (60871.0d0/5040.0d0)*v2) &
         + v2*((1727.0d0/1260.0d0)*v2)
    b0 = v5*((53959.0d0/2520.0d0)*v5 - (649501.0d0/5040.0d0)*v4 + (252941.0d0/1680.0d0)*v3 - (411487.0d0/5040.0d0)*v2 + (86329.0d0/5040.0d0)*v1) &
         + v4*((1020563.0d0/5040.0d0)*v4 - (68391.0d0/140.0d0)*v3 + (679229.0d0/2520.0d0)*v2 - (288007.0d0/5040.0d0)*v1) &
         + v3*((507131.0d0/1680.0d0)*v3 - (142033.0d0/420.0d0)*v2 + (121621.0d0/1680.0d0)*v1) &
         + v2*((482963.0d0/5040.0d0)*v2 - (208501.0d0/5040.0d0)*v1) &
         + v1*((11329.0d0/2520.0d0)*v1)
    tau = abs(1.0d0*b0 + 2.0d0*b1 - 6.0d0*b2 + 2.0d0*b3 + 1.0d0*b4)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    c3 = b3 + eps
    c4 = b4 + eps
    e01 = c0*c1
    e012 = e01*c2
    e34 = c3*c4
    e234 = c2*e34
    invd = 1.0d0 / (e01*e234)
    tinv = tau * invd
    r0 = tinv * (c1*e234)
    r1 = tinv * (c0*e234)
    r2 = tinv * (e01*e34)
    r3 = tinv * (e012*c4)
    r4 = tinv * (e012*c3)
    a0 = dd0 * (1.0d0 + r0*r0)
    a1 = dd1 * (1.0d0 + r1*r1)
    a2 = dd2 * (1.0d0 + r2*r2)
    a3 = dd3 * (1.0d0 + r3*r3)
    a4 = dd4 * (1.0d0 + r4*r4)
    vf = (a0*q0 + a1*q1 + a2*q2 + a3*q3 + a4*q4) / (60.0d0*(a0 + a1 + a2 + a3 + a4))
  end function weno9z_right

  !> 10-point-window sibling: al/ar at the face between a(5) and a(6).
  pure attributes(device) subroutine delta10_weno(a, al, ar)
    real(8), intent(in), contiguous :: a(10)
    real(8), intent(out)            :: al, ar
    al = weno9z_left (a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9))
    ar = weno9z_right(a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), a(10))
  end subroutine delta10_weno

  !> 8-point-window sibling of delta6_weno: al/ar at the face between a(4) and
  !> a(5). al's 7-point window is centred on a(4), ar's on a(5).
  pure attributes(device) subroutine delta8_weno(a, al, ar)
    real(8), intent(in), contiguous :: a(8)
    real(8), intent(out)            :: al, ar
    al = weno7z_left (a(1), a(2), a(3), a(4), a(5), a(6), a(7))
    ar = weno7z_right(a(2), a(3), a(4), a(5), a(6), a(7), a(8))
  end subroutine delta8_weno

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
