!> WENO5-Z microbenchmark kernels.
!>
!> NOTE -- do NOT put `value` on the real(8) dummies of the attributes(device)
!> routines below. nvfortran 24.7 miscompiles a device routine with
!> `intent(in), value` dummies when it is called from another device routine
!> that also has `value` dummies: at -O2/-fast every result comes back NaN,
!> while -O0 is correct and -Kieee / -Mnofprelaxed / -Mnoflushz / -Mnovect do
!> not help. One level of nesting (kernel -> device routine) is fine, which is
!> why the weight_poly32_* modes always worked while every var_* mode returned
!> NaN -- including var_fp64_seq, which contains no FP32 at all. That NaN was
!> previously misattributed to FP32 overflow in the WENO-Z tau5/(beta+eps)
!> ratio; the ratio_cap guard in weights5_32_left is still correct and still needed,
!> but it was not the cause. Plain `intent(in)` is semantically equivalent here
!> (these are read-only scalars) and generates correct code.
!>
!> DIVISION BUDGET -- kept in step with src/calc_weno.f90 (see
!> report/rtx4060_weno_division_reduction.md). `div.rn.f64` is ~20 SASS
!> instructions (MUFU.RCP64H + Newton DFMAs + a CALL to the slow path), and
!> `-fast` does NOT fold `/6.0d0` into a multiply. The textbook split form used
!> here had 18 divisions per weno5_64_pair: 6 per weights5_64_left (3 ratios + 3
!> normalisations) and 3 per poly64_*. Now 4:
!>
!>   weights5_64_left  6 -> 2   batched ratio inversion + one reciprocal for 1/s
!>   weights5_32_left  6 -> 4   normalisation batched only; the ratio batch would
!>                       underflow (c_i = eps = 1e-20 on a plateau, so
!>                       c0*c1*c2 = 1e-60 flushes to 0 in FP32 -> 1/0 -> NaN),
!>                       exactly as in the fltflt twin in the solver
!>   poly64_*   3 -> 0   the 1/6 folded into the polynomial coefficients
!>   poly32_*   3 -> 0   ditto
!>
!> The weights/poly SPLIT is deliberately preserved -- it is the thing this
!> benchmark exists to measure -- so the 1/6 is folded into the coefficients
!> rather than into the final normalisation the way the solver does it. Folding
!> it into the weights instead would move FP64 work across the very boundary
!> being measured.
module weno_micro_kernels
  use cudafor
  implicit none

  integer, parameter :: face_threads = 128
  integer, parameter :: block_threads = 2 * face_threads
  integer, parameter :: mode_fp64 = 1
  integer, parameter :: mode_fp32 = 2
  integer, parameter :: mode_rho64 = 3
  integer, parameter :: mode_u64 = 4
  integer, parameter :: mode_p64 = 5

contains
  attributes(device) subroutine weights5_64_left(v1, v2, v3, v4, v5, w0, w1, w2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(8), intent(out) :: w0, w1, w2
    real(8) :: b0, b1, b2, c0, c1, c2, c01, c12, c02, a0, a1, a2, tau5
    real(8) :: invd, tinv, r0, r1, r2, invs
    real(8), parameter :: eps = 1.0d-20
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    ! One reciprocal for all three ratios: tau5/c_i = tau5 * prod_{j/=i} c_j / D.
    ! Safe in FP64 only -- D >= 1e-60, ~250 decades above min normal.
    c01 = c0*c1
    c12 = c1*c2
    c02 = c0*c2
    invd = 1.0d0 / (c01*c2)
    tinv = tau5 * invd
    r0 = tinv * c12
    r1 = tinv * c02
    r2 = tinv * c01
    a0 = 0.1d0 * (1.0d0 + r0*r0)
    a1 = 0.6d0 * (1.0d0 + r1*r1)
    a2 = 0.3d0 * (1.0d0 + r2*r2)
    invs = 1.0d0 / (a0 + a1 + a2)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
  end subroutine weights5_64_left

  ! Right-biased twin: identical betas, MIRRORED optimal weights. The mirror
  ! about the face maps cell i-k to i+1+k, so the near-face candidate takes
  ! 3/10 and the one reaching furthest takes 1/10 -- d reversed. Using the
  ! left-biased d here (as this file did until 2026-08-17) silently makes v^+
  ! THIRD order instead of fifth; instruction counts are identical, so it does
  ! not show up in timings. report/check_weno_order.py's micro_split_test
  ! is what catches it.
  attributes(device) subroutine weights5_64_right(v1, v2, v3, v4, v5, w0, w1, w2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(8), intent(out) :: w0, w1, w2
    real(8) :: b0, b1, b2, c0, c1, c2, c01, c12, c02, a0, a1, a2, tau5
    real(8) :: invd, tinv, r0, r1, r2, invs
    real(8), parameter :: eps = 1.0d-20
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    c0 = b0 + eps
    c1 = b1 + eps
    c2 = b2 + eps
    ! One reciprocal for all three ratios: tau5/c_i = tau5 * prod_{j/=i} c_j / D.
    ! Safe in FP64 only -- D >= 1e-60, ~250 decades above min normal.
    c01 = c0*c1
    c12 = c1*c2
    c02 = c0*c2
    invd = 1.0d0 / (c01*c2)
    tinv = tau5 * invd
    r0 = tinv * c12
    r1 = tinv * c02
    r2 = tinv * c01
    a0 = 0.3d0 * (1.0d0 + r0*r0)
    a1 = 0.6d0 * (1.0d0 + r1*r1)
    a2 = 0.1d0 * (1.0d0 + r2*r2)
    invs = 1.0d0 / (a0 + a1 + a2)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
  end subroutine weights5_64_right

  attributes(device) subroutine weights5_32_left(v1, v2, v3, v4, v5, w0, w1, w2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(4), intent(out) :: w0, w1, w2
    real(4) :: x1, x2, x3, x4, x5
    real(4) :: b0, b1, b2, a0, a1, a2, tau5, r0, r1, r2, invs
    real(4), parameter :: eps = 1.0e-20
    real(4), parameter :: ratio_cap = 1.0e9
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    b0 = (13.0_4/12.0_4)*(x1 - 2.0_4*x2 + x3)**2 + 0.25_4*(x1 - 4.0_4*x2 + 3.0_4*x3)**2
    b1 = (13.0_4/12.0_4)*(x2 - 2.0_4*x3 + x4)**2 + 0.25_4*(x2 - x4)**2
    b2 = (13.0_4/12.0_4)*(x3 - 2.0_4*x4 + x5)**2 + 0.25_4*(3.0_4*x3 - 4.0_4*x4 + x5)**2
    tau5 = abs(b0 - b2)
    ! NOT batched: in FP32, c0*c1*c2 underflows to 0 on a smooth plateau where
    ! every b_i is 0 (eps^3 = 1e-60), giving 1/0 = Inf and then 0*Inf = NaN.
    r0 = min(tau5/(b0+eps), ratio_cap)
    r1 = min(tau5/(b1+eps), ratio_cap)
    r2 = min(tau5/(b2+eps), ratio_cap)
    a0 = 0.1_4 * (1.0_4 + r0*r0)
    a1 = 0.6_4 * (1.0_4 + r1*r1)
    a2 = 0.3_4 * (1.0_4 + r2*r2)
    ! The normalisation IS safe to batch: a_i >= d_i so s >= 1 always, and
    ! s <= 1.8e18 with ratio_cap=1e9, well inside FP32 range.
    invs = 1.0_4 / (a0 + a1 + a2)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
  end subroutine weights5_32_left

  ! Right-biased twin: identical betas, MIRRORED optimal weights. The mirror
  ! about the face maps cell i-k to i+1+k, so the near-face candidate takes
  ! 3/10 and the one reaching furthest takes 1/10 -- d reversed. Using the
  ! left-biased d here (as this file did until 2026-08-17) silently makes v^+
  ! THIRD order instead of fifth; instruction counts are identical, so it does
  ! not show up in timings. report/check_weno_order.py's micro_split_test
  ! is what catches it.
  attributes(device) subroutine weights5_32_right(v1, v2, v3, v4, v5, w0, w1, w2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(4), intent(out) :: w0, w1, w2
    real(4) :: x1, x2, x3, x4, x5
    real(4) :: b0, b1, b2, a0, a1, a2, tau5, r0, r1, r2, invs
    real(4), parameter :: eps = 1.0e-20
    real(4), parameter :: ratio_cap = 1.0e9
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    b0 = (13.0_4/12.0_4)*(x1 - 2.0_4*x2 + x3)**2 + 0.25_4*(x1 - 4.0_4*x2 + 3.0_4*x3)**2
    b1 = (13.0_4/12.0_4)*(x2 - 2.0_4*x3 + x4)**2 + 0.25_4*(x2 - x4)**2
    b2 = (13.0_4/12.0_4)*(x3 - 2.0_4*x4 + x5)**2 + 0.25_4*(3.0_4*x3 - 4.0_4*x4 + x5)**2
    tau5 = abs(b0 - b2)
    ! NOT batched: in FP32, c0*c1*c2 underflows to 0 on a smooth plateau where
    ! every b_i is 0 (eps^3 = 1e-60), giving 1/0 = Inf and then 0*Inf = NaN.
    r0 = min(tau5/(b0+eps), ratio_cap)
    r1 = min(tau5/(b1+eps), ratio_cap)
    r2 = min(tau5/(b2+eps), ratio_cap)
    a0 = 0.3_4 * (1.0_4 + r0*r0)
    a1 = 0.6_4 * (1.0_4 + r1*r1)
    a2 = 0.1_4 * (1.0_4 + r2*r2)
    ! The normalisation IS safe to batch: a_i >= d_i so s >= 1 always, and
    ! s <= 1.8e18 with ratio_cap=1e9, well inside FP32 range.
    invs = 1.0_4 / (a0 + a1 + a2)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
  end subroutine weights5_32_right

  ! The 1/6 is folded into the coefficients, so each candidate polynomial is a
  ! 3-term FMA chain with no division. Written as named parameters rather than
  ! decimal literals so the FP32 twins below get correctly rounded values at
  ! real(4) instead of truncated digits.
  attributes(device) subroutine poly5_64_left(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(8), intent(out) :: p0, p1, p2
    real(8), parameter :: c1 = 1.0d0/6.0d0, c2 = 2.0d0/6.0d0, c5 = 5.0d0/6.0d0
    real(8), parameter :: c7 = 7.0d0/6.0d0, c11 = 11.0d0/6.0d0
    p0 =  c2*v1 - c7*v2 + c11*v3
    p1 = -c1*v2 + c5*v3 +  c2*v4
    p2 =  c2*v3 + c5*v4 -  c1*v5
  end subroutine poly5_64_left

  attributes(device) subroutine poly5_64_right(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(8), intent(out) :: p0, p1, p2
    real(8), parameter :: c1 = 1.0d0/6.0d0, c2 = 2.0d0/6.0d0, c5 = 5.0d0/6.0d0
    real(8), parameter :: c7 = 7.0d0/6.0d0, c11 = 11.0d0/6.0d0
    p0 = -c1*v1 + c5*v2 +  c2*v3
    p1 =  c2*v2 + c5*v3 -  c1*v4
    p2 = c11*v3 - c7*v4 +  c2*v5
  end subroutine poly5_64_right

  attributes(device) subroutine poly5_32_left(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(4), intent(out) :: p0, p1, p2
    real(4) :: x1, x2, x3, x4, x5
    real(4), parameter :: c1 = 1.0_4/6.0_4, c2 = 2.0_4/6.0_4, c5 = 5.0_4/6.0_4
    real(4), parameter :: c7 = 7.0_4/6.0_4, c11 = 11.0_4/6.0_4
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    p0 =  c2*x1 - c7*x2 + c11*x3
    p1 = -c1*x2 + c5*x3 +  c2*x4
    p2 =  c2*x3 + c5*x4 -  c1*x5
  end subroutine poly5_32_left

  attributes(device) subroutine poly5_32_right(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(4), intent(out) :: p0, p1, p2
    real(4) :: x1, x2, x3, x4, x5
    real(4), parameter :: c1 = 1.0_4/6.0_4, c2 = 2.0_4/6.0_4, c5 = 5.0_4/6.0_4
    real(4), parameter :: c7 = 7.0_4/6.0_4, c11 = 11.0_4/6.0_4
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    p0 = -c1*x1 + c5*x2 +  c2*x3
    p1 =  c2*x2 + c5*x3 -  c1*x4
    p2 = c11*x3 - c7*x4 +  c2*x5
  end subroutine poly5_32_right

  attributes(device) subroutine weno5_64_pair(a0, a1, a2, a3, a4, a5, ql, qr)
    real(8), intent(in) :: a0, a1, a2, a3, a4, a5
    real(8), intent(out) :: ql, qr
    real(8) :: w0, w1, w2, p0, p1, p2
    call weights5_64_left(a0, a1, a2, a3, a4, w0, w1, w2)
    call poly5_64_left(a0, a1, a2, a3, a4, p0, p1, p2)
    ql = w0*p0 + w1*p1 + w2*p2
    call weights5_64_right(a1, a2, a3, a4, a5, w0, w1, w2)
    call poly5_64_right(a1, a2, a3, a4, a5, p0, p1, p2)
    qr = w0*p0 + w1*p1 + w2*p2
  end subroutine weno5_64_pair

  attributes(device) subroutine weno5_32_pair(a0, a1, a2, a3, a4, a5, ql, qr)
    real(8), intent(in) :: a0, a1, a2, a3, a4, a5
    real(8), intent(out) :: ql, qr
    real(4) :: w0, w1, w2, p0, p1, p2
    call weights5_32_left(a0, a1, a2, a3, a4, w0, w1, w2)
    call poly5_32_left(a0, a1, a2, a3, a4, p0, p1, p2)
    ql = real(w0*p0 + w1*p1 + w2*p2, 8)
    call weights5_32_right(a1, a2, a3, a4, a5, w0, w1, w2)
    call poly5_32_right(a1, a2, a3, a4, a5, p0, p1, p2)
    qr = real(w0*p0 + w1*p1 + w2*p2, 8)
  end subroutine weno5_32_pair


  ! ==================== WENO7-Z (r=4) and WENO9-Z (r=5) ====================
  ! GENERATED by report/check_weno_order.py --gen-micro 4 / --gen-micro 5.
  ! Do not hand edit; regenerate. Every coefficient -- candidates, optimal
  ! weights, the beta quadratic forms and the tau combination -- is solved from
  ! its defining condition in exact rational arithmetic there, and the split
  ! form is order-verified on BOTH biases by micro_split_test().
  !
  ! WENO5-Z above stays hand-written on purpose: its betas use the textbook
  ! sum-of-two-squares form, which is a WENO5 special case (the published
  ! WENO7/9 betas are the dense quadratic forms emitted below). Regenerating
  ! r=3 would swap 2 squares for 6 products and change WENO5's instruction
  ! count, invalidating every prior WENO5 timing.

  attributes(device) subroutine weights7_64_left(v1, v2, v3, v4, v5, v6, v7, w0, w1, w2, w3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(8), intent(out) :: w0, w1, w2, w3
    real(8) :: b0, b1, b2, b3, tau, invs
    real(8) :: a0, a1, a2, a3, r0, r1, r2, r3
    real(8) :: c0, c1, c2, c3, invd, tinv, e01, e23
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=1.0d0/35.0d0, dd1=12.0d0/35.0d0, dd2=18.0d0/35.0d0, dd3=4.0d0/35.0d0
    ! betas scaled by 240; the common factor cancels in tau/(beta+eps)
    b0 = v1*(547.0d0*v1 - 3882.0d0*v2 + 4642.0d0*v3 - 1854.0d0*v4) &
         + v2*(7043.0d0*v2 - 17246.0d0*v3 + 7042.0d0*v4) &
         + v3*(11003.0d0*v3 - 9402.0d0*v4) &
         + v4*(2107.0d0*v4)
    b1 = v2*(267.0d0*v2 - 1642.0d0*v3 + 1602.0d0*v4 - 494.0d0*v5) &
         + v3*(2843.0d0*v3 - 5966.0d0*v4 + 1922.0d0*v5) &
         + v4*(3443.0d0*v4 - 2522.0d0*v5) &
         + v5*(547.0d0*v5)
    b2 = v3*(547.0d0*v3 - 2522.0d0*v4 + 1922.0d0*v5 - 494.0d0*v6) &
         + v4*(3443.0d0*v4 - 5966.0d0*v5 + 1602.0d0*v6) &
         + v5*(2843.0d0*v5 - 1642.0d0*v6) &
         + v6*(267.0d0*v6)
    b3 = v4*(2107.0d0*v4 - 9402.0d0*v5 + 7042.0d0*v6 - 1854.0d0*v7) &
         + v5*(11003.0d0*v5 - 17246.0d0*v6 + 4642.0d0*v7) &
         + v6*(7043.0d0*v6 - 3882.0d0*v7) &
         + v7*(547.0d0*v7)
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
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0d0 / (a0 + a1 + a2 + a3)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
  end subroutine weights7_64_left

  attributes(device) subroutine weights7_64_right(v1, v2, v3, v4, v5, v6, v7, w0, w1, w2, w3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(8), intent(out) :: w0, w1, w2, w3
    real(8) :: b0, b1, b2, b3, tau, invs
    real(8) :: a0, a1, a2, a3, r0, r1, r2, r3
    real(8) :: c0, c1, c2, c3, invd, tinv, e01, e23
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=4.0d0/35.0d0, dd1=18.0d0/35.0d0, dd2=12.0d0/35.0d0, dd3=1.0d0/35.0d0
    ! betas scaled by 240; the common factor cancels in tau/(beta+eps)
    b3 = v7*(547.0d0*v7 - 3882.0d0*v6 + 4642.0d0*v5 - 1854.0d0*v4) &
         + v6*(7043.0d0*v6 - 17246.0d0*v5 + 7042.0d0*v4) &
         + v5*(11003.0d0*v5 - 9402.0d0*v4) &
         + v4*(2107.0d0*v4)
    b2 = v6*(267.0d0*v6 - 1642.0d0*v5 + 1602.0d0*v4 - 494.0d0*v3) &
         + v5*(2843.0d0*v5 - 5966.0d0*v4 + 1922.0d0*v3) &
         + v4*(3443.0d0*v4 - 2522.0d0*v3) &
         + v3*(547.0d0*v3)
    b1 = v5*(547.0d0*v5 - 2522.0d0*v4 + 1922.0d0*v3 - 494.0d0*v2) &
         + v4*(3443.0d0*v4 - 5966.0d0*v3 + 1602.0d0*v2) &
         + v3*(2843.0d0*v3 - 1642.0d0*v2) &
         + v2*(267.0d0*v2)
    b0 = v4*(2107.0d0*v4 - 9402.0d0*v3 + 7042.0d0*v2 - 1854.0d0*v1) &
         + v3*(11003.0d0*v3 - 17246.0d0*v2 + 4642.0d0*v1) &
         + v2*(7043.0d0*v2 - 3882.0d0*v1) &
         + v1*(547.0d0*v1)
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
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0d0 / (a0 + a1 + a2 + a3)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
  end subroutine weights7_64_right

  attributes(device) subroutine weights7_32_left(v1, v2, v3, v4, v5, v6, v7, w0, w1, w2, w3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(4), intent(out) :: w0, w1, w2, w3
    real(4) :: x1, x2, x3, x4, x5, x6, x7
    real(4) :: b0, b1, b2, b3, tau, invs
    real(4) :: a0, a1, a2, a3, r0, r1, r2, r3
    real(4), parameter :: eps = 1.0e-20
    real(4), parameter :: ratio_cap = 1.0e9
    real(4), parameter :: dd0=1.0_4/35.0_4, dd1=12.0_4/35.0_4, dd2=18.0_4/35.0_4, dd3=4.0_4/35.0_4
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    ! betas scaled by 240; the common factor cancels in tau/(beta+eps)
    b0 = x1*(547.0_4*x1 - 3882.0_4*x2 + 4642.0_4*x3 - 1854.0_4*x4) &
         + x2*(7043.0_4*x2 - 17246.0_4*x3 + 7042.0_4*x4) &
         + x3*(11003.0_4*x3 - 9402.0_4*x4) &
         + x4*(2107.0_4*x4)
    b1 = x2*(267.0_4*x2 - 1642.0_4*x3 + 1602.0_4*x4 - 494.0_4*x5) &
         + x3*(2843.0_4*x3 - 5966.0_4*x4 + 1922.0_4*x5) &
         + x4*(3443.0_4*x4 - 2522.0_4*x5) &
         + x5*(547.0_4*x5)
    b2 = x3*(547.0_4*x3 - 2522.0_4*x4 + 1922.0_4*x5 - 494.0_4*x6) &
         + x4*(3443.0_4*x4 - 5966.0_4*x5 + 1602.0_4*x6) &
         + x5*(2843.0_4*x5 - 1642.0_4*x6) &
         + x6*(267.0_4*x6)
    b3 = x4*(2107.0_4*x4 - 9402.0_4*x5 + 7042.0_4*x6 - 1854.0_4*x7) &
         + x5*(11003.0_4*x5 - 17246.0_4*x6 + 4642.0_4*x7) &
         + x6*(7043.0_4*x6 - 3882.0_4*x7) &
         + x7*(547.0_4*x7)
    tau = abs(1.0_4*b0 + 3.0_4*b1 - 3.0_4*b2 - 1.0_4*b3)
    ! NOT batched: the product of the c_i underflows to 0 in FP32 on a
    ! smooth plateau (eps**r), giving 1/0 = Inf and then 0*Inf = NaN.
    r0 = min(tau/(b0+eps), ratio_cap)
    r1 = min(tau/(b1+eps), ratio_cap)
    r2 = min(tau/(b2+eps), ratio_cap)
    r3 = min(tau/(b3+eps), ratio_cap)
    a0 = dd0 * (1.0_4 + r0*r0)
    a1 = dd1 * (1.0_4 + r1*r1)
    a2 = dd2 * (1.0_4 + r2*r2)
    a3 = dd3 * (1.0_4 + r3*r3)
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0_4 / (a0 + a1 + a2 + a3)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
  end subroutine weights7_32_left

  attributes(device) subroutine weights7_32_right(v1, v2, v3, v4, v5, v6, v7, w0, w1, w2, w3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(4), intent(out) :: w0, w1, w2, w3
    real(4) :: x1, x2, x3, x4, x5, x6, x7
    real(4) :: b0, b1, b2, b3, tau, invs
    real(4) :: a0, a1, a2, a3, r0, r1, r2, r3
    real(4), parameter :: eps = 1.0e-20
    real(4), parameter :: ratio_cap = 1.0e9
    real(4), parameter :: dd0=4.0_4/35.0_4, dd1=18.0_4/35.0_4, dd2=12.0_4/35.0_4, dd3=1.0_4/35.0_4
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    ! betas scaled by 240; the common factor cancels in tau/(beta+eps)
    b3 = x7*(547.0_4*x7 - 3882.0_4*x6 + 4642.0_4*x5 - 1854.0_4*x4) &
         + x6*(7043.0_4*x6 - 17246.0_4*x5 + 7042.0_4*x4) &
         + x5*(11003.0_4*x5 - 9402.0_4*x4) &
         + x4*(2107.0_4*x4)
    b2 = x6*(267.0_4*x6 - 1642.0_4*x5 + 1602.0_4*x4 - 494.0_4*x3) &
         + x5*(2843.0_4*x5 - 5966.0_4*x4 + 1922.0_4*x3) &
         + x4*(3443.0_4*x4 - 2522.0_4*x3) &
         + x3*(547.0_4*x3)
    b1 = x5*(547.0_4*x5 - 2522.0_4*x4 + 1922.0_4*x3 - 494.0_4*x2) &
         + x4*(3443.0_4*x4 - 5966.0_4*x3 + 1602.0_4*x2) &
         + x3*(2843.0_4*x3 - 1642.0_4*x2) &
         + x2*(267.0_4*x2)
    b0 = x4*(2107.0_4*x4 - 9402.0_4*x3 + 7042.0_4*x2 - 1854.0_4*x1) &
         + x3*(11003.0_4*x3 - 17246.0_4*x2 + 4642.0_4*x1) &
         + x2*(7043.0_4*x2 - 3882.0_4*x1) &
         + x1*(547.0_4*x1)
    tau = abs(- 1.0_4*b0 - 3.0_4*b1 + 3.0_4*b2 + 1.0_4*b3)
    ! NOT batched: the product of the c_i underflows to 0 in FP32 on a
    ! smooth plateau (eps**r), giving 1/0 = Inf and then 0*Inf = NaN.
    r0 = min(tau/(b0+eps), ratio_cap)
    r1 = min(tau/(b1+eps), ratio_cap)
    r2 = min(tau/(b2+eps), ratio_cap)
    r3 = min(tau/(b3+eps), ratio_cap)
    a0 = dd0 * (1.0_4 + r0*r0)
    a1 = dd1 * (1.0_4 + r1*r1)
    a2 = dd2 * (1.0_4 + r2*r2)
    a3 = dd3 * (1.0_4 + r3*r3)
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0_4 / (a0 + a1 + a2 + a3)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
  end subroutine weights7_32_right

  attributes(device) subroutine poly7_64_left(v1, v2, v3, v4, v5, v6, v7, p0, p1, p2, p3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(8), intent(out) :: p0, p1, p2, p3
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p0 = - (1.0d0/4.0d0)*v1 + (13.0d0/12.0d0)*v2 - (23.0d0/12.0d0)*v3 + (25.0d0/12.0d0)*v4
    p1 = (1.0d0/12.0d0)*v2 - (5.0d0/12.0d0)*v3 + (13.0d0/12.0d0)*v4 + (1.0d0/4.0d0)*v5
    p2 = - (1.0d0/12.0d0)*v3 + (7.0d0/12.0d0)*v4 + (7.0d0/12.0d0)*v5 - (1.0d0/12.0d0)*v6
    p3 = (1.0d0/4.0d0)*v4 + (13.0d0/12.0d0)*v5 - (5.0d0/12.0d0)*v6 + (1.0d0/12.0d0)*v7
  end subroutine poly7_64_left

  attributes(device) subroutine poly7_64_right(v1, v2, v3, v4, v5, v6, v7, p0, p1, p2, p3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(8), intent(out) :: p0, p1, p2, p3
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p3 = - (1.0d0/4.0d0)*v7 + (13.0d0/12.0d0)*v6 - (23.0d0/12.0d0)*v5 + (25.0d0/12.0d0)*v4
    p2 = (1.0d0/12.0d0)*v6 - (5.0d0/12.0d0)*v5 + (13.0d0/12.0d0)*v4 + (1.0d0/4.0d0)*v3
    p1 = - (1.0d0/12.0d0)*v5 + (7.0d0/12.0d0)*v4 + (7.0d0/12.0d0)*v3 - (1.0d0/12.0d0)*v2
    p0 = (1.0d0/4.0d0)*v4 + (13.0d0/12.0d0)*v3 - (5.0d0/12.0d0)*v2 + (1.0d0/12.0d0)*v1
  end subroutine poly7_64_right

  attributes(device) subroutine poly7_32_left(v1, v2, v3, v4, v5, v6, v7, p0, p1, p2, p3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(4), intent(out) :: p0, p1, p2, p3
    real(4) :: x1, x2, x3, x4, x5, x6, x7
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p0 = - (1.0_4/4.0_4)*x1 + (13.0_4/12.0_4)*x2 - (23.0_4/12.0_4)*x3 + (25.0_4/12.0_4)*x4
    p1 = (1.0_4/12.0_4)*x2 - (5.0_4/12.0_4)*x3 + (13.0_4/12.0_4)*x4 + (1.0_4/4.0_4)*x5
    p2 = - (1.0_4/12.0_4)*x3 + (7.0_4/12.0_4)*x4 + (7.0_4/12.0_4)*x5 - (1.0_4/12.0_4)*x6
    p3 = (1.0_4/4.0_4)*x4 + (13.0_4/12.0_4)*x5 - (5.0_4/12.0_4)*x6 + (1.0_4/12.0_4)*x7
  end subroutine poly7_32_left

  attributes(device) subroutine poly7_32_right(v1, v2, v3, v4, v5, v6, v7, p0, p1, p2, p3)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7
    real(4), intent(out) :: p0, p1, p2, p3
    real(4) :: x1, x2, x3, x4, x5, x6, x7
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p3 = - (1.0_4/4.0_4)*x7 + (13.0_4/12.0_4)*x6 - (23.0_4/12.0_4)*x5 + (25.0_4/12.0_4)*x4
    p2 = (1.0_4/12.0_4)*x6 - (5.0_4/12.0_4)*x5 + (13.0_4/12.0_4)*x4 + (1.0_4/4.0_4)*x3
    p1 = - (1.0_4/12.0_4)*x5 + (7.0_4/12.0_4)*x4 + (7.0_4/12.0_4)*x3 - (1.0_4/12.0_4)*x2
    p0 = (1.0_4/4.0_4)*x4 + (13.0_4/12.0_4)*x3 - (5.0_4/12.0_4)*x2 + (1.0_4/12.0_4)*x1
  end subroutine poly7_32_right

  attributes(device) subroutine weno7_64_pair(a1, a2, a3, a4, a5, a6, a7, a8, ql, qr)
    real(8), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8
    real(8), intent(out) :: ql, qr
    real(8) :: w0, w1, w2, w3, p0, p1, p2, p3
    call weights7_64_left(a1, a2, a3, a4, a5, a6, a7, w0, w1, w2, w3)
    call poly7_64_left(a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3)
    ql = w0*p0 + w1*p1 + w2*p2 + w3*p3
    call weights7_64_right(a2, a3, a4, a5, a6, a7, a8, w0, w1, w2, w3)
    call poly7_64_right(a2, a3, a4, a5, a6, a7, a8, p0, p1, p2, p3)
    qr = w0*p0 + w1*p1 + w2*p2 + w3*p3
  end subroutine weno7_64_pair

  attributes(device) subroutine weno7_32_pair(a1, a2, a3, a4, a5, a6, a7, a8, ql, qr)
    real(8), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8
    real(8), intent(out) :: ql, qr
    real(4) :: w0, w1, w2, w3, p0, p1, p2, p3
    call weights7_32_left(a1, a2, a3, a4, a5, a6, a7, w0, w1, w2, w3)
    call poly7_32_left(a1, a2, a3, a4, a5, a6, a7, p0, p1, p2, p3)
    ql = real(w0*p0 + w1*p1 + w2*p2 + w3*p3, 8)
    call weights7_32_right(a2, a3, a4, a5, a6, a7, a8, w0, w1, w2, w3)
    call poly7_32_right(a2, a3, a4, a5, a6, a7, a8, p0, p1, p2, p3)
    qr = real(w0*p0 + w1*p1 + w2*p2 + w3*p3, 8)
  end subroutine weno7_32_pair

  attributes(device) subroutine weights9_64_left(v1, v2, v3, v4, v5, v6, v7, v8, v9, w0, w1, w2, w3, w4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: w0, w1, w2, w3, w4
    real(8) :: b0, b1, b2, b3, b4, tau, invs
    real(8) :: a0, a1, a2, a3, a4, r0, r1, r2, r3, r4
    real(8) :: c0, c1, c2, c3, c4, invd, tinv, e01, e012, e34, e234
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=1.0d0/126.0d0, dd1=10.0d0/63.0d0, dd2=10.0d0/21.0d0, dd3=20.0d0/63.0d0, dd4=5.0d0/126.0d0
    ! betas scaled by 10080; the common factor cancels in tau/(beta+eps)
    b0 = v1*(45316.0d0*v1 - 417002.0d0*v2 + 729726.0d0*v3 - 576014.0d0*v4 + 172658.0d0*v5) &
         + v2*(965926.0d0*v2 - 3408792.0d0*v3 + 2716916.0d0*v4 - 822974.0d0*v5) &
         + v3*(3042786.0d0*v3 - 4924152.0d0*v4 + 1517646.0d0*v5) &
         + v4*(2041126.0d0*v4 - 1299002.0d0*v5) &
         + v5*(215836.0d0*v5)
    b1 = v2*(13816.0d0*v2 - 121742.0d0*v3 + 198426.0d0*v4 - 140474.0d0*v5 + 36158.0d0*v6) &
         + v3*(277126.0d0*v3 - 929952.0d0*v4 + 674036.0d0*v5 - 176594.0d0*v6) &
         + v4*(812586.0d0*v4 - 1223952.0d0*v5 + 330306.0d0*v6) &
         + v5*(485446.0d0*v5 - 280502.0d0*v6) &
         + v6*(45316.0d0*v6)
    b2 = v3*(13816.0d0*v3 - 102002.0d0*v4 + 135846.0d0*v5 - 77894.0d0*v6 + 16418.0d0*v7) &
         + v4*(209926.0d0*v4 - 598152.0d0*v5 + 358196.0d0*v6 - 77894.0d0*v7) &
         + v5*(462306.0d0*v5 - 598152.0d0*v6 + 135846.0d0*v7) &
         + v6*(209926.0d0*v6 - 102002.0d0*v7) &
         + v7*(13816.0d0*v7)
    b3 = v4*(45316.0d0*v4 - 280502.0d0*v5 + 330306.0d0*v6 - 176594.0d0*v7 + 36158.0d0*v8) &
         + v5*(485446.0d0*v5 - 1223952.0d0*v6 + 674036.0d0*v7 - 140474.0d0*v8) &
         + v6*(812586.0d0*v6 - 929952.0d0*v7 + 198426.0d0*v8) &
         + v7*(277126.0d0*v7 - 121742.0d0*v8) &
         + v8*(13816.0d0*v8)
    b4 = v5*(215836.0d0*v5 - 1299002.0d0*v6 + 1517646.0d0*v7 - 822974.0d0*v8 + 172658.0d0*v9) &
         + v6*(2041126.0d0*v6 - 4924152.0d0*v7 + 2716916.0d0*v8 - 576014.0d0*v9) &
         + v7*(3042786.0d0*v7 - 3408792.0d0*v8 + 729726.0d0*v9) &
         + v8*(965926.0d0*v8 - 417002.0d0*v9) &
         + v9*(45316.0d0*v9)
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
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0d0 / (a0 + a1 + a2 + a3 + a4)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
    w4 = a4 * invs
  end subroutine weights9_64_left

  attributes(device) subroutine weights9_64_right(v1, v2, v3, v4, v5, v6, v7, v8, v9, w0, w1, w2, w3, w4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: w0, w1, w2, w3, w4
    real(8) :: b0, b1, b2, b3, b4, tau, invs
    real(8) :: a0, a1, a2, a3, a4, r0, r1, r2, r3, r4
    real(8) :: c0, c1, c2, c3, c4, invd, tinv, e01, e012, e34, e234
    real(8), parameter :: eps = 1.0d-20
    real(8), parameter :: dd0=5.0d0/126.0d0, dd1=20.0d0/63.0d0, dd2=10.0d0/21.0d0, dd3=10.0d0/63.0d0, dd4=1.0d0/126.0d0
    ! betas scaled by 10080; the common factor cancels in tau/(beta+eps)
    b4 = v9*(45316.0d0*v9 - 417002.0d0*v8 + 729726.0d0*v7 - 576014.0d0*v6 + 172658.0d0*v5) &
         + v8*(965926.0d0*v8 - 3408792.0d0*v7 + 2716916.0d0*v6 - 822974.0d0*v5) &
         + v7*(3042786.0d0*v7 - 4924152.0d0*v6 + 1517646.0d0*v5) &
         + v6*(2041126.0d0*v6 - 1299002.0d0*v5) &
         + v5*(215836.0d0*v5)
    b3 = v8*(13816.0d0*v8 - 121742.0d0*v7 + 198426.0d0*v6 - 140474.0d0*v5 + 36158.0d0*v4) &
         + v7*(277126.0d0*v7 - 929952.0d0*v6 + 674036.0d0*v5 - 176594.0d0*v4) &
         + v6*(812586.0d0*v6 - 1223952.0d0*v5 + 330306.0d0*v4) &
         + v5*(485446.0d0*v5 - 280502.0d0*v4) &
         + v4*(45316.0d0*v4)
    b2 = v7*(13816.0d0*v7 - 102002.0d0*v6 + 135846.0d0*v5 - 77894.0d0*v4 + 16418.0d0*v3) &
         + v6*(209926.0d0*v6 - 598152.0d0*v5 + 358196.0d0*v4 - 77894.0d0*v3) &
         + v5*(462306.0d0*v5 - 598152.0d0*v4 + 135846.0d0*v3) &
         + v4*(209926.0d0*v4 - 102002.0d0*v3) &
         + v3*(13816.0d0*v3)
    b1 = v6*(45316.0d0*v6 - 280502.0d0*v5 + 330306.0d0*v4 - 176594.0d0*v3 + 36158.0d0*v2) &
         + v5*(485446.0d0*v5 - 1223952.0d0*v4 + 674036.0d0*v3 - 140474.0d0*v2) &
         + v4*(812586.0d0*v4 - 929952.0d0*v3 + 198426.0d0*v2) &
         + v3*(277126.0d0*v3 - 121742.0d0*v2) &
         + v2*(13816.0d0*v2)
    b0 = v5*(215836.0d0*v5 - 1299002.0d0*v4 + 1517646.0d0*v3 - 822974.0d0*v2 + 172658.0d0*v1) &
         + v4*(2041126.0d0*v4 - 4924152.0d0*v3 + 2716916.0d0*v2 - 576014.0d0*v1) &
         + v3*(3042786.0d0*v3 - 3408792.0d0*v2 + 729726.0d0*v1) &
         + v2*(965926.0d0*v2 - 417002.0d0*v1) &
         + v1*(45316.0d0*v1)
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
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0d0 / (a0 + a1 + a2 + a3 + a4)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
    w4 = a4 * invs
  end subroutine weights9_64_right

  attributes(device) subroutine weights9_32_left(v1, v2, v3, v4, v5, v6, v7, v8, v9, w0, w1, w2, w3, w4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(4), intent(out) :: w0, w1, w2, w3, w4
    real(4) :: x1, x2, x3, x4, x5, x6, x7, x8, x9
    real(4) :: b0, b1, b2, b3, b4, tau, invs
    real(4) :: a0, a1, a2, a3, a4, r0, r1, r2, r3, r4
    real(4), parameter :: eps = 1.0e-20
    real(4), parameter :: ratio_cap = 1.0e9
    real(4), parameter :: dd0=1.0_4/126.0_4, dd1=10.0_4/63.0_4, dd2=10.0_4/21.0_4, dd3=20.0_4/63.0_4, dd4=5.0_4/126.0_4
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    x8 = real(v8,4)
    x9 = real(v9,4)
    ! betas scaled by 10080; the common factor cancels in tau/(beta+eps)
    b0 = x1*(45316.0_4*x1 - 417002.0_4*x2 + 729726.0_4*x3 - 576014.0_4*x4 + 172658.0_4*x5) &
         + x2*(965926.0_4*x2 - 3408792.0_4*x3 + 2716916.0_4*x4 - 822974.0_4*x5) &
         + x3*(3042786.0_4*x3 - 4924152.0_4*x4 + 1517646.0_4*x5) &
         + x4*(2041126.0_4*x4 - 1299002.0_4*x5) &
         + x5*(215836.0_4*x5)
    b1 = x2*(13816.0_4*x2 - 121742.0_4*x3 + 198426.0_4*x4 - 140474.0_4*x5 + 36158.0_4*x6) &
         + x3*(277126.0_4*x3 - 929952.0_4*x4 + 674036.0_4*x5 - 176594.0_4*x6) &
         + x4*(812586.0_4*x4 - 1223952.0_4*x5 + 330306.0_4*x6) &
         + x5*(485446.0_4*x5 - 280502.0_4*x6) &
         + x6*(45316.0_4*x6)
    b2 = x3*(13816.0_4*x3 - 102002.0_4*x4 + 135846.0_4*x5 - 77894.0_4*x6 + 16418.0_4*x7) &
         + x4*(209926.0_4*x4 - 598152.0_4*x5 + 358196.0_4*x6 - 77894.0_4*x7) &
         + x5*(462306.0_4*x5 - 598152.0_4*x6 + 135846.0_4*x7) &
         + x6*(209926.0_4*x6 - 102002.0_4*x7) &
         + x7*(13816.0_4*x7)
    b3 = x4*(45316.0_4*x4 - 280502.0_4*x5 + 330306.0_4*x6 - 176594.0_4*x7 + 36158.0_4*x8) &
         + x5*(485446.0_4*x5 - 1223952.0_4*x6 + 674036.0_4*x7 - 140474.0_4*x8) &
         + x6*(812586.0_4*x6 - 929952.0_4*x7 + 198426.0_4*x8) &
         + x7*(277126.0_4*x7 - 121742.0_4*x8) &
         + x8*(13816.0_4*x8)
    b4 = x5*(215836.0_4*x5 - 1299002.0_4*x6 + 1517646.0_4*x7 - 822974.0_4*x8 + 172658.0_4*x9) &
         + x6*(2041126.0_4*x6 - 4924152.0_4*x7 + 2716916.0_4*x8 - 576014.0_4*x9) &
         + x7*(3042786.0_4*x7 - 3408792.0_4*x8 + 729726.0_4*x9) &
         + x8*(965926.0_4*x8 - 417002.0_4*x9) &
         + x9*(45316.0_4*x9)
    tau = abs(1.0_4*b0 + 2.0_4*b1 - 6.0_4*b2 + 2.0_4*b3 + 1.0_4*b4)
    ! NOT batched: the product of the c_i underflows to 0 in FP32 on a
    ! smooth plateau (eps**r), giving 1/0 = Inf and then 0*Inf = NaN.
    r0 = min(tau/(b0+eps), ratio_cap)
    r1 = min(tau/(b1+eps), ratio_cap)
    r2 = min(tau/(b2+eps), ratio_cap)
    r3 = min(tau/(b3+eps), ratio_cap)
    r4 = min(tau/(b4+eps), ratio_cap)
    a0 = dd0 * (1.0_4 + r0*r0)
    a1 = dd1 * (1.0_4 + r1*r1)
    a2 = dd2 * (1.0_4 + r2*r2)
    a3 = dd3 * (1.0_4 + r3*r3)
    a4 = dd4 * (1.0_4 + r4*r4)
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0_4 / (a0 + a1 + a2 + a3 + a4)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
    w4 = a4 * invs
  end subroutine weights9_32_left

  attributes(device) subroutine weights9_32_right(v1, v2, v3, v4, v5, v6, v7, v8, v9, w0, w1, w2, w3, w4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(4), intent(out) :: w0, w1, w2, w3, w4
    real(4) :: x1, x2, x3, x4, x5, x6, x7, x8, x9
    real(4) :: b0, b1, b2, b3, b4, tau, invs
    real(4) :: a0, a1, a2, a3, a4, r0, r1, r2, r3, r4
    real(4), parameter :: eps = 1.0e-20
    real(4), parameter :: ratio_cap = 1.0e9
    real(4), parameter :: dd0=5.0_4/126.0_4, dd1=20.0_4/63.0_4, dd2=10.0_4/21.0_4, dd3=10.0_4/63.0_4, dd4=1.0_4/126.0_4
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    x8 = real(v8,4)
    x9 = real(v9,4)
    ! betas scaled by 10080; the common factor cancels in tau/(beta+eps)
    b4 = x9*(45316.0_4*x9 - 417002.0_4*x8 + 729726.0_4*x7 - 576014.0_4*x6 + 172658.0_4*x5) &
         + x8*(965926.0_4*x8 - 3408792.0_4*x7 + 2716916.0_4*x6 - 822974.0_4*x5) &
         + x7*(3042786.0_4*x7 - 4924152.0_4*x6 + 1517646.0_4*x5) &
         + x6*(2041126.0_4*x6 - 1299002.0_4*x5) &
         + x5*(215836.0_4*x5)
    b3 = x8*(13816.0_4*x8 - 121742.0_4*x7 + 198426.0_4*x6 - 140474.0_4*x5 + 36158.0_4*x4) &
         + x7*(277126.0_4*x7 - 929952.0_4*x6 + 674036.0_4*x5 - 176594.0_4*x4) &
         + x6*(812586.0_4*x6 - 1223952.0_4*x5 + 330306.0_4*x4) &
         + x5*(485446.0_4*x5 - 280502.0_4*x4) &
         + x4*(45316.0_4*x4)
    b2 = x7*(13816.0_4*x7 - 102002.0_4*x6 + 135846.0_4*x5 - 77894.0_4*x4 + 16418.0_4*x3) &
         + x6*(209926.0_4*x6 - 598152.0_4*x5 + 358196.0_4*x4 - 77894.0_4*x3) &
         + x5*(462306.0_4*x5 - 598152.0_4*x4 + 135846.0_4*x3) &
         + x4*(209926.0_4*x4 - 102002.0_4*x3) &
         + x3*(13816.0_4*x3)
    b1 = x6*(45316.0_4*x6 - 280502.0_4*x5 + 330306.0_4*x4 - 176594.0_4*x3 + 36158.0_4*x2) &
         + x5*(485446.0_4*x5 - 1223952.0_4*x4 + 674036.0_4*x3 - 140474.0_4*x2) &
         + x4*(812586.0_4*x4 - 929952.0_4*x3 + 198426.0_4*x2) &
         + x3*(277126.0_4*x3 - 121742.0_4*x2) &
         + x2*(13816.0_4*x2)
    b0 = x5*(215836.0_4*x5 - 1299002.0_4*x4 + 1517646.0_4*x3 - 822974.0_4*x2 + 172658.0_4*x1) &
         + x4*(2041126.0_4*x4 - 4924152.0_4*x3 + 2716916.0_4*x2 - 576014.0_4*x1) &
         + x3*(3042786.0_4*x3 - 3408792.0_4*x2 + 729726.0_4*x1) &
         + x2*(965926.0_4*x2 - 417002.0_4*x1) &
         + x1*(45316.0_4*x1)
    tau = abs(1.0_4*b0 + 2.0_4*b1 - 6.0_4*b2 + 2.0_4*b3 + 1.0_4*b4)
    ! NOT batched: the product of the c_i underflows to 0 in FP32 on a
    ! smooth plateau (eps**r), giving 1/0 = Inf and then 0*Inf = NaN.
    r0 = min(tau/(b0+eps), ratio_cap)
    r1 = min(tau/(b1+eps), ratio_cap)
    r2 = min(tau/(b2+eps), ratio_cap)
    r3 = min(tau/(b3+eps), ratio_cap)
    r4 = min(tau/(b4+eps), ratio_cap)
    a0 = dd0 * (1.0_4 + r0*r0)
    a1 = dd1 * (1.0_4 + r1*r1)
    a2 = dd2 * (1.0_4 + r2*r2)
    a3 = dd3 * (1.0_4 + r3*r3)
    a4 = dd4 * (1.0_4 + r4*r4)
    ! normalisation IS safe to batch: a_i >= d_i so the sum is >= 1
    invs = 1.0_4 / (a0 + a1 + a2 + a3 + a4)
    w0 = a0 * invs
    w1 = a1 * invs
    w2 = a2 * invs
    w3 = a3 * invs
    w4 = a4 * invs
  end subroutine weights9_32_right

  attributes(device) subroutine poly9_64_left(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p0 = (1.0d0/5.0d0)*v1 - (21.0d0/20.0d0)*v2 + (137.0d0/60.0d0)*v3 - (163.0d0/60.0d0)*v4 + (137.0d0/60.0d0)*v5
    p1 = - (1.0d0/20.0d0)*v2 + (17.0d0/60.0d0)*v3 - (43.0d0/60.0d0)*v4 + (77.0d0/60.0d0)*v5 + (1.0d0/5.0d0)*v6
    p2 = (1.0d0/30.0d0)*v3 - (13.0d0/60.0d0)*v4 + (47.0d0/60.0d0)*v5 + (9.0d0/20.0d0)*v6 - (1.0d0/20.0d0)*v7
    p3 = - (1.0d0/20.0d0)*v4 + (9.0d0/20.0d0)*v5 + (47.0d0/60.0d0)*v6 - (13.0d0/60.0d0)*v7 + (1.0d0/30.0d0)*v8
    p4 = (1.0d0/5.0d0)*v5 + (77.0d0/60.0d0)*v6 - (43.0d0/60.0d0)*v7 + (17.0d0/60.0d0)*v8 - (1.0d0/20.0d0)*v9
  end subroutine poly9_64_left

  attributes(device) subroutine poly9_64_right(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p4 = (1.0d0/5.0d0)*v9 - (21.0d0/20.0d0)*v8 + (137.0d0/60.0d0)*v7 - (163.0d0/60.0d0)*v6 + (137.0d0/60.0d0)*v5
    p3 = - (1.0d0/20.0d0)*v8 + (17.0d0/60.0d0)*v7 - (43.0d0/60.0d0)*v6 + (77.0d0/60.0d0)*v5 + (1.0d0/5.0d0)*v4
    p2 = (1.0d0/30.0d0)*v7 - (13.0d0/60.0d0)*v6 + (47.0d0/60.0d0)*v5 + (9.0d0/20.0d0)*v4 - (1.0d0/20.0d0)*v3
    p1 = - (1.0d0/20.0d0)*v6 + (9.0d0/20.0d0)*v5 + (47.0d0/60.0d0)*v4 - (13.0d0/60.0d0)*v3 + (1.0d0/30.0d0)*v2
    p0 = (1.0d0/5.0d0)*v5 + (77.0d0/60.0d0)*v4 - (43.0d0/60.0d0)*v3 + (17.0d0/60.0d0)*v2 - (1.0d0/20.0d0)*v1
  end subroutine poly9_64_right

  attributes(device) subroutine poly9_32_left(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(4), intent(out) :: p0, p1, p2, p3, p4
    real(4) :: x1, x2, x3, x4, x5, x6, x7, x8, x9
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    x8 = real(v8,4)
    x9 = real(v9,4)
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p0 = (1.0_4/5.0_4)*x1 - (21.0_4/20.0_4)*x2 + (137.0_4/60.0_4)*x3 - (163.0_4/60.0_4)*x4 + (137.0_4/60.0_4)*x5
    p1 = - (1.0_4/20.0_4)*x2 + (17.0_4/60.0_4)*x3 - (43.0_4/60.0_4)*x4 + (77.0_4/60.0_4)*x5 + (1.0_4/5.0_4)*x6
    p2 = (1.0_4/30.0_4)*x3 - (13.0_4/60.0_4)*x4 + (47.0_4/60.0_4)*x5 + (9.0_4/20.0_4)*x6 - (1.0_4/20.0_4)*x7
    p3 = - (1.0_4/20.0_4)*x4 + (9.0_4/20.0_4)*x5 + (47.0_4/60.0_4)*x6 - (13.0_4/60.0_4)*x7 + (1.0_4/30.0_4)*x8
    p4 = (1.0_4/5.0_4)*x5 + (77.0_4/60.0_4)*x6 - (43.0_4/60.0_4)*x7 + (17.0_4/60.0_4)*x8 - (1.0_4/20.0_4)*x9
  end subroutine poly9_32_left

  attributes(device) subroutine poly9_32_right(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(4), intent(out) :: p0, p1, p2, p3, p4
    real(4) :: x1, x2, x3, x4, x5, x6, x7, x8, x9
    x1 = real(v1,4)
    x2 = real(v2,4)
    x3 = real(v3,4)
    x4 = real(v4,4)
    x5 = real(v5,4)
    x6 = real(v6,4)
    x7 = real(v7,4)
    x8 = real(v8,4)
    x9 = real(v9,4)
    ! the 1/den is folded into the coefficients, so each candidate is a
    ! division-free FMA chain and the weights/poly boundary stays clean
    p4 = (1.0_4/5.0_4)*x9 - (21.0_4/20.0_4)*x8 + (137.0_4/60.0_4)*x7 - (163.0_4/60.0_4)*x6 + (137.0_4/60.0_4)*x5
    p3 = - (1.0_4/20.0_4)*x8 + (17.0_4/60.0_4)*x7 - (43.0_4/60.0_4)*x6 + (77.0_4/60.0_4)*x5 + (1.0_4/5.0_4)*x4
    p2 = (1.0_4/30.0_4)*x7 - (13.0_4/60.0_4)*x6 + (47.0_4/60.0_4)*x5 + (9.0_4/20.0_4)*x4 - (1.0_4/20.0_4)*x3
    p1 = - (1.0_4/20.0_4)*x6 + (9.0_4/20.0_4)*x5 + (47.0_4/60.0_4)*x4 - (13.0_4/60.0_4)*x3 + (1.0_4/30.0_4)*x2
    p0 = (1.0_4/5.0_4)*x5 + (77.0_4/60.0_4)*x4 - (43.0_4/60.0_4)*x3 + (17.0_4/60.0_4)*x2 - (1.0_4/20.0_4)*x1
  end subroutine poly9_32_right

  attributes(device) subroutine weno9_64_pair(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, ql, qr)
    real(8), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
    real(8), intent(out) :: ql, qr
    real(8) :: w0, w1, w2, w3, w4, p0, p1, p2, p3, p4
    call weights9_64_left(a1, a2, a3, a4, a5, a6, a7, a8, a9, w0, w1, w2, w3, w4)
    call poly9_64_left(a1, a2, a3, a4, a5, a6, a7, a8, a9, p0, p1, p2, p3, p4)
    ql = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
    call weights9_64_right(a2, a3, a4, a5, a6, a7, a8, a9, a10, w0, w1, w2, w3, w4)
    call poly9_64_right(a2, a3, a4, a5, a6, a7, a8, a9, a10, p0, p1, p2, p3, p4)
    qr = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
  end subroutine weno9_64_pair

  attributes(device) subroutine weno9_32_pair(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, ql, qr)
    real(8), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
    real(8), intent(out) :: ql, qr
    real(4) :: w0, w1, w2, w3, w4, p0, p1, p2, p3, p4
    call weights9_32_left(a1, a2, a3, a4, a5, a6, a7, a8, a9, w0, w1, w2, w3, w4)
    call poly9_32_left(a1, a2, a3, a4, a5, a6, a7, a8, a9, p0, p1, p2, p3, p4)
    ql = real(w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4, 8)
    call weights9_32_right(a2, a3, a4, a5, a6, a7, a8, a9, a10, w0, w1, w2, w3, w4)
    call poly9_32_right(a2, a3, a4, a5, a6, a7, a8, a9, a10, p0, p1, p2, p3, p4)
    qr = real(w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4, 8)
  end subroutine weno9_32_pair

  attributes(global) subroutine init_input(n, x)
    integer, intent(in), value :: n
    real(8), intent(out), device :: x(n,3)
    integer :: i
    real(8) :: z, stepv
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n) return
    z = dble(i) / dble(max(n,1))
    if (i > n/2) then
      stepv = 1.d0
    else
      stepv = 0.d0
    endif
    x(i,1) = 1.0d0 + 0.03d0*z + 0.15d0*stepv
    x(i,2) = 0.2d0 + 0.02d0*z - 0.04d0*stepv
    x(i,3) = 1.0d0 - 0.01d0*z + 0.30d0*stepv
  end subroutine init_input

  attributes(global) subroutine weno_var_seq(n, nrepeat, mode, x, out)
    integer, intent(in), value :: n, nrepeat, mode
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: ql(3), qr(3), acc
    logical :: fp64_field
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-5) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        fp64_field = (mode == mode_fp64) .or. (mode == mode_rho64 .and. f == 1) .or. &
                     (mode == mode_u64 .and. f == 2) .or. (mode == mode_p64 .and. f == 3)
        if (fp64_field) then
          call weno5_64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), ql(f), qr(f))
        else
          call weno5_32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), ql(f), qr(f))
        endif
      enddo
      acc = acc + ql(1) + qr(1) + ql(2) + qr(2) + ql(3) + qr(3)
    enddo
    out(i,1) = ql(1) + 1.d-30*acc; out(i,2) = qr(1)
    out(i,3) = ql(2); out(i,4) = qr(2)
    out(i,5) = ql(3); out(i,6) = qr(3)
  end subroutine weno_var_seq

  attributes(global) subroutine weno_var_warp(n, nrepeat, mode, x, out)
    integer, intent(in), value :: n, nrepeat, mode
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f
    real(8), shared :: q(face_threads,6)
    logical :: fp64_field
    real(8) :: acc
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-5) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            fp64_field = (mode == mode_fp64) .or. (mode == mode_rho64 .and. f == 1) .or. &
                         (mode == mode_u64 .and. f == 2) .or. (mode == mode_p64 .and. f == 3)
            if (fp64_field) call weno5_64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                                             q(idx,2*f-1), q(idx,2*f))
          enddo
        else
          do f = 1, 3
            fp64_field = (mode == mode_fp64) .or. (mode == mode_rho64 .and. f == 1) .or. &
                         (mode == mode_u64 .and. f == 2) .or. (mode == mode_p64 .and. f == 3)
            if (.not. fp64_field) call weno5_32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                                                   q(idx,2*f-1), q(idx,2*f))
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-5) then
      acc = q(idx,1) + q(idx,2) + q(idx,3) + q(idx,4) + q(idx,5) + q(idx,6)
      out(i,1) = q(idx,1) + 1.d-30*acc; out(i,2) = q(idx,2)
      out(i,3) = q(idx,3); out(i,4) = q(idx,4)
      out(i,5) = q(idx,5); out(i,6) = q(idx,6)
    endif
  end subroutine weno_var_warp

  attributes(global) subroutine weno_weight_poly_seq(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, q(6), acc
    real(4) :: p0, p1, p2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-5) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
        call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
        q(2*f-1) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8)
        call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
        call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
        q(2*f) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_weight_poly_seq

  attributes(global) subroutine weno_weight_poly_warp(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,18)
    real(4), shared :: psh(face_threads,18)
    real(8) :: w0, w1, w2, acc
    real(4) :: p0, p1, p2
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-5) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            b = 6*(f-1)
            call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
            call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
            psh(idx,b+4) = p0; psh(idx,b+5) = p1; psh(idx,b+6) = p2
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8)
        w0 = wsh(idx,b+4); w1 = wsh(idx,b+5); w2 = wsh(idx,b+6)
        out(i,2*f) = w0*real(psh(idx,b+4),8) + w1*real(psh(idx,b+5),8) + w2*real(psh(idx,b+6),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_warp

  attributes(global) subroutine weno_weight_poly_serial_warp(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,18)
    real(4), shared :: psh(face_threads,18)
    real(8) :: w0, w1, w2, acc
    real(4) :: p0, p1, p2
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-5) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            b = 6*(f-1)
            call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
            call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
            call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
            psh(idx,b+4) = p0; psh(idx,b+5) = p1; psh(idx,b+6) = p2
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8)
        w0 = wsh(idx,b+4); w1 = wsh(idx,b+5); w2 = wsh(idx,b+6)
        out(i,2*f) = w0*real(psh(idx,b+4),8) + w1*real(psh(idx,b+5),8) + w2*real(psh(idx,b+6),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_serial_warp

  attributes(global) subroutine weno_weight_poly_warp_oncebar(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,18)
    real(4), shared :: psh(face_threads,18)
    real(8) :: w0, w1, w2, acc
    real(4) :: p0, p1, p2
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (i <= n-5) then
        if (it <= face_threads) then
          do f = 1, 3
            b = 6*(f-1)
            call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
            call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
            psh(idx,b+4) = p0; psh(idx,b+5) = p1; psh(idx,b+6) = p2
          enddo
        endif
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8)
        w0 = wsh(idx,b+4); w1 = wsh(idx,b+5); w2 = wsh(idx,b+6)
        out(i,2*f) = w0*real(psh(idx,b+4),8) + w1*real(psh(idx,b+5),8) + w2*real(psh(idx,b+6),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_warp_oncebar

  attributes(global) subroutine weno_weight_poly_serial_oncebar(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,18)
    real(4), shared :: psh(face_threads,18)
    real(8) :: w0, w1, w2, acc
    real(4) :: p0, p1, p2
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (it <= face_threads .and. i <= n-5) then
        do f = 1, 3
          b = 6*(f-1)
          call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                         wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
          call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                         wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
          psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
          call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
          psh(idx,b+4) = p0; psh(idx,b+5) = p1; psh(idx,b+6) = p2
        enddo
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8)
        w0 = wsh(idx,b+4); w1 = wsh(idx,b+5); w2 = wsh(idx,b+6)
        out(i,2*f) = w0*real(psh(idx,b+4),8) + w1*real(psh(idx,b+5),8) + w2*real(psh(idx,b+6),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_serial_oncebar

  attributes(global) subroutine weno_weight_poly_wsmem_warp(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,18)
    real(4) :: p(18), p0, p1, p2
    real(8) :: w0, w1, w2, acc
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    p = 0.0_4
    do k = 1, nrepeat
      if (i <= n-5) then
        if (it <= face_threads) then
          do f = 1, 3
            b = 6*(f-1)
            call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            p(b+1) = p0; p(b+2) = p1; p(b+3) = p2
            call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
            p(b+4) = p0; p(b+5) = p1; p(b+6) = p2
          enddo
        endif
      endif
    enddo
    call syncthreads()
    if (it > face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3)
        out(i,2*f-1) = w0*real(p(b+1),8) + w1*real(p(b+2),8) + w2*real(p(b+3),8)
        w0 = wsh(idx,b+4); w1 = wsh(idx,b+5); w2 = wsh(idx,b+6)
        out(i,2*f) = w0*real(p(b+4),8) + w1*real(p(b+5),8) + w2*real(p(b+6),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_wsmem_warp

  attributes(global) subroutine weno_weight_poly_wsmem_serial(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,18)
    real(4) :: p(18), p0, p1, p2
    real(8) :: w0, w1, w2, acc
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    p = 0.0_4
    do k = 1, nrepeat
      if (it <= face_threads .and. i <= n-5) then
        do f = 1, 3
          b = 6*(f-1)
          call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                         wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
          call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                         wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
          p(b+1) = p0; p(b+2) = p1; p(b+3) = p2
          call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
          p(b+4) = p0; p(b+5) = p1; p(b+6) = p2
        enddo
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3)
        out(i,2*f-1) = w0*real(p(b+1),8) + w1*real(p(b+2),8) + w2*real(p(b+3),8)
        w0 = wsh(idx,b+4); w1 = wsh(idx,b+5); w2 = wsh(idx,b+6)
        out(i,2*f) = w0*real(p(b+4),8) + w1*real(p(b+5),8) + w2*real(p(b+6),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_wsmem_serial

  attributes(global) subroutine weno_weight_poly_wsmem_tile2_warp(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b, t
    real(8), shared :: wsh(face_threads,18,2)
    real(4) :: p(18,2), p0, p1, p2
    real(8) :: w0, w1, w2, q(6), acc
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    p = 0.0_4
    do k = 1, nrepeat
      do t = 1, 2
        i = (blockIdx%x-1)*face_threads*2 + idx + (t-1)*face_threads
        if (i <= n-5) then
          if (it <= face_threads) then
            do f = 1, 3
              b = 6*(f-1)
              call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                             wsh(idx,b+1,t), wsh(idx,b+2,t), wsh(idx,b+3,t))
              call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                             wsh(idx,b+4,t), wsh(idx,b+5,t), wsh(idx,b+6,t))
            enddo
          else
            do f = 1, 3
              b = 6*(f-1)
              call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
              p(b+1,t) = p0; p(b+2,t) = p1; p(b+3,t) = p2
              call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
              p(b+4,t) = p0; p(b+5,t) = p1; p(b+6,t) = p2
            enddo
          endif
        endif
      enddo
    enddo
    call syncthreads()
    if (it > face_threads) then
      do t = 1, 2
        i = (blockIdx%x-1)*face_threads*2 + idx + (t-1)*face_threads
        if (i <= n-5) then
          do f = 1, 3
            b = 6*(f-1)
            w0 = wsh(idx,b+1,t); w1 = wsh(idx,b+2,t); w2 = wsh(idx,b+3,t)
            q(2*f-1) = w0*real(p(b+1,t),8) + w1*real(p(b+2,t),8) + w2*real(p(b+3,t),8)
            w0 = wsh(idx,b+4,t); w1 = wsh(idx,b+5,t); w2 = wsh(idx,b+6,t)
            q(2*f) = w0*real(p(b+4,t),8) + w1*real(p(b+5,t),8) + w2*real(p(b+6,t),8)
          enddo
          acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
          out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
          out(i,3) = q(3); out(i,4) = q(4)
          out(i,5) = q(5); out(i,6) = q(6)
        endif
      enddo
    endif
  end subroutine weno_weight_poly_wsmem_tile2_warp

  attributes(global) subroutine weno_weight_poly_wsmem_tile2_serial(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b, t
    real(8), shared :: wsh(face_threads,18,2)
    real(4) :: p(18,2), p0, p1, p2
    real(8) :: w0, w1, w2, q(6), acc
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    p = 0.0_4
    do k = 1, nrepeat
      if (it <= face_threads) then
        do t = 1, 2
          i = (blockIdx%x-1)*face_threads*2 + idx + (t-1)*face_threads
          if (i <= n-5) then
            do f = 1, 3
              b = 6*(f-1)
              call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                             wsh(idx,b+1,t), wsh(idx,b+2,t), wsh(idx,b+3,t))
              call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                             wsh(idx,b+4,t), wsh(idx,b+5,t), wsh(idx,b+6,t))
              call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
              p(b+1,t) = p0; p(b+2,t) = p1; p(b+3,t) = p2
              call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
              p(b+4,t) = p0; p(b+5,t) = p1; p(b+6,t) = p2
            enddo
          endif
        enddo
      endif
    enddo
    call syncthreads()
    if (it <= face_threads) then
      do t = 1, 2
        i = (blockIdx%x-1)*face_threads*2 + idx + (t-1)*face_threads
        if (i <= n-5) then
          do f = 1, 3
            b = 6*(f-1)
            w0 = wsh(idx,b+1,t); w1 = wsh(idx,b+2,t); w2 = wsh(idx,b+3,t)
            q(2*f-1) = w0*real(p(b+1,t),8) + w1*real(p(b+2,t),8) + w2*real(p(b+3,t),8)
            w0 = wsh(idx,b+4,t); w1 = wsh(idx,b+5,t); w2 = wsh(idx,b+6,t)
            q(2*f) = w0*real(p(b+4,t),8) + w1*real(p(b+5,t),8) + w2*real(p(b+6,t),8)
          enddo
          acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
          out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
          out(i,3) = q(3); out(i,4) = q(4)
          out(i,5) = q(5); out(i,6) = q(6)
        endif
      enddo
    endif
  end subroutine weno_weight_poly_wsmem_tile2_serial

  attributes(global) subroutine weno_weight_poly_halfwarp_serial(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, lane, warp, i, k, f
    real(8) :: w0, w1, w2, q(6), acc
    real(4) :: p0, p1, p2
    it = threadIdx%x - 1
    lane = mod(it, 32)
    warp = it / 32
    i = (blockIdx%x-1) * ((blockDim%x/32) * 16) + warp*16 + lane + 1
    if (lane >= 16 .or. i > n-5) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
        call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
        q(2*f-1) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8)
        call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
        call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
        q(2*f) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_weight_poly_halfwarp_serial

  attributes(global) subroutine weno_weight_poly_halfwarp_shfl(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, lane, face_lane, warp, i, k, f
    real(8) :: w0, w1, w2, q(6), acc
    real(4) :: p0, p1, p2, pp0, pp1, pp2
    logical :: lower
    it = threadIdx%x - 1
    lane = mod(it, 32)
    face_lane = mod(lane, 16)
    warp = it / 32
    i = (blockIdx%x-1) * ((blockDim%x/32) * 16) + warp*16 + face_lane + 1
    lower = lane < 16
    q = 0.d0
    acc = 0.d0
    do k = 1, nrepeat
      if (i <= n-5) then
        do f = 1, 3
          p0 = 0.0_4; p1 = 0.0_4; p2 = 0.0_4
          if (lower) then
            call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
          else
            call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
          endif
          pp0 = __shfl_xor(p0, 16)
          pp1 = __shfl_xor(p1, 16)
          pp2 = __shfl_xor(p2, 16)
          if (lower) q(2*f-1) = w0*real(pp0,8) + w1*real(pp1,8) + w2*real(pp2,8)

          p0 = 0.0_4; p1 = 0.0_4; p2 = 0.0_4
          if (lower) then
            call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
          else
            call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
          endif
          pp0 = __shfl_xor(p0, 16)
          pp1 = __shfl_xor(p1, 16)
          pp2 = __shfl_xor(p2, 16)
          if (lower) q(2*f) = w0*real(pp0,8) + w1*real(pp1,8) + w2*real(pp2,8)
        enddo
        if (lower) acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      endif
    enddo
    if (lower .and. i <= n-5) then
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_weight_poly_halfwarp_shfl
  ! ---- GENERATED by scratchpad/genkernels.py: WENO7/9 copies of the
  ! ---- four core modes. Structure is transcribed from the WENO5
  ! ---- originals above; only the arity differs.

  attributes(global) subroutine weno_var_seq7(n, nrepeat, mode, x, out)
    integer, intent(in), value :: n, nrepeat, mode
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: ql(3), qr(3), acc
    logical :: fp64_field
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-7) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        fp64_field = (mode == mode_fp64) .or. (mode == mode_rho64 .and. f == 1) .or. &
                     (mode == mode_u64 .and. f == 2) .or. (mode == mode_p64 .and. f == 3)
        if (fp64_field) then
          call weno7_64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), ql(f), qr(f))
        else
          call weno7_32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), ql(f), qr(f))
        endif
      enddo
      acc = acc + ql(1) + qr(1) + ql(2) + qr(2) + ql(3) + qr(3)
    enddo
    out(i,1) = ql(1) + 1.d-30*acc; out(i,2) = qr(1)
    out(i,3) = ql(2); out(i,4) = qr(2)
    out(i,5) = ql(3); out(i,6) = qr(3)
  end subroutine weno_var_seq7

  attributes(global) subroutine weno_weight_poly_seq7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, q(6), acc
    real(4) :: p0, p1, p2, p3
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-7) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), w0, w1, w2, w3)
        call poly7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p0, p1, p2, p3)
        q(2*f-1) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8) + w3*real(p3,8)
        call weights7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), w0, w1, w2, w3)
        call poly7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p0, p1, p2, p3)
        q(2*f) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8) + w3*real(p3,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_weight_poly_seq7

  attributes(global) subroutine weno_weight_poly_warp7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,24)
    real(4), shared :: psh(face_threads,24)
    real(8) :: w0, w1, w2, w3, acc
    real(4) :: p0, p1, p2, p3
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-7) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            b = 8*(f-1)
            call weights7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3), wsh(idx,b+4))
            call weights7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), &
                           wsh(idx,b+5), wsh(idx,b+6), wsh(idx,b+7), wsh(idx,b+8))
          enddo
        else
          do f = 1, 3
            b = 8*(f-1)
            call poly7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p0, p1, p2, p3)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2; psh(idx,b+4) = p3
            call poly7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p0, p1, p2, p3)
            psh(idx,b+5) = p0; psh(idx,b+6) = p1; psh(idx,b+7) = p2; psh(idx,b+8) = p3
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-7) then
      do f = 1, 3
        b = 8*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3); w3 = wsh(idx,b+4)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8) + w3*real(psh(idx,b+4),8)
        w0 = wsh(idx,b+5); w1 = wsh(idx,b+6); w2 = wsh(idx,b+7); w3 = wsh(idx,b+8)
        out(i,2*f) = w0*real(psh(idx,b+5),8) + w1*real(psh(idx,b+6),8) + w2*real(psh(idx,b+7),8) + w3*real(psh(idx,b+8),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_warp7

  attributes(global) subroutine weno_weight_poly_serial_warp7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,24)
    real(4), shared :: psh(face_threads,24)
    real(8) :: w0, w1, w2, w3, acc
    real(4) :: p0, p1, p2, p3
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-7) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            b = 8*(f-1)
            call weights7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3), wsh(idx,b+4))
            call weights7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), &
                           wsh(idx,b+5), wsh(idx,b+6), wsh(idx,b+7), wsh(idx,b+8))
            call poly7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p0, p1, p2, p3)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2; psh(idx,b+4) = p3
            call poly7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p0, p1, p2, p3)
            psh(idx,b+5) = p0; psh(idx,b+6) = p1; psh(idx,b+7) = p2; psh(idx,b+8) = p3
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-7) then
      do f = 1, 3
        b = 8*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3); w3 = wsh(idx,b+4)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8) + w3*real(psh(idx,b+4),8)
        w0 = wsh(idx,b+5); w1 = wsh(idx,b+6); w2 = wsh(idx,b+7); w3 = wsh(idx,b+8)
        out(i,2*f) = w0*real(psh(idx,b+5),8) + w1*real(psh(idx,b+6),8) + w2*real(psh(idx,b+7),8) + w3*real(psh(idx,b+8),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_serial_warp7

  attributes(global) subroutine weno_var_seq9(n, nrepeat, mode, x, out)
    integer, intent(in), value :: n, nrepeat, mode
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: ql(3), qr(3), acc
    logical :: fp64_field
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        fp64_field = (mode == mode_fp64) .or. (mode == mode_rho64 .and. f == 1) .or. &
                     (mode == mode_u64 .and. f == 2) .or. (mode == mode_p64 .and. f == 3)
        if (fp64_field) then
          call weno9_64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), ql(f), qr(f))
        else
          call weno9_32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), ql(f), qr(f))
        endif
      enddo
      acc = acc + ql(1) + qr(1) + ql(2) + qr(2) + ql(3) + qr(3)
    enddo
    out(i,1) = ql(1) + 1.d-30*acc; out(i,2) = qr(1)
    out(i,3) = ql(2); out(i,4) = qr(2)
    out(i,5) = ql(3); out(i,6) = qr(3)
  end subroutine weno_var_seq9

  attributes(global) subroutine weno_weight_poly_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, w4, q(6), acc
    real(4) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
        call poly9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8) + w3*real(p3,8) + w4*real(p4,8)
        call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
        call poly9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8) + w3*real(p3,8) + w4*real(p4,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_weight_poly_seq9

  attributes(global) subroutine weno_weight_poly_warp9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,30)
    real(4), shared :: psh(face_threads,30)
    real(8) :: w0, w1, w2, w3, w4, acc
    real(4) :: p0, p1, p2, p3, p4
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-9) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            b = 10*(f-1)
            call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3), wsh(idx,b+4), wsh(idx,b+5))
            call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), &
                           wsh(idx,b+6), wsh(idx,b+7), wsh(idx,b+8), wsh(idx,b+9), wsh(idx,b+10))
          enddo
        else
          do f = 1, 3
            b = 10*(f-1)
            call poly9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2; psh(idx,b+4) = p3; psh(idx,b+5) = p4
            call poly9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
            psh(idx,b+6) = p0; psh(idx,b+7) = p1; psh(idx,b+8) = p2; psh(idx,b+9) = p3; psh(idx,b+10) = p4
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-9) then
      do f = 1, 3
        b = 10*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3); w3 = wsh(idx,b+4); w4 = wsh(idx,b+5)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8) + w3*real(psh(idx,b+4),8) + w4*real(psh(idx,b+5),8)
        w0 = wsh(idx,b+6); w1 = wsh(idx,b+7); w2 = wsh(idx,b+8); w3 = wsh(idx,b+9); w4 = wsh(idx,b+10)
        out(i,2*f) = w0*real(psh(idx,b+6),8) + w1*real(psh(idx,b+7),8) + w2*real(psh(idx,b+8),8) + w3*real(psh(idx,b+9),8) + w4*real(psh(idx,b+10),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_warp9

  attributes(global) subroutine weno_weight_poly_serial_warp9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(8), shared :: wsh(face_threads,30)
    real(4), shared :: psh(face_threads,30)
    real(8) :: w0, w1, w2, w3, w4, acc
    real(4) :: p0, p1, p2, p3, p4
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    if (i <= n-9) then
      do k = 1, nrepeat
        if (it <= face_threads) then
          do f = 1, 3
            b = 10*(f-1)
            call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3), wsh(idx,b+4), wsh(idx,b+5))
            call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), &
                           wsh(idx,b+6), wsh(idx,b+7), wsh(idx,b+8), wsh(idx,b+9), wsh(idx,b+10))
            call poly9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2; psh(idx,b+4) = p3; psh(idx,b+5) = p4
            call poly9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
            psh(idx,b+6) = p0; psh(idx,b+7) = p1; psh(idx,b+8) = p2; psh(idx,b+9) = p3; psh(idx,b+10) = p4
          enddo
        endif
        call syncthreads()
      enddo
    endif
    if (it <= face_threads .and. i <= n-9) then
      do f = 1, 3
        b = 10*(f-1)
        w0 = wsh(idx,b+1); w1 = wsh(idx,b+2); w2 = wsh(idx,b+3); w3 = wsh(idx,b+4); w4 = wsh(idx,b+5)
        out(i,2*f-1) = w0*real(psh(idx,b+1),8) + w1*real(psh(idx,b+2),8) + w2*real(psh(idx,b+3),8) + w3*real(psh(idx,b+4),8) + w4*real(psh(idx,b+5),8)
        w0 = wsh(idx,b+6); w1 = wsh(idx,b+7); w2 = wsh(idx,b+8); w3 = wsh(idx,b+9); w4 = wsh(idx,b+10)
        out(i,2*f) = w0*real(psh(idx,b+6),8) + w1*real(psh(idx,b+7),8) + w2*real(psh(idx,b+8),8) + w3*real(psh(idx,b+9),8) + w4*real(psh(idx,b+10),8)
      enddo
      acc = out(i,1) + out(i,2) + out(i,3) + out(i,4) + out(i,5) + out(i,6)
      out(i,1) = out(i,1) + 1.d-30*acc
    endif
  end subroutine weno_weight_poly_serial_warp9
end module weno_micro_kernels

program weno_micro
  use cudafor
  use weno_micro_kernels
  implicit none
  integer :: n, nrepeat, mode, argn, ierr, nlaunch
  character(len=128) :: mode_name, arg
  real(8), allocatable, device :: x(:,:), out(:,:)
  real(8), allocatable :: h(:,:)
  type(dim3) :: b128, b256, g128, gface, gface2, ghalf

  n = 4194304
  nrepeat = 1
  nlaunch = 10
  mode_name = 'weight_poly32_seq'
  argn = command_argument_count()
  if (argn >= 1) call get_command_argument(1, mode_name)
  if (argn >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) n
  endif
  if (argn >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) nrepeat
  endif
  ! Timed launches. This repeats the LAUNCH, not the in-kernel `nrepeat` loop:
  ! nrepeat must stay 1, because at nrepeat>1 the modes that write only private
  ! registers get their loop hoisted by nvfortran while the modes that touch
  ! shared memory or a barrier do not, which silently destroys the comparison.
  if (argn >= 4) then
    call get_command_argument(4, arg)
    read(arg, *) nlaunch
  endif

  allocate(x(n,3), out(n,6), h(8,6))
  ! Poison the output so a kernel that never runs cannot masquerade as a
  ! successful run producing zeros -- this is what hid cudaErrorInvalidPtx.
  out = -1.d0
  b128 = dim3(128,1,1)
  b256 = dim3(block_threads,1,1)
  g128 = dim3((n + 127)/128,1,1)
  gface = dim3((n + face_threads - 1)/face_threads,1,1)
  gface2 = dim3((n + 2*face_threads - 1)/(2*face_threads),1,1)
  ghalf = dim3((n + 63)/64,1,1)
  call init_input<<<g128,b128>>>(n, x)
  call check_launch('init_input')

  ! Warm-up launch (JIT, caches, clock ramp), then timed launches. The harness
  ! times itself with CUDA events rather than relying on ncu/nsys: both
  ! profilers hung repeatedly on this machine, and an external profiler is a
  ! poor dependency for the one number this benchmark exists to produce.
  call launch_selected()
  call check_launch(trim(mode_name))

  block
    type(cudaEvent) :: ev0, ev1
    real(4) :: ms, ms_min, ms_sum
    integer :: il, ie
    ie = cudaEventCreate(ev0)
    ie = cudaEventCreate(ev1)
    ms_min = huge(ms_min)
    ms_sum = 0.0
    do il = 1, nlaunch
      ie = cudaEventRecord(ev0, 0)
      call launch_selected()
      ie = cudaEventRecord(ev1, 0)
      ie = cudaEventSynchronize(ev1)
      ie = cudaEventElapsedTime(ms, ev0, ev1)
      ms_min = min(ms_min, ms)
      ms_sum = ms_sum + ms
    enddo
    ie = cudaEventDestroy(ev0)
    ie = cudaEventDestroy(ev1)
    call check_launch(trim(mode_name))
    ! Clocks cannot be locked on every box, so min is the figure to compare;
    ! a mean far above it means the run throttled.
    print '(a,f14.3)', 'time_min_us=', ms_min * 1.0e3
    print '(a,f14.3)', 'time_avg_us=', (ms_sum / real(nlaunch)) * 1.0e3
    print '(a,i0)',    'nlaunch=', nlaunch
  end block


  ! Whole-array checksum. The old version summed out(1:8,:) -- 8 faces out of
  ! millions -- which could not see a NaN or a wrong value anywhere else, and
  ! is how the FP32 WENO overflow went unnoticed.
  block
    real(8), allocatable :: hall(:,:)
    real(8) :: csum
    integer :: nbad
    allocate(hall(n,6))
    hall = out
    csum = sum(hall(1:n-5,:))
    nbad = count(.not. ieee_is_finite_all(hall(1:n-5,:)))
    h = hall(1:8,:)
    print '(a)', 'mode=' // trim(mode_name)
    print '(a,i0,a,i0)', 'n=', n, ' nrepeat=', nrepeat
    print '(a,es16.8)', 'checksum=', sum(h)
    print '(a,es24.16)', 'checksum_all=', csum
    print '(a,i0)', 'nonfinite=', nbad
    if (nbad > 0) then
      print *, 'FAIL: non-finite values in output'
      error stop 4
    endif
    deallocate(hall)
  end block
  deallocate(x, out, h)

contains
  !> Launch whichever kernel `mode_name` selects. Split out of the main body so
  !> the CUDA-event timing loop can call it repeatedly; every variable it needs
  !> is host-associated from the program.
  subroutine launch_selected()
    select case (trim(mode_name))
    case ('var_fp64_seq')
      mode = mode_fp64
      call weno_var_seq<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_fp64_warp')
      mode = mode_fp64
      call weno_var_warp<<<gface,b256>>>(n, nrepeat, mode, x, out)
    case ('var_fp32_seq')
      mode = mode_fp32
      call weno_var_seq<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_fp32_warp')
      mode = mode_fp32
      call weno_var_warp<<<gface,b256>>>(n, nrepeat, mode, x, out)
    case ('var_rho64_seq')
      mode = mode_rho64
      call weno_var_seq<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_rho64_warp')
      mode = mode_rho64
      call weno_var_warp<<<gface,b256>>>(n, nrepeat, mode, x, out)
    case ('var_u64_seq')
      mode = mode_u64
      call weno_var_seq<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_u64_warp')
      mode = mode_u64
      call weno_var_warp<<<gface,b256>>>(n, nrepeat, mode, x, out)
    case ('var_p64_seq')
      mode = mode_p64
      call weno_var_seq<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_p64_warp')
      mode = mode_p64
      call weno_var_warp<<<gface,b256>>>(n, nrepeat, mode, x, out)
    case ('weight_poly32_seq')
      call weno_weight_poly_seq<<<g128,b128>>>(n, nrepeat, x, out)
    case ('weight_poly32_serial')
      call weno_weight_poly_serial_warp<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_warp')
      call weno_weight_poly_warp<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_serial_oncebar')
      call weno_weight_poly_serial_oncebar<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_warp_oncebar')
      call weno_weight_poly_warp_oncebar<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_wsmem_serial')
      call weno_weight_poly_wsmem_serial<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_wsmem_warp')
      call weno_weight_poly_wsmem_warp<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_wsmem_tile2_serial')
      call weno_weight_poly_wsmem_tile2_serial<<<gface2,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_wsmem_tile2_warp')
      call weno_weight_poly_wsmem_tile2_warp<<<gface2,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_halfwarp_serial')
      call weno_weight_poly_halfwarp_serial<<<ghalf,b128>>>(n, nrepeat, x, out)
    case ('weight_poly32_halfwarp_shfl')
      call weno_weight_poly_halfwarp_shfl<<<ghalf,b128>>>(n, nrepeat, x, out)
    case ('var_fp64_seq7')
      mode = mode_fp64
      call weno_var_seq7<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_fp32_seq7')
      mode = mode_fp32
      call weno_var_seq7<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('weight_poly32_seq7')
      call weno_weight_poly_seq7<<<g128,b128>>>(n, nrepeat, x, out)
    case ('weight_poly32_serial7')
      call weno_weight_poly_serial_warp7<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_warp7')
      call weno_weight_poly_warp7<<<gface,b256>>>(n, nrepeat, x, out)
    case ('var_fp64_seq9')
      mode = mode_fp64
      call weno_var_seq9<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('var_fp32_seq9')
      mode = mode_fp32
      call weno_var_seq9<<<g128,b128>>>(n, nrepeat, mode, x, out)
    case ('weight_poly32_seq9')
      call weno_weight_poly_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('weight_poly32_serial9')
      call weno_weight_poly_serial_warp9<<<gface,b256>>>(n, nrepeat, x, out)
    case ('weight_poly32_warp9')
      call weno_weight_poly_warp9<<<gface,b256>>>(n, nrepeat, x, out)
    case default
      print *, 'bad mode: ', trim(mode_name)
      error stop 2
    end select
  end subroutine launch_selected

  !> Catch launch-time failures. cudaDeviceSynchronize alone returns 0 when the
  !> launch itself was rejected (e.g. cudaErrorInvalidPtx=218 from a failed JIT),
  !> so nothing is enqueued and the untouched output buffer gets reported as a
  !> successful result.
  subroutine check_launch(what)
    character(len=*), intent(in) :: what
    integer :: elaunch, esync
    elaunch = cudaGetLastError()
    if (elaunch /= 0) then
      print *, 'LAUNCH FAILED for ', what, ' err=', elaunch, ' ', &
               trim(cudaGetErrorString(elaunch))
      error stop 3
    endif
    esync = cudaDeviceSynchronize()
    if (esync /= 0) then
      print *, 'SYNC FAILED for ', what, ' err=', esync, ' ', &
               trim(cudaGetErrorString(esync))
      error stop 3
    endif
  end subroutine check_launch

  !> elementwise finite test as an array (no ieee_arithmetic dependency on
  !> device-adjacent host arrays).
  function ieee_is_finite_all(a) result(m)
    real(8), intent(in) :: a(:,:)
    logical :: m(size(a,1), size(a,2))
    m = (a == a) .and. (abs(a) <= huge(1.d0))
  end function ieee_is_finite_all
end program weno_micro
