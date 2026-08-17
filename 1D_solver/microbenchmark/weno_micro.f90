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
!> DOUBLE-FLOAT (DF) ARM -- the fltflt library from 1D_solver/src is included
!> below, not linked, because it is written to be textually included so
!> nvfortran's ordinary inliner can reach the bodies without crossing a module
!> boundary (rationale in src/fltflt.f90). Only src/fltflt.f90 is compiled; the
!> four include files carry the interfaces and bodies. CMakeLists.txt adds
!> ../src/fltflt.f90 to the target and ../src was already on the include path.
!>
!> Verified 2026-08-17 that the error-free transforms survive this benchmark's
!> exact `-fast -Mfma ... maxregcount:128` flag set on cc89:
!> report/df_accuracy_test.f90 reports DF-vs-FP64 3.689e-14 (gate is < 1e-11)
!> and plain FP32 1.316e-06. Re-run that probe after any flag change -- FMA
!> contraction silently collapses DF to ~1e-6 and every DF timing below would
!> then be measuring the wrong arithmetic.
!>
!> TWO COSTS THE FP32 MODES ABOVE DO NOT PAY, both measured by the modes in
!> family G:
!>   1. op multiplier -- add_ff_ff is 20 flops, mul_ff_ff ~9, fltflt_fma ~30,
!>      dot4_ff ~38. On A100 (FP32:FP64 = 2:1) the break-even multiplier is
!>      ~11.6, so the library operators sit on the wrong side of it and the
!>      hand-written relaxed arm (poly*_dfr_*) exists to get under it.
!>   2. boundary conversions issue on the FP64 PIPE -- fltflt_init(real(8)) is
!>      one DADD plus two F2F, and both F2F directions use the FP64 unit. That
!>      is ~3 FP64-pipe ops per value per crossing, which is why DF over a
!>      block as small as the candidate polynomials cannot pay for itself.
module weno_micro_kernels
  use cudafor
  !> `only:` is required, not stylistic. A bare `use cudadevice` re-exports
  !> cudadevice's DEVICE-side cudaGetLastError/cudaGetErrorString through this
  !> module into `program weno_micro`, where they shadow the host generics that
  !> check_launch() needs -- the build then fails with "Illegal call from host
  !> code to device subprogram __pgi_get_last_error". These five are everything
  !> fltflt references (it calls the _rn intrinsics explicitly so that -Mfma
  !> cannot contract its error-free transforms).
  use cudadevice, only: __fadd_rn, __fmaf_rn, __shfl_down, __shfl_xor
  use fltflt
  implicit none
  include 'fltflt_operator_interfaces.f90'
  include 'fltflt_subroutines_interfaces.f90'

  integer, parameter :: face_threads = 128
  integer, parameter :: block_threads = 2 * face_threads
  integer, parameter :: mode_fp64 = 1
  integer, parameter :: mode_fp32 = 2
  integer, parameter :: mode_rho64 = 3
  integer, parameter :: mode_u64 = 4
  integer, parameter :: mode_p64 = 5

contains
  include 'fltflt_operator.f90'
  include 'fltflt_subroutines.f90'

  !> Widen a double-float back to real(8). The fltflt library deliberately does
  !> not provide this; every consumer defines its own. Identical to
  !> ff_to_r8 in 1D_solver/src/calc_slau_kernel.f90.fypp.
  pure attributes(device) function ff_to_r8(a) result(r)
    type(fltflt), intent(in) :: a
    real(8) :: r
    r = real(a%hi, 8) + real(a%lo, 8)
  end function ff_to_r8

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

  ! ------------------------------------------------------------------------
  ! DOUBLE-FLOAT (DF) candidate polynomials -- two arms, family G
  ! ------------------------------------------------------------------------
  ! Both arms scale every candidate by 60, the common denominator of the WENO9
  ! coefficients, so the multipliers become EXACT real(4) integers
  ! (12,63,137,163,3,17,43,77,2,13,47,27 -- all < 2**24). That is the same trick
  ! check_weno_order.py --gen-df uses for the solver twins, and it matters twice:
  ! a coefficient with no lo word needs no cross-term FMA, and the common factor
  ! is undone once at the end instead of five times.
  !
  ! The coefficient rows below are transcribed term-by-term from poly9_64_left /
  ! poly9_64_right above rather than mirrored by hand. Hand-mirroring the
  ! right-biased set is exactly what made v+ third order in both weights64 here
  ! and weno5z_right in the solver, with an identical instruction count and a
  ! checksum that could not see it -- see check_weno_order.py.

  !> Relaxed-DF 5-term dot: (c1*v1 + ... + c5*v5) with each v given as a hi/lo
  !> real(4) pair and each c an EXACT real(4) integer. Returns a renormalised
  !> hi/lo pair, ~2**-48 relative.
  !>
  !> 50 FP32 flops: 15 for the five error-free products, 4 to pre-sum the error
  !> stream, 28 for the four TwoSum accumulations with their errors folded in,
  !> 3 for the final FastTwoSum. Against 5 FP64 FMAs that is C = 10, versus
  !> C ~= 13 for the fltflt route (5*mul_ff_r4 + fltflt_add5 = 65 flops), which
  !> is the whole point of carrying two arms: A100 break-even is C ~= 11.6.
  !>
  !> Every add/sub is an explicit __fadd_rn. Plain `+`/`-` would let -Mfma
  !> contract the TwoSum residuals into FMAs and silently collapse this to FP32
  !> accuracy -- the same reason fltflt_operator.f90 spells them out. Note
  !> subtraction is written __fadd_rn(a, -b) because negation is exact and
  !> __fsub_rn is NOT exposed to CUDA Fortran (fltflt_operator.f90:11-16).
  attributes(device) subroutine dfr_dot5(c1, c2, c3, c4, c5, &
                                         h1, l1, h2, l2, h3, l3, h4, l4, h5, l5, rhi, rlo)
    real(4), intent(in) :: c1, c2, c3, c4, c5
    real(4), intent(in) :: h1, l1, h2, l2, h3, l3, h4, l4, h5, l5
    real(4), intent(out) :: rhi, rlo
    real(4) :: p1, p2, p3, p4, p5, e1, e2, e3, e4, e5
    real(4) :: s, t, d, v
    ! error-free products: p + e = c*v exactly, no c-lo term needed
    p1 = c1*h1; e1 = __fmaf_rn(c1, h1, -p1); e1 = __fmaf_rn(c1, l1, e1)
    p2 = c2*h2; e2 = __fmaf_rn(c2, h2, -p2); e2 = __fmaf_rn(c2, l2, e2)
    p3 = c3*h3; e3 = __fmaf_rn(c3, h3, -p3); e3 = __fmaf_rn(c3, l3, e3)
    p4 = c4*h4; e4 = __fmaf_rn(c4, h4, -p4); e4 = __fmaf_rn(c4, l4, e4)
    p5 = c5*h5; e5 = __fmaf_rn(c5, h5, -p5); e5 = __fmaf_rn(c5, l5, e5)
    t = __fadd_rn(__fadd_rn(__fadd_rn(e1, e2), __fadd_rn(e3, e4)), e5)
    ! TwoSum chain on the hi stream; each residual joins the error stream
    s = __fadd_rn(p1, p2)
    v = __fadd_rn(s, -p1)
    t = __fadd_rn(t, __fadd_rn(__fadd_rn(p1, -__fadd_rn(s, -v)), __fadd_rn(p2, -v)))
    d = __fadd_rn(s, p3)
    v = __fadd_rn(d, -s)
    t = __fadd_rn(t, __fadd_rn(__fadd_rn(s, -__fadd_rn(d, -v)), __fadd_rn(p3, -v)))
    s = d
    d = __fadd_rn(s, p4)
    v = __fadd_rn(d, -s)
    t = __fadd_rn(t, __fadd_rn(__fadd_rn(s, -__fadd_rn(d, -v)), __fadd_rn(p4, -v)))
    s = d
    d = __fadd_rn(s, p5)
    v = __fadd_rn(d, -s)
    t = __fadd_rn(t, __fadd_rn(__fadd_rn(s, -__fadd_rn(d, -v)), __fadd_rn(p5, -v)))
    s = d
    ! FastTwoSum renormalisation
    rhi = __fadd_rn(s, t)
    rlo = __fadd_rn(t, -__fadd_rn(rhi, -s))
  end subroutine dfr_dot5

  attributes(device) subroutine poly9_dfr_left(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    real(8), parameter :: inv60 = 1.0d0/60.0d0
    real(4) :: h1, h2, h3, h4, h5, h6, h7, h8, h9
    real(4) :: l1, l2, l3, l4, l5, l6, l7, l8, l9
    real(4) :: shi, slo
    ! the FP64-pipe boundary tax: one DADD + two F2F per value, and per
    ! 1D_solver/CLAUDE.md both F2F directions issue on the FP64 unit
    h1 = real(v1,4); l1 = real(v1 - real(h1,8), 4)
    h2 = real(v2,4); l2 = real(v2 - real(h2,8), 4)
    h3 = real(v3,4); l3 = real(v3 - real(h3,8), 4)
    h4 = real(v4,4); l4 = real(v4 - real(h4,8), 4)
    h5 = real(v5,4); l5 = real(v5 - real(h5,8), 4)
    h6 = real(v6,4); l6 = real(v6 - real(h6,8), 4)
    h7 = real(v7,4); l7 = real(v7 - real(h7,8), 4)
    h8 = real(v8,4); l8 = real(v8 - real(h8,8), 4)
    h9 = real(v9,4); l9 = real(v9 - real(h9,8), 4)
    call dfr_dot5(12.0_4, -63.0_4, 137.0_4, -163.0_4, 137.0_4, &
                  h1,l1, h2,l2, h3,l3, h4,l4, h5,l5, shi, slo)
    p0 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 17.0_4, -43.0_4, 77.0_4, 12.0_4, &
                  h2,l2, h3,l3, h4,l4, h5,l5, h6,l6, shi, slo)
    p1 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(2.0_4, -13.0_4, 47.0_4, 27.0_4, -3.0_4, &
                  h3,l3, h4,l4, h5,l5, h6,l6, h7,l7, shi, slo)
    p2 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 27.0_4, 47.0_4, -13.0_4, 2.0_4, &
                  h4,l4, h5,l5, h6,l6, h7,l7, h8,l8, shi, slo)
    p3 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(12.0_4, 77.0_4, -43.0_4, 17.0_4, -3.0_4, &
                  h5,l5, h6,l6, h7,l7, h8,l8, h9,l9, shi, slo)
    p4 = (real(shi,8) + real(slo,8)) * inv60
  end subroutine poly9_dfr_left

  attributes(device) subroutine poly9_dfr_right(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    real(8), parameter :: inv60 = 1.0d0/60.0d0
    real(4) :: h1, h2, h3, h4, h5, h6, h7, h8, h9
    real(4) :: l1, l2, l3, l4, l5, l6, l7, l8, l9
    real(4) :: shi, slo
    h1 = real(v1,4); l1 = real(v1 - real(h1,8), 4)
    h2 = real(v2,4); l2 = real(v2 - real(h2,8), 4)
    h3 = real(v3,4); l3 = real(v3 - real(h3,8), 4)
    h4 = real(v4,4); l4 = real(v4 - real(h4,8), 4)
    h5 = real(v5,4); l5 = real(v5 - real(h5,8), 4)
    h6 = real(v6,4); l6 = real(v6 - real(h6,8), 4)
    h7 = real(v7,4); l7 = real(v7 - real(h7,8), 4)
    h8 = real(v8,4); l8 = real(v8 - real(h8,8), 4)
    h9 = real(v9,4); l9 = real(v9 - real(h9,8), 4)
    ! mirrors poly9_64_right term for term: p4 first, on v9..v5
    call dfr_dot5(12.0_4, -63.0_4, 137.0_4, -163.0_4, 137.0_4, &
                  h9,l9, h8,l8, h7,l7, h6,l6, h5,l5, shi, slo)
    p4 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 17.0_4, -43.0_4, 77.0_4, 12.0_4, &
                  h8,l8, h7,l7, h6,l6, h5,l5, h4,l4, shi, slo)
    p3 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(2.0_4, -13.0_4, 47.0_4, 27.0_4, -3.0_4, &
                  h7,l7, h6,l6, h5,l5, h4,l4, h3,l3, shi, slo)
    p2 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 27.0_4, 47.0_4, -13.0_4, 2.0_4, &
                  h6,l6, h5,l5, h4,l4, h3,l3, h2,l2, shi, slo)
    p1 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(12.0_4, 77.0_4, -43.0_4, 17.0_4, -3.0_4, &
                  h5,l5, h4,l4, h3,l3, h2,l2, h1,l1, shi, slo)
    p0 = (real(shi,8) + real(slo,8)) * inv60
  end subroutine poly9_dfr_right

  !> ***THIS ARM IS NOT DOUBLE-FLOAT UNDER THIS BENCHMARK'S FLAGS.*** Keep it
  !> only as the diagnostic that demonstrates why, and never quote its timing as
  !> a DF datapoint -- it is doing FP32 work at FP32 accuracy.
  !>
  !> Intended as the faithful arm: the same candidates through the fltflt
  !> library operators (mul_ff_r4 = 6 flops, fltflt_add5 = 35 => 65 per
  !> candidate, C ~= 13), so the delta against poly9_dfr_* would measure what
  !> the abstraction costs. Measured instead, at -fast -Mfma -gpu=...,lto on
  !> cc89 against poly9_64_left on a smooth non-linear field:
  !>
  !>   poly9_df_left  (this)        2.1e-08 left / 1.5e-08 right   <- FP32 class
  !>   poly9_dfr_left (relaxed)     3.9e-14 left / 3.0e-14 right   <- genuine DF
  !>
  !> and the SASS shows only 14 FP32 ops per candidate where the algorithm needs
  !> ~65: the lo-word arithmetic was eliminated, not merely reordered.
  !>
  !> Bisected: the SAME source text placed in a small standalone module (with
  !> the same flags, the same includes and the same fltflt.f90 object) gives
  !> 3.1e-14, and stays correct through five overlapping candidates and the
  !> trailing FP64 scale. It only degrades inside this large module, where
  !> fltflt_init gets inlined across the LTO boundary and its exact split
  !>     lo = real(a - real(hi,8), 4)
  !> is folded to zero -- the compiler treats the real(8)->real(4)->real(8)
  !> round trip as the identity. poly9_dfr_* is immune because it writes that
  !> same split inline in the consuming routine rather than calling into the
  !> library, which is the only difference that survived the bisect.
  !>
  !> Consequence beyond this benchmark: every DF twin in the solver
  !> (weno5z_left_df, weno7z/9z_*_df, delta4_df, delta6_df) opens with a block
  !> of fltflt_init calls, so the same fold can silently demote them too. That
  !> is worth checking directly -- report/rtx4060_weno_division_reduction.md
  !> quotes DF-vs-FP64 agreement of 1.0e-09, but that is the Q.dat print floor
  !> and cannot distinguish 1e-13 from 1e-8.
  attributes(device) subroutine poly9_df_left(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    real(8), parameter :: inv60 = 1.0d0/60.0d0
    type(fltflt) :: x1, x2, x3, x4, x5, x6, x7, x8, x9
    x1 = fltflt_init(v1); x2 = fltflt_init(v2); x3 = fltflt_init(v3)
    x4 = fltflt_init(v4); x5 = fltflt_init(v5); x6 = fltflt_init(v6)
    x7 = fltflt_init(v7); x8 = fltflt_init(v8); x9 = fltflt_init(v9)
    p0 = ff_to_r8(fltflt_add5(x1 * 12.0_4, x2 * (-63.0_4), x3 * 137.0_4, &
                              x4 * (-163.0_4), x5 * 137.0_4)) * inv60
    p1 = ff_to_r8(fltflt_add5(x2 * (-3.0_4), x3 * 17.0_4, x4 * (-43.0_4), &
                              x5 * 77.0_4, x6 * 12.0_4)) * inv60
    p2 = ff_to_r8(fltflt_add5(x3 * 2.0_4, x4 * (-13.0_4), x5 * 47.0_4, &
                              x6 * 27.0_4, x7 * (-3.0_4))) * inv60
    p3 = ff_to_r8(fltflt_add5(x4 * (-3.0_4), x5 * 27.0_4, x6 * 47.0_4, &
                              x7 * (-13.0_4), x8 * 2.0_4)) * inv60
    p4 = ff_to_r8(fltflt_add5(x5 * 12.0_4, x6 * 77.0_4, x7 * (-43.0_4), &
                              x8 * 17.0_4, x9 * (-3.0_4))) * inv60
  end subroutine poly9_df_left

  attributes(device) subroutine poly9_df_right(v1, v2, v3, v4, v5, v6, v7, v8, v9, p0, p1, p2, p3, p4)
    real(8), intent(in) :: v1, v2, v3, v4, v5, v6, v7, v8, v9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    real(8), parameter :: inv60 = 1.0d0/60.0d0
    type(fltflt) :: x1, x2, x3, x4, x5, x6, x7, x8, x9
    x1 = fltflt_init(v1); x2 = fltflt_init(v2); x3 = fltflt_init(v3)
    x4 = fltflt_init(v4); x5 = fltflt_init(v5); x6 = fltflt_init(v6)
    x7 = fltflt_init(v7); x8 = fltflt_init(v8); x9 = fltflt_init(v9)
    ! mirrors poly9_64_right term for term: p4 first, on v9..v5
    p4 = ff_to_r8(fltflt_add5(x9 * 12.0_4, x8 * (-63.0_4), x7 * 137.0_4, &
                              x6 * (-163.0_4), x5 * 137.0_4)) * inv60
    p3 = ff_to_r8(fltflt_add5(x8 * (-3.0_4), x7 * 17.0_4, x6 * (-43.0_4), &
                              x5 * 77.0_4, x4 * 12.0_4)) * inv60
    p2 = ff_to_r8(fltflt_add5(x7 * 2.0_4, x6 * (-13.0_4), x5 * 47.0_4, &
                              x4 * 27.0_4, x3 * (-3.0_4))) * inv60
    p1 = ff_to_r8(fltflt_add5(x6 * (-3.0_4), x5 * 27.0_4, x4 * 47.0_4, &
                              x3 * (-13.0_4), x2 * 2.0_4)) * inv60
    p0 = ff_to_r8(fltflt_add5(x5 * 12.0_4, x4 * 77.0_4, x3 * (-43.0_4), &
                              x2 * 17.0_4, x1 * (-3.0_4))) * inv60
  end subroutine poly9_df_right

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

  ! ==========================================================================
  ! Balanced-split families (A100 pipe-balance experiments, 2026-08-17).
  !
  ! The families above answered "which half goes to FP32" with the FP64 pipe
  ! left as the bottleneck in every mode: weight_poly32_* moves only the
  ! polynomials (~12% of the per-face FP64-pipe instructions), var_fp32_*
  ! moves everything, leaving no FP64 co-issue partner at all. On A100
  ! (FP64:FP32 = 1:2) pipe-time balance wants N_FP32 = 2*N_FP64, and none of
  ! the modes above sit anywhere near that point.
  !
  !   w32_poly64_*   the REVERSE split: WENO-Z weights (betas, tau, ratios,
  !                  normalisation -- ~88% of the arithmetic) in FP32,
  !                  candidate polynomials + final combine in FP64. Order-
  !                  verified by report/check_weno_order.py --w32-split: the
  !                  design order survives down to an ~2e-8 L1 floor because
  !                  weight errors multiply O(h^r) candidate differences.
  !   w32mix*_       balance tuning: the first `nsel` of the 6 var-bias weight
  !                  calls stay FP64, the rest go FP32; polys always FP64.
  !   var{3,4,5}_*   variable-granular split: k of nv variables fully FP64,
  !                  the rest fully FP32 (nv = 5 models the 3D solver's
  !                  [rho,u,v,w,p]). The warp form needs NO shared memory and
  !                  NO barrier -- each role owns complete reconstructions and
  !                  stores its own columns of `out` directly.
  !   w64_only_*     ablation: the weights64 stream alone, no polynomials.
  !   poly32_only_*  ablation: the candidate polynomials alone, no weights.
  !                  Together these let us separate
  !                  (a) "the polynomial half is intrinsically small" from
  !                  (b) "the polynomial half is hidden under the FP64 stream".
  !                  If t(weight_poly32_seqX) - t(w64_only_seqX) << t(poly32_only_seqX),
  !                  the FP32 half is not merely light -- it is being hidden by
  !                  single-thread instruction-level overlap.
  !
  ! The serial/warp pair here fixes two known unfairnesses of the older pair:
  ! only the real(4) weights cross shared memory (9216/12288/15360 B at
  ! WENO5/7/9 vs 27648..46080 B), and the checksum acc is formed from
  ! registers instead of re-reading out() from global memory. The single
  ! barrier sits outside both the repeat loop and the tail guard, copying the
  ! _oncebar template (the in-guard barrier of the non-oncebar kernels is a
  ! divergent-barrier hazard in tail blocks).
  ! ==========================================================================

  !> nv-column twin of init_input for the var4_*/var5_* modes. Columns 1..3
  !> are bit-identical to init_input; 4 and 5 are the same piecewise-linear
  !> ramp+step class (every WENO candidate is exact on linear data, so the
  !> NaN-guard semantics of checksum_all are preserved).
  attributes(global) subroutine init_input_nv(n, nv, x)
    integer, intent(in), value :: n, nv
    real(8), intent(out), device :: x(n,nv)
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
    if (nv >= 4) x(i,4) = 0.5d0 + 0.015d0*z - 0.06d0*stepv
    if (nv >= 5) x(i,5) = 0.8d0 - 0.025d0*z + 0.10d0*stepv
  end subroutine init_input_nv

  ! ---------------- family A: w32_poly64_seq{,7,9} ----------------

  attributes(global) subroutine weno_w32_poly64_seq(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(4) :: w0, w1, w2
    real(8) :: p0, p1, p2, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-5) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
        call poly5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
        q(2*f-1) = real(w0,8)*p0 + real(w1,8)*p1 + real(w2,8)*p2
        call weights5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
        call poly5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
        q(2*f) = real(w0,8)*p0 + real(w1,8)*p1 + real(w2,8)*p2
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w32_poly64_seq

  attributes(global) subroutine weno_w32_poly64_seq7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(4) :: w0, w1, w2, w3
    real(8) :: p0, p1, p2, p3, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-7) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), w0, w1, w2, w3)
        call poly7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p0, p1, p2, p3)
        q(2*f-1) = real(w0,8)*p0 + real(w1,8)*p1 + real(w2,8)*p2 + real(w3,8)*p3
        call weights7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), w0, w1, w2, w3)
        call poly7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p0, p1, p2, p3)
        q(2*f) = real(w0,8)*p0 + real(w1,8)*p1 + real(w2,8)*p2 + real(w3,8)*p3
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w32_poly64_seq7

  attributes(global) subroutine weno_w32_poly64_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(4) :: w0, w1, w2, w3, w4
    real(8) :: p0, p1, p2, p3, p4, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
        call poly9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = real(w0,8)*p0 + real(w1,8)*p1 + real(w2,8)*p2 + real(w3,8)*p3 + real(w4,8)*p4
        call weights9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
        call poly9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = real(w0,8)*p0 + real(w1,8)*p1 + real(w2,8)*p2 + real(w3,8)*p3 + real(w4,8)*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w32_poly64_seq9

  ! ---------------- family B: w32_poly64_{serial,warp}{,7,9} ----------------
  ! FP32 role (upper warps): weights32 -> real(4) shared. FP64 role (lower
  ! warps): poly64 in registers, then combine + register acc after ONE barrier.

  attributes(global) subroutine weno_w32_poly64_warp(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(4), shared :: wsh(face_threads,18)
    real(8) :: p(18), q(6), acc
    real(4) :: w0, w1, w2
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (i <= n-5) then
        if (it <= face_threads) then
          do f = 1, 3
            b = 6*(f-1)
            call poly5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p(b+1), p(b+2), p(b+3))
            call poly5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p(b+4), p(b+5), p(b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call weights5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
            wsh(idx,b+1) = w0; wsh(idx,b+2) = w1; wsh(idx,b+3) = w2
            call weights5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
            wsh(idx,b+4) = w0; wsh(idx,b+5) = w1; wsh(idx,b+6) = w2
          enddo
        endif
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        q(2*f-1) = real(wsh(idx,b+1),8)*p(b+1) + real(wsh(idx,b+2),8)*p(b+2) + real(wsh(idx,b+3),8)*p(b+3)
        q(2*f)   = real(wsh(idx,b+4),8)*p(b+4) + real(wsh(idx,b+5),8)*p(b+5) + real(wsh(idx,b+6),8)*p(b+6)
      enddo
      acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_w32_poly64_warp

  attributes(global) subroutine weno_w32_poly64_serial(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(4), shared :: wsh(face_threads,18)
    real(8) :: p(18), q(6), acc
    real(4) :: w0, w1, w2
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (it <= face_threads .and. i <= n-5) then
        do f = 1, 3
          b = 6*(f-1)
          call poly5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p(b+1), p(b+2), p(b+3))
          call poly5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p(b+4), p(b+5), p(b+6))
          call weights5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
          wsh(idx,b+1) = w0; wsh(idx,b+2) = w1; wsh(idx,b+3) = w2
          call weights5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
          wsh(idx,b+4) = w0; wsh(idx,b+5) = w1; wsh(idx,b+6) = w2
        enddo
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-5) then
      do f = 1, 3
        b = 6*(f-1)
        q(2*f-1) = real(wsh(idx,b+1),8)*p(b+1) + real(wsh(idx,b+2),8)*p(b+2) + real(wsh(idx,b+3),8)*p(b+3)
        q(2*f)   = real(wsh(idx,b+4),8)*p(b+4) + real(wsh(idx,b+5),8)*p(b+5) + real(wsh(idx,b+6),8)*p(b+6)
      enddo
      acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_w32_poly64_serial

  attributes(global) subroutine weno_w32_poly64_warp7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(4), shared :: wsh(face_threads,24)
    real(8) :: p(24), q(6), acc
    real(4) :: w0, w1, w2, w3
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (i <= n-7) then
        if (it <= face_threads) then
          do f = 1, 3
            b = 8*(f-1)
            call poly7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p(b+1), p(b+2), p(b+3), p(b+4))
            call poly7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p(b+5), p(b+6), p(b+7), p(b+8))
          enddo
        else
          do f = 1, 3
            b = 8*(f-1)
            call weights7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), w0, w1, w2, w3)
            wsh(idx,b+1) = w0; wsh(idx,b+2) = w1; wsh(idx,b+3) = w2; wsh(idx,b+4) = w3
            call weights7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), w0, w1, w2, w3)
            wsh(idx,b+5) = w0; wsh(idx,b+6) = w1; wsh(idx,b+7) = w2; wsh(idx,b+8) = w3
          enddo
        endif
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-7) then
      do f = 1, 3
        b = 8*(f-1)
        q(2*f-1) = real(wsh(idx,b+1),8)*p(b+1) + real(wsh(idx,b+2),8)*p(b+2) + real(wsh(idx,b+3),8)*p(b+3) + real(wsh(idx,b+4),8)*p(b+4)
        q(2*f)   = real(wsh(idx,b+5),8)*p(b+5) + real(wsh(idx,b+6),8)*p(b+6) + real(wsh(idx,b+7),8)*p(b+7) + real(wsh(idx,b+8),8)*p(b+8)
      enddo
      acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_w32_poly64_warp7

  attributes(global) subroutine weno_w32_poly64_serial7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(4), shared :: wsh(face_threads,24)
    real(8) :: p(24), q(6), acc
    real(4) :: w0, w1, w2, w3
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (it <= face_threads .and. i <= n-7) then
        do f = 1, 3
          b = 8*(f-1)
          call poly7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p(b+1), p(b+2), p(b+3), p(b+4))
          call poly7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p(b+5), p(b+6), p(b+7), p(b+8))
          call weights7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), w0, w1, w2, w3)
          wsh(idx,b+1) = w0; wsh(idx,b+2) = w1; wsh(idx,b+3) = w2; wsh(idx,b+4) = w3
          call weights7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), w0, w1, w2, w3)
          wsh(idx,b+5) = w0; wsh(idx,b+6) = w1; wsh(idx,b+7) = w2; wsh(idx,b+8) = w3
        enddo
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-7) then
      do f = 1, 3
        b = 8*(f-1)
        q(2*f-1) = real(wsh(idx,b+1),8)*p(b+1) + real(wsh(idx,b+2),8)*p(b+2) + real(wsh(idx,b+3),8)*p(b+3) + real(wsh(idx,b+4),8)*p(b+4)
        q(2*f)   = real(wsh(idx,b+5),8)*p(b+5) + real(wsh(idx,b+6),8)*p(b+6) + real(wsh(idx,b+7),8)*p(b+7) + real(wsh(idx,b+8),8)*p(b+8)
      enddo
      acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_w32_poly64_serial7

  attributes(global) subroutine weno_w32_poly64_warp9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(4), shared :: wsh(face_threads,30)
    real(8) :: p(30), q(6), acc
    real(4) :: w0, w1, w2, w3, w4
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (i <= n-9) then
        if (it <= face_threads) then
          do f = 1, 3
            b = 10*(f-1)
            call poly9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p(b+1), p(b+2), p(b+3), p(b+4), p(b+5))
            call poly9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p(b+6), p(b+7), p(b+8), p(b+9), p(b+10))
          enddo
        else
          do f = 1, 3
            b = 10*(f-1)
            call weights9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
            wsh(idx,b+1) = w0; wsh(idx,b+2) = w1; wsh(idx,b+3) = w2; wsh(idx,b+4) = w3; wsh(idx,b+5) = w4
            call weights9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
            wsh(idx,b+6) = w0; wsh(idx,b+7) = w1; wsh(idx,b+8) = w2; wsh(idx,b+9) = w3; wsh(idx,b+10) = w4
          enddo
        endif
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-9) then
      do f = 1, 3
        b = 10*(f-1)
        q(2*f-1) = real(wsh(idx,b+1),8)*p(b+1) + real(wsh(idx,b+2),8)*p(b+2) + real(wsh(idx,b+3),8)*p(b+3) + real(wsh(idx,b+4),8)*p(b+4) + real(wsh(idx,b+5),8)*p(b+5)
        q(2*f)   = real(wsh(idx,b+6),8)*p(b+6) + real(wsh(idx,b+7),8)*p(b+7) + real(wsh(idx,b+8),8)*p(b+8) + real(wsh(idx,b+9),8)*p(b+9) + real(wsh(idx,b+10),8)*p(b+10)
      enddo
      acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_w32_poly64_warp9

  attributes(global) subroutine weno_w32_poly64_serial9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: it, idx, i, k, f, b
    real(4), shared :: wsh(face_threads,30)
    real(8) :: p(30), q(6), acc
    real(4) :: w0, w1, w2, w3, w4
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    do k = 1, nrepeat
      if (it <= face_threads .and. i <= n-9) then
        do f = 1, 3
          b = 10*(f-1)
          call poly9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p(b+1), p(b+2), p(b+3), p(b+4), p(b+5))
          call poly9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p(b+6), p(b+7), p(b+8), p(b+9), p(b+10))
          call weights9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
          wsh(idx,b+1) = w0; wsh(idx,b+2) = w1; wsh(idx,b+3) = w2; wsh(idx,b+4) = w3; wsh(idx,b+5) = w4
          call weights9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
          wsh(idx,b+6) = w0; wsh(idx,b+7) = w1; wsh(idx,b+8) = w2; wsh(idx,b+9) = w3; wsh(idx,b+10) = w4
        enddo
      endif
    enddo
    call syncthreads()
    if (it <= face_threads .and. i <= n-9) then
      do f = 1, 3
        b = 10*(f-1)
        q(2*f-1) = real(wsh(idx,b+1),8)*p(b+1) + real(wsh(idx,b+2),8)*p(b+2) + real(wsh(idx,b+3),8)*p(b+3) + real(wsh(idx,b+4),8)*p(b+4) + real(wsh(idx,b+5),8)*p(b+5)
        q(2*f)   = real(wsh(idx,b+6),8)*p(b+6) + real(wsh(idx,b+7),8)*p(b+7) + real(wsh(idx,b+8),8)*p(b+8) + real(wsh(idx,b+9),8)*p(b+9) + real(wsh(idx,b+10),8)*p(b+10)
      enddo
      acc = q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
      out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
      out(i,3) = q(3); out(i,4) = q(4)
      out(i,5) = q(5); out(i,6) = q(6)
    endif
  end subroutine weno_w32_poly64_serial9

  ! ---------------- family C: w32mix{h,1}_poly64_seq9 ----------------
  ! The first `nsel` of the 6 var-bias slots (1=rho left, 2=rho right, ...)
  ! keep FP64 weights; the rest use FP32 weights. Polys + combine always FP64.
  ! Branches depend only on (f, nsel), so they are warp-uniform.

  attributes(global) subroutine weno_w32mix_poly64_seq9(n, nrepeat, nsel, x, out)
    integer, intent(in), value :: n, nrepeat, nsel
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, w4
    real(4) :: s0, s1, s2, s3, s4
    real(8) :: p0, p1, p2, p3, p4, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        if (2*f-1 <= nsel) then
          call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
          q(2*f-1) = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
        else
          call weights9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), s0, s1, s2, s3, s4)
          q(2*f-1) = real(s0,8)*p0 + real(s1,8)*p1 + real(s2,8)*p2 + real(s3,8)*p3 + real(s4,8)*p4
        endif
        call poly9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        if (2*f <= nsel) then
          call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
          q(2*f) = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
        else
          call weights9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), s0, s1, s2, s3, s4)
          q(2*f) = real(s0,8)*p0 + real(s1,8)*p1 + real(s2,8)*p2 + real(s3,8)*p3 + real(s4,8)*p4
        endif
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w32mix_poly64_seq9

  ! ---------------- family D: var{3,4,5} splits, WENO9 ----------------
  ! Variables 1..kk are reconstructed fully in FP64, kk+1..nv fully in FP32.
  ! Scalar ql/qr with stores inside the loop (no runtime-indexed local array,
  ! which would spill to local memory under the runtime trip count); the first
  ! variable of each role carries the 1.d-30*acc anti-DCE term, written last.

  attributes(global) subroutine weno_varsplit_seq9(n, nrepeat, nv, kk, x, out)
    integer, intent(in), value :: n, nrepeat, nv, kk
    real(8), intent(in), device :: x(n,nv)
    real(8), intent(out), device :: out(n,2*nv)
    integer :: i, k, f
    real(8) :: ql, qr, ql1, acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    ql1 = 0.d0
    do k = 1, nrepeat
      do f = 1, nv
        if (f <= kk) then
          call weno9_64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), ql, qr)
        else
          call weno9_32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), ql, qr)
        endif
        acc = acc + ql + qr
        if (f == 1) ql1 = ql
        if (f > 1) out(i,2*f-1) = ql
        out(i,2*f) = qr
      enddo
    enddo
    out(i,1) = ql1 + 1.d-30*acc
  end subroutine weno_varsplit_seq9

  attributes(global) subroutine weno_varsplit_warp9(n, nrepeat, nv, kk, x, out)
    integer, intent(in), value :: n, nrepeat, nv, kk
    real(8), intent(in), device :: x(n,nv)
    real(8), intent(out), device :: out(n,2*nv)
    integer :: it, idx, i, k, f
    real(8) :: ql, qr, qfirst, acc
    it = threadIdx%x
    idx = mod(it-1, face_threads) + 1
    i = (blockIdx%x-1)*face_threads + idx
    ! no shared memory and no barrier anywhere, so the early return is safe
    if (i > n-9) return
    acc = 0.d0
    qfirst = 0.d0
    if (it <= face_threads) then
      do k = 1, nrepeat
        do f = 1, kk
          call weno9_64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), ql, qr)
          acc = acc + ql + qr
          if (f == 1) qfirst = ql
          if (f > 1) out(i,2*f-1) = ql
          out(i,2*f) = qr
        enddo
      enddo
      if (kk >= 1) out(i,1) = qfirst + 1.d-30*acc
    else
      do k = 1, nrepeat
        do f = kk+1, nv
          call weno9_32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), ql, qr)
          acc = acc + ql + qr
          if (f == kk+1) qfirst = ql
          if (f > kk+1) out(i,2*f-1) = ql
          out(i,2*f) = qr
        enddo
      enddo
      if (kk+1 <= nv) out(i,2*kk+1) = qfirst + 1.d-30*acc
    endif
  end subroutine weno_varsplit_warp9

  ! ---------------- family E: w64_only_seq{,7,9} ----------------
  ! The FP64 weights stream ALONE (no polynomials): the ablation partner of
  ! weight_poly32_seqX. If t(weight_poly32_seqX) - t(w64_only_seqX) is only
  ! the combine+F2F delta, the FP32 polynomial half is hidden under the FP64
  ! stream (ILP co-issue); a serialized FP32 half would add its full pipe time
  ! on top. The output is a fixed combination of the weights, so checksums are
  ! NOT comparable with any other family.

  attributes(global) subroutine weno_w64_only_seq(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-5) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights5_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
        q(2*f-1) = w0 + 2.d0*w1 + 3.d0*w2
        call weights5_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
        q(2*f) = w0 + 2.d0*w1 + 3.d0*w2
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w64_only_seq

  attributes(global) subroutine weno_w64_only_seq7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-7) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights7_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), w0, w1, w2, w3)
        q(2*f-1) = w0 + 2.d0*w1 + 3.d0*w2 + 4.d0*w3
        call weights7_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), w0, w1, w2, w3)
        q(2*f) = w0 + 2.d0*w1 + 3.d0*w2 + 4.d0*w3
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w64_only_seq7

  attributes(global) subroutine weno_w64_only_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, w4, q(6), acc
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
        q(2*f-1) = w0 + 2.d0*w1 + 3.d0*w2 + 4.d0*w3 + 5.d0*w4
        call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
        q(2*f) = w0 + 2.d0*w1 + 3.d0*w2 + 4.d0*w3 + 5.d0*w4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w64_only_seq9

  ! ---------------- family F: poly32_only_seq{,7,9} ----------------
  ! Candidate polynomials ALONE (no weights): the partner ablation to
  ! w64_only_seqX. The output is a fixed linear combination of the polynomial
  ! candidates so checksums are NOT comparable with the mixed/full families.

  attributes(global) subroutine weno_poly32_only_seq(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(4) :: p0, p1, p2
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-5) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly5_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
        q(2*f-1) = real(p0,8) + 2.d0*real(p1,8) + 3.d0*real(p2,8)
        call poly5_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
        q(2*f) = real(p0,8) + 2.d0*real(p1,8) + 3.d0*real(p2,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_poly32_only_seq

  attributes(global) subroutine weno_poly32_only_seq7(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(4) :: p0, p1, p2, p3
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-7) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly7_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), p0, p1, p2, p3)
        q(2*f-1) = real(p0,8) + 2.d0*real(p1,8) + 3.d0*real(p2,8) + 4.d0*real(p3,8)
        call poly7_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), p0, p1, p2, p3)
        q(2*f) = real(p0,8) + 2.d0*real(p1,8) + 3.d0*real(p2,8) + 4.d0*real(p3,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_poly32_only_seq7

  attributes(global) subroutine weno_poly32_only_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(4) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly9_32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = real(p0,8) + 2.d0*real(p1,8) + 3.d0*real(p2,8) + 4.d0*real(p3,8) + 5.d0*real(p4,8)
        call poly9_32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = real(p0,8) + 2.d0*real(p1,8) + 3.d0*real(p2,8) + 4.d0*real(p3,8) + 5.d0*real(p4,8)
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_poly32_only_seq9

  ! ------------------------------------------------------------------------
  ! family G: double-float (DF) hybrid -- ablations first
  ! ------------------------------------------------------------------------
  ! These two mirror weno_poly32_only_seq9 exactly, including its weight-free
  ! integer combine, so that
  !     C  =  fp32_instructions(polyXX_only_seq9) / (FP64 ops it replaces)
  ! and the FP64-pipe conversion tax both fall straight out of
  ! analyze_weno_static.py. The ordering of the sweep matters: nothing about the
  ! hybrid modes is worth building until C is a measured number, because A100
  ! break-even sits at C ~= 11.6.
  !
  ! As with every *_only_* ablation, the output is a fixed combination rather
  ! than a real reconstruction, so checksum_all is NOT comparable to any other
  ! family -- it is a launch/NaN guard here and nothing more. Correctness of the
  ! DF polynomials is gated by check_weno_order.py, not by this number.
  !
  ! Reading the static counts: analyze_weno_static.py classifies MUFU as FP32,
  ! so an FP64 division's MUFU.RCP64H lands in the fp32 column. These two
  ! kernels contain no FP64 division, so their fp32 column is clean.

  !> FP64 reference for the polynomial ablation. Needed twice: it is the
  !> denominator for C (the FP64 work the DF arms replace), and it is the
  !> correctness reference the DF arms must reproduce to ~1e-13, since
  !> poly32_only_seq9 is itself only FP32-accurate and cannot serve as one.
  attributes(global) subroutine weno_poly64_only_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(8) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
        call poly9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_poly64_only_seq9

  attributes(global) subroutine weno_polydf_only_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(8) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly9_df_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
        call poly9_df_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_polydf_only_seq9

  !> The hypothesis mode: FP64 weights, relaxed-DF candidate polynomials and
  !> combine. Same layout as weno_weight_poly_seq9 (128 threads, no shared
  !> memory, no barrier) so the only difference from that mode and from
  !> weno_var_seq9 is which arithmetic the polynomial half uses.
  !>
  !> Measured on this GPU with analyze_weno_static.py, per-thread SASS:
  !>   poly64_only_seq9    213 FP64,    0 FP32,   0 F2F
  !>   polydfr_only_seq9   140 FP64, 1404 FP32, 150 F2F
  !> so C = 1404/213 = 6.6, comfortably under A100's 11.6 break-even -- BUT the
  !> FP64-pipe cost does not fall. 140 DADD + 150 F2F = 290 FP64-pipe ops
  !> replace 213, because both F2F directions issue on the FP64 unit and the
  !> per-value split is one DADD. Moving the polynomials to DF therefore ADDS
  !> ~36% FP64-pipe work to the block it was supposed to unload. The candidate
  !> polynomials are simply too small a block to amortise a DF boundary; the
  !> weights (973 FP64 ops behind the same ~10-value boundary) are not.
  attributes(global) subroutine weno_w64_polydfr_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, w4, q(6), acc
    real(8) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
        call poly9_dfr_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
        call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
        call poly9_dfr_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w64_polydfr_seq9

  !> Splits an FP64 field into the real(4) hi/lo pair the DF polynomials want.
  !> Run once outside the timed region, which is the whole point: it moves the
  !> per-value split OFF the measured kernel's FP64 pipe. Same 8 bytes per value
  !> as the real(8) array, so this is a layout change, not a precision change.
  attributes(global) subroutine init_input_df(n, x, xhi, xlo)
    integer, intent(in), value :: n
    real(8), intent(in), device :: x(n,3)
    real(4), intent(out), device :: xhi(n,3), xlo(n,3)
    integer :: i, f
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n) return
    do f = 1, 3
      xhi(i,f) = real(x(i,f),4)
      xlo(i,f) = real(x(i,f) - real(xhi(i,f),8), 4)
    enddo
  end subroutine init_input_df

  !> Relaxed-DF 5-term dot straight from a pre-split hi/lo field: identical
  !> arithmetic to dfr_dot5, but the caller never pays the FP64-pipe split.
  attributes(device) subroutine poly9_dfr_pre(h1,l1,h2,l2,h3,l3,h4,l4,h5,l5,h6,l6,h7,l7,h8,l8,h9,l9, &
                                              p0, p1, p2, p3, p4)
    real(4), intent(in) :: h1,l1,h2,l2,h3,l3,h4,l4,h5,l5,h6,l6,h7,l7,h8,l8,h9,l9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    real(8), parameter :: inv60 = 1.0d0/60.0d0
    real(4) :: shi, slo
    call dfr_dot5(12.0_4, -63.0_4, 137.0_4, -163.0_4, 137.0_4, h1,l1, h2,l2, h3,l3, h4,l4, h5,l5, shi, slo)
    p0 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 17.0_4, -43.0_4, 77.0_4, 12.0_4, h2,l2, h3,l3, h4,l4, h5,l5, h6,l6, shi, slo)
    p1 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(2.0_4, -13.0_4, 47.0_4, 27.0_4, -3.0_4, h3,l3, h4,l4, h5,l5, h6,l6, h7,l7, shi, slo)
    p2 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 27.0_4, 47.0_4, -13.0_4, 2.0_4, h4,l4, h5,l5, h6,l6, h7,l7, h8,l8, shi, slo)
    p3 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(12.0_4, 77.0_4, -43.0_4, 17.0_4, -3.0_4, h5,l5, h6,l6, h7,l7, h8,l8, h9,l9, shi, slo)
    p4 = (real(shi,8) + real(slo,8)) * inv60
  end subroutine poly9_dfr_pre

  attributes(device) subroutine poly9_dfr_pre_r(h1,l1,h2,l2,h3,l3,h4,l4,h5,l5,h6,l6,h7,l7,h8,l8,h9,l9, &
                                                p0, p1, p2, p3, p4)
    real(4), intent(in) :: h1,l1,h2,l2,h3,l3,h4,l4,h5,l5,h6,l6,h7,l7,h8,l8,h9,l9
    real(8), intent(out) :: p0, p1, p2, p3, p4
    real(8), parameter :: inv60 = 1.0d0/60.0d0
    real(4) :: shi, slo
    ! mirrors poly9_64_right term for term: p4 first, on v9..v5
    call dfr_dot5(12.0_4, -63.0_4, 137.0_4, -163.0_4, 137.0_4, h9,l9, h8,l8, h7,l7, h6,l6, h5,l5, shi, slo)
    p4 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 17.0_4, -43.0_4, 77.0_4, 12.0_4, h8,l8, h7,l7, h6,l6, h5,l5, h4,l4, shi, slo)
    p3 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(2.0_4, -13.0_4, 47.0_4, 27.0_4, -3.0_4, h7,l7, h6,l6, h5,l5, h4,l4, h3,l3, shi, slo)
    p2 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(-3.0_4, 27.0_4, 47.0_4, -13.0_4, 2.0_4, h6,l6, h5,l5, h4,l4, h3,l3, h2,l2, shi, slo)
    p1 = (real(shi,8) + real(slo,8)) * inv60
    call dfr_dot5(12.0_4, 77.0_4, -43.0_4, 17.0_4, -3.0_4, h5,l5, h4,l4, h3,l3, h2,l2, h1,l1, shi, slo)
    p0 = (real(shi,8) + real(slo,8)) * inv60
  end subroutine poly9_dfr_pre_r

  !> Step 4 of the DF plan: the same hybrid as weno_w64_polydfr_seq9 with the
  !> input-split tax removed, isolating whether that tax (rather than the DF op
  !> multiplier C = 6.6) is what blocks the hypothesis.
  !>
  !> Costs +50% read traffic: the FP64 weights half still needs x(n,3) while the
  !> DF half reads xhi/xlo, so 302 MB -> 453 MB at nx=4194304, nvar=3. That is
  !> affordable precisely where this experiment matters -- at WENO9 the kernels
  !> run near 310 GB/s, ~16% of the A100's 1935 GB/s peak, so they are
  !> arithmetic-bound with bandwidth to spare. It would NOT be affordable at
  !> WENO5, which sits within 1.6-1.9x of the ~1497 GB/s memory floor.
  attributes(global) subroutine weno_w64_polydfr_dfin_seq9(n, nrepeat, x, xhi, xlo, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(4), intent(in), device :: xhi(n,3), xlo(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: w0, w1, w2, w3, w4, q(6), acc
    real(8) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call weights9_64_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), w0, w1, w2, w3, w4)
        call poly9_dfr_pre(xhi(i,f),xlo(i,f), xhi(i+1,f),xlo(i+1,f), xhi(i+2,f),xlo(i+2,f), &
                           xhi(i+3,f),xlo(i+3,f), xhi(i+4,f),xlo(i+4,f), xhi(i+5,f),xlo(i+5,f), &
                           xhi(i+6,f),xlo(i+6,f), xhi(i+7,f),xlo(i+7,f), xhi(i+8,f),xlo(i+8,f), &
                           p0, p1, p2, p3, p4)
        q(2*f-1) = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
        call weights9_64_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), w0, w1, w2, w3, w4)
        call poly9_dfr_pre_r(xhi(i+1,f),xlo(i+1,f), xhi(i+2,f),xlo(i+2,f), xhi(i+3,f),xlo(i+3,f), &
                             xhi(i+4,f),xlo(i+4,f), xhi(i+5,f),xlo(i+5,f), xhi(i+6,f),xlo(i+6,f), &
                             xhi(i+7,f),xlo(i+7,f), xhi(i+8,f),xlo(i+8,f), xhi(i+9,f),xlo(i+9,f), &
                             p0, p1, p2, p3, p4)
        q(2*f) = w0*p0 + w1*p1 + w2*p2 + w3*p3 + w4*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_w64_polydfr_dfin_seq9

  !> Pre-split twin of weno_polydfr_only_seq9: the polynomial-only ablation with
  !> no input-split tax, so t(polydfr_only) - t(polydfr_dfin_only) prices the tax
  !> on its own.
  attributes(global) subroutine weno_polydfr_dfin_only_seq9(n, nrepeat, xhi, xlo, out)
    integer, intent(in), value :: n, nrepeat
    real(4), intent(in), device :: xhi(n,3), xlo(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(8) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly9_dfr_pre(xhi(i,f),xlo(i,f), xhi(i+1,f),xlo(i+1,f), xhi(i+2,f),xlo(i+2,f), &
                           xhi(i+3,f),xlo(i+3,f), xhi(i+4,f),xlo(i+4,f), xhi(i+5,f),xlo(i+5,f), &
                           xhi(i+6,f),xlo(i+6,f), xhi(i+7,f),xlo(i+7,f), xhi(i+8,f),xlo(i+8,f), &
                           p0, p1, p2, p3, p4)
        q(2*f-1) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
        call poly9_dfr_pre_r(xhi(i+1,f),xlo(i+1,f), xhi(i+2,f),xlo(i+2,f), xhi(i+3,f),xlo(i+3,f), &
                             xhi(i+4,f),xlo(i+4,f), xhi(i+5,f),xlo(i+5,f), xhi(i+6,f),xlo(i+6,f), &
                             xhi(i+7,f),xlo(i+7,f), xhi(i+8,f),xlo(i+8,f), xhi(i+9,f),xlo(i+9,f), &
                             p0, p1, p2, p3, p4)
        q(2*f) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_polydfr_dfin_only_seq9

  attributes(global) subroutine weno_polydfr_only_seq9(n, nrepeat, x, out)
    integer, intent(in), value :: n, nrepeat
    real(8), intent(in), device :: x(n,3)
    real(8), intent(out), device :: out(n,6)
    integer :: i, k, f
    real(8) :: q(6), acc
    real(8) :: p0, p1, p2, p3, p4
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n-9) return
    acc = 0.d0
    do k = 1, nrepeat
      do f = 1, 3
        call poly9_dfr_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), p0, p1, p2, p3, p4)
        q(2*f-1) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
        call poly9_dfr_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), x(i+6,f), x(i+7,f), x(i+8,f), x(i+9,f), p0, p1, p2, p3, p4)
        q(2*f) = p0 + 2.d0*p1 + 3.d0*p2 + 4.d0*p3 + 5.d0*p4
      enddo
      acc = acc + q(1) + q(2) + q(3) + q(4) + q(5) + q(6)
    enddo
    out(i,1) = q(1) + 1.d-30*acc; out(i,2) = q(2)
    out(i,3) = q(3); out(i,4) = q(4)
    out(i,5) = q(5); out(i,6) = q(6)
  end subroutine weno_polydfr_only_seq9
end module weno_micro_kernels

program weno_micro
  use cudafor
  use weno_micro_kernels
  implicit none
  integer :: n, nrepeat, mode, argn, ierr, nlaunch, nvar
  character(len=128) :: mode_name, arg
  real(8), allocatable, device :: x(:,:), out(:,:)
  ! Pre-split DF-format copy of the input, allocated only for the *dfin* modes
  ! so every other mode keeps its original memory footprint and stays
  ! bit-identical and time-comparable to earlier runs.
  real(4), allocatable, device :: xhi(:,:), xlo(:,:)
  logical :: want_dfin
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

  ! var4_*/var5_* modes carry 4/5 physical variables (rho,u,v,p / rho,u,v,w,p);
  ! everything else keeps the original 3-variable layout bit-identically.
  nvar = 3
  if (len_trim(mode_name) >= 4) then
    if (mode_name(1:4) == 'var4') nvar = 4
    if (mode_name(1:4) == 'var5') nvar = 5
  endif

  allocate(x(n,nvar), out(n,2*nvar), h(8,2*nvar))
  ! Poison the output so a kernel that never runs cannot masquerade as a
  ! successful run producing zeros -- this is what hid cudaErrorInvalidPtx.
  out = -1.d0
  b128 = dim3(128,1,1)
  b256 = dim3(block_threads,1,1)
  g128 = dim3((n + 127)/128,1,1)
  gface = dim3((n + face_threads - 1)/face_threads,1,1)
  gface2 = dim3((n + 2*face_threads - 1)/(2*face_threads),1,1)
  ghalf = dim3((n + 63)/64,1,1)
  if (nvar == 3) then
    call init_input<<<g128,b128>>>(n, x)
  else
    call init_input_nv<<<g128,b128>>>(n, nvar, x)
  endif
  call check_launch('init_input')

  ! The DF-format copy is built ONCE here, outside the timed region: that is
  ! exactly what the *dfin* modes are testing, since it moves the per-value
  ! real(8)->hi/lo split off the measured kernel's FP64 pipe.
  want_dfin = index(mode_name, 'dfin') > 0
  if (want_dfin) then
    allocate(xhi(n,3), xlo(n,3))
    call init_input_df<<<g128,b128>>>(n, x, xhi, xlo)
    call check_launch('init_input_df')
  endif

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
    allocate(hall(n,2*nvar))
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
    case ('w32_poly64_seq')
      call weno_w32_poly64_seq<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w32_poly64_serial')
      call weno_w32_poly64_serial<<<gface,b256>>>(n, nrepeat, x, out)
    case ('w32_poly64_warp')
      call weno_w32_poly64_warp<<<gface,b256>>>(n, nrepeat, x, out)
    case ('w32_poly64_seq7')
      call weno_w32_poly64_seq7<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w32_poly64_serial7')
      call weno_w32_poly64_serial7<<<gface,b256>>>(n, nrepeat, x, out)
    case ('w32_poly64_warp7')
      call weno_w32_poly64_warp7<<<gface,b256>>>(n, nrepeat, x, out)
    case ('w32_poly64_seq9')
      call weno_w32_poly64_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w32_poly64_serial9')
      call weno_w32_poly64_serial9<<<gface,b256>>>(n, nrepeat, x, out)
    case ('w32_poly64_warp9')
      call weno_w32_poly64_warp9<<<gface,b256>>>(n, nrepeat, x, out)
    case ('w32mixh_poly64_seq9')
      call weno_w32mix_poly64_seq9<<<g128,b128>>>(n, nrepeat, 1, x, out)
    case ('w32mix1_poly64_seq9')
      call weno_w32mix_poly64_seq9<<<g128,b128>>>(n, nrepeat, 2, x, out)
    case ('w64_only_seq')
      call weno_w64_only_seq<<<g128,b128>>>(n, nrepeat, x, out)
    case ('poly32_only_seq')
      call weno_poly32_only_seq<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w64_only_seq7')
      call weno_w64_only_seq7<<<g128,b128>>>(n, nrepeat, x, out)
    case ('poly32_only_seq7')
      call weno_poly32_only_seq7<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w64_only_seq9')
      call weno_w64_only_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('poly32_only_seq9')
      call weno_poly32_only_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('var3_fp64_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 3, 3, x, out)
    case ('var3_fp32_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 3, 0, x, out)
    case ('var3_k1_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 3, 1, x, out)
    case ('var3_k1_warp9')
      call weno_varsplit_warp9<<<gface,b256>>>(n, nrepeat, 3, 1, x, out)
    case ('var4_fp64_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 4, 4, x, out)
    case ('var4_fp32_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 4, 0, x, out)
    case ('var4_k1_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 4, 1, x, out)
    case ('var4_k2_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 4, 2, x, out)
    case ('var4_k1_warp9')
      call weno_varsplit_warp9<<<gface,b256>>>(n, nrepeat, 4, 1, x, out)
    case ('var5_fp64_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 5, 5, x, out)
    case ('var5_fp32_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 5, 0, x, out)
    case ('var5_k1_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 5, 1, x, out)
    case ('var5_k2_seq9')
      call weno_varsplit_seq9<<<g128,b128>>>(n, nrepeat, 5, 2, x, out)
    case ('var5_k2_warp9')
      call weno_varsplit_warp9<<<gface,b256>>>(n, nrepeat, 5, 2, x, out)
    ! family G: double-float hybrid. Ablations pair with w64_only_seq9 the same
    ! way poly32_only_seq9 does; see the family header in the module.
    case ('poly64_only_seq9')
      call weno_poly64_only_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('polydf_only_seq9')
      call weno_polydf_only_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('polydfr_only_seq9')
      call weno_polydfr_only_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w64_polydfr_seq9')
      call weno_w64_polydfr_seq9<<<g128,b128>>>(n, nrepeat, x, out)
    case ('w64_polydfr_dfin_seq9')
      call weno_w64_polydfr_dfin_seq9<<<g128,b128>>>(n, nrepeat, x, xhi, xlo, out)
    case ('polydfr_dfin_only_seq9')
      call weno_polydfr_dfin_only_seq9<<<g128,b128>>>(n, nrepeat, xhi, xlo, out)
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
