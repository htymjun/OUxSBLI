! fltflt operator bodies: EFT backend + arithmetic/comparison overloads named
! in fltflt_operator_interfaces.f90. `include`d after `contains` — textually local to
! the including module so nvfortran's ordinary inliner can inline these into
! the caller without crossing a module boundary. Needs `use cudadevice` in the
! including module (fltflt_two_prod_fma calls __fmaf_rn).

! ================================================================
! EFT Backend (private) — names match fltflt.h counterparts
! ================================================================

! fltflt_two_sum: Knuth 1969 exact split, 6 flops, no precondition.
! Written with __fadd_rn (negation is exact, so a-b == __fadd_rn(a,-b);
! __fsub_rn is not exposed to CUDA Fortran) rather than plain +/-: at -O1 and
! above nvfortran contracts `a - (r%hi - v)` style subexpressions into an FMA,
! which destroys the error-free transform and silently degrades the whole
! double-float type to about FP32 accuracy (measured: 7.2e-8 instead of 1.7e-16).
! -Mnofma also fixes it but would strip DFMA from the native FP64 paths.
attributes(device) function fltflt_two_sum(a, b) result(r)
  real(4), intent(in) :: a, b
  type(fltflt) :: r
  real(4) :: v
  r%hi = __fadd_rn(a, b)
  v    = __fadd_rn(r%hi, -a)
  r%lo = __fadd_rn(__fadd_rn(a, -__fadd_rn(r%hi, -v)), __fadd_rn(b, -v))
end function fltflt_two_sum

! fltflt_fast_two_sum: Dekker 1971, 3 flops, requires |a| >= |b|.
! Intrinsics for the same reason as fltflt_two_sum above.
attributes(device) function fltflt_fast_two_sum(a, b) result(r)
  real(4), intent(in) :: a, b
  type(fltflt) :: r
  r%hi = __fadd_rn(a, b)
  r%lo = __fadd_rn(b, -__fadd_rn(r%hi, -a))
end function fltflt_fast_two_sum

! fltflt_two_prod_fma: exact product via FMA. 1 MUL + 1 FMAF.F32.
! Uses __fmaf_rn (cudadevice intrinsic) instead of a*b+c wrapper: the
! wrapper allows the GPU optimizer to CSE a*b=hi and fold hi-hi to 0.
! __fmaf_rn is opaque and maps directly to hardware FMAF.F32.
attributes(device) function fltflt_two_prod_fma(a, b) result(r)
  real(4), intent(in) :: a, b
  type(fltflt) :: r
  r%hi = a * b
  r%lo = __fmaf_rn(a, b, -r%hi)
end function fltflt_two_prod_fma

! ================================================================
! Negation
! ================================================================

pure attributes(device) function neg_ff(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c
  c%hi = -a%hi
  c%lo = -a%lo
end function neg_ff

! ================================================================
! Addition — FPAN algorithm (Zhang & Aiken SC'25, Figure 2)
!
! ff+ff: critical-path depth 10 ops vs 13 for Thall, same 20 flops.
! Two TwoSums on hi paths run in parallel; their error words join
! the lo accumulation before the final FastTwoSum.
! ================================================================

pure attributes(device) function add_ff_ff(a, b) result(c)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: c, s, t, q
  real(4) :: st_lo
  s     = fltflt_two_sum(a%hi, b%hi)
  t     = fltflt_two_sum(a%lo, b%lo)
  q     = fltflt_fast_two_sum(s%hi, t%hi)
  st_lo = s%lo + t%lo
  c     = fltflt_fast_two_sum(q%hi, st_lo + q%lo)
end function add_ff_ff

! ff + scalar: b%lo = 0 inlined, saves one TwoSum. ~9 flops.
pure attributes(device) function add_ff_r4(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: c, s
  s    = fltflt_two_sum(a%hi, b)
  s%lo = s%lo + a%lo
  c    = fltflt_fast_two_sum(s%hi, s%lo)
end function add_ff_r4

pure attributes(device) function add_r4_ff(a, b) result(c)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = add_ff_r4(b, a)
end function add_r4_ff

pure attributes(device) function add_ff_r8(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(8),      intent(in) :: b
  type(fltflt) :: c, b_ff
  b_ff%hi = real(b, 4)
  b_ff%lo = real(b - real(b_ff%hi, 8), 4)
  c = add_ff_ff(a, b_ff)
end function add_ff_r8

pure attributes(device) function add_r8_ff(a, b) result(c)
  real(8),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = add_ff_r8(b, a)
end function add_r8_ff

! ================================================================
! Subtraction
! ================================================================

pure attributes(device) function sub_ff_ff(a, b) result(c)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: c
  c = add_ff_ff(a, neg_ff(b))
end function sub_ff_ff

pure attributes(device) function sub_ff_r4(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: c
  c = add_ff_r4(a, -b)
end function sub_ff_r4

pure attributes(device) function sub_r4_ff(a, b) result(c)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = add_ff_r4(neg_ff(b), a)
end function sub_r4_ff

pure attributes(device) function sub_ff_r8(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(8),      intent(in) :: b
  type(fltflt) :: c
  c = add_ff_r8(a, -b)
end function sub_ff_r8

pure attributes(device) function sub_r8_ff(a, b) result(c)
  real(8),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = add_ff_r8(neg_ff(b), a)
end function sub_r8_ff

! ================================================================
! Multiplication (Thall / Hida: TwoProdFMA + cross terms)
! ================================================================

! ff*ff: TwoProdFMA + cross terms. ~8 flops + 1 FMA.
attributes(device) function mul_ff_ff(a, b) result(c)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: c, p
  p    = fltflt_two_prod_fma(a%hi, b%hi)
  p%lo = __fmaf_rn(a%hi, b%lo, __fmaf_rn(a%lo, b%hi, p%lo))
  c    = fltflt_fast_two_sum(p%hi, p%lo)
end function mul_ff_ff

! ff*scalar: b%lo = 0 inlined. ~5 flops + 1 FMA.
attributes(device) function mul_ff_r4(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: c, p
  p    = fltflt_two_prod_fma(a%hi, b)
  p%lo = __fmaf_rn(a%lo, b, p%lo)
  c    = fltflt_fast_two_sum(p%hi, p%lo)
end function mul_ff_r4

attributes(device) function mul_r4_ff(a, b) result(c)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = mul_ff_r4(b, a)
end function mul_r4_ff

attributes(device) function mul_ff_r8(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(8),      intent(in) :: b
  type(fltflt) :: c, b_ff
  b_ff%hi = real(b, 4)
  b_ff%lo = real(b - real(b_ff%hi, 8), 4)
  c = mul_ff_ff(a, b_ff)
end function mul_ff_r8

attributes(device) function mul_r8_ff(a, b) result(c)
  real(8),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = mul_ff_r8(b, a)
end function mul_r8_ff

! ================================================================
! Division (Dekker 1971, one Newton-step refinement)
!
! q1 = a%hi / b%hi                  first quotient approximation
! (p%hi, p%lo) = TwoProdFMA(q1, b%hi)  exact q1*b%hi
! s = a%hi - p%hi                   residual high part
! e = s - p%lo + a%lo - q1*b%lo     full residual ≈ a - q1*b
! q2 = e / b%hi                     Newton correction
! ================================================================

attributes(device) function div_ff_ff(a, b) result(c)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: c, p
  real(4) :: q1, q2, s, e
  q1   = a%hi / b%hi
  p    = fltflt_two_prod_fma(q1, b%hi)
  s    = a%hi - p%hi
  e    = __fmaf_rn(-q1, b%lo, s - p%lo + a%lo)
  q2   = e / b%hi
  c    = fltflt_fast_two_sum(q1, q2)
end function div_ff_ff

! ff/scalar: b%lo = 0 inlined.
attributes(device) function div_ff_r4(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: c, p
  real(4) :: q1, q2, s, e
  q1   = a%hi / b
  p    = fltflt_two_prod_fma(q1, b)
  s    = a%hi - p%hi
  e    = s - p%lo + a%lo
  q2   = e / b
  c    = fltflt_fast_two_sum(q1, q2)
end function div_ff_r4

attributes(device) function div_r4_ff(a, b) result(c)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = div_ff_ff(fltflt_init(a), b)
end function div_r4_ff

attributes(device) function div_ff_r8(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(8),      intent(in) :: b
  type(fltflt) :: c, b_ff
  b_ff%hi = real(b, 4)
  b_ff%lo = real(b - real(b_ff%hi, 8), 4)
  c = div_ff_ff(a, b_ff)
end function div_ff_r8

attributes(device) function div_r8_ff(a, b) result(c)
  real(8),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c, a_ff
  a_ff%hi = real(a, 4)
  a_ff%lo = real(a - real(a_ff%hi, 8), 4)
  c = div_ff_ff(a_ff, b)
end function div_r8_ff

! ================================================================
! Comparison operators (match fltflt.h semantics exactly)
! For normalized fltflt, sign(value) == sign(hi).
! ================================================================

pure attributes(device) function eq_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  logical :: r
  r = a%hi == b%hi .and. a%lo == b%lo
end function eq_ff_ff

pure attributes(device) function eq_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  logical :: r
  r = a%hi == b .and. a%lo == 0.0
end function eq_ff_r4

pure attributes(device) function eq_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  logical :: r
  r = eq_ff_r4(b, a)
end function eq_r4_ff

! ne_*: avoid .not. — in nvfortran device code .true.=+1, so bitwise
! .not.(+1) = -2 (non-zero = "true"), making .not. .true. still "true".
pure attributes(device) function ne_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  logical :: r
  r = a%hi /= b%hi .or. a%lo /= b%lo
end function ne_ff_ff

pure attributes(device) function ne_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  logical :: r
  r = a%hi /= b .or. a%lo /= 0.0
end function ne_ff_r4

pure attributes(device) function ne_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  logical :: r
  r = a /= b%hi .or. b%lo /= 0.0
end function ne_r4_ff

pure attributes(device) function lt_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  logical :: r
  r = a%hi < b%hi .or. (a%hi == b%hi .and. a%lo < b%lo)
end function lt_ff_ff

pure attributes(device) function lt_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  logical :: r
  r = a%hi < b .or. (a%hi == b .and. a%lo < 0.0)
end function lt_ff_r4

pure attributes(device) function lt_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  logical :: r
  r = a < b%hi .or. (a == b%hi .and. b%lo > 0.0)
end function lt_r4_ff

pure attributes(device) function gt_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  logical :: r
  r = a%hi > b%hi .or. (a%hi == b%hi .and. a%lo > b%lo)
end function gt_ff_ff

pure attributes(device) function gt_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  logical :: r
  r = a%hi > b .or. (a%hi == b .and. a%lo > 0.0)
end function gt_ff_r4

pure attributes(device) function gt_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  logical :: r
  r = a > b%hi .or. (a == b%hi .and. b%lo < 0.0)
end function gt_r4_ff

pure attributes(device) function le_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  logical :: r
  r = lt_ff_ff(a, b) .or. eq_ff_ff(a, b)
end function le_ff_ff

pure attributes(device) function le_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  logical :: r
  r = lt_ff_r4(a, b) .or. eq_ff_r4(a, b)
end function le_ff_r4

pure attributes(device) function le_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  logical :: r
  r = lt_r4_ff(a, b) .or. eq_r4_ff(a, b)
end function le_r4_ff

pure attributes(device) function ge_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  logical :: r
  r = gt_ff_ff(a, b) .or. eq_ff_ff(a, b)
end function ge_ff_ff

pure attributes(device) function ge_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  logical :: r
  r = gt_ff_r4(a, b) .or. eq_ff_r4(a, b)
end function ge_ff_r4

pure attributes(device) function ge_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  logical :: r
  r = gt_r4_ff(a, b) .or. eq_r4_ff(a, b)
end function ge_r4_ff
