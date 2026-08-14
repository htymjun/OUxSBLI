! fltflt "additional subroutine" bodies (no C++ fltflt.h equivalent for most
! of these): fma, dot products, sqrt, rounding, min/max/clamp, etc. `include`d
! after `contains`, alongside fltflt_operator.f90 (these bodies use the `+ - *`
! operators and the EFT backend, so a module including this file must also
! include fltflt_operator.f90). Needs `use cudadevice` in the including module
! (fltflt_shfl_down/xor call __shfl_down/__shfl_xor).

! ================================================================
! Warp shuffle — NOT pure (hardware warp communication)
! Both components shuffled independently.
! Uses __shfl_down/__shfl_xor (nvfortran 25.x mask-implicit API).
! ================================================================

attributes(device) function fltflt_shfl_down(val, delta) result(r)
  type(fltflt), intent(in) :: val
  integer(4),   intent(in), value :: delta
  type(fltflt) :: r
  r%hi = __shfl_down(val%hi, delta)
  r%lo = __shfl_down(val%lo, delta)
end function fltflt_shfl_down

attributes(device) function fltflt_shfl_xor(val, lane_mask) result(r)
  type(fltflt), intent(in) :: val
  integer(4),   intent(in), value :: lane_mask
  type(fltflt) :: r
  r%hi = __shfl_xor(val%hi, lane_mask)
  r%lo = __shfl_xor(val%lo, lane_mask)
end function fltflt_shfl_xor

! ================================================================
! fltflt_add_same_sign: 11 flops; only valid when a and b share
! the same sign (Handbook of FP Arithmetic, Muller et al., Alg 14.1).
! Saves 9 flops vs fltflt_add_ff_ff for sign-guaranteed inputs.
! ================================================================

pure attributes(device) function add_same_sign_ff_ff(a, b) result(c)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: c, r
  real(4) :: s
  r = fltflt_two_sum(a%hi, b%hi)
  s = (r%lo + b%lo) + a%lo
  c = fltflt_fast_two_sum(r%hi, s)
end function add_same_sign_ff_ff

! ff + r4 variant (b%lo = 0).
pure attributes(device) function add_same_sign_ff_r4(a, b) result(c)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: c, r
  real(4) :: s
  r = fltflt_two_sum(a%hi, b)
  s = r%lo + a%lo
  c = fltflt_fast_two_sum(r%hi, s)
end function add_same_sign_ff_r4

! r4 + ff variant: commutative, delegates to ff+r4.
pure attributes(device) function add_same_sign_r4_ff(a, b) result(c)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: c
  c = add_same_sign_ff_r4(b, a)
end function add_same_sign_r4_ff

! ================================================================
! fltflt_fma: a*b + c (7 overloads)
!
! More accurate and efficient than fltflt_add(fltflt_mul(a,b), c)
! because it skips one normalization step and preserves the exact
! product from TwoProdFMA all the way into the accumulation.
! ================================================================

! Full ff*ff+ff.
attributes(device) function fma_ff_ff_ff(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b, c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a%hi, b%hi)
  p%lo = __fmaf_rn(a%hi, b%lo, p%lo)
  p%lo = __fmaf_rn(a%lo, b%hi, p%lo)
  s    = fltflt_two_sum(p%hi, c%hi)
  s%lo = s%lo + p%lo
  s    = fltflt_fast_two_sum(s%hi, s%lo)
  s%lo = s%lo + c%lo
  s%lo = __fmaf_rn(a%lo, b%lo, s%lo)
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_ff_ff_ff

! ff*ff + r4: c%lo = 0, skip a%lo*b%lo term.
attributes(device) function fma_ff_ff_r4(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b
  real(4),      intent(in) :: c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a%hi, b%hi)
  p%lo = __fmaf_rn(a%hi, b%lo, p%lo)
  p%lo = __fmaf_rn(a%lo, b%hi, p%lo)
  s    = fltflt_two_sum(p%hi, c)
  s%lo = s%lo + p%lo
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_ff_ff_r4

! r4*ff + ff: a is scalar, no a%lo cross terms.
attributes(device) function fma_r4_ff_ff(a, b, c) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b, c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a, b%hi)
  p%lo = __fmaf_rn(a, b%lo, p%lo)
  s    = fltflt_two_sum(p%hi, c%hi)
  s%lo = s%lo + p%lo
  s    = fltflt_fast_two_sum(s%hi, s%lo)
  s%lo = s%lo + c%lo
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_r4_ff_ff

! ff*r4 + ff: b is scalar, commutative with fma_r4_ff_ff.
attributes(device) function fma_ff_r4_ff(a, b, c) result(r)
  type(fltflt), intent(in) :: a, c
  real(4),      intent(in) :: b
  type(fltflt) :: r
  r = fma_r4_ff_ff(b, a, c)
end function fma_ff_r4_ff

! ff*r4 + r4: both b and c are scalars.
attributes(device) function fma_ff_r4_r4(a, b, c) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b, c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a%hi, b)
  p%lo = __fmaf_rn(a%lo, b, p%lo)
  s    = fltflt_two_sum(p%hi, c)
  s%lo = s%lo + p%lo
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_ff_r4_r4

! r4*ff + r4: commutative with fma_ff_r4_r4.
attributes(device) function fma_r4_ff_r4(a, b, c) result(r)
  real(4),      intent(in) :: a, c
  type(fltflt), intent(in) :: b
  type(fltflt) :: r
  r = fma_ff_r4_r4(b, a, c)
end function fma_r4_ff_r4

! r4*r4 + ff: both a and b are scalars.
attributes(device) function fma_r4_r4_ff(a, b, c) result(r)
  real(4),      intent(in) :: a, b
  type(fltflt), intent(in) :: c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a, b)
  s    = fltflt_two_sum(p%hi, c%hi)
  s%lo = s%lo + p%lo
  s    = fltflt_fast_two_sum(s%hi, s%lo)
  s%lo = s%lo + c%lo
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_r4_r4_ff

! ================================================================
! fltflt_fma_approx: a*b + c omitting a%lo*b%lo term (~1 ULP less
! accurate than fltflt_fma, but saves one FMA in hot loops).
! Only the two ff*ff overloads differ from fltflt_fma; the rest
! are identical and delegate to fltflt_fma implementations.
! ================================================================

! ff*ff + ff (approx): no a%lo*b%lo correction.
attributes(device) function fma_approx_ff_ff_ff(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b, c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a%hi, b%hi)
  p%lo = __fmaf_rn(a%hi, b%lo, p%lo)
  p%lo = __fmaf_rn(a%lo, b%hi, p%lo)
  s    = fltflt_two_sum(p%hi, c%hi)
  s%lo = s%lo + p%lo
  s    = fltflt_fast_two_sum(s%hi, s%lo)
  s%lo = s%lo + c%lo
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_approx_ff_ff_ff

! ff*ff + r4 (approx): no a%lo*b%lo correction.
attributes(device) function fma_approx_ff_ff_r4(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b
  real(4),      intent(in) :: c
  type(fltflt) :: r, p, s
  p    = fltflt_two_prod_fma(a%hi, b%hi)
  p%lo = __fmaf_rn(a%hi, b%lo, p%lo)
  p%lo = __fmaf_rn(a%lo, b%hi, p%lo)
  s    = fltflt_two_sum(p%hi, c)
  s%lo = s%lo + p%lo
  r    = fltflt_fast_two_sum(s%hi, s%lo)
end function fma_approx_ff_ff_r4

! Remaining overloads: identical to fltflt_fma (scalar args have no lo term).
attributes(device) function fma_approx_r4_ff_ff(a, b, c) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b, c
  type(fltflt) :: r
  r = fma_r4_ff_ff(a, b, c)
end function fma_approx_r4_ff_ff

attributes(device) function fma_approx_ff_r4_ff(a, b, c) result(r)
  type(fltflt), intent(in) :: a, c
  real(4),      intent(in) :: b
  type(fltflt) :: r
  r = fma_ff_r4_ff(a, b, c)
end function fma_approx_ff_r4_ff

attributes(device) function fma_approx_ff_r4_r4(a, b, c) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b, c
  type(fltflt) :: r
  r = fma_ff_r4_r4(a, b, c)
end function fma_approx_ff_r4_r4

attributes(device) function fma_approx_r4_ff_r4(a, b, c) result(r)
  real(4),      intent(in) :: a, c
  type(fltflt), intent(in) :: b
  type(fltflt) :: r
  r = fma_r4_ff_r4(a, b, c)
end function fma_approx_r4_ff_r4

attributes(device) function fma_approx_r4_r4_ff(a, b, c) result(r)
  real(4),      intent(in) :: a, b
  type(fltflt), intent(in) :: c
  type(fltflt) :: r
  r = fma_r4_r4_ff(a, b, c)
end function fma_approx_r4_r4_ff

! ================================================================
! Multi-operand addition (Fortran-specific; no C++ equivalent)
! ================================================================

! fltflt_add3: 2 TwoSums + lo accumulation. ~19 flops vs 28 sequential.
pure attributes(device) function fltflt_add3(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b, c
  type(fltflt) :: r, s1, s2
  real(4) :: lo
  s1 = fltflt_two_sum(a%hi, b%hi)
  s2 = fltflt_two_sum(s1%hi, c%hi)
  lo = (a%lo + b%lo + c%lo) + (s1%lo + s2%lo)
  r  = fltflt_fast_two_sum(s2%hi, lo)
end function fltflt_add3

! fltflt_add4: balanced tree of 3 TwoSums. ~27 flops vs 42 sequential.
pure attributes(device) function fltflt_add4(a, b, c, d) result(r)
  type(fltflt), intent(in) :: a, b, c, d
  type(fltflt) :: r, s1, s2, s3
  real(4) :: lo
  s1 = fltflt_two_sum(a%hi, b%hi)
  s2 = fltflt_two_sum(c%hi, d%hi)
  s3 = fltflt_two_sum(s1%hi, s2%hi)
  lo = ((a%lo + b%lo) + (c%lo + d%lo)) + (s1%lo + s2%lo + s3%lo)
  r  = fltflt_fast_two_sum(s3%hi, lo)
end function fltflt_add4

! fltflt_add5: 4 TwoSums (balanced + 1). ~35 flops vs 56 sequential.
pure attributes(device) function fltflt_add5(a, b, c, d, e) result(r)
  type(fltflt), intent(in) :: a, b, c, d, e
  type(fltflt) :: r, s1, s2, s3, s4
  real(4) :: lo
  s1 = fltflt_two_sum(a%hi, b%hi)
  s2 = fltflt_two_sum(c%hi, d%hi)
  s3 = fltflt_two_sum(s1%hi, s2%hi)
  s4 = fltflt_two_sum(s3%hi, e%hi)
  lo = ((a%lo + b%lo) + (c%lo + d%lo) + e%lo) + (s1%lo + s2%lo + s3%lo + s4%lo)
  r  = fltflt_fast_two_sum(s4%hi, lo)
end function fltflt_add5

! ================================================================
! Squaring and reciprocal (Fortran-specific; no C++ equivalent)
! ================================================================

! fltflt_square: symmetric cross-term; saves one FMA vs mul_ff_ff. ~6 ops + 2 FMAs.
attributes(device) function fltflt_square(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c, p
  p    = fltflt_two_prod_fma(a%hi, a%hi)
  p%lo = __fmaf_rn(a%hi + a%hi, a%lo, p%lo)
  c    = fltflt_fast_two_sum(p%hi, p%lo)
end function fltflt_square

! fltflt_recip: 1/a with numerator {1,0} inlined; saves 1 add vs div_ff_ff.
attributes(device) function fltflt_recip(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c, p
  real(4) :: q1, q2, s, e
  q1 = 1.0 / a%hi
  p  = fltflt_two_prod_fma(q1, a%hi)
  s  = 1.0 - p%hi
  e  = __fmaf_rn(-q1, a%lo, s - p%lo)
  q2 = e / a%hi
  c  = fltflt_fast_two_sum(q1, q2)
end function fltflt_recip

! ================================================================
! Compensated dot products (Ogita-Rump-Oishi 2005)
! All inputs are real(4).  TwoProdFMA gives the exact product;
! accumulating all error words before the final fltflt_fast_two_sum
! yields a result faithful to the true dot product.
! ================================================================

! ---- dot2: a*b + c*d ----

! dot2_r4: real(4) inputs. ~8 flops + 2 FMAs.
attributes(device) function dot2_r4(a, b, c, d) result(r)
  real(4), intent(in) :: a, b, c, d
  type(fltflt) :: r, p1, p2, s
  p1 = fltflt_two_prod_fma(a, b)
  p2 = fltflt_two_prod_fma(c, d)
  s  = fltflt_two_sum(p1%hi, p2%hi)
  r  = fltflt_fast_two_sum(s%hi, s%lo + p1%lo + p2%lo)
end function dot2_r4

! dot2_r8: real(8) inputs — each converted to fltflt, then delegates to dot2_ff.
attributes(device) function dot2_r8(a, b, c, d) result(r)
  real(8), intent(in) :: a, b, c, d
  type(fltflt) :: r
  r = dot2_ff(fltflt_init(a), fltflt_init(b), fltflt_init(c), fltflt_init(d))
end function dot2_r8

! dot2_ff: fltflt inputs. TwoProd(hi,hi) + 2 cross FMAs per product. ~12 flops + 6 FMAs.
! a%lo*b%lo is omitted — beyond ~48-bit budget, same convention as mul_ff_ff.
attributes(device) function dot2_ff(a, b, c, d) result(r)
  type(fltflt), intent(in) :: a, b, c, d
  type(fltflt) :: r, p1, p2, s
  p1    = fltflt_two_prod_fma(a%hi, b%hi)
  p2    = fltflt_two_prod_fma(c%hi, d%hi)
  p1%lo = __fmaf_rn(a%hi, b%lo, p1%lo)
  p1%lo = __fmaf_rn(a%lo, b%hi, p1%lo)
  p2%lo = __fmaf_rn(c%hi, d%lo, p2%lo)
  p2%lo = __fmaf_rn(c%lo, d%hi, p2%lo)
  s     = fltflt_two_sum(p1%hi, p2%hi)
  r     = fltflt_fast_two_sum(s%hi, s%lo + p1%lo + p2%lo)
end function dot2_ff

! ---- dot3: a*b + c*d + e*f ----

! dot3_r4: real(4) inputs. ~13 flops + 3 FMAs.
attributes(device) function dot3_r4(a, b, c, d, e, f) result(r)
  real(4), intent(in) :: a, b, c, d, e, f
  type(fltflt) :: r, p1, p2, p3, s1, s2
  real(4) :: lo
  p1 = fltflt_two_prod_fma(a, b)
  p2 = fltflt_two_prod_fma(c, d)
  p3 = fltflt_two_prod_fma(e, f)
  s1 = fltflt_two_sum(p1%hi, p2%hi)
  s2 = fltflt_two_sum(s1%hi, p3%hi)
  lo = (s1%lo + s2%lo) + (p1%lo + p2%lo + p3%lo)
  r  = fltflt_fast_two_sum(s2%hi, lo)
end function dot3_r4

! dot3_r8: real(8) inputs.
attributes(device) function dot3_r8(a, b, c, d, e, f) result(r)
  real(8), intent(in) :: a, b, c, d, e, f
  type(fltflt) :: r
  r = dot3_ff(fltflt_init(a), fltflt_init(b), fltflt_init(c), &
              fltflt_init(d), fltflt_init(e), fltflt_init(f))
end function dot3_r8

! dot3_ff: fltflt inputs. ~19 flops + 9 FMAs.
attributes(device) function dot3_ff(a, b, c, d, e, f) result(r)
  type(fltflt), intent(in) :: a, b, c, d, e, f
  type(fltflt) :: r, p1, p2, p3, s1, s2
  real(4) :: lo
  p1    = fltflt_two_prod_fma(a%hi, b%hi)
  p2    = fltflt_two_prod_fma(c%hi, d%hi)
  p3    = fltflt_two_prod_fma(e%hi, f%hi)
  p1%lo = __fmaf_rn(a%hi, b%lo, p1%lo)
  p1%lo = __fmaf_rn(a%lo, b%hi, p1%lo)
  p2%lo = __fmaf_rn(c%hi, d%lo, p2%lo)
  p2%lo = __fmaf_rn(c%lo, d%hi, p2%lo)
  p3%lo = __fmaf_rn(e%hi, f%lo, p3%lo)
  p3%lo = __fmaf_rn(e%lo, f%hi, p3%lo)
  s1    = fltflt_two_sum(p1%hi, p2%hi)
  s2    = fltflt_two_sum(s1%hi, p3%hi)
  lo    = (s1%lo + s2%lo) + (p1%lo + p2%lo + p3%lo)
  r     = fltflt_fast_two_sum(s2%hi, lo)
end function dot3_ff

! ---- dot4: a*b + c*d + e*f + g*h ----

! dot4_r4: real(4) inputs. ~23 flops + 4 FMAs.
attributes(device) function dot4_r4(a, b, c, d, e, f, g, h) result(r)
  real(4), intent(in) :: a, b, c, d, e, f, g, h
  type(fltflt) :: r, p1, p2, p3, p4, s1, s2, s3
  real(4) :: lo
  p1 = fltflt_two_prod_fma(a, b)
  p2 = fltflt_two_prod_fma(c, d)
  p3 = fltflt_two_prod_fma(e, f)
  p4 = fltflt_two_prod_fma(g, h)
  s1 = fltflt_two_sum(p1%hi, p2%hi)
  s2 = fltflt_two_sum(s1%hi, p3%hi)
  s3 = fltflt_two_sum(s2%hi, p4%hi)
  lo = ((s1%lo + s2%lo) + s3%lo) + (p1%lo + p2%lo + p3%lo + p4%lo)
  r  = fltflt_fast_two_sum(s3%hi, lo)
end function dot4_r4

! dot4_r8: real(8) inputs.
attributes(device) function dot4_r8(a, b, c, d, e, f, g, h) result(r)
  real(8), intent(in) :: a, b, c, d, e, f, g, h
  type(fltflt) :: r
  r = dot4_ff(fltflt_init(a), fltflt_init(b), fltflt_init(c), fltflt_init(d), &
              fltflt_init(e), fltflt_init(f), fltflt_init(g), fltflt_init(h))
end function dot4_r8

! dot4_ff: fltflt inputs. ~26 flops + 12 FMAs.
attributes(device) function dot4_ff(a, b, c, d, e, f, g, h) result(r)
  type(fltflt), intent(in) :: a, b, c, d, e, f, g, h
  type(fltflt) :: r, p1, p2, p3, p4, s1, s2, s3
  real(4) :: lo
  p1    = fltflt_two_prod_fma(a%hi, b%hi)
  p2    = fltflt_two_prod_fma(c%hi, d%hi)
  p3    = fltflt_two_prod_fma(e%hi, f%hi)
  p4    = fltflt_two_prod_fma(g%hi, h%hi)
  p1%lo = __fmaf_rn(a%hi, b%lo, p1%lo)
  p1%lo = __fmaf_rn(a%lo, b%hi, p1%lo)
  p2%lo = __fmaf_rn(c%hi, d%lo, p2%lo)
  p2%lo = __fmaf_rn(c%lo, d%hi, p2%lo)
  p3%lo = __fmaf_rn(e%hi, f%lo, p3%lo)
  p3%lo = __fmaf_rn(e%lo, f%hi, p3%lo)
  p4%lo = __fmaf_rn(g%hi, h%lo, p4%lo)
  p4%lo = __fmaf_rn(g%lo, h%hi, p4%lo)
  s1    = fltflt_two_sum(p1%hi, p2%hi)
  s2    = fltflt_two_sum(s1%hi, p3%hi)
  s3    = fltflt_two_sum(s2%hi, p4%hi)
  lo    = ((s1%lo + s2%lo) + s3%lo) + (p1%lo + p2%lo + p3%lo + p4%lo)
  r     = fltflt_fast_two_sum(s3%hi, lo)
end function dot4_ff

! ---- dot6: a*b + c*d + e*f + g*h + i*j + k*l ----

! dot6_ff: fltflt inputs only (no r4/r8 use site needs those overloads).
! Same balanced TwoProdFMA + TwoSum-tree + single fast_two_sum shape as dot4_ff,
! extended to 6 products. ~39 flops + 18 FMAs.
attributes(device) function dot6_ff(a, b, c, d, e, f, g, h, i, j, k, l) result(r)
  type(fltflt), intent(in) :: a, b, c, d, e, f, g, h, i, j, k, l
  type(fltflt) :: r, p1, p2, p3, p4, p5, p6, s1, s2, s3, s4, s5
  real(4) :: lo
  p1    = fltflt_two_prod_fma(a%hi, b%hi)
  p2    = fltflt_two_prod_fma(c%hi, d%hi)
  p3    = fltflt_two_prod_fma(e%hi, f%hi)
  p4    = fltflt_two_prod_fma(g%hi, h%hi)
  p5    = fltflt_two_prod_fma(i%hi, j%hi)
  p6    = fltflt_two_prod_fma(k%hi, l%hi)
  p1%lo = __fmaf_rn(a%hi, b%lo, p1%lo)
  p1%lo = __fmaf_rn(a%lo, b%hi, p1%lo)
  p2%lo = __fmaf_rn(c%hi, d%lo, p2%lo)
  p2%lo = __fmaf_rn(c%lo, d%hi, p2%lo)
  p3%lo = __fmaf_rn(e%hi, f%lo, p3%lo)
  p3%lo = __fmaf_rn(e%lo, f%hi, p3%lo)
  p4%lo = __fmaf_rn(g%hi, h%lo, p4%lo)
  p4%lo = __fmaf_rn(g%lo, h%hi, p4%lo)
  p5%lo = __fmaf_rn(i%hi, j%lo, p5%lo)
  p5%lo = __fmaf_rn(i%lo, j%hi, p5%lo)
  p6%lo = __fmaf_rn(k%hi, l%lo, p6%lo)
  p6%lo = __fmaf_rn(k%lo, l%hi, p6%lo)
  s1    = fltflt_two_sum(p1%hi, p2%hi)
  s2    = fltflt_two_sum(s1%hi, p3%hi)
  s3    = fltflt_two_sum(s2%hi, p4%hi)
  s4    = fltflt_two_sum(s3%hi, p5%hi)
  s5    = fltflt_two_sum(s4%hi, p6%hi)
  lo    = (((s1%lo + s2%lo) + s3%lo) + (s4%lo + s5%lo)) + &
          ((p1%lo + p2%lo) + (p3%lo + p4%lo) + (p5%lo + p6%lo))
  r     = fltflt_fast_two_sum(s5%hi, lo)
end function dot6_ff

! ---- add_mul: (a + b) * c ----
! TwoSum(a%hi,b%hi) folds a%lo/b%lo in, then the combined sum is multiplied
! via the same TwoProdFMA + cross-term shape as mul_ff_ff -- skips fully
! normalizing (a+b) into its own fltflt before the multiply. Uses __fmaf_rn
! (not a bare fma()/x*y+z expression): a naive fma() here risks the optimizer
! recognizing sh*c%hi was already computed as p%hi and folding the "-p%hi"
! term to a CSE'd zero, silently destroying the rounding-error term -- the
! exact same hazard fltflt_two_prod_fma's header comment warns about.
attributes(device) function fltflt_add_mul(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b, c
  type(fltflt) :: r, s, p
  real(4) :: sl
  s    = fltflt_two_sum(a%hi, b%hi)
  sl   = s%lo + a%lo + b%lo
  p    = fltflt_two_prod_fma(s%hi, c%hi)
  p%lo = __fmaf_rn(s%hi, c%lo, p%lo)
  p%lo = __fmaf_rn(sl,   c%hi, p%lo)
  r    = fltflt_fast_two_sum(p%hi, p%lo)
end function fltflt_add_mul

! ---- sub_mul: (a - b) * c ----
! Same shape as fltflt_add_mul, TwoSum(a%hi, -b%hi) in place of TwoSum(a%hi, b%hi).
attributes(device) function fltflt_sub_mul(a, b, c) result(r)
  type(fltflt), intent(in) :: a, b, c
  type(fltflt) :: r, s, p
  real(4) :: sl
  s    = fltflt_two_sum(a%hi, -b%hi)
  sl   = s%lo + a%lo - b%lo
  p    = fltflt_two_prod_fma(s%hi, c%hi)
  p%lo = __fmaf_rn(s%hi, c%lo, p%lo)
  p%lo = __fmaf_rn(sl,   c%hi, p%lo)
  r    = fltflt_fast_two_sum(p%hi, p%lo)
end function fltflt_sub_mul

! ---- add_add_mul: (a + b) * (c + d) ----
! Both sums stay unnormalized (raw TwoSum hi/lo) through the multiply --
! same mul_ff_ff cross-term shape, operating on the two combined sums
! instead of two fully-normalized fltflt operands. Skips two FastTwoSum
! normalizations compared to add_ff_ff(a,b) * add_ff_ff(c,d) via mul_ff_ff.
attributes(device) function fltflt_add_add_mul(a, b, c, d) result(r)
  type(fltflt), intent(in) :: a, b, c, d
  type(fltflt) :: r, s1, s2, p
  real(4) :: s1l, s2l
  s1    = fltflt_two_sum(a%hi, b%hi)
  s1l   = s1%lo + a%lo + b%lo
  s2    = fltflt_two_sum(c%hi, d%hi)
  s2l   = s2%lo + c%lo + d%lo
  p     = fltflt_two_prod_fma(s1%hi, s2%hi)
  p%lo  = __fmaf_rn(s1%hi, s2l,  p%lo)
  p%lo  = __fmaf_rn(s1l,   s2%hi, p%lo)
  r     = fltflt_fast_two_sum(p%hi, p%lo)
end function fltflt_add_add_mul

! ================================================================
! Absolute value
! ================================================================

! fltflt_abs: branchless; multiplies both components by sign(hi). 2 flops.
! Assumes normalized input (sign(value) == sign(hi)).
pure attributes(device) function fltflt_abs(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c
  real(4) :: s
  s    = sign(1.0, a%hi)
  c%hi = a%hi * s
  c%lo = a%lo * s
end function fltflt_abs

! A + abs(A): 2*A if A >= 0, else 0. Branchless via a 0/2 mask (matches
! fltflt_abs's own branchless-sign-mask convention). Assumes normalized
! input (sign(value) == sign(hi)).
attributes(device) pure function fltflt_add_abs(a) result(r)
  type(fltflt), intent(in) :: a
  type(fltflt) :: r
  real(4) :: mask
  mask = merge(2.0, 0.0, a%hi >= 0.0)
  r%hi = a%hi * mask
  r%lo = a%lo * mask
end function fltflt_add_abs

! A - abs(A): 0 if A >= 0, else 2*A.
attributes(device) pure function fltflt_sub_abs(a) result(r)
  type(fltflt), intent(in) :: a
  type(fltflt) :: r
  real(4) :: mask
  mask = merge(2.0, 0.0, a%hi < 0.0)
  r%hi = a%hi * mask
  r%lo = a%lo * mask
end function fltflt_sub_abs

! ================================================================
! Square root
! ================================================================

! fltflt_sqrt: one Newton step after real(4) sqrt. ~9 flops + 1 FMA.
! q1 = sqrt(a%hi)              initial approximation
! (p%hi, p%lo) = TwoProdFMA(q1, q1)  exact q1^2
! s = a - q1^2                 residual in extended precision
! q2 = s / (2*q1)              Newton correction
attributes(device) function fltflt_sqrt(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c, p
  real(4) :: q1, s, q2
  if (a%hi == 0.0) then
    c%hi = 0.0;  c%lo = 0.0;  return
  end if
  q1 = sqrt(a%hi)
  p  = fltflt_two_prod_fma(q1, q1)
  s  = (a%hi - p%hi - p%lo) + a%lo
  q2 = s / (q1 + q1)
  c  = fltflt_fast_two_sum(q1, q2)
end function fltflt_sqrt

! fltflt_sqrt_fast: ~7 flops via rsqrt initial estimate.
! Slightly less precise than fltflt_sqrt (~45 bits vs ~48 bits)
! but ~5x faster on throughput-bound kernels.
! The -fast flag maps 1.0/sqrt(x) to hardware rsqrtf instruction.
attributes(device) function fltflt_sqrt_fast(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c
  real(4) :: xn, yn, residual, correction
  if (a%hi == 0.0) then
    c%hi = 0.0
    c%lo = 0.0
    return
  end if
  xn         = 1.0 / sqrt(a%hi)        ! rsqrt approx: -fast maps to rsqrtf
  yn         = a%hi * xn               ! approximate sqrt
  residual   = __fmaf_rn(-yn, yn, a%hi) + a%lo  ! a - yn^2 in extended precision
  correction = (xn * 0.5) * residual   ! Newton correction term
  c          = fltflt_fast_two_sum(yn, correction)
end function fltflt_sqrt_fast

! ================================================================
! fltflt_norm3d: sqrt(dx^2 + dy^2 + dz^2)
!
! Three exact squares accumulated with one normalization, then
! fltflt_sqrt_fast. Cross-terms (2*hi*lo) from each square are
! folded into the lo accumulation via FMA. ~39 flops.
! ================================================================

attributes(device) function fltflt_norm3d(dx, dy, dz) result(c)
  type(fltflt), intent(in) :: dx, dy, dz
  type(fltflt) :: c, px, py, pz, s, t, sum_sq
  real(4) :: lo
  px  = fltflt_two_prod_fma(dx%hi, dx%hi)
  py  = fltflt_two_prod_fma(dy%hi, dy%hi)
  pz  = fltflt_two_prod_fma(dz%hi, dz%hi)
  s   = fltflt_two_sum(px%hi, py%hi)
  t   = fltflt_two_sum(s%hi, pz%hi)
  lo  = (t%lo + s%lo) + (px%lo + py%lo + pz%lo)
  lo  = __fmaf_rn(dx%hi + dx%hi, dx%lo, lo)
  lo  = __fmaf_rn(dy%hi + dy%hi, dy%lo, lo)
  lo  = __fmaf_rn(dz%hi + dz%hi, dz%lo, lo)
  sum_sq = fltflt_fast_two_sum(t%hi, lo)
  c   = fltflt_sqrt_fast(sum_sq)
end function fltflt_norm3d

! ================================================================
! Rounding functions (branchful but correct)
! Uses Fortran intrinsics: aint (truncate toward zero), anint
! (round nearest, ties away from zero), floor, mod, abs, sign.
! ================================================================

! fltflt_round_toward_zero: largest-magnitude integer with |result| <= |a|.
pure attributes(device) function fltflt_round_toward_zero(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c
  real(4) :: hi_trunc, lo_trunc
  logical :: opp_sign
  if (abs(a%hi) < 8388608.0) then  ! 2^23: a%hi not yet an integer
    hi_trunc = aint(a%hi)
    opp_sign = (a%hi > 0.0 .and. a%lo < 0.0) .or. &
               (a%hi < 0.0 .and. a%lo > 0.0)
    if (hi_trunc == a%hi .and. opp_sign) then
      c%hi = a%hi + merge(-1.0, 1.0, a%hi > 0.0)
      c%lo = 0.0
    else
      c%hi = hi_trunc
      c%lo = 0.0
    end if
  else  ! a%hi is already an integer; refine with lo
    lo_trunc = aint(a%lo)
    opp_sign = (a%hi > 0.0 .and. a%lo < 0.0) .or. &
               (a%hi < 0.0 .and. a%lo > 0.0)
    if (lo_trunc /= a%lo .and. opp_sign) then
      lo_trunc = lo_trunc + merge(-1.0, 1.0, a%hi > 0.0)
    end if
    c = fltflt_fast_two_sum(a%hi, lo_trunc)
  end if
end function fltflt_round_toward_zero

! fltflt_floor: largest integer not greater than a.
pure attributes(device) function fltflt_floor(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c
  real(4) :: hi_floor
  if (abs(a%hi) < 8388608.0) then  ! 2^23
    hi_floor = real(floor(a%hi), 4)
    if (hi_floor == a%hi .and. a%lo < 0.0) then
      c%hi = a%hi - 1.0
      c%lo = 0.0
    else
      c%hi = hi_floor
      c%lo = 0.0
    end if
  else  ! a%hi is already an integer; refine with lo
    c = fltflt_fast_two_sum(a%hi, real(floor(a%lo), 4))
  end if
end function fltflt_floor

! fltflt_round_to_nearest: round to nearest integer, ties to even.
pure attributes(device) function fltflt_round_to_nearest(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c
  real(4) :: candidate, err, r_lo, frac
  if (abs(a%hi) < 8388608.0) then  ! 2^23
    candidate = anint(a%hi)
    err       = a%hi - candidate
    if (abs(err) /= 0.5) then
      c%hi = candidate
      c%lo = 0.0
    else if (a%lo == 0.0) then
      ! Exact tie at 0.5: round to even
      if (mod(candidate, 2.0) == 0.0) then
        c%hi = candidate
      else
        c%hi = candidate + sign(1.0, err)
      end if
      c%lo = 0.0
    else if ((err > 0.0 .and. a%lo > 0.0) .or. &
             (err < 0.0 .and. a%lo < 0.0)) then
      c%hi = candidate + sign(1.0, err)
      c%lo = 0.0
    else
      c%hi = candidate
      c%lo = 0.0
    end if
  else  ! a%hi is already an integer; round the lo part
    r_lo = anint(a%lo)
    frac = a%lo - r_lo
    if (abs(frac) > 0.5) then
      r_lo = r_lo + sign(1.0, frac)
    else if (abs(frac) == 0.5) then
      if (mod(a%hi + r_lo, 2.0) /= 0.0) then
        r_lo = r_lo + sign(1.0, frac)
      end if
    end if
    c = fltflt_fast_two_sum(a%hi, r_lo)
  end if
end function fltflt_round_to_nearest

! ================================================================
! fltflt_fmod: a - trunc(a/b)*b  (floating-point remainder)
! Returns {NaN, NaN} if b is zero.
! ================================================================

! ff/ff overload.
attributes(device) function fmod_ff_ff(a_in, b_in) result(c)
  type(fltflt), intent(in) :: a_in, b_in
  type(fltflt) :: c
  type(fltflt) :: a_abs, b_abs, q, trunc_q, neg_tq, rem
  real(4) :: sgn
  if (b_in%hi == 0.0 .and. b_in%lo == 0.0) then
    c%hi = b_in%hi / b_in%hi  ! NaN via hardware (runtime 0/0, not a folded literal constant)
    c%lo = c%hi
    return
  end if
  sgn   = sign(1.0, a_in%hi)
  a_abs = fltflt_abs(a_in)
  b_abs = fltflt_abs(b_in)
  q       = div_ff_ff(a_abs, b_abs)
  trunc_q = fltflt_round_toward_zero(q)
  neg_tq%hi = -trunc_q%hi
  neg_tq%lo = -trunc_q%lo
  rem = fma_ff_ff_ff(neg_tq, b_abs, a_abs)  ! a_abs - trunc_q * b_abs
  do while (rem >= b_abs)
    rem = sub_ff_ff(rem, b_abs)
  end do
  do while (lt_ff_r4(rem, 0.0))
    rem = add_ff_ff(rem, b_abs)
  end do
  c%hi = sgn * rem%hi
  c%lo = sgn * rem%lo
end function fmod_ff_ff

! ff/r4 overload.
attributes(device) function fmod_ff_r4(a_in, b_in) result(c)
  type(fltflt), intent(in) :: a_in
  real(4),      intent(in) :: b_in
  type(fltflt) :: c
  type(fltflt) :: a_abs, b_abs, q, trunc_q, neg_tq, rem
  real(4) :: sgn, b_pos
  if (b_in == 0.0) then
    c%hi = b_in / b_in  ! NaN via hardware (runtime 0/0, not a folded literal constant)
    c%lo = c%hi
    return
  end if
  sgn   = sign(1.0, a_in%hi)
  a_abs = fltflt_abs(a_in)
  b_pos = abs(b_in)
  b_abs = fltflt_init(b_pos)
  q       = div_ff_r4(a_abs, b_pos)
  trunc_q = fltflt_round_toward_zero(q)
  neg_tq%hi = -trunc_q%hi
  neg_tq%lo = -trunc_q%lo
  rem = fma_ff_r4_ff(neg_tq, b_pos, a_abs)  ! a_abs - trunc_q * b_pos
  do while (rem >= b_abs)
    rem = sub_ff_ff(rem, b_abs)
  end do
  do while (lt_ff_r4(rem, 0.0))
    rem = add_ff_ff(rem, b_abs)
  end do
  c%hi = sgn * rem%hi
  c%lo = sgn * rem%lo
end function fmod_ff_r4

! ================================================================
! min / max / clamp
! ================================================================

pure attributes(device) function fltflt_min_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: r
  if (a <= b) then;  r = a;  else;  r = b;  end if
end function fltflt_min_ff_ff

pure attributes(device) function fltflt_min_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: r, bf
  bf = fltflt_init(b)
  if (a <= bf) then;  r = a;  else;  r = bf;  end if
end function fltflt_min_ff_r4

pure attributes(device) function fltflt_min_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: r, af
  af = fltflt_init(a)
  if (af <= b) then;  r = af;  else;  r = b;  end if
end function fltflt_min_r4_ff

pure attributes(device) function fltflt_max_ff_ff(a, b) result(r)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: r
  if (a >= b) then;  r = a;  else;  r = b;  end if
end function fltflt_max_ff_ff

pure attributes(device) function fltflt_max_ff_r4(a, b) result(r)
  type(fltflt), intent(in) :: a
  real(4),      intent(in) :: b
  type(fltflt) :: r, bf
  bf = fltflt_init(b)
  if (a >= bf) then;  r = a;  else;  r = bf;  end if
end function fltflt_max_ff_r4

pure attributes(device) function fltflt_max_r4_ff(a, b) result(r)
  real(4),      intent(in) :: a
  type(fltflt), intent(in) :: b
  type(fltflt) :: r, af
  af = fltflt_init(a)
  if (af >= b) then;  r = af;  else;  r = b;  end if
end function fltflt_max_r4_ff

! fltflt_clamp: clamp a to [lo, hi].
pure attributes(device) function fltflt_clamp(a, lo, hi) result(r)
  type(fltflt), intent(in) :: a, lo, hi
  type(fltflt) :: r
  r = fltflt_min_ff_ff(fltflt_max_ff_ff(a, lo), hi)
end function fltflt_clamp

! ================================================================
! ceil
! ================================================================

! fltflt_ceil: smallest integer not less than a. Uses identity ceil(x) = -floor(-x).
pure attributes(device) function fltflt_ceil(a) result(c)
  type(fltflt), intent(in) :: a
  type(fltflt) :: c, neg
  neg%hi = -a%hi
  neg%lo = -a%lo
  c = fltflt_floor(neg)
  c%hi = -c%hi
  c%lo = -c%lo
end function fltflt_ceil

! ================================================================
! hypot
! ================================================================

! fltflt_hypot: sqrt(a^2 + b^2) for real(4) inputs. Uses exact dot2 for the sum of squares.
attributes(device) function fltflt_hypot(a, b) result(c)
  real(4), intent(in) :: a, b
  type(fltflt) :: c
  c = fltflt_sqrt(fltflt_dot2(a, a, b, b))
end function fltflt_hypot

! ================================================================
! cross product
! ================================================================

! fltflt_cross3d: a x b = (ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx).
! Uses fltflt_dot2 for exact cancellation in each component. 3 x (8 flops + 2 FMAs).
attributes(device) subroutine fltflt_cross3d(cx, cy, cz, ax, ay, az, bx, by, bz)
  type(fltflt), intent(out) :: cx, cy, cz
  real(4),      intent(in)  :: ax, ay, az, bx, by, bz
  cx = fltflt_dot2( ay, bz, -az, by)
  cy = fltflt_dot2( az, bx, -ax, bz)
  cz = fltflt_dot2( ax, by, -ay, bx)
end subroutine fltflt_cross3d

! ================================================================
! integer power, sign, lerp, warp reduction
! ================================================================

! fltflt_pow_int: a^n for non-negative integer n via square-and-multiply.
attributes(device) function fltflt_pow_int(a, n) result(r)
  type(fltflt), intent(in) :: a
  integer,      intent(in) :: n
  type(fltflt) :: r, base
  integer :: m
  r%hi = 1.0;  r%lo = 0.0
  base = a;    m = n
  do while (m > 0)
    if (mod(m, 2) == 1) r = r * base
    base = fltflt_square(base)
    m = m / 2
  end do
end function fltflt_pow_int

! fltflt_sign: copysign — magnitude of a with sign of b.
pure attributes(device) function fltflt_sign(a, b) result(c)
  type(fltflt), intent(in) :: a, b
  type(fltflt) :: c, abs_a
  real(4) :: s
  abs_a = fltflt_abs(a)
  s     = sign(1.0, b%hi)
  c%hi  = abs_a%hi * s
  c%lo  = abs_a%lo * s
end function fltflt_sign

! fltflt_lerp: linear interpolation a + t*(b-a) for real(4) t in [0,1].
attributes(device) function fltflt_lerp(a, b, t) result(c)
  type(fltflt), intent(in) :: a, b
  real(4),      intent(in) :: t
  type(fltflt) :: c
  c = a + t * (b - a)
end function fltflt_lerp

! fltflt_warp_reduce_sum: XOR-butterfly all-reduce across 32 lanes.
! After the call every lane holds the sum of all 32 input values.
attributes(device) function fltflt_warp_reduce_sum(val) result(r)
  type(fltflt), intent(in) :: val
  type(fltflt) :: r
  integer :: mask
  r    = val
  mask = 16
  do while (mask > 0)
    r    = r + fltflt_shfl_xor(r, mask)
    mask = mask / 2
  end do
end function fltflt_warp_reduce_sum
