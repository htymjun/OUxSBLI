! Generic-name declarations for the fltflt "additional subroutine" bodies in
! fltflt_subroutines.f90. `include`d into a module's specification part
! (after `implicit none`, before `contains`) — every module procedure named
! below must also be brought in via `include 'fltflt_subroutines.f90'` after
! `contains`. Requires `fltflt_operator_interfaces.f90` / 'fltflt_operator.f90'
! to also be included in the same module (these bodies use `+ - *` and the
! EFT backend).

! fltflt_add_same_sign: 11 flops; valid only when a and b share the same sign.
interface fltflt_add_same_sign
  module procedure add_same_sign_ff_ff, add_same_sign_ff_r4, add_same_sign_r4_ff
end interface fltflt_add_same_sign

! fltflt_fma: a*b + c — more accurate and efficient than fltflt_add(fltflt_mul(a,b),c).
interface fltflt_fma
  module procedure fma_ff_ff_ff, fma_ff_ff_r4, fma_r4_ff_ff, fma_ff_r4_ff, &
                   fma_ff_r4_r4, fma_r4_ff_r4, fma_r4_r4_ff
end interface fltflt_fma

! fltflt_fma_approx: a*b + c omitting the a%lo*b%lo term (~1 ULP less precise,
! but faster in throughput-bound kernels).
interface fltflt_fma_approx
  module procedure fma_approx_ff_ff_ff, fma_approx_ff_ff_r4, &
                   fma_approx_r4_ff_ff, fma_approx_ff_r4_ff, &
                   fma_approx_ff_r4_r4, fma_approx_r4_ff_r4, &
                   fma_approx_r4_r4_ff
end interface fltflt_fma_approx

! fltflt_fmod: floating-point remainder a - trunc(a/b)*b (2 overloads: ff/ff and ff/r4).
interface fltflt_fmod
  module procedure fmod_ff_ff, fmod_ff_r4
end interface fltflt_fmod

! fltflt_dot2/3/4: compensated dot products. 3 overloads each (r4 / r8 / fltflt).
interface fltflt_dot2
  module procedure dot2_r4, dot2_r8, dot2_ff
end interface fltflt_dot2

interface fltflt_dot3
  module procedure dot3_r4, dot3_r8, dot3_ff
end interface fltflt_dot3

interface fltflt_dot4
  module procedure dot4_r4, dot4_r8, dot4_ff
end interface fltflt_dot4

! fltflt_dot6: a*b+c*d+e*f+g*h+i*j+k*l, fltflt operands only (no r4/r8 use site
! needs those overloads -- see calc_keep_3d_ff.f90's KEEP6_ff).
interface fltflt_dot6
  module procedure dot6_ff
end interface fltflt_dot6

! fltflt_min / fltflt_max: branchless on GPU (predicated select). 3 overloads each.
interface fltflt_min
  module procedure fltflt_min_ff_ff, fltflt_min_ff_r4, fltflt_min_r4_ff
end interface fltflt_min

interface fltflt_max
  module procedure fltflt_max_ff_ff, fltflt_max_ff_r4, fltflt_max_r4_ff
end interface fltflt_max
