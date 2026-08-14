! Generic-name declarations for the fltflt operator bodies in
! fltflt_operator.f90. `include`d into a module's specification part (after
! `implicit none`, before `contains`) — every module procedure named below
! must also be brought in via `include 'fltflt_operator.f90'` after `contains`.

interface operator(+)
  module procedure add_ff_ff, add_ff_r4, add_r4_ff, add_ff_r8, add_r8_ff
end interface

interface operator(-)
  module procedure neg_ff, sub_ff_ff, sub_ff_r4, sub_r4_ff, sub_ff_r8, sub_r8_ff
end interface

interface operator(*)
  module procedure mul_ff_ff, mul_ff_r4, mul_r4_ff, mul_ff_r8, mul_r8_ff
end interface

interface operator(/)
  module procedure div_ff_ff, div_ff_r4, div_r4_ff, div_ff_r8, div_r8_ff
end interface

interface operator(==)
  module procedure eq_ff_ff, eq_ff_r4, eq_r4_ff
end interface

interface operator(/=)
  module procedure ne_ff_ff, ne_ff_r4, ne_r4_ff
end interface

interface operator(<)
  module procedure lt_ff_ff, lt_ff_r4, lt_r4_ff
end interface

interface operator(>)
  module procedure gt_ff_ff, gt_ff_r4, gt_r4_ff
end interface

interface operator(<=)
  module procedure le_ff_ff, le_ff_r4, le_r4_ff
end interface

interface operator(>=)
  module procedure ge_ff_ff, ge_ff_r4, ge_r4_ff
end interface
