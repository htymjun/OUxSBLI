!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CUDA Fortran Implementation of fltflt
!
! Copyright (c) 2026, Jun Hatayama
! This is a Fortran translation/derivative work of the original C++ fltflt.h.
!
! Original C++ implementation:
! Copyright (c) 2026, NVIDIA Corporation
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Float-float (FF) precision arithmetic — type definition + constructors only.
!
! Represents a value as an unevaluated sum (hi + lo) of two real(4) numbers,
! giving ~48 mantissa bits — equivalent to real(8) precision.
!
! This module carries only the `fltflt` type itself and its `fltflt_init`
! constructors: it has to stay a real, `use`-able module because `type(fltflt)`
! is declared as module-level `constant`/`shared`/`device` data across many
! separately-compiled files (mod_constant.f90.fypp, calc_keep_kernel_internal_ff,
! calc_slau_kernel_internal_ff, ...), which all need to agree on the identical
! type. `fltflt_init` runs on the host once at startup (see
! init_fltflt_constants / init_keep_ff_constants), so there's no hot-path
! inlining concern for it either.
!
! Everything else that used to live here — arithmetic/comparison operators and
! the "additional subroutines" (fltflt_fma, fltflt_dot2/3/4, fltflt_sqrt, ...)
! — has moved to fltflt_operator_interfaces.f90 / fltflt_operator.f90 /
! fltflt_subroutines_interfaces.f90 / fltflt_subroutines.f90, which are
! `include`d directly into each kernel module
! that needs them (calc_keep_kernel_internal_ff.f90.fypp,
! calc_slau_kernel_internal_ff.f90.fypp, calc_hybrid_kernel_ff.f90.fypp,
! calc_hybrid_kernel_internal_ff.f90.fypp, calc_muscl.f90.fypp, calc_hybrid_ff.f90)
! instead of `use fltflt`d from here. That keeps the hot arithmetic as local,
! same-translation-unit module procedures so nvfortran's ordinary inliner can
! inline them into the caller, rather than relying on cross-module inlining
! through `-Minline=name:...` (3D_solver_fltflt/CMakeLists.txt).
!
! Companion C++ header: fltflt.h  (NVIDIA MatX, BSD-3)
! References:
!   Thall 2006 "Extended-Precision Floating-Point Numbers for GPU Computation"
!   Zhang & Aiken SC'25 "High-Performance Branch-Free Algorithms for
!     Extended-Precision Floating-Point Arithmetic" (FPAN addition)
!   Ogita, Rump, Oishi 2005 "Accurate Sum and Dot Product" (compensated dot)
module fltflt
  use cudafor
  use cudadevice
  implicit none
  private

  ! ================================================================
  ! The fltflt type: an unevaluated sum hi + lo where |lo| <= 0.5*ulp(hi).
  ! No default initializers — nvfortran crashes with default-init on device types.
  ! Always use fltflt_init() or assign %hi/%lo explicitly before first use.
  ! ================================================================
  type, public :: fltflt
    real(4) :: hi, lo
  end type fltflt

  ! ================================================================
  ! Constructors
  ! ================================================================
  interface fltflt_init
    module procedure init_from_r4
    module procedure init_from_r8
  end interface fltflt_init
  public :: fltflt_init

contains

  pure attributes(device, host) function init_from_r4(a) result(r)
    real(4), intent(in) :: a
    type(fltflt) :: r
    r%hi = a
    r%lo = 0.0
  end function init_from_r4

  ! Splits a double exactly into two non-overlapping single-precision halves.
  pure attributes(device, host) function init_from_r8(a) result(r)
    real(8), intent(in) :: a
    type(fltflt) :: r
    r%hi = real(a, 4)
    r%lo = real(a - real(r%hi, 8), 4)
  end function init_from_r8

end module fltflt
