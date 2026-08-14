! Does the recovered fltflt library survive the 1D solver's -Mfma?
!
! Replicates the accuracy probe from ~/fp64fp32/df_concurrent.cu: iterate
! x <- 0.99*x + 0.013 a thousand times over many series and compare DF and
! plain FP32 against an FP64 reference. FMA contraction would destroy the
! error-free transforms that double-float depends on, and that shows up here
! as DF error collapsing from ~1e-13 to ~1e-6.
module df_test
  use cudafor
  use cudadevice
  use fltflt
  implicit none
  include 'fltflt_operator_interfaces.f90'
  include 'fltflt_subroutines_interfaces.f90'
contains
  include 'fltflt_operator.f90'
  include 'fltflt_subroutines.f90'

  attributes(global) subroutine k_df(n, nit, out)
    integer, value :: n, nit
    real(8), device :: out(n)
    integer :: i, k
    type(fltflt) :: x, a, b
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n) return
    x = fltflt_init(0.5d0 + 1.0d-4 * dble(i))
    a = fltflt_init(0.99d0)
    b = fltflt_init(0.013d0)
    do k = 1, nit
      x = a * x + b
    enddo
    out(i) = real(x%hi, 8) + real(x%lo, 8)
  end subroutine k_df

  attributes(global) subroutine k_f32(n, nit, out)
    integer, value :: n, nit
    real(8), device :: out(n)
    integer :: i, k
    real(4) :: x
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n) return
    x = real(0.5d0 + 1.0d-4 * dble(i), 4)
    do k = 1, nit
      x = 0.99_4 * x + 0.013_4
    enddo
    out(i) = real(x, 8)
  end subroutine k_f32

  attributes(global) subroutine k_f64(n, nit, out)
    integer, value :: n, nit
    real(8), device :: out(n)
    integer :: i, k
    real(8) :: x
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    if (i > n) return
    x = 0.5d0 + 1.0d-4 * dble(i)
    do k = 1, nit
      x = 0.99d0 * x + 0.013d0
    enddo
    out(i) = x
  end subroutine k_f64
end module df_test

program main
  use cudafor
  use df_test
  implicit none
  integer, parameter :: n = 4096, nit = 1000
  real(8), allocatable, device :: d_df(:), d_32(:), d_64(:)
  real(8) :: h_df(n), h_32(n), h_64(n), e_df, e_32
  type(dim3) :: g, b

  allocate(d_df(n), d_32(n), d_64(n))
  b = dim3(128,1,1); g = dim3((n+127)/128,1,1)
  call k_df <<<g,b>>>(n, nit, d_df)
  call k_f32<<<g,b>>>(n, nit, d_32)
  call k_f64<<<g,b>>>(n, nit, d_64)
  h_df = d_df; h_32 = d_32; h_64 = d_64

  e_df = maxval(abs(h_df - h_64) / max(abs(h_64), 1.d-300))
  e_32 = maxval(abs(h_32 - h_64) / max(abs(h_64), 1.d-300))
  print "(a,es10.3,a)", "  double-float vs FP64 : ", e_df, "   (CUDA-C reference: 1.250e-13)"
  print "(a,es10.3,a)", "  plain FP32   vs FP64 : ", e_32, "   (CUDA-C reference: 5.524e-06)"
  if (e_df < 1.d-11) then
    print *, " -> PASS: error-free transforms survived the build flags"
  else
    print *, " -> FAIL: DF accuracy collapsed; FMA contraction likely broke the EFTs"
  endif
  deallocate(d_df, d_32, d_64)
end program main
