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
!> ratio; the ratio_cap guard in weights32 is still correct and still needed,
!> but it was not the cause. Plain `intent(in)` is semantically equivalent here
!> (these are read-only scalars) and generates correct code.
!>
!> DIVISION BUDGET -- kept in step with src/calc_weno.f90 (see
!> report/rtx4060_weno_division_reduction.md). `div.rn.f64` is ~20 SASS
!> instructions (MUFU.RCP64H + Newton DFMAs + a CALL to the slow path), and
!> `-fast` does NOT fold `/6.0d0` into a multiply. The textbook split form used
!> here had 18 divisions per weno64_pair: 6 per weights64 (3 ratios + 3
!> normalisations) and 3 per poly64_*. Now 4:
!>
!>   weights64  6 -> 2   batched ratio inversion + one reciprocal for 1/s
!>   weights32  6 -> 4   normalisation batched only; the ratio batch would
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
  attributes(device) subroutine weights64(v1, v2, v3, v4, v5, w0, w1, w2)
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
  end subroutine weights64

  attributes(device) subroutine weights32(v1, v2, v3, v4, v5, w0, w1, w2)
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
  end subroutine weights32

  ! The 1/6 is folded into the coefficients, so each candidate polynomial is a
  ! 3-term FMA chain with no division. Written as named parameters rather than
  ! decimal literals so the FP32 twins below get correctly rounded values at
  ! real(4) instead of truncated digits.
  attributes(device) subroutine poly64_left(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(8), intent(out) :: p0, p1, p2
    real(8), parameter :: c1 = 1.0d0/6.0d0, c2 = 2.0d0/6.0d0, c5 = 5.0d0/6.0d0
    real(8), parameter :: c7 = 7.0d0/6.0d0, c11 = 11.0d0/6.0d0
    p0 =  c2*v1 - c7*v2 + c11*v3
    p1 = -c1*v2 + c5*v3 +  c2*v4
    p2 =  c2*v3 + c5*v4 -  c1*v5
  end subroutine poly64_left

  attributes(device) subroutine poly64_right(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in) :: v1, v2, v3, v4, v5
    real(8), intent(out) :: p0, p1, p2
    real(8), parameter :: c1 = 1.0d0/6.0d0, c2 = 2.0d0/6.0d0, c5 = 5.0d0/6.0d0
    real(8), parameter :: c7 = 7.0d0/6.0d0, c11 = 11.0d0/6.0d0
    p0 = -c1*v1 + c5*v2 +  c2*v3
    p1 =  c2*v2 + c5*v3 -  c1*v4
    p2 = c11*v3 - c7*v4 +  c2*v5
  end subroutine poly64_right

  attributes(device) subroutine poly32_left(v1, v2, v3, v4, v5, p0, p1, p2)
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
  end subroutine poly32_left

  attributes(device) subroutine poly32_right(v1, v2, v3, v4, v5, p0, p1, p2)
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
  end subroutine poly32_right

  attributes(device) subroutine weno64_pair(a0, a1, a2, a3, a4, a5, ql, qr)
    real(8), intent(in) :: a0, a1, a2, a3, a4, a5
    real(8), intent(out) :: ql, qr
    real(8) :: w0, w1, w2, p0, p1, p2
    call weights64(a0, a1, a2, a3, a4, w0, w1, w2)
    call poly64_left(a0, a1, a2, a3, a4, p0, p1, p2)
    ql = w0*p0 + w1*p1 + w2*p2
    call weights64(a1, a2, a3, a4, a5, w0, w1, w2)
    call poly64_right(a1, a2, a3, a4, a5, p0, p1, p2)
    qr = w0*p0 + w1*p1 + w2*p2
  end subroutine weno64_pair

  attributes(device) subroutine weno32_pair(a0, a1, a2, a3, a4, a5, ql, qr)
    real(8), intent(in) :: a0, a1, a2, a3, a4, a5
    real(8), intent(out) :: ql, qr
    real(4) :: w0, w1, w2, p0, p1, p2
    call weights32(a0, a1, a2, a3, a4, w0, w1, w2)
    call poly32_left(a0, a1, a2, a3, a4, p0, p1, p2)
    ql = real(w0*p0 + w1*p1 + w2*p2, 8)
    call weights32(a1, a2, a3, a4, a5, w0, w1, w2)
    call poly32_right(a1, a2, a3, a4, a5, p0, p1, p2)
    qr = real(w0*p0 + w1*p1 + w2*p2, 8)
  end subroutine weno32_pair

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
          call weno64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), ql(f), qr(f))
        else
          call weno32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), ql(f), qr(f))
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
            if (fp64_field) call weno64_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                                             q(idx,2*f-1), q(idx,2*f))
          enddo
        else
          do f = 1, 3
            fp64_field = (mode == mode_fp64) .or. (mode == mode_rho64 .and. f == 1) .or. &
                         (mode == mode_u64 .and. f == 2) .or. (mode == mode_p64 .and. f == 3)
            if (.not. fp64_field) call weno32_pair(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
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
        call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
        call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
        q(2*f-1) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8)
        call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
        call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
            call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
            call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
            call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
            call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
            call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
            call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
            call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
          call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                         wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
          call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                         wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
          psh(idx,b+1) = p0; psh(idx,b+2) = p1; psh(idx,b+3) = p2
          call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
            call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                           wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
            call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                           wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          enddo
        else
          do f = 1, 3
            b = 6*(f-1)
            call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
            p(b+1) = p0; p(b+2) = p1; p(b+3) = p2
            call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
          call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                         wsh(idx,b+1), wsh(idx,b+2), wsh(idx,b+3))
          call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                         wsh(idx,b+4), wsh(idx,b+5), wsh(idx,b+6))
          call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
          p(b+1) = p0; p(b+2) = p1; p(b+3) = p2
          call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
              call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                             wsh(idx,b+1,t), wsh(idx,b+2,t), wsh(idx,b+3,t))
              call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                             wsh(idx,b+4,t), wsh(idx,b+5,t), wsh(idx,b+6,t))
            enddo
          else
            do f = 1, 3
              b = 6*(f-1)
              call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
              p(b+1,t) = p0; p(b+2,t) = p1; p(b+3,t) = p2
              call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
              call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), &
                             wsh(idx,b+1,t), wsh(idx,b+2,t), wsh(idx,b+3,t))
              call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), &
                             wsh(idx,b+4,t), wsh(idx,b+5,t), wsh(idx,b+6,t))
              call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
              p(b+1,t) = p0; p(b+2,t) = p1; p(b+3,t) = p2
              call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
        call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
        call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
        q(2*f-1) = w0*real(p0,8) + w1*real(p1,8) + w2*real(p2,8)
        call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
        call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
            call weights64(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), w0, w1, w2)
          else
            call poly32_left(x(i,f), x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), p0, p1, p2)
          endif
          pp0 = __shfl_xor(p0, 16)
          pp1 = __shfl_xor(p1, 16)
          pp2 = __shfl_xor(p2, 16)
          if (lower) q(2*f-1) = w0*real(pp0,8) + w1*real(pp1,8) + w2*real(pp2,8)

          p0 = 0.0_4; p1 = 0.0_4; p2 = 0.0_4
          if (lower) then
            call weights64(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), w0, w1, w2)
          else
            call poly32_right(x(i+1,f), x(i+2,f), x(i+3,f), x(i+4,f), x(i+5,f), p0, p1, p2)
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
