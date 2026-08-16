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
    real(8), intent(in), value :: v1, v2, v3, v4, v5
    real(8), intent(out) :: w0, w1, w2
    real(8) :: b0, b1, b2, a0, a1, a2, s, tau5
    real(8), parameter :: eps = 1.0d-20
    b0 = (13.0d0/12.0d0)*(v1 - 2.0d0*v2 + v3)**2 + 0.25d0*(v1 - 4.0d0*v2 + 3.0d0*v3)**2
    b1 = (13.0d0/12.0d0)*(v2 - 2.0d0*v3 + v4)**2 + 0.25d0*(v2 - v4)**2
    b2 = (13.0d0/12.0d0)*(v3 - 2.0d0*v4 + v5)**2 + 0.25d0*(3.0d0*v3 - 4.0d0*v4 + v5)**2
    tau5 = abs(b0 - b2)
    a0 = 0.1d0 * (1.0d0 + (tau5/(b0+eps))**2)
    a1 = 0.6d0 * (1.0d0 + (tau5/(b1+eps))**2)
    a2 = 0.3d0 * (1.0d0 + (tau5/(b2+eps))**2)
    s = a0 + a1 + a2
    w0 = a0 / s
    w1 = a1 / s
    w2 = a2 / s
  end subroutine weights64

  attributes(device) subroutine weights32(v1, v2, v3, v4, v5, w0, w1, w2)
    real(8), intent(in), value :: v1, v2, v3, v4, v5
    real(4), intent(out) :: w0, w1, w2
    real(4) :: x1, x2, x3, x4, x5
    real(4) :: b0, b1, b2, a0, a1, a2, s, tau5
    real(4), parameter :: eps = 1.0e-20
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    b0 = (13.0_4/12.0_4)*(x1 - 2.0_4*x2 + x3)**2 + 0.25_4*(x1 - 4.0_4*x2 + 3.0_4*x3)**2
    b1 = (13.0_4/12.0_4)*(x2 - 2.0_4*x3 + x4)**2 + 0.25_4*(x2 - x4)**2
    b2 = (13.0_4/12.0_4)*(x3 - 2.0_4*x4 + x5)**2 + 0.25_4*(3.0_4*x3 - 4.0_4*x4 + x5)**2
    tau5 = abs(b0 - b2)
    a0 = 0.1_4 * (1.0_4 + (tau5/(b0+eps))**2)
    a1 = 0.6_4 * (1.0_4 + (tau5/(b1+eps))**2)
    a2 = 0.3_4 * (1.0_4 + (tau5/(b2+eps))**2)
    s = a0 + a1 + a2
    w0 = a0 / s
    w1 = a1 / s
    w2 = a2 / s
  end subroutine weights32

  attributes(device) subroutine poly64_left(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in), value :: v1, v2, v3, v4, v5
    real(8), intent(out) :: p0, p1, p2
    p0 = ( 2.0d0*v1 - 7.0d0*v2 + 11.0d0*v3) / 6.0d0
    p1 = (-1.0d0*v2 + 5.0d0*v3 +  2.0d0*v4) / 6.0d0
    p2 = ( 2.0d0*v3 + 5.0d0*v4 -  1.0d0*v5) / 6.0d0
  end subroutine poly64_left

  attributes(device) subroutine poly64_right(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in), value :: v1, v2, v3, v4, v5
    real(8), intent(out) :: p0, p1, p2
    p0 = (-1.0d0*v1 + 5.0d0*v2 +  2.0d0*v3) / 6.0d0
    p1 = ( 2.0d0*v2 + 5.0d0*v3 -  1.0d0*v4) / 6.0d0
    p2 = (11.0d0*v3 - 7.0d0*v4 +  2.0d0*v5) / 6.0d0
  end subroutine poly64_right

  attributes(device) subroutine poly32_left(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in), value :: v1, v2, v3, v4, v5
    real(4), intent(out) :: p0, p1, p2
    real(4) :: x1, x2, x3, x4, x5
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    p0 = ( 2.0_4*x1 - 7.0_4*x2 + 11.0_4*x3) / 6.0_4
    p1 = (-1.0_4*x2 + 5.0_4*x3 +  2.0_4*x4) / 6.0_4
    p2 = ( 2.0_4*x3 + 5.0_4*x4 -  1.0_4*x5) / 6.0_4
  end subroutine poly32_left

  attributes(device) subroutine poly32_right(v1, v2, v3, v4, v5, p0, p1, p2)
    real(8), intent(in), value :: v1, v2, v3, v4, v5
    real(4), intent(out) :: p0, p1, p2
    real(4) :: x1, x2, x3, x4, x5
    x1 = real(v1,4); x2 = real(v2,4); x3 = real(v3,4)
    x4 = real(v4,4); x5 = real(v5,4)
    p0 = (-1.0_4*x1 + 5.0_4*x2 +  2.0_4*x3) / 6.0_4
    p1 = ( 2.0_4*x2 + 5.0_4*x3 -  1.0_4*x4) / 6.0_4
    p2 = (11.0_4*x3 - 7.0_4*x4 +  2.0_4*x5) / 6.0_4
  end subroutine poly32_right

  attributes(device) subroutine weno64_pair(a0, a1, a2, a3, a4, a5, ql, qr)
    real(8), intent(in), value :: a0, a1, a2, a3, a4, a5
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
    real(8), intent(in), value :: a0, a1, a2, a3, a4, a5
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
end module weno_micro_kernels

program weno_micro
  use cudafor
  use weno_micro_kernels
  implicit none
  integer :: n, nrepeat, mode, argn, ierr
  character(len=128) :: mode_name, arg
  real(8), allocatable, device :: x(:,:), out(:,:)
  real(8), allocatable :: h(:,:)
  type(dim3) :: b128, b256, g128, gface

  n = 4194304
  nrepeat = 1
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

  allocate(x(n,3), out(n,6), h(8,6))
  b128 = dim3(128,1,1)
  b256 = dim3(block_threads,1,1)
  g128 = dim3((n + 127)/128,1,1)
  gface = dim3((n + face_threads - 1)/face_threads,1,1)
  call init_input<<<g128,b128>>>(n, x)
  ierr = cudaDeviceSynchronize()

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
  case default
    print *, 'bad mode: ', trim(mode_name)
    error stop 2
  end select
  ierr = cudaDeviceSynchronize()
  if (ierr /= 0) then
    print *, 'cudaDeviceSynchronize failed: ', ierr
    error stop 3
  endif

  h = out(1:8,:)
  print '(a)', 'mode=' // trim(mode_name)
  print '(a,i0,a,i0)', 'n=', n, ' nrepeat=', nrepeat
  print '(a,es16.8)', 'checksum=', sum(h)
  deallocate(x, out, h)
end program weno_micro
