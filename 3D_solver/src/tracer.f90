!> Lagrangian tracer particles advected with the resolved velocity (TRACER=True in config.fypp)
!> Used by SWLBLI to show that the fluid which carried the initial disturbance of
!> convert_2d_to_3d.f90 leaves the domain (and the separation bubble).
!> Particles are seeded at the grid nodes i_LE <= i <= nx-1, 2 <= j (y <= trace_ymax),
!> 4 <= k <= nz-3 (strides trace_si, trace_sj, trace_sk) and moved with Heun's method:
!> tracer_predict (u^n, before the RK step) and tracer_correct (u^(n+1), after it).
!> The velocity QJ_2..4 / QJ_1 is interpolated trilinearly; z is periodic.
!> Requirements: one compute rank (mpirun -n 2), COMMZ=False, RK=3, RESCALE=False.
!>
!> Files (stream, little endian):
!>   data/data_trace/P######.dat  every trace_nstep steps (###### = step / trace_nstep):
!>                     int32 n, int64 step, real64 t, real32 x(n), y(n), z(n),
!>                     int32 status(n), real32 t_exit(n)
!>                     status: 0 in the domain, 1 left through the outlet, 2 the top, 3 the inlet
!>   data/tracer.d     step, t, number of particles with status 0, 1, 2, 3
!>   recal/P00001.dat  restart: int32 n, int64 step, real64 t, real64 x, y, z, int32 status,
!>                     real64 t_exit (read when trace_restart = .true.)
module tracer
  use cudafor
  use mpi
  use mod_globals, only : dt, trace_restart, trace_ymax, trace_si, trace_sj, trace_sk, trace_nstep, i_LE
  implicit none
  private
  public init_tracer, tracer_predict, tracer_correct, print_tracer, finalize_tracer
  integer, parameter :: nthreads_p = 256
  integer, save :: n_p
  integer(8), save :: step_p
  real(8), save :: z0_p, dz_p, lz_p, xmin_p, xmax_p, ymax_p
  real(8), allocatable, device, save :: x_d(:), y_d(:)
  real(8), allocatable, device, save :: xp(:), yp(:), zp(:), xq(:), yq(:), zq(:), up(:), vp(:), wp(:), t_exit(:)
  integer, allocatable, device, save :: stat_p(:)
contains
  !> index of the cell [a(i), a(i+1)] containing s (clamped to 1..n-1) and the linear weight
  attributes(device) subroutine locate(n, a, s, i, f)
    integer, intent(in), value  :: n
    real(8), intent(in)         :: a(n)
    real(8), intent(in), value  :: s
    integer, intent(out)        :: i
    real(8), intent(out)        :: f
    integer lo, hi, mid
    lo = 1; hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (a(mid) <= s) then
        lo = mid
      else
        hi = mid
      endif
    enddo
    i = lo
    f = min(max((s - a(i)) / (a(i+1) - a(i)), 0.d0), 1.d0)
  end subroutine locate


  !> trilinear interpolation of u, v, w = QJ_2..4 / QJ_1 at (xs, ys, zs)
  attributes(device) subroutine interp_vel(nx, ny, nz, z0, dz, xs, ys, zs, Q_1, Q_2, Q_3, Q_4, u, v, w)
    integer, intent(in), value  :: nx, ny, nz
    real(8), intent(in), value  :: z0, dz, xs, ys, zs
    real(8), intent(in)         :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
    real(8), intent(out)        :: u, v, w
    integer i, j, k, a, b, c, ia, jb, kc
    real(8) fx, fy, fz, sz, wt, r
    call locate(nx, x_d, xs, i, fx)
    call locate(ny, y_d, ys, j, fy)
    ! unique planes 4..nz-3 cover [z0, z0 + (nz-6)*dz); plane nz-2 is the periodic image of 4
    sz = modulo((zs - z0) / dz, dble(nz - 6))
    k  = min(int(sz), nz - 7) + 4
    fz = sz - dble(k - 4)
    u = 0.d0; v = 0.d0; w = 0.d0
    do c = 0, 1
      kc = k + c
      do b = 0, 1
        jb = j + b
        do a = 0, 1
          ia = i + a
          wt = (dble(1-a) + dble(2*a-1) * fx) * (dble(1-b) + dble(2*b-1) * fy) * (dble(1-c) + dble(2*c-1) * fz)
          r  = wt / Q_1(ia,jb,kc)
          u  = u + r * Q_2(ia,jb,kc)
          v  = v + r * Q_3(ia,jb,kc)
          w  = w + r * Q_4(ia,jb,kc)
    enddo;enddo;enddo
  end subroutine interp_vel


  attributes(global) subroutine k_predict(n, nx, ny, nz, z0, dz, dt_, Q_1, Q_2, Q_3, Q_4)
    integer, intent(in), value  :: n, nx, ny, nz
    real(8), intent(in), value  :: z0, dz, dt_
    real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
    integer p
    real(8) u, v, w
    p = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (p > n) return
    if (stat_p(p) /= 0) return
    call interp_vel(nx, ny, nz, z0, dz, xp(p), yp(p), zp(p), Q_1, Q_2, Q_3, Q_4, u, v, w)
    up(p) = u; vp(p) = v; wp(p) = w
    xq(p) = xp(p) + dt_ * u
    yq(p) = abs(yp(p) + dt_ * v)
    zq(p) = zp(p) + dt_ * w
  end subroutine k_predict


  attributes(global) subroutine k_correct(n, nx, ny, nz, z0, dz, lz, xmin, xmax, ymax, dt_, t_new, Q_1, Q_2, Q_3, Q_4)
    integer, intent(in), value  :: n, nx, ny, nz
    real(8), intent(in), value  :: z0, dz, lz, xmin, xmax, ymax, dt_, t_new
    real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
    integer p
    real(8) u, v, w, xn, yn, zn
    p = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (p > n) return
    if (stat_p(p) /= 0) return
    call interp_vel(nx, ny, nz, z0, dz, xq(p), yq(p), zq(p), Q_1, Q_2, Q_3, Q_4, u, v, w)
    xn = xp(p) + 0.5d0 * dt_ * (up(p) + u)
    yn = abs(yp(p) + 0.5d0 * dt_ * (vp(p) + v))   ! reflect at the wall y = 0
    zn = z0 + modulo(zp(p) + 0.5d0 * dt_ * (wp(p) + w) - z0, lz)
    xp(p) = xn; yp(p) = yn; zp(p) = zn
    if (xn >= xmax) then
      stat_p(p) = 1; t_exit(p) = t_new
    else if (yn >= ymax) then
      stat_p(p) = 2; t_exit(p) = t_new
    else if (xn < xmin) then
      stat_p(p) = 3; t_exit(p) = t_new
    endif
  end subroutine k_correct


  !> seed the particles (or read recal/P00001.dat) on the compute rank and write the first file
  subroutine init_tracer(myrank, nranks, nx, ny, nz, x, y, z)
    integer, intent(in) :: myrank, nranks, nx, ny, nz
    real(8), intent(in) :: x(nx), y(ny), z(nz)
    real(8), allocatable :: h_x(:), h_y(:), h_z(:), h_t(:)
    integer, allocatable :: h_s(:)
    integer i, j, k, p, u, ios, ierr, nj
    real(8) t_p
    if (nranks /= 2) then
      print *, "tracer: only one compute rank is supported (mpirun -n 2), nranks =", nranks
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
    allocate(x_d(nx), y_d(ny))
    x_d = x; y_d = y
    z0_p = z(4); dz_p = z(5) - z(4); lz_p = dble(nz - 6) * dz_p
    xmin_p = x(1); xmax_p = x(nx); ymax_p = y(ny)
    call execute_command_line('mkdir -p data/data_trace', wait=.true., exitstat=ierr)

    if (trace_restart) then
      open(newunit=u, file="recal/P00001.dat", form="unformatted", access="stream", status="old", &
           action="read", iostat=ios)
      if (ios /= 0) then
        print *, "tracer: cannot open recal/P00001.dat (trace_restart = .true.)"
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      endif
      read(u) n_p, step_p, t_p
      allocate(h_x(n_p), h_y(n_p), h_z(n_p), h_s(n_p), h_t(n_p))
      read(u) h_x, h_y, h_z, h_s, h_t
      close(u)
      print *, "tracer: restarted", n_p, " particles at step", step_p, " t =", t_p
    else
      nj = 0
      do j = 2, ny, trace_sj
        if (y(j) <= trace_ymax) nj = nj + 1
      enddo
      n_p = ((nx - 1 - i_LE) / trace_si + 1) * nj * ((nz - 7) / trace_sk + 1)
      allocate(h_x(n_p), h_y(n_p), h_z(n_p), h_s(n_p), h_t(n_p))
      p = 0
      do k = 4, nz-3, trace_sk
        do j = 2, ny, trace_sj
          if (y(j) > trace_ymax) cycle
          do i = i_LE, nx-1, trace_si
            p = p + 1
            h_x(p) = x(i); h_y(p) = y(j); h_z(p) = z(k)
      enddo;enddo;enddo
      h_s = 0; h_t = 0.d0
      step_p = 0
      print *, "tracer: seeded", n_p, " particles"
    endif

    allocate(xp(n_p), yp(n_p), zp(n_p), xq(n_p), yq(n_p), zq(n_p), up(n_p), vp(n_p), wp(n_p), &
             t_exit(n_p), stat_p(n_p))
    xp = h_x; yp = h_y; zp = h_z; stat_p = h_s; t_exit = h_t
    deallocate(h_x, h_y, h_z, h_s, h_t)
    if (.not. trace_restart) call print_tracer()
  end subroutine init_tracer


  subroutine tracer_predict(nx, ny, nz, Q_1, Q_2, Q_3, Q_4)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
    call k_predict<<<(n_p + nthreads_p - 1) / nthreads_p, nthreads_p>>>(n_p, nx, ny, nz, z0_p, dz_p, dt, &
                                                                        Q_1, Q_2, Q_3, Q_4)
  end subroutine tracer_predict


  !> corrector step; writes data/data_trace/P######.dat every trace_nstep steps
  subroutine tracer_correct(nx, ny, nz, Q_1, Q_2, Q_3, Q_4)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in), device :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz)
    step_p = step_p + 1
    call k_correct<<<(n_p + nthreads_p - 1) / nthreads_p, nthreads_p>>>(n_p, nx, ny, nz, z0_p, dz_p, lz_p, &
      xmin_p, xmax_p, ymax_p, dt, dble(step_p) * dt, Q_1, Q_2, Q_3, Q_4)
    if (mod(step_p, int(trace_nstep, 8)) == 0) call print_tracer()
  end subroutine tracer_correct


  subroutine print_tracer()
    real(8), allocatable :: h(:)
    integer, allocatable :: h_s(:)
    integer u, cnt(0:3), s
    character(len=40) filename
    allocate(h(n_p), h_s(n_p))
    h_s = stat_p
    write(filename, "(a, i6.6, a)") "data/data_trace/P", int(step_p / trace_nstep), ".dat"
    open(newunit=u, file=filename, status="replace", action="write", form="unformatted", access="stream", &
         convert="LITTLE_ENDIAN")
    write(u) n_p, step_p, dble(step_p) * dt
    h = xp; write(u) real(h, 4)
    h = yp; write(u) real(h, 4)
    h = zp; write(u) real(h, 4)
    write(u) h_s
    h = t_exit; write(u) real(h, 4)
    close(u)
    do s = 0, 3
      cnt(s) = count(h_s == s)
    enddo
    open(newunit=u, file="data/tracer.d", position="append")
    write(u, "(i12, es16.8, 4i12)") step_p, dble(step_p) * dt, cnt
    close(u)
    deallocate(h, h_s)
  end subroutine print_tracer


  !> write the restart file recal/P00001.dat and free the device memory
  subroutine finalize_tracer()
    real(8), allocatable :: h_x(:), h_y(:), h_z(:), h_t(:)
    integer, allocatable :: h_s(:)
    integer u
    allocate(h_x(n_p), h_y(n_p), h_z(n_p), h_s(n_p), h_t(n_p))
    h_x = xp; h_y = yp; h_z = zp; h_s = stat_p; h_t = t_exit
    open(newunit=u, file="recal/P00001.dat", status="replace", action="write", form="unformatted", access="stream")
    write(u) n_p, step_p, dble(step_p) * dt
    write(u) h_x, h_y, h_z, h_s, h_t
    close(u)
    print *, "tracer: wrote recal/P00001.dat at step", step_p
    deallocate(h_x, h_y, h_z, h_s, h_t)
    deallocate(x_d, y_d, xp, yp, zp, xq, yq, zq, up, vp, wp, t_exit, stat_p)
  end subroutine finalize_tracer
end module tracer
