program main
  use, intrinsic :: iso_fortran_env
  use cudafor
  use mpi
  use mod_globals, only : nx, threadsE, threadsEv, threads, blocksE, blocksEv, blocks
  use set
  use calc_time_dev
  use print_1d
  implicit none
  integer mygpu, s, m
  real(8) t_start, t_end
  real(8), allocatable :: x(:), Q(:,:)
  ! MPI (single-rank; kept for structural consistency with the 2D/3D solvers)
  integer nranks, myrank, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  mygpu = 0

  ! blocksE/blocksEv/blocks are set by RungeKutta: they depend on the ghost-cell
  ! depth ng, which is a fypp quantity and so is only available there.

  allocate(Q(nx,3), x(nx))
  call set_grid(myrank, nx, x)
  call set_init(myrank, nx, x, Q)

  call cpu_time(t_start)
  call RungeKutta(myrank, mygpu, nx, Q)
  call cpu_time(t_end)

  if (t_end - t_start <= 60.d0) then
    s = int(t_end - t_start)
    print *, "calculation time:", s, " [sec]"
  else
    m = int(t_end - t_start) / 60
    s = int(t_end - t_start) - 60 * m
    print *, "calculation time:", m, " [min] ", s, " [sec]"
  endif

  ! Q.dat is ASCII, ~70 bytes/cell. At the grid sizes used for throughput
  ! benchmarking (nx of a few million) writing it costs more than the run and
  ! nothing reads it -- the correctness checkers all work at nx=4096.
  if (nx <= 1048576) then
    call write_Q_dat(nx, x, Q)
  else
    print *, "nx =", nx, "> 2^20: skipping Q.dat (benchmark mode)"
  endif

  deallocate(Q, x)
  call MPI_FINALIZE(ierr)
end program main
