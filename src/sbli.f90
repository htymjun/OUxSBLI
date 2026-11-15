!> Main program of the two-domain SWTBLI case (4 MPI ranks): ranks 0/1 run the
!> recycling-rescaling TBL generator, ranks 2/3 the SBLI domain. The single-domain
!> counterpart is 3D_solver/src/main.f90; the shared steps are in main_common.
program main
  use cudafor
  use mpi
  use mod_globals,  only : nx1, nx2, ny1, ny2, nz1, nz2, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use mod_constant, only : id_recal
  use set_coordinate
  use calc_time_dev
  use main_common
  implicit none
  integer mygpu, nx, ny, nz
  real(8) Lx, Ly, Lz, t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:,:), Q(:,:,:,:)
  ! MPI
  integer nranks, myrank, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  if (nranks /= 4) then
    if (myrank == 0) print *, "sbli.f90: the two-domain case needs exactly 4 MPI ranks, got", nranks
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  endif
  if (myrank <= 1) then
    ! boundary layer (recycling-rescaling generator)
    nx = nx1
    ny = ny1
    nz = nz1
    Lx = Lx1
    Ly = Ly1
    Lz = Lz1
  elseif (2 <= myrank) then
    ! boundary layer + oblique shock
    nx = nx2
    ny = ny2
    nz = nz2
    Lx = Lx2
    Ly = Ly2
    Lz = Lz2
  endif
  ! each compute/IO pair shares one device, so with two visible GPUs the TBL pair
  ! runs on device 0 and the SBLI pair on device 1; with one GPU both share it.
  call select_gpu(myrank, mygpu)
  call execute_command_line('mkdir -p recal data/1 data/3', wait=.true., exitstat=ierr)

  call set_block_3D(nx, ny, nz, threads, threadsE, threadsEv, threadsF, threadsFv, threadsG, threadsGv, &
                    blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv)
  allocate(Q(nx,ny,nz,5), x(nx), dx(nx), y(ny), dy(ny), z(nz), dz(nz), Jacobian(nx,ny))

  ! set grid information (the I/O ranks take the grid of their compute rank)
  if (mod(myrank,2) == 0) then
    call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, 0.9d0 * Lx1, x, y, z, dx, dy, dz)
    if (id_recal) then
      call read_restart(myrank, nx, ny, nz, Q)
    else
      print *, "myrank is ", myrank, "set initial condition"
      call set_init(myrank, nx, ny, nz, x, y, z, Q)
    endif
  else
    call set_grid(myrank-1, nx, ny, nz, Lx, Ly, Lz, 0.9d0 * Lx1, x, y, z, dx, dy, dz)
  endif
  call set_Jacobian_xy3(nx, ny, nz, dx, dy, dz, Jacobian)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call cpu_time(t_start)
  call RungeKutta(myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)

  if (mod(myrank,2) == 0) call report_and_save(myrank, nx, ny, nz, Jacobian, Q, t_end - t_start)

  deallocate(Q, x, dx, y, dy, z, dz, Jacobian)
  call MPI_FINALIZE(ierr)
end program main
