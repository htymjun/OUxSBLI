program main
  use mpi
  use cudafor
  use mod_globals, only : nx, ny, nz, Lx, Ly, Lz, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use mod_constant, only : id_recal
  use set
  use set_coordinate
  use calc_time_dev
  use main_common
  implicit none
  integer mygpu
  real(8) t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:,:), Q(:,:,:,:)
  ! MPI
  integer nranks, myrank, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call select_gpu(myrank, mygpu)

  call execute_command_line('mkdir -p recal', wait=.true., exitstat=ierr)

  print *, "my rank is", myrank
  if (mod(myrank,2) == 0) then
    call set_block_3D(nx, ny, nz, threads, threadsE, threadsEv, threadsF, threadsFv, threadsG, threadsGv, &
                      blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv)
  endif
  allocate(Q(nx,ny,nz,5), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))
  call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
  call set_Jacobian_xy3(nx, ny, nz, dx, dy, dz, Jacobian)

  if (mod(myrank,2) == 0) then
    if (id_recal) then
      call read_restart(myrank, nx, ny, nz, Q)
    else
      write(*,*) "set initial condition"
      call set_init(myrank, nx, ny, nz, x, y, z, Q)
    endif
  endif

  call cpu_time(t_start)
  call RungeKutta(myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)

  if (mod(myrank,2) == 0) call report_and_save(myrank, nx, ny, nz, Jacobian, Q, t_end - t_start)

  deallocate(Q, x, dx, y, dy, z, dz, Jacobian)
  call MPI_FINALIZE(ierr)
end program main
