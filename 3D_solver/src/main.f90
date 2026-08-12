program main
  use, intrinsic :: iso_fortran_env
  use mpi
  use mod_globals, only : dimension, nx, ny, nz, Lx, Ly, Lz, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use mod_constant, only : id_recal
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j, l, m, s, mygpu, ios, errorcode
  real(8) t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:,:), Q(:,:,:,:)
  character(len=8) header
  character(len=40) filename
  logical is_sequential
  ! MPI
  integer nranks, myrank, ierr, ireq, istat(MPI_STATUS_SIZE)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  mygpu = myrank / 2
  
  call execute_command_line('mkdir -p recal', wait=.true., exitstat=ierr)

  print *, "my rank is", myrank
  if (mod(myrank,2) == 0) then
    call set_block_3D(nx, ny, nz, threads, threadsE, threadsEv, threadsF, threadsFv, threadsG, threadsGv, &
                      blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv)
  endif
  allocate(Q(nx,ny,nz,dimension+2), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))
  call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
  call set_Jacobian_xy3(nx, ny, nz, dx, dy, dz, Jacobian)

  if (mod(myrank,2) == 0) then
    if (id_recal) then
      write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
      open(10, file=filename, action="read", form="unformatted", access="sequential", status="old", iostat=ios)
      if (ios /= 0) then
        print *, "Error opening file."
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
        stop
      endif
      read(10, iostat=ios) header
      close(10)
      is_sequential = (ios == 0 .and. header == 'SEQFMT01')
      if (is_sequential) then
        open(10, file=filename, action="read", form="unformatted", access="sequential", status="old")
        read(10) header
        read(10) Q
        print *, "myrank is ", myrank, "simulation has been restarted. access is sequential"
      else
        open(10, file=filename, action="read", form="unformatted", access="stream", status="old")
        rewind(10)
        read(10) Q
        print *, "myrank is ", myrank, "simulation has been restarted. access is stream"
      endif
      close(10)
    else
      write(*,*) "set initial condition"
      call set_init(myrank, nx, ny, nz, x, y, z, Q)
    endif
  endif

  call cpu_time(t_start)
  call RungeKutta(myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)

  if (mod(myrank,2) == 0) then
    ! calculation time
    if (t_end - t_start <= 60.d0) then
      s = int(t_end - t_start)
      print *, "calculation time:", s, " [sec]"
    else
      m = int(t_end - t_start) / 60
      s = int(t_end - t_start) - 60 * m
      print *, "calculation time:", m, " [min] ", s, " [sec]"
    endif
    ! save data
    do l = 1, nz
      do j = 1, ny
        do m = 1, dimension+2
          do i = 1, nx
            Q(i,j,l,m) = Jacobian(i,j) * Q(i,j,l,m)
    enddo;enddo;enddo;enddo
    call cpu_time(t_start)
    
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="stream")
    write(10) Q
    close(10)
    call cpu_time(t_end)
    s = t_end - t_start
    print *, "output time:", s, " [sec]"
  endif

  deallocate(Q, x, dx, y, dy, z, dz, Jacobian)
  call MPI_FINALIZE(ierr)
end program main
