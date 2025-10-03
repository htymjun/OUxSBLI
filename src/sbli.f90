program main
  use, intrinsic :: iso_fortran_env
  use cudafor
  use mpi
  use mod_globals,  only : id_RungeKutta, id_recal, nx1, nx2, ny1, ny2, nz1, nz2, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, &
  & mygpu1, mygpu2, blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j, l, mygpu, m, s, nx, ny, nz, ios
  real(8) Lx, Ly, Lz, t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:,:), Q(:,:,:,:)
  character(len=8) header
  character(len=40) filename
  logical is_sequential
  ! MPI
  integer nranks, myrank, ierr, ireq, istat(MPI_STATUS_SIZE)

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  if (myrank <= 1) then
    ! boundary layer
    nx = nx1
    ny = ny1
    nz = nz1
    Lx = Lx1
    Ly = Ly1
    Lz = Lz1
    mygpu = mygpu1
  elseif (2 <= myrank) then
    ! boundary layer + oblique shock
    nx = nx2
    ny = ny2
    nz = nz2
    Lx = Lx2
    Ly = Ly2
    Lz = Lz2
    mygpu = mygpu2
  endif
  call set_block(nx, ny, nz, threads, threadsE, threadsEv, threadsF, threadsFv, threadsG, threadsGv, &
                 blocks, blocksE, blocksEv, blocksF, blocksFv, blocksG, blocksGv)
  allocate(Q(5,nx,ny,nz), x(nx), dx(nx), y(ny), dy(ny), z(nz), dz(nz), Jacobian(nx,ny))

  ! set grid information
  if (mod(myrank,2) == 0) then
    call set_grid(myrank, nx, ny, nz, Lx, Ly, Lz, 0.9d0 * Lx1, x, y, z, dx, dy, dz)
    if (kind(id_recal) == 4) then
      write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
      open(10, file=filename, action="read", form="unformatted", access="sequential", status="old", iostat=ios)
      if (ios /= 0) then
        print *, "Error opening file."
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
    elseif (kind(id_recal) == 2) then
      print *, "myrank is ", myrank, "set initial condition"
      call set_init(myrank, nx, ny, nz, x, y, z, Q)
    else
      print *, "wrong paramater was found"
    endif
  else
    call set_grid(myrank-1, nx, ny, nz, Lx, Ly, Lz, 0.9d0 * Lx1, x, y, z, dx, dy, dz)
  endif
  call set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call cpu_time(t_start)
  call RungeKutta(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
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
        do i = 1, nx
          do m = 1, 5
            Q(m,i,j,l) = Jacobian(i,j) * Q(m,i,j,l)
    enddo;enddo;enddo;enddo
    call cpu_time(t_start)
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10,file=filename,status="replace",action="write",form="unformatted",access="sequential")
    header = 'SEQFMT01'
    write(10) header
    write(10) Q
    close(10)
    call cpu_time(t_end)
    s = t_end - t_start
    print *, "output time:", s, " [sec]"
  endif

  deallocate(Q,x,dx,y,dy,z,dz,Jacobian)
  call MPI_FINALIZE(ierr)
end program main

