!> Steps shared by the single-domain main (3D_solver/src/main.f90) and the two-domain
!> main of SWTBLI (src/sbli.f90): GPU selection, restart input, timing and restart output.
module main_common
  use cudafor
  use mpi
  implicit none
contains
  !> Pick the GPU from this rank's position WITHIN ITS NODE. A global
  !> myrank/2 only addresses the right device on a single node: spread over
  !> several nodes the higher ranks ask for device numbers their node does not
  !> have, and cudaSetDevice does not reject that (it returns success and the
  !> failure surfaces later, or not at all). Ranks are paired compute/IO, so
  !> two consecutive node-local ranks share one device.
  subroutine select_gpu(myrank, mygpu)
    integer, intent(in)  :: myrank
    integer, intent(out) :: mygpu
    integer comm_node, rank_node, ndevices, cstat, ierr
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_node, ierr)
    call MPI_COMM_RANK(comm_node, rank_node, ierr)
    call MPI_COMM_FREE(comm_node, ierr)
    cstat = cudaGetDeviceCount(ndevices)
    if (cstat /= 0 .or. ndevices < 1) then
      print *, "rank", myrank, ": no CUDA device visible (cudaGetDeviceCount status", cstat, ")"
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
    mygpu = mod(rank_node / 2, ndevices)
  end subroutine select_gpu


  !> Read recal/Q<myrank/2+1>.dat: a raw stream of Q, or a sequential file that
  !> starts with the header 'SEQFMT01'.
  subroutine read_restart(myrank, nx, ny, nz, Q)
    integer, intent(in)  :: myrank, nx, ny, nz
    real(8), intent(out) :: Q(nx,ny,nz,5)
    character(len=8) header
    character(len=40) filename
    integer ios, ierr
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10, file=filename, action="read", form="unformatted", access="sequential", status="old", iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file."
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
    read(10, iostat=ios) header
    close(10)
    if (ios == 0 .and. header == 'SEQFMT01') then
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
  end subroutine read_restart


  !> Print the run time t_calc [s] and write the restart file recal/Q<myrank/2+1>.dat.
  !> Q comes back from RungeKutta as Q/J and is turned into the physical state in place.
  subroutine report_and_save(myrank, nx, ny, nz, Jacobian, Q, t_calc)
    integer, intent(in)    :: myrank, nx, ny, nz
    real(8), intent(in)    :: Jacobian(nx,ny), t_calc
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    character(len=40) filename
    integer i, j, l, m, s
    real(8) t_start, t_end
    if (t_calc <= 60.d0) then
      s = int(t_calc)
      print *, "calculation time:", s, " [sec]"
    else
      m = int(t_calc) / 60
      s = int(t_calc) - 60 * m
      print *, "calculation time:", m, " [min] ", s, " [sec]"
    endif
    do l = 1, nz
      do j = 1, ny
        do m = 1, 5
          do i = 1, nx
            Q(i,j,l,m) = Jacobian(i,j) * Q(i,j,l,m)
    enddo;enddo;enddo;enddo
    call cpu_time(t_start)
    write(filename, "(a, i5.5, a)") "recal/Q", int(myrank/2+1), ".dat"
    open(10, file=filename, status="replace", action="write", form="unformatted", access="stream")
    write(10) Q
    close(10)
    call cpu_time(t_end)
    s = t_end - t_start
    print *, "output time:", s, " [sec]"
  end subroutine report_and_save
end module main_common
