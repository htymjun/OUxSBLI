module cpu_gpu_mpi
  use mod_precision
  use cudafor
  use mpi
  implicit none
  interface CPUGPU_MPI_SEND
    module procedure CPU_MPI_SEND, GPU_MPI_SEND
  end interface

  interface CPUGPU_MPI_RECV
    module procedure CPU_MPI_RECV, GPU_MPI_RECV
  end interface
  
  interface CPUGPU_MPI_ISEND
    module procedure CPU_MPI_ISEND, GPU_MPI_ISEND
  end interface

  interface CPUGPU_MPI_IRECV
    module procedure CPU_MPI_IRECV, GPU_MPI_IRECV
  end interface
contains
  subroutine CPU_MPI_SEND(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(2), intent(in)      :: id_gpu_mpi
    real(kd_arr), intent(in), device :: buf(count)
    integer, intent(in), value  :: count, dest, tag, comm
    integer, intent(inout)      :: ireq, ierr
    real(kd_arr), allocatable :: buf_cpu(:)
    allocate(buf_cpu(count))
    buf_cpu = buf
    call MPI_SEND(buf_cpu, count, MPI_KD_ARR, dest, tag, comm, ierr)
    deallocate(buf_cpu)
    ireq = MPI_REQUEST_NULL
  end subroutine CPU_MPI_SEND


  subroutine GPU_MPI_SEND(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(4), intent(in)      :: id_gpu_mpi
    real(kd_arr), intent(in), device :: buf(count)
    integer, intent(in), value  :: count, dest, tag, comm
    integer, intent(inout)      :: ireq, ierr
    call MPI_SEND(buf, count, MPI_KD_ARR, dest, tag, comm, ierr)
    ireq = MPI_REQUEST_NULL
  end subroutine GPU_MPI_SEND


  subroutine CPU_MPI_RECV(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(2), intent(in)       :: id_gpu_mpi
    real(kd_arr), intent(out), device :: buf(count)
    integer, intent(in), value   :: count, dest, tag, comm
    integer, intent(inout)       :: ireq, ierr
    integer istat(MPI_STATUS_SIZE)
    real(kd_arr), allocatable :: buf_cpu(:)
    allocate(buf_cpu(count))
    call MPI_RECV(buf_cpu, count, MPI_KD_ARR, dest, tag, comm, istat, ierr)
    buf = buf_cpu
    deallocate(buf_cpu)
    ireq = MPI_REQUEST_NULL
  end subroutine CPU_MPI_RECV


  subroutine GPU_MPI_RECV(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(4), intent(in)       :: id_gpu_mpi
    real(kd_arr), intent(out), device :: buf(count)
    integer, intent(in), value   :: count, dest, tag, comm
    integer, intent(inout)       :: ireq, ierr
    integer istat(MPI_STATUS_SIZE)
    call MPI_RECV(buf, count, MPI_KD_ARR, dest, tag, comm, istat, ierr)
    ireq = MPI_REQUEST_NULL
  end subroutine GPU_MPI_RECV
  

  ! NOTE: CPU path must stage through a host buffer. MPI_WAIT is called here
  ! before deallocating buf_cpu to avoid freeing memory while the send is
  ! in flight (use-after-free). This makes the CPU variant synchronous, unlike
  ! the GPU variant which returns immediately. Caller's ireq is set to
  ! MPI_REQUEST_NULL. Will be revisited when cpu_gpu_mpi is refactored.
  subroutine CPU_MPI_ISEND(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(2), intent(in)      :: id_gpu_mpi
    real(kd_arr), intent(in), device :: buf(count)
    integer, intent(in), value  :: count, dest, tag, comm
    integer, intent(inout)      :: ireq, ierr
    real(kd_arr), allocatable :: buf_cpu(:)
    integer stat(MPI_STATUS_SIZE)
    allocate(buf_cpu(count))
    buf_cpu = buf
    call MPI_ISEND(buf_cpu, count, MPI_KD_ARR, dest, tag, comm, ireq, ierr)
    call MPI_WAIT(ireq, stat, ierr)
    deallocate(buf_cpu)
  end subroutine CPU_MPI_ISEND


  subroutine GPU_MPI_ISEND(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(4), intent(in)      :: id_gpu_mpi
    real(kd_arr), intent(in), device :: buf(count)
    integer, intent(in), value  :: count, dest, tag, comm
    integer, intent(inout)      :: ireq, ierr
    call MPI_ISEND(buf, count, MPI_KD_ARR, dest, tag, comm, ireq, ierr)
  end subroutine GPU_MPI_ISEND


  subroutine CPU_MPI_IRECV(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(2), intent(in)       :: id_gpu_mpi
    real(kd_arr), intent(out), device :: buf(count)
    integer, intent(in), value   :: count, dest, tag, comm
    integer, intent(inout)       :: ireq, ierr
    integer stat(MPI_STATUS_SIZE)
    real(kd_arr), allocatable :: buf_cpu(:)
    allocate(buf_cpu(count))
    call MPI_IRECV(buf_cpu, count, MPI_KD_ARR, dest, tag, comm, ireq, ierr)
    call MPI_WAIT(ireq, stat, ierr)
    buf = buf_cpu
    deallocate(buf_cpu)
  end subroutine CPU_MPI_IRECV


  subroutine GPU_MPI_IRECV(id_gpu_mpi, buf, count, dest, tag, comm, ireq, ierr)
    integer(4), intent(in)       :: id_gpu_mpi
    real(kd_arr), intent(out), device :: buf(count)
    integer, intent(in), value   :: count, dest, tag, comm
    integer, intent(inout)       :: ireq, ierr
    call MPI_IRECV(buf, count, MPI_KD_ARR, dest, tag, comm, ireq, ierr)
  end subroutine GPU_MPI_IRECV
end module cpu_gpu_mpi

