module calc_time_dev
  use cudafor
  use mpi
  use nvtx
  use mod_globals, only : id_rescale, id_exchange, nt, np, nre2, rerank
  use calc_flux_base
  use calc_steps
  use calc_rescale
  use calc_para
  use calc_rand
  use set
  use print
  implicit none
  interface RungeKutta
    module procedure RungeKutta_3rd, RungeKutta_4th, Gauss_RungeKutta
  end interface
contains
  subroutine check_gpu(mygpu)
    integer, intent(in) :: mygpu
    integer ilen, stat
    type(cudaDeviceProp) prop
    stat = cudaSetDevice(mygpu)
    stat = cudaGetDeviceProperties(prop, mygpu)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
  end subroutine check_gpu


  subroutine allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
    use mod_globals, only : id_visc
    integer, intent(in)                       :: myrank, nx, ny, nz
    real(8), intent(out), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    real(8), intent(out), allocatable, device :: ruvwp(:,:,:,:), T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), intent(out), allocatable, device :: E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    integer ierr
    allocate(ruvwp(5,nx,ny,nz), E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), stat=ierr)
    allocate(dx(nx-1), dy(ny-1), dz(nz-1), xix(nx-1), etay(ny-1), zetaz(nz-1), Jacobian(nx,ny), stat=ierr)
    if (kind(id_visc) == 2) then
      allocate(T(1,1,1), mu(1,1,1), mut(1,1,1), qc2(1,1,1), stat=ierr)
    elseif (kind(id_visc) == 4) then
      allocate(T(nx,ny,nz), mu(nx,ny,nz), mut(1,1,1), qc2(1,1,1), stat=ierr)
    elseif (kind(id_visc) == 8) then
      allocate(T(nx,ny,nz), mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz), stat=ierr)
    endif
    if (ierr /= 0) then
      print *, "myrank is ", myrank, " memory allocation failed", ierr
    else
      print *, "myrank is ", myrank, " memory allocation has completed"
    endif
  end subroutine allocate_device_mem


  subroutine pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                      dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
    use mod_globals, only : id_accuracy
    integer, intent(in)    :: nx, ny, nz, myrank, nranks
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1), z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    integer, intent(out)   :: overlap
    real(8), intent(out), device :: dx(nx-1), dy(ny-1), dz(nz-1), xix(nx-1), etay(ny-1), zetaz(nz-1), Jacobian(nx,ny)
    real(8), intent(out), device :: QJ(5,nx,ny,nz)
    real(4), intent(inout)       :: ke0, entropy0
    real(8) xix_cpu(nx-1), etay_cpu(ny-1), zetaz_cpu(nz-1)
    real(4) rho1d(nx*ny*nz), p1d(nx*ny*nz), v1d(nx*ny*nz*3)
    integer i, j, k, l, ierr
    ! set Q / Jacobian
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            Q(l,i,j,k) = Q(l,i,j,k) / Jacobian_cpu(i,j)
    enddo;enddo;enddo;enddo
    ! copy on GPU
    xix_cpu   = 1.d0 / dx_cpu
    etay_cpu  = 1.d0 / dy_cpu
    zetaz_cpu = 1.d0 / dz_cpu
    dx       = dx_cpu
    dy       = dy_cpu
    dz       = dz_cpu
    xix      = xix_cpu
    etay     = etay_cpu
    zetaz    = zetaz_cpu
    Jacobian = Jacobian_cpu
    QJ = Q
    ! for multi GPU
    if (kind(id_accuracy) == 8) then
      overlap = 3
    elseif (kind(id_accuracy) == 4) then
      overlap = 2
    else
      overlap = 1
    endif
    call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho1d, p1d, v1d)
    call print_vtk(0, nx, ny, nz, myrank+1, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
    call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
  end subroutine pre_calc


  subroutine pre_rescale(myrank, ny, nz, Qre, Qm, Qm_cpu)
    integer, intent(in)                         :: myrank, ny, nz
    real(8), intent(inout), allocatable, device :: Qre(:), Qm(:)
    real(8), intent(inout), allocatable         :: Qm_cpu(:)
    integer stat, ilen, ierr
    type(cudaDeviceProp) prop
    if (mod(myrank,2) == 0) then
      allocate(Qre(ny*(nz-6)*5), Qm(ny*5), stat=ierr)
      if (ierr /= 0) then
        print *, "myrank is ", myrank, " memory allocation failed (Qm)", ierr
      else
        print *, "myrank is ", myrank, " memory allocation has completed (Qm)"
      endif
    elseif (myrank == rerank+1) then
      stat = cudaSetDevice(0)
      stat = cudaGetDeviceProperties(prop, 0)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", 0, ") calculates rescaling"
      allocate(Qm_cpu(ny*5))
    endif
  end subroutine pre_rescale


  subroutine RungeKutta_3rd(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=2), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! rescal_cpu!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step, flag_re = 0
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJ2(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! Landau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(8), allocatable, device :: seed(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    print *, "rank", myrank, " has found ", ndevices, " GPU devices"
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJ2(5,nx,ny,nz), stat=ierr)
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      if (kind(id_LL) == 4) then
        allocate(seed(nx,ny,nz))
        call init_seed(nx, ny, nz, seed)
      endif
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    if (kind(id_rescale) == 4) then
      call pre_rescale(myrank, ny, nz, Qre, Qm, Qm_cpu)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    do t2 = 1, np
      do t1 = 1, nt
        step = np * (t2-1) + t1
        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(1, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          endif
          call nvtxStartRange("calc flux", 1)
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          !print *, "myrank is ", myrank, " calc EFG"
          call nvtxEndRange
          call nvtxStartRange("calc step", 2)
          call calc_step1(nx, ny, nz, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          !print *, "myrank is ", myrank, " calc step"
          call nvtxEndRange
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call nvtxStartRange("exchange", 3)
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ2)
            call nvtxEndRange
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
            call nvtxStartRange("set bc", 4)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)
          endif
          !print *, "myrank is ", myrank, " set bc"
          call nvtxEndRange
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call nvtxStartRange("calc rescale", 5)
          call rescale_recv_send(1, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
          call nvtxEndRange
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(2, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          call calc_step2_3(nx, ny, nz, 0.75d0, 0.25d0, 0.25d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJ2)
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ2)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ2)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(3, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ2, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ2, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          call calc_step2_3(nx, ny, nz, 2.d0, 1.d0, 2.d0, 3.d0, dx, dy, dz, E, F, G, QJ2, QJ)
          if (ndevices >= 2 .and. kind(id_exchange) == 4) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJ2, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(seed)
      endif
    endif
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        deallocate(Qre, Qm)
      elseif (myrank == rerank+1) then
        deallocate(Qm_cpu)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_3rd


  subroutine RungeKutta_4th(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=4), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step, flag_re = 0
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! Landau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(8), allocatable, device :: seed(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), Rs(5,nx-2,ny-2,nz-2))
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      Rs = 0.d0
      if (kind(id_LL) == 4) then
        allocate(seed(nx,ny,nz))
        call init_seed(nx, ny, nz, seed)
      endif
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    if (kind(id_rescale) == 4) then
      call pre_rescale(myrank, ny, nz, Qre, Qm, Qm_cpu)
    endif
    
    do t2 = 1, np
      do t1 = 1, nt
        step = np * (t2-1) + t1
        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(1, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          call calc_step(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(1, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(2, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          call calc_step(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(2, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(3, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          call calc_step(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJs)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(3, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif

        if (mod(myrank,2) == 0) then
          if (kind(id_rescale) == 4) then
            call step_rescale(4, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJs, Qm, Qre)
          endif
          if (kind(id_LL) == 4) then
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G, seed)
          else
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          endif
          call calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
          if (ndevices >= 2) then
            call exchange(id_rescale, myrank, nranks, overlap, nx, ny, nz, QJ)
          endif
          if (kind(id_rescale) == 4) then
            call wait_rescale(myrank, ireq, ireq2, istat, istat2)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
          else
            call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
          endif
        elseif (myrank == rerank+1 .and. kind(id_rescale) == 4) then
          call rescale_recv_send(4, flag_re, nx, ny, nz, step, y, Jacobian_cpu, Qm_cpu)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, Rs, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(seed)
      endif
    endif
    if (kind(id_rescale) == 4) then
      if (myrank == 0 .or. myrank == rerank) then
        deallocate(Qre, Qm)
      elseif (myrank == rerank+1) then
        deallocate(Qm_cpu)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine RungeKutta_4th
 

  subroutine Gauss_RungeKutta(id_RungeKutta, myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    use mod_globals, only : id_visc
    integer(kind=8), intent(in) :: id_RungeKutta
    integer, intent(in)         :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)         :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)         :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)         :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, itr, max_itr, t1, t2, overlap, ierr, nranks, ndevices, stat, ireq, ireqs(2)
    integer istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2)
    ! rescale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: step, flag_re = 0
    real(8) :: c1, c2, a11, a12, a21, a22, b1, b2, err, tol = 1.d-16
    real(8), allocatable, device :: Qre(:), Qm(:)
    real(8), allocatable, pinned :: Qm_cpu(:)
    ! GPU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, device :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable, device :: T(:,:,:), mu(:,:,:), mut(:,:,:), qc2(:,:,:)
    real(8), allocatable, device :: R1(:,:,:,:), R2(:,:,:,:), R1_new(:,:,:,:), R2_new(:,:,:,:)
    real(8), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! Landau !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(8), allocatable, device :: seed(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    ! count GPU
    stat = cudaGetDeviceCount(ndevices)
    if (myrank == 0) then
      print '(2x, i2, a)', ndevices, " GPU devices are found"
    endif
    if (mod(myrank,2) == 0) then
      call check_gpu(mygpu)
      call allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, mu, mut, qc2, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), R1(5,nx-2,ny-2,nz-2), R2(5,nx-2,ny-2,nz-2))
      allocate(R1_new(5,nx-2,ny-2,nz-2), R2_new(5,nx-2,ny-2,nz-2))
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      if (kind(id_LL) == 4) then
        allocate(seed(nx,ny,nz))
        call init_seed(nx, ny, nz, seed)
      endif
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif
    if (kind(id_rescale) == 4) then
      call pre_rescale(myrank, ny, nz, Qre, Qm, Qm_cpu)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "myrank is ", myrank, " start Runge-Kutta"
    c1  = 0.5d0 - sqrt(3.d0) / 6.d0
    c2  = 0.5d0 - sqrt(3.d0) / 6.d0
    a11 = 0.25d0
    a12 = 0.25d0 - sqrt(3.d0) / 6.d0
    a21 = 0.25d0 + sqrt(3.d0) / 6.d0
    a22 = 0.25d0
    b1  = 0.5d0
    b2  = 0.5d0
    max_itr = 100
    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          ! calc R1
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step1(nx, ny, nz, c1, dx, dy, dz, E, F, G, QJ, QJs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R1)
          ! calc R2
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_step1(nx, ny, nz, c2, dx, dy, dz, E, F, G, QJ, QJs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
          call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
          call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R2)
          do itr = 1, max_itr
            ! calc R1
            call calc_Gauss_step(nx, ny, nz, a11, a12, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
            call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R1_new)
            ! calc R2
            call calc_Gauss_step(nx, ny, nz, a21, a22, R1, R2, QJ, QJs)
            call set_bc(myrank, nx, ny, nz, Jacobian, QJs)
            call calc_EFG(id_visc, nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, mu, mut, qc2, E, F, G)
            call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R2_new)
            call calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
            if (err < tol) exit
            R1 = R1_new
            R2 = R2_new
          enddo
          if (err < tol) then
            print *, "Converged at itr=", itr
          else
            print *, "Didn't Converged error=", err
          endif
          call calc_Gauss_step_Q(nx, ny, nz, b1, b2, R1, R2, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ, Qre)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, mu, mut, qc2, QJ, QJs, R1, R2, R1_new, R2_new, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
      if (kind(id_LL) == 4) then
        deallocate(seed)
      endif
    endif
    print *, "myrank is ", myrank, " deallocate GPU memory"
  end subroutine Gauss_RungeKutta
end module calc_time_dev

