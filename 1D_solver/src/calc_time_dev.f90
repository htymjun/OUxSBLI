module calc_time_dev
  use cudafor
  use mpi
  use mod_globals, only : gamma, Pr, R_gas, nt, np, dt, dx, dtdx, blocks, threads
  use set
  use print
  implicit none
contains
  attributes(global) subroutine calc_E(nx, rho, u, p, E)
    integer, intent(in), value                 :: nx
    real(8), intent(in), dimension(nx), device :: rho, u, p
    real(8), intent(out), device               :: E(nx-1,3)
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    ! KEEP scheme
    E(i,1) = 0.25d0 * (rho(i) + rho(i+1)) * (u(i) + u(i+1))
    E(i,2) = E(i,1) * 0.5d0 * (u(i) + u(i+1)) + 0.5d0 * (p(i) + p(i+1))
    E(i,3) = 0.5d0 * E(i,1) * u(i) * u(i+1) &
           + E(i,1) * 0.5d0 * (p(i) / rho(i) + p(i+1) / rho(i+1)) / (gamma - 1.d0) &
           + 0.5d0 * (u(i) * p(i+1) + u(i+1) * p(i))
  end subroutine calc_E

  
  attributes(device) function Sutherland(T) result(mu)
    real(8), intent(in), value :: T
    real(8) :: mu, C = 1.461d-6, S = 110.3d0
    mu = C * T**1.5d0 / (T + S)
  end function Sutherland

  
  attributes(global) subroutine calc_Ev(nx, rho, u, p, Ev)
    integer, intent(in), value                 :: nx
    real(8), intent(in), dimension(nx), device :: rho, u, p
    real(8), intent(out), device               :: Ev(nx-1,3)
    real(8) :: mu, Cp = gamma * R_gas / (gamma - 1.d0)
    real(8), device :: T(2)
    real(8) txx, utxx, kTx
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    T(:) = p(i:i+1) / (R_gas * rho(i:i+1))
    mu   = 0.5d0 * (Sutherland(T(1)) + Sutherland(T(2)))
    txx  = 4.d0 * mu * (-u(i) + u(i+1)) / (3.d0 * dx)
    utxx = 0.5d0 * (u(i) + u(i+1)) * txx
    kTx  = Cp * mu * (-T(1) + T(2)) / (dx * Pr)
    Ev(i,2) = txx
    Ev(i,3) = utxx + kTx
  end subroutine calc_Ev


  attributes(global) subroutine calc_Ev_LL(nx, rho, u, p, Z, Ev)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: rho(nx), u(nx), p(nx)
    real(8), intent(in), device  :: Z(2,nx)            
    real(8), intent(out), device :: Ev(nx-1,3)
    real(8) :: mu, Cp = gamma * R_gas / (gamma - 1.d0)
    real(8) :: kb = 1.380650d-23
    real(8), device :: T(2)
    real(8) txx, utxx, kappa, kTx, s, q
    integer i
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    T(:)  = p(i:i+1) / (R_gas * rho(i:i+1))
    mu    = 0.5d0 * (Sutherland(T(1)) + Sutherland(T(2)))
    txx   = 4.d0 * mu * (-u(i) + u(i+1)) / (3.d0 * dx)
    utxx  = 0.5d0 * (u(i) + u(i+1)) * txx
    kappa = Cp * mu / Pr
    kTx   = kappa * (-T(1) + T(2)) / dx
    s     = sqrt(4.d0 * kb * mu * (T(1) + T(2)) / (3.d0 * dt * dx)) * 0.5d0 * (Z(1,i) + Z(1,i+1))
    q     = sqrt(kb * kappa * (T(1)**2 + T(2)**2) / (dt * dx)) * 0.5d0 * (Z(2,i) + Z(2,i+1))
    Ev(i,2) = txx + sqrt(2.d0) * s
    Ev(i,3) = utxx + kTx + sqrt(2.d0) * (q + 0.5d0 * (u(i) + u(i+1)) * s)
  end subroutine calc_Ev_LL


  subroutine calc_R(nx, Q, Z, R)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(nx,3)
    real(8), intent(inout), device :: Z(2,nx)
    real(8), intent(out), device   :: R(nx-2,3)
    real(8), dimension(nx), device :: rho, u, p
    integer stat, i, j
    real(8), dimension(nx-1,3), device :: E, Ev
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, nx
      rho(i) = Q(i,1)
      u(i)   = Q(i,2) / rho(i)
      p(i)   = (gamma - 1.d0) * (Q(i,3) - 0.5d0 * rho(i) * u(i)**2)
    enddo
    call calc_E<<<blocks,threads,0>>>(nx, rho, u, p, E)
    if (kind(id_LL) == 2) then
      call calc_Ev<<<blocks,threads,1>>>(nx, rho, u, p, Ev)
    else
      call calc_Ev_LL<<<blocks,threads,1>>>(nx, rho, u, p, Z, Ev)
    endif
    stat = cudaDeviceSynchronize() 
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, 3
      do i = 1, nx-2
        R(i,j) = dtdx * (-E(i,j) + Ev(i,j) + E(i+1,j) - Ev(i+1,j))
    enddo;enddo
  end subroutine calc_R


  subroutine calc_step1(nx, Q, R, Q2)
    integer, intent(in), value   :: nx
    real(8), intent(in), device  :: Q(nx,3)
    real(8), intent(in), device  :: R(nx-2,3)
    real(8), intent(out), device :: Q2(nx,3)
    integer i, j
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, 3
      do i = 2, nx-1
        Q2(i,j) = Q(i,j) - R(i-1,j)
    enddo;enddo
  end subroutine calc_step1


  subroutine calc_step2(nx, Q, R, Q2)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q(nx,3)
    real(8), intent(in), device    :: R(nx-2,3)
    real(8), intent(inout), device :: Q2(nx,3)
    integer i, j
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, 3
      do i = 2, nx-1
        Q2(i,j) = 0.25d0 * (3.d0 * Q(i,j) + Q2(i,j) - R(i-1,j))
    enddo;enddo
  end subroutine calc_step2
  

  subroutine calc_step3(nx, Q2, R, Q)
    integer, intent(in), value     :: nx
    real(8), intent(in), device    :: Q2(nx,3)
    real(8), intent(in), device    :: R(nx-2,3)
    real(8), intent(inout), device :: Q(nx,3)
    integer i, j
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, 3
      do i = 2, nx-1
        Q(i,j) = (Q(i,j) + 2.d0 * Q2(i,j) - 2.d0 * R(i-1,j)) / 3.d0
    enddo;enddo
  end subroutine calc_step3
  

  subroutine RungeKutta(myrank, nx, x, Q_cpu)
    integer, intent(in)    :: myrank, nx
    real(8), intent(in)    :: x(nx)
    real(8), intent(inout) :: Q_cpu(nx,3)
    integer t1, t2, ndevices, ilen, ierr, stat, ireq, istat(MPI_STATUS_SIZE)
    type(cudaDeviceProp)         :: prop
    real(8), allocatable, device :: Q(:,:), Q2(:,:), Z(:,:), R(:,:)

    if (myrank == 0) then
      stat = cudaGetDeviceCount(ndevices)
      print '(2x, i2, a)', ndevices, " GPU devices are found"
      stat = cudaSetDevice(myrank)
      stat = cudaGetDeviceProperties(prop,myrank)
      ilen = verify(prop%name, ' ', .true.)
      print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", myrank, ") is available"

      allocate(Q(nx,3), Q2(nx,3), Z(2,nx), R(nx-2,3))

      call print_1d(0, nx, real(x), real(Q_cpu))
      Q = Q_cpu
    endif

    do t2 = 1, np
      if (myrank == 0) then
        do t1 = 1, nt
          call calc_R(nx, Q, Z, R)
          call calc_step1(nx, Q, R, Q2)
          call set_bc(nx, Q2)

          call calc_R(nx, Q2, Z, R)
          call calc_step2(nx, Q, R, Q2)
          call set_bc(nx, Q2)

          call calc_R(nx, Q2, Z, R)
          call calc_step3(nx, Q2, R, Q)
          call set_bc(nx, Q)
        enddo
      endif

      ! send and recv device arrays
      if (myrank == 0) then
        Q_cpu = Q
        call MPI_SEND(Q_cpu, nx*3, MPI_REAL8, 1, 0, MPI_COMM_WORLD, ierr) 
      elseif (myrank == 1) then
        call MPI_RECV(Q_cpu, nx*3, MPI_REAL8, 0, 0, MPI_COMM_WORLD, istat, ierr)
        call print_1d(t2, nx, real(x), real(Q_cpu))
      endif
    enddo
    
    if (myrank == 0) then
      deallocate(Q, Q2, Z, R)
    endif
  end subroutine RungeKutta
end module calc_time_dev

