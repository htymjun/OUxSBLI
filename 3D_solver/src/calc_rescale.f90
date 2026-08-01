module calc_rescale
  use cudafor
  use mpi
  use mod_globals, only : nre1, nre2, rerank, nt, np, dt, gamma , R, Pr, u0, rho0, p0, M0, blt, start_rescale
  use mod_constant, only : Cp, gamma_1, over_gamma_1, mu0_T0_S_over_T0_2_3, over_T0, id_gpumpi, id_recal
  use cpu_gpu_mpi
contains
  subroutine calc_mean(step, ireq, flag_re, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qm_1, Qm_2, Qm_3, Qm_4, Qm_5)
    integer, intent(inout)         :: step, ireq
    integer, intent(in)            :: flag_re, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(in), device    :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout), device :: Qm_1(ny), Qm_2(ny), Qm_3(ny), Qm_4(ny), Qm_5(ny)
    real(8) Q1, Q2, Q3, Q4, Q5, rhoinv, Jacobian_tmp, volinv, step1, step2
    logical arrived
    integer i, k, istat, ierr
    volinv = 1.d0 / dble((nre2 - nre1 + 1) * (nz - 6))
    call MPI_TEST(ireq, arrived, MPI_STATUS_IGNORE, ierr)
    if (flag_re == 0) then
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        Q1 = 0.d0; Q2 = 0.d0; Q3 = 0.d0; Q4 = 0.d0; Q5 = 0.d0
        do k = 4, nz-3
          do i = nre1, nre2
            Jacobian_tmp = Jacobian(i,j)
            rhoinv = 1.d0 / QJ_1(i,j,k)
            Q1 = Q1 + QJ_1(i,j,k) * Jacobian_tmp
            Q2 = Q2 + QJ_2(i,j,k) * rhoinv
            Q3 = Q3 + QJ_3(i,j,k) * rhoinv
            Q4 = Q4 + QJ_4(i,j,k) * rhoinv
            Q5 = Q5 + gamma_1 * Jacobian_tmp * (QJ_5(i,j,k) &
                    - 0.5d0 * (QJ_2(i,j,k)**2 + QJ_3(i,j,k)**2 + QJ_4(i,j,k)**2) * rhoinv)
        enddo;enddo
        Qm_1(j) = Q1 * volinv
        Qm_2(j) = Q2 * volinv
        Qm_3(j) = Q3 * volinv
        Qm_4(j) = Q4 * volinv
        Qm_5(j) = Q5 * volinv
      enddo
    else
      step1 = dble(step-1); step2 = 1.d0 / dble(step)
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        Q1 = 0.d0; Q2 = 0.d0; Q3 = 0.d0; Q4 = 0.d0; Q5 = 0.d0
        do k = 4, nz-3
          do i = nre1, nre2
            Jacobian_tmp = Jacobian(i,j)
            rhoinv = 1.d0 / QJ_1(i,j,k)
            Q1 = Q1 + QJ_1(i,j,k) * Jacobian_tmp
            Q2 = Q2 + QJ_2(i,j,k) * rhoinv
            Q3 = Q3 + QJ_3(i,j,k) * rhoinv
            Q4 = Q4 + QJ_4(i,j,k) * rhoinv
            Q5 = Q5 + gamma_1 * Jacobian_tmp * (QJ_5(i,j,k) &
                    - 0.5d0 * (QJ_2(i,j,k)**2 + QJ_3(i,j,k)**2 + QJ_4(i,j,k)**2) * rhoinv)
        enddo;enddo
        Qm_1(j) = (step1 * Qm_1(j) + Q1 * volinv) * step2
        Qm_2(j) = (step1 * Qm_2(j) + Q2 * volinv) * step2
        Qm_3(j) = (step1 * Qm_3(j) + Q3 * volinv) * step2
        Qm_4(j) = (step1 * Qm_4(j) + Q4 * volinv) * step2
        Qm_5(j) = (step1 * Qm_5(j) + Q5 * volinv) * step2
      enddo
      step = step + 1
    endif
    istat = cudaDeviceSynchronize() ! This line is necessary for next MPI_SEND
  end subroutine calc_mean


  subroutine copy(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in)          :: nx, ny, nz
    real(8), intent(in), device  :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(out), device :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    integer j, k, offset
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz-6
      do j = 1, ny
        offset = ny*(k-1) + j
        Qre_1(offset) = QJ_1(nre2,j,k+3)
        Qre_2(offset) = QJ_2(nre2,j,k+3)
        Qre_3(offset) = QJ_3(nre2,j,k+3)
        Qre_4(offset) = QJ_4(nre2,j,k+3)
        Qre_5(offset) = QJ_5(nre2,j,k+3)
    enddo;enddo
  end subroutine copy


  subroutine step_rescale(num, myrank, nx, ny, nz, step, flag_re, flag_req, ireq, ireq2, Jacobian, &
                           QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qm_1, Qm_2, Qm_3, Qm_4, Qm_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(in)            :: num, myrank, nx, ny, nz
    integer, intent(inout)         :: step, flag_re, flag_req, ireq, ireq2(2)
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(in), device    :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout), device :: Qm_1(ny), Qm_2(ny), Qm_3(ny), Qm_4(ny), Qm_5(ny)
    real(8), intent(inout), device :: Qre_1(ny*(nz-6)), Qre_2(ny*(nz-6)), Qre_3(ny*(nz-6)), Qre_4(ny*(nz-6)), Qre_5(ny*(nz-6))
    real(8), device :: Qre_flat(ny*(nz-6)*5), Qm_flat(ny*5)
    integer ierr, j, n
    n = ny*(nz-6)
    if (myrank == rerank) then
      call copy(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
      Qre_flat(1:n)       = Qre_1
      Qre_flat(n+1:2*n)   = Qre_2
      Qre_flat(2*n+1:3*n) = Qre_3
      Qre_flat(3*n+1:4*n) = Qre_4
      Qre_flat(4*n+1:5*n) = Qre_5
      call CPUGPU_MPI_SEND(id_gpumpi, Qre_flat, 5*ny*(nz-6), rerank+1, 0, MPI_COMM_WORLD, ireq2(1), ierr)
      if (num == 1 .and. id_recal == 0) then
        call calc_mean(step, flag_req, flag_re, nx, ny, nz, Jacobian, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Qm_1, Qm_2, Qm_3, Qm_4, Qm_5)
        Qm_flat(1:ny)        = Qm_1
        Qm_flat(ny+1:2*ny)   = Qm_2
        Qm_flat(2*ny+1:3*ny) = Qm_3
        Qm_flat(3*ny+1:4*ny) = Qm_4
        Qm_flat(4*ny+1:5*ny) = Qm_5
        call CPUGPU_MPI_SEND(id_gpumpi, Qm_flat, 5*ny, rerank+1, 1, MPI_COMM_WORLD, ireq2(2), ierr)
      endif
    endif
    if (myrank == 0) then
      call CPUGPU_MPI_RECV(id_gpumpi, Qre_flat, 5*ny*(nz-6), rerank+1, 0, MPI_COMM_WORLD, ireq, ierr)
      Qre_1 = Qre_flat(1:n)
      Qre_2 = Qre_flat(n+1:2*n)
      Qre_3 = Qre_flat(2*n+1:3*n)
      Qre_4 = Qre_flat(3*n+1:4*n)
      Qre_5 = Qre_flat(4*n+1:5*n)
    endif
  end subroutine step_rescale


  subroutine wait_rescale(myrank, ireq, ireq2, istat, istat2)
    integer, intent(in)    :: myrank
    integer, intent(inout) :: ireq, ireq2(2), istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    integer ierr
    if (myrank == rerank) then
      call MPI_WAITALL(2, ireq2, istat2, ierr)
    endif
    if (myrank == 0) then
      call MPI_WAIT(ireq, istat, ierr)
    endif
  end subroutine wait_rescale


  subroutine rescale_recv_send(num, flag_re, nx, ny, nz, step, y, Jacobian, Qm_cpu_1, Qm_cpu_2, Qm_cpu_3, Qm_cpu_4, Qm_cpu_5)
    integer, intent(in)    :: num
    integer, intent(inout) :: flag_re
    integer, intent(in)    :: nx, ny, nz, step
    real(8), intent(in)    :: y(ny), Jacobian(nx,ny)
    real(8), intent(inout) :: Qm_cpu_1(ny), Qm_cpu_2(ny), Qm_cpu_3(ny), Qm_cpu_4(ny), Qm_cpu_5(ny)
    real(8)         :: Qre_cpu(ny*(nz-6)*5), Qm_cpu_flat(ny*5), bltre
    real(8)         :: Qre_cpu_1(ny*(nz-6)), Qre_cpu_2(ny*(nz-6)), Qre_cpu_3(ny*(nz-6)), Qre_cpu_4(ny*(nz-6)), Qre_cpu_5(ny*(nz-6))
    real(8), device ::     Qre(ny*(nz-6)*5), Qm(ny*5)
    integer stat, errorcode, ierr, ireq, ireqs(2), flag_req, n
    integer istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2), j
    real(8) t
    character(len=40) filename
    write(filename, "(a)") "data/rescaling.d"
    n = ny*(nz-6)

    if (num == 1 .and. id_recal == 0) then
      call CPUGPU_MPI_RECV(id_gpumpi, Qre, 5*ny*(nz-6), rerank, 0, MPI_COMM_WORLD, ireqs(1), ierr)
      call CPUGPU_MPI_RECV(id_gpumpi, Qm,  5*ny,        rerank, 1, MPI_COMM_WORLD, ireqs(2), ierr)
      call MPI_WAITALL(2, ireqs, istats, ierr)
      stat = cudaDeviceSynchronize()
      stat = cudaMemcpy(Qm_cpu_flat, Qm, 5*ny, cudaMemcpyDeviceToHost)
      if (stat /= cudaSuccess) then
        print *, "Qm  cudaMemcpy failed:", trim(cudaGetErrorString(stat))
      endif
      Qm_cpu_1 = Qm_cpu_flat(1:ny)
      Qm_cpu_2 = Qm_cpu_flat(ny+1:2*ny)
      Qm_cpu_3 = Qm_cpu_flat(2*ny+1:3*ny)
      Qm_cpu_4 = Qm_cpu_flat(3*ny+1:4*ny)
      Qm_cpu_5 = Qm_cpu_flat(4*ny+1:5*ny)
    else
      call CPUGPU_MPI_RECV(id_gpumpi, Qre, 5*ny*(nz-6), rerank, 0, MPI_COMM_WORLD, ireq, ierr)
      stat = cudaDeviceSynchronize()
    endif

    stat = cudaMemcpy(Qre_cpu, Qre, 5*ny*(nz-6), cudaMemcpyDeviceToHost)
    if (stat /= cudaSuccess) then
      print *, "Qre cudaMemcpy failed:", trim(cudaGetErrorString(stat))
    endif
    stat = cudaDeviceSynchronize()
    Qre_cpu_1 = Qre_cpu(1:n)
    Qre_cpu_2 = Qre_cpu(n+1:2*n)
    Qre_cpu_3 = Qre_cpu(2*n+1:3*n)
    Qre_cpu_4 = Qre_cpu(3*n+1:4*n)
    Qre_cpu_5 = Qre_cpu(4*n+1:5*n)
    call set_rescale(flag_re, step, nx, ny, nz-6, y, Jacobian, Qm_cpu_1, Qm_cpu_2, Qm_cpu_3, Qm_cpu_4, Qm_cpu_5, bltre, &
                      Qre_cpu_1, Qre_cpu_2, Qre_cpu_3, Qre_cpu_4, Qre_cpu_5)
    Qre_cpu(1:n)       = Qre_cpu_1
    Qre_cpu(n+1:2*n)   = Qre_cpu_2
    Qre_cpu(2*n+1:3*n) = Qre_cpu_3
    Qre_cpu(3*n+1:4*n) = Qre_cpu_4
    Qre_cpu(4*n+1:5*n) = Qre_cpu_5
    stat = cudaMemcpy(Qre, Qre_cpu, 5*ny*(nz-6), cudaMemcpyHostToDevice)
    call CPUGPU_MPI_SEND(id_gpumpi, Qre, 5*ny*(nz-6), 0, 0, MPI_COMM_WORLD, ireq, ierr)
    if (flag_re == 1) then
      call MPI_ISEND(flag_re, 1, MPI_INTEGER, rerank, 1001, MPI_COMM_WORLD, flag_req, ierr)
      call MPI_WAIT(flag_req, istat, ierr)
      print *, "Start calculating time average for rescale"
    endif
    if (num == 1) then
      if (bltre == 0.d0) then
        print *, "Invalid boundary layer thickness was detected"
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      endif
      t = nt * step * dt
      if (flag_re >= 1 .and. step >= start_rescale) then
        open(10, file=filename, position="append")
        write(10, "(2e12.4, a)") t, bltre/blt, " rescale"
        close(10)
      else
        open(10, file=filename, position="append")
        write(10, "(2e12.4, a)") t, bltre/blt, " cyclic"
        close(10)
      endif
    endif
  end subroutine rescale_recv_send


  subroutine write_Qm(ny, step, y, Qm_1, Qm_2, Qm_3, Qm_4, Qm_5)
    integer, intent(in) :: ny, step
    real(8), intent(in) :: y(ny), Qm_1(ny), Qm_2(ny), Qm_3(ny), Qm_4(ny), Qm_5(ny)
    character(len=40) filename
    integer j
    write(filename, "(a, i5.5, a)") "recal/Qm", int(step), ".d"
    open(10, file=filename, status="replace", action="write")
    write(10, "(a)") "y    rho   u     v     w     p"
    do j = 1, ny
      write(10, "(6e12.4)") y(j), Qm_1(j), Qm_2(j), Qm_3(j), Qm_4(j), Qm_5(j)
    enddo
    close(10)
    if (step == np) then
      write(filename, "(a)") "recal/Qm.dat"
      open(10, file=filename, status="replace", action="write", form="unformatted", access="stream")
      write(10) Qm_1, Qm_2, Qm_3, Qm_4, Qm_5
      close(10)
    endif
  end subroutine write_Qm


  subroutine set_rescale(flag_re, step, nx, ny, nz, y, Jacobian, Qm_1, Qm_2, Qm_3, Qm_4, Qm_5, bltre, &
                          Qre_1, Qre_2, Qre_3, Qre_4, Qre_5)
    integer, intent(inout) :: flag_re
    integer, intent(in)    :: step, nx, ny, nz! nz-6
    real(8), intent(in)    :: y(ny), Jacobian(nx,ny)
    real(8), intent(in)    :: Qm_1(ny), Qm_2(ny), Qm_3(ny), Qm_4(ny), Qm_5(ny)
    real(8), intent(out)   :: bltre
    real(8), intent(inout) :: Qre_1(ny,nz), Qre_2(ny,nz), Qre_3(ny,nz), Qre_4(ny,nz), Qre_5(ny,nz) ! Q / J
    integer i, j, jup, jdown, jj, k, kh, ierr, errorcode, flag
    integer, dimension(ny) :: jj_y, jj_e
    real(8) t, u99, dudy, taure, utre, utin, beta, mu, nu, ady, ade, one_ady, one_ade, one_weight
    ! mean properties at rescaling plane
    real(8), dimension(ny)    :: Um, Vm, Wm, rhom, Tm, pm
    ! fluctuating properties at rescaling plane
    real(8), dimension(ny,nz) :: ufre, vfre, wfre, Tfre, pfre
    real(8), dimension(ny)    :: ypre, ypin, etre, etin
    ! fluctuating properties at both inner and outer region
    real(8), dimension(ny,nz) :: ufin, vfin, wfin, Tfin, pfin
    real(8), dimension(ny,nz) :: ufout, vfout, wfout, Tfout, pfout
    ! mean properties at both inner and outer region
    real(8), dimension(ny)    :: Umin, Vmin, pmin, Tmin, Umout, Vmout, pmout, Tmout
    ! weighting function
    real(8), dimension(ny)    :: weight
    ! properties at rescaling plane
    real(8) ure, vre, wre, rhore, Tre, pre
    ! rescaled properties at inlet
    real(8) uin, vin, win, rhoin, Tin, pin
    ! cache
    real(8) :: u_tmp, v_tmp, p_tmp, T_tmp, weight_tmp, Jacobian_tmp, over_rhore, utin_nu, utre_nu
    real(8) :: over_blt = 1.d0 / blt
    ! for exception
    real(8) :: blt_min = 0.5d0 * blt, blt_max = 2.d0 * blt
    do k = 1, nz
      do j = 1, ny
        Qre_1(j,k) = Qre_1(j,k) * Jacobian(nre2,j)
        Qre_2(j,k) = Qre_2(j,k) * Jacobian(nre2,j)
        Qre_3(j,k) = Qre_3(j,k) * Jacobian(nre2,j)
        Qre_4(j,k) = Qre_4(j,k) * Jacobian(nre2,j)
        Qre_5(j,k) = Qre_5(j,k) * Jacobian(nre2,j)
    enddo;enddo

    do j = 1, ny
      rhom(j)  = Qm_1(j)
      u_tmp    = Qm_2(j)
      v_tmp    = Qm_3(j)
      Wm(j)    = Qm_4(j)
      p_tmp    = Qm_5(j)
      T_tmp    = p_tmp / (R * rhom(j))
      Um(j)    = u_tmp
      Umin(j)  = u_tmp
      Umout(j) = u_tmp
      Vm(j)    = v_tmp
      Vmin(j)  = v_tmp
      Vmout(j) = v_tmp
      pm(j)    = p_tmp
      pmin(j)  = p_tmp
      pmout(j) = p_tmp
      Tm(j)    = T_tmp
      Tmin(j)  = T_tmp
      Tmout(j) = T_tmp
    enddo

    ! free stream
    jup = -1
    do j = 1, ny
      if (y(j) >= 2.d0 * blt) then
        jup = j
      endif
    enddo
    if (jup < 0) then
      print *, "Ly is not enough"
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
    endif
    u99 = 0.99d0 * sum(Um(jup:ny)) / size(Um(jup:ny))
    if (abs(u99 - 0.99d0 * u0) > 0.01d0 * u0) then
      print *, "Free stream streamwise velocity is not constant"
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
    endif

    ! blt up
    jup   = -1
    bltre = 0.d0
    do j = 2, ny
      if (Um(j) >= u99) then
        dudy  = (u99 - Um(j-1)) / (-Um(j-1) + Um(j) + 1.d-15)
        bltre = y(j-1) + (-y(j-1) + y(j)) * dudy
        jup   = j
        exit
      endif
    enddo

    flag = 1
    if (jup < 0) then
      print *, "No 99% doundary layer thickness"
      call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
    endif
    if (bltre < blt_min .or. bltre > blt_max) flag = 0

    if (bltre > blt) then
      flag_re = flag_re + 1
    endif

    if (bltre > blt .and. flag_re >= 1 .and. step >= start_rescale .and. flag == 1) then
      ufin(:,:)  = 0.d0; vfin(:,:)  = 0.d0; wfin(:,:)  = 0.d0; Tfin(:,:)  = 0.d0; pfin(:,:)  = 0.d0
      ufout(:,:) = 0.d0; vfout(:,:) = 0.d0; wfout(:,:) = 0.d0; Tfout(:,:) = 0.d0; pfout(:,:) = 0.d0
      do k = 1, nz
        do j = 1, ny
          rhore = Qre_1(j,k)
          over_rhore = 1.d0 / rhore
          ure   = Qre_2(j,k) * over_rhore
          vre   = Qre_3(j,k) * over_rhore
          wre   = Qre_4(j,k) * over_rhore
          pre   = gamma_1 * (Qre_5(j,k) - 0.5d0 * rhore * (ure**2 + vre**2 + wre**2))
          Tre   = pre / (rhore * R)
          ufre(j,k) = ure - Um(j)
          vfre(j,k) = vre - Vm(j)
          wfre(j,k) = wre - Wm(j)
          Tfre(j,k) = Tre - Tm(j)
          pfre(j,k) = pre - pm(j)
      enddo;enddo

      ! friction velocity
      mu    = mu0_T0_S_over_T0_2_3 / (Tm(1) + 111.d0) * Tm(1)**1.5
      nu    = mu / rhom(1)
      taure = mu * abs(-Um(1) + Um(2)) / (-y(1) + y(2))
      utre  = sqrt(taure / rhom(1))
      beta  = (bltre * over_blt)**0.1
      utin  = beta * utre

      utin_nu = utin / nu
      utre_nu = utre / nu
      do j = 1, ny
        ypin(j) = y(j) * utin_nu
        ypre(j) = y(j) * utre_nu
        etin(j) = y(j) * over_blt
        etre(j) = y(j) / bltre
        weight(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / (0.6d0 * etin(j) + 0.2d0)) / tanh(4.d0)))
      enddo

      jj_y(:) = 0
      do j = 1, ny
        do jj = 2, ny
          if (ypre(jj) >= ypin(j)) then
            ady     = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            one_ady = 1.d0 - ady
            ! mean
            Umin(j) = beta * (one_ady * Um(jj-1) + ady * Um(jj))!Um(jj-1) + ady * (-Um(jj-1) + Um(jj))
            Vmin(j) =         one_ady * Vm(jj-1) + ady * Vm(jj) !Vm(jj-1) + ady * (-Vm(jj-1) + Vm(jj))
            Tmin(j) =         one_ady * Tm(jj-1) + ady * Tm(jj) !Tm(jj-1) + ady * (-Tm(jj-1) + Tm(jj))
            pmin(j) =         one_ady * pm(jj-1) + ady * pm(jj) !pm(jj-1) + ady * (-pm(jj-1) + pm(jj))
            jj_y(j) = jj
            exit
          endif
      enddo;enddo

      jj_e(:) = 0
      do j = 1, ny
        do jj = 2, ny
          if (etre(jj) >= etin(j)) then
            ade     = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            one_ade = 1.d0 - ade
            ! mean
            Umout(j) = beta * (one_ade * Um(jj-1) + ade * Um(jj)) + (1.d0 - beta) * u0
            Vmout(j) =         one_ade * Vm(jj-1) + ade * Vm(jj) !Vm(jj-1) + ade * (-Vm(jj-1) + Vm(jj))
            Tmout(j) =         one_ade * Tm(jj-1) + ade * Tm(jj) !Tm(jj-1) + ade * (-Tm(jj-1) + Tm(jj))
            pmout(j) =         one_ade * pm(jj-1) + ade * pm(jj) !pm(jj-1) + ade * (-pm(jj-1) + pm(jj))
            jj_e(j)  = jj
            exit
          endif
      enddo;enddo

      do k = 1, nz
        do j = 2, ny
          jj = jj_y(j)
          if (jj > 1) then
            ady     = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            one_ady = 1.d0 - ady
            ufin(j,k) = beta * (one_ady * ufre(jj-1,k) + ady * ufre(jj,k))!ufre(jj-1,k) + ady * (-ufre(jj-1,k) + ufre(jj,k))
            vfin(j,k) = beta * (one_ady * vfre(jj-1,k) + ady * vfre(jj,k))!vfre(jj-1,k) + ady * (-vfre(jj-1,k) + vfre(jj,k))
            wfin(j,k) = beta * (one_ady * wfre(jj-1,k) + ady * wfre(jj,k))!wfre(jj-1,k) + ady * (-wfre(jj-1,k) + wfre(jj,k))
            Tfin(j,k) =         one_ady * Tfre(jj-1,k) + ady * Tfre(jj,k) !Tfre(jj-1,k) + ady * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfin(j,k) =         one_ady * pfre(jj-1,k) + ady * pfre(jj,k) !pfre(jj-1,k) + ady * (-pfre(jj-1,k) + pfre(jj,k))
          endif
          jj = jj_e(j)
          if (jj > 1) then
            ade     = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            one_ade = 1.d0 - ade
            ufout(j,k) = beta * (one_ade * ufre(jj-1,k) + ade * ufre(jj,k))!ufre(jj-1,k) + ade * (-ufre(jj-1,k) + ufre(jj,k))
            vfout(j,k) = beta * (one_ade * vfre(jj-1,k) + ade * vfre(jj,k))!vfre(jj-1,k) + ade * (-vfre(jj-1,k) + vfre(jj,k))
            wfout(j,k) = beta * (one_ade * wfre(jj-1,k) + ade * wfre(jj,k))!wfre(jj-1,k) + ade * (-wfre(jj-1,k) + wfre(jj,k))
            Tfout(j,k) =         one_ade * Tfre(jj-1,k) + ade * Tfre(jj,k) !Tfre(jj-1,k) + ade * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfout(j,k) =         one_ade * pfre(jj-1,k) + ade * pfre(jj,k) !pfre(jj-1,k) + ade * (-pfre(jj-1,k) + pfre(jj,k))
          endif
      enddo;enddo

      ! re-introducing
      do k = 1, nz
        kh = mod(k+nz/2,nz) + 1
        do j = 1, ny
          weight_tmp   = weight(j)
          one_weight   = 1.d0 - weight_tmp
          Jacobian_tmp = 1.d0 / Jacobian(nre2,j)
          uin = (Umin(j) + ufin(j,kh)) * one_weight + (Umout(j) + ufout(j,kh)) * weight_tmp
          vin = (Vmin(j) + vfin(j,kh)) * one_weight + (Vmout(j) + vfout(j,kh)) * weight_tmp
          win =            wfin(j,kh)  * one_weight +             wfout(j,kh)  * weight_tmp
          Tin = (Tmin(j) + Tfin(j,kh)) * one_weight + (Tmout(j) + Tfout(j,kh)) * weight_tmp
          pin = (pmin(j) + pfin(j,kh)) * one_weight + (pmout(j) + pfout(j,kh)) * weight_tmp
          rhoin = pin / (R * Tin)
          Qre_1(j,k) = rhoin * Jacobian_tmp
          Qre_2(j,k) = rhoin * uin * Jacobian_tmp
          Qre_3(j,k) = rhoin * vin * Jacobian_tmp
          Qre_4(j,k) = rhoin * win * Jacobian_tmp
          Qre_5(j,k) = (pin * over_gamma_1 + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)) * Jacobian_tmp
      enddo;enddo
    else
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = 1, nz
        do j = 1, ny
          Qre_1(j,k) = Qre_1(j,k) / Jacobian(nre2,j)
          Qre_2(j,k) = Qre_2(j,k) / Jacobian(nre2,j)
          Qre_3(j,k) = Qre_3(j,k) / Jacobian(nre2,j)
          Qre_4(j,k) = Qre_4(j,k) / Jacobian(nre2,j)
          Qre_5(j,k) = Qre_5(j,k) / Jacobian(nre2,j)
      enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale
