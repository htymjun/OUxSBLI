module calc_rescale
  use cudafor
  use mpi
  use mod_globals, only : nre1, nre2, rerank, nt, dt, gamma , R, Pr, u0, rho0, p0, M0, blt, start_rescale
  use mod_constant, only : Cp, gamma_1, over_gamma_1, mu0_T0_S, over_T0
contains
  subroutine calc_mean(step, flag_re, nx, ny, nz, Jacobian, QJ, Qm)
    integer, intent(in)            :: step, flag_re, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(inout), device :: Qm(ny*5)
    real(8) Q1, Q2, Q3, Q4, Q5, rhoinv, Jacobian_tmp, volinv, step1, step2
    integer i, k
    volinv = 1.d0 / dble((nre2 - nre1 + 1) * (nz - 6))
    if (flag_re == 0) then
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        Q1 = 0.d0; Q2 = 0.d0; Q3 = 0.d0; Q4 = 0.d0; Q5 = 0.d0
        do k = 4, nz-3
          do i = nre1, nre2
            Jacobian_tmp = Jacobian(i,j)
            rhoinv = 1.d0 / QJ(1,i,j,k)
            Q1 = Q1 + QJ(1,i,j,k) * Jacobian_tmp
            Q2 = Q2 + QJ(2,i,j,k) * rhoinv
            Q3 = Q3 + QJ(3,i,j,k) * rhoinv
            Q4 = Q4 + QJ(4,i,j,k) * rhoinv
            Q5 = Q5 + gamma_1 * Jacobian_tmp * (QJ(5,i,j,k) &
                    - 0.5d0 * (QJ(2,i,j,k)**2 + QJ(3,i,j,k)**2 + QJ(4,i,j,k)**2) * rhoinv)
        enddo;enddo
        Qm(5*(j-1)+1) = Q1 * volinv
        Qm(5*(j-1)+2) = Q2 * volinv
        Qm(5*(j-1)+3) = Q3 * volinv
        Qm(5*(j-1)+4) = Q4 * volinv
        Qm(5*(j-1)+5) = Q5 * volinv
      enddo
    else
      step1 = dble(step-1); step2 = 1.d0 / dble(step)
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        Q1 = 0.d0; Q2 = 0.d0; Q3 = 0.d0; Q4 = 0.d0; Q5 = 0.d0
        do k = 4, nz-3
          do i = nre1, nre2
            Jacobian_tmp = Jacobian(i,j)
            rhoinv = 1.d0 / QJ(1,i,j,k)
            Q1 = Q1 + QJ(1,i,j,k) * Jacobian_tmp
            Q2 = Q2 + QJ(2,i,j,k) * rhoinv
            Q3 = Q3 + QJ(3,i,j,k) * rhoinv
            Q4 = Q4 + QJ(4,i,j,k) * rhoinv
            Q5 = Q5 + gamma_1 * Jacobian_tmp * (QJ(5,i,j,k) &
                    - 0.5d0 * (QJ(2,i,j,k)**2 + QJ(3,i,j,k)**2 + QJ(4,i,j,k)**2) * rhoinv)
        enddo;enddo
        Qm(5*(j-1)+1) = (step1 * Qm(5*(j-1)+1) + Q1 * volinv) * step2
        Qm(5*(j-1)+2) = (step1 * Qm(5*(j-1)+2) + Q2 * volinv) * step2
        Qm(5*(j-1)+3) = (step1 * Qm(5*(j-1)+3) + Q3 * volinv) * step2
        Qm(5*(j-1)+4) = (step1 * Qm(5*(j-1)+4) + Q4 * volinv) * step2
        Qm(5*(j-1)+5) = (step1 * Qm(5*(j-1)+5) + Q5 * volinv) * step2
      enddo
    endif
  end subroutine calc_mean

  subroutine copy(nx, ny, nz, QJ, Qre)
    integer, intent(in)          :: nx, ny, nz
    real(8), intent(in), device  :: QJ(5,nx,ny,nz)
    real(8), intent(out), device :: Qre(ny*(nz-6)*5)
    integer j, k, l, j_offset, k_offset
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz-6
      k_offset = ny * 5 * (k-1)
      do j = 1, ny
        j_offset = 5 * (j-1)
        do l = 1, 5
          Qre(k_offset+j_offset+l) = QJ(l,nre2,j,k+3)
    enddo;enddo;enddo
  end subroutine copy
   
  subroutine step_rescale(num, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, Jacobian, QJ, Qm, Qre)
    integer, intent(in)            :: num, myrank, step, nx, ny, nz
    integer, intent(inout)         :: flag_re, ireq, ireq2(2)
    real(8), intent(in), device    :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(inout), device :: Qm(ny*5), Qre(ny*(nz-6)*5)
    real(8) Qm_cpu(ny*5)
    integer ierr, j
    if (myrank == rerank) then
      call copy(nx, ny, nz, QJ, Qre)
      call MPI_ISEND(Qre, 5*ny*(nz-6), MPI_REAL8, rerank+1, 0, MPI_COMM_WORLD, ireq2(1), ierr)
      if (num == 1) then
        call calc_mean(step, flag_re, nx, ny, nz, Jacobian, QJ, Qm)
        Qm_cpu = Qm
        !do j = 1, ny
        !  print *, "send j=", j, "Q", Qm_cpu(5*(j-1)+1), Qm_cpu(5*(j-1)+2)
        !enddo
        call MPI_ISEND(Qm, 5*ny, MPI_REAL8, rerank+1, 1, MPI_COMM_WORLD, ireq2(2), ierr)
      endif
      !print *, "myrank=", myrank, "send Qre"
    endif
    if (myrank == 0) then
      call MPI_IRECV(Qre, 5*ny*(nz-6), MPI_REAL8, rerank+1, 0, MPI_COMM_WORLD, ireq, ierr)
    endif
  end subroutine step_rescale

  subroutine wait_rescale(myrank, ireq, ireq2, istat, istat2)
    integer, intent(in)    :: myrank
    integer, intent(inout) :: ireq, ireq2(2), istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    integer ierr
    if (myrank == rerank) then
      call MPI_WAITALL(2, ireq2, istat2, ierr)
      !print *, "myrank=", myrank, "WAITALL send Qm, Qre"
    endif
    if (myrank == 0) then
      call MPI_WAIT(ireq, istat, ierr)
      !print *, "myrank=", myrank, "WAIT recv Qre"
    endif
  end subroutine wait_rescale

  subroutine rescale_recv_send(num, flag_re, nx, ny, nz, step, y, Jacobian, Qm_cpu)
    integer, intent(in)    :: num
    integer, intent(inout) :: flag_re
    integer, intent(in)    :: nx, ny, nz, step
    real(8), intent(in)    :: y(ny), Jacobian(nx,ny)
    real(8), intent(inout) :: Qm_cpu(ny*5)
    real(8)         :: Qre_cpu(ny*(nz-6)*5), bltre
    real(8), device ::     Qre(ny*(nz-6)*5), Qm(ny*5)
    integer stat, errorcode, ierr, ireq, ireqs(2), istat(MPI_STATUS_SIZE), istats(MPI_STATUS_SIZE,2), j
    real(8) t
    character(len=40) filename
    write(filename, "(a)") "data/rescaling.d"

    if (num == 1) then
      call MPI_IRECV(Qre, 5*ny*(nz-6), MPI_REAL8, rerank, 0, MPI_COMM_WORLD, ireqs(1), ierr)
      call MPI_IRECV(Qm,  5*ny,        MPI_REAL8, rerank, 1, MPI_COMM_WORLD, ireqs(2), ierr)
      call MPI_WAITALL(2, ireqs, istats, ierr)
      stat = cudaDeviceSynchronize()
      stat = cudaMemcpy(Qm_cpu,  Qm,  5*ny,        cudaMemcpyDeviceToHost)
      if (stat /= cudaSuccess) then
        print *, "Qm  cudaMemcpy failed:", trim(cudaGetErrorString(stat))
      endif
    else
      call MPI_IRECV(Qre, 5*ny*(nz-6), MPI_REAL8, rerank, 0, MPI_COMM_WORLD, ireq, ierr)
      call MPI_WAIT(ireq, istat, ierr)
      stat = cudaDeviceSynchronize()
    endif

    stat = cudaMemcpy(Qre_cpu, Qre, 5*ny*(nz-6), cudaMemcpyDeviceToHost)
    if (stat /= cudaSuccess) then
      print *, "Qre cudaMemcpy failed:", trim(cudaGetErrorString(stat))
    endif
    stat = cudaDeviceSynchronize()
    call set_rescale(flag_re, step, nx, ny, nz-6, y, Jacobian, Qm_cpu, bltre, Qre_cpu)
    stat = cudaMemcpy(Qre, Qre_cpu, 5*ny*(nz-6), cudaMemcpyHostToDevice)
    call MPI_ISEND(Qre, 5*ny*(nz-6), MPI_REAL8, 0, 0, MPI_COMM_WORLD, ireq, ierr)
    call MPI_WAIT(ireq, istat, ierr)
    if (flag_re == 1) then
      call MPI_BCAST(flag_re, 1, MPI_INTEGER, rerank+1, MPI_COMM_WORLD, ierr)
    endif
    if (num == 1) then
      if (bltre == 0.d0) then
        print *, "Invalid boundary layer thickness was detected"
        call MPI_ABORT(MPI_COMM_WORLD, errorcode, ierr)
      endif
      t = nt * step * dt
      if (flag_re >= 1 .and. step >= start_rescale) then
        open(10, file=filename, position="append")
        write(10, "(2e12.4, a)") t*1d3, bltre, "rescale"
        close(10)
      else
        open(10, file=filename, position="append")
        write(10, "(2e12.4, a)") t*1d3, bltre, "cyclic"
        close(10)
      endif
    endif
  end subroutine rescale_recv_send

  subroutine set_rescale(flag_re, step, nx, ny, nz, y, Jacobian, Qm, bltre, Qre)
    integer, intent(inout) :: flag_re
    integer, intent(in)    :: step, nx, ny, nz! nz-6
    real(8), intent(in)    :: y(ny), Jacobian(nx,ny), Qm(ny*5)
    real(8), intent(out)   :: bltre
    real(8), intent(inout) :: Qre(ny*nz*5) ! Q / J
    integer i, j, jj, k, kh, l, j_offset, k_offset, ierr
    integer, dimension(ny) :: jj_y, jj_e
    real(8) t, dudy, taure, utre, utin, beta, mu, nu, ady, ade 
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
    real(8) :: u_tmp, v_tmp, p_tmp, T_tmp, weight_tmp, Jacobian_tmp, over_rhore
    do k = 1, nz
      k_offset = ny * 5 * (k-1)
      do j = 1, ny
        j_offset = 5 * (j-1)
        do l = 1, 5
          i = k_offset + j_offset + l
          Qre(i) = Qre(i) * Jacobian(nre2,j)
    enddo;enddo;enddo

    do j = 1, ny
      rhom(j)  = Qm(5*(j-1)+1)
      u_tmp    = Qm(5*(j-1)+2)
      v_tmp    = Qm(5*(j-1)+3)
      Wm(j)    = Qm(5*(j-1)+4)
      p_tmp    = Qm(5*(j-1)+5)
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

    ! check boundary layer thickness at rescaling plane
    bltre = 0.d0
    do j = 2, ny
      if (Um(j) >= 0.99d0 * u0) then
        dudy  = (Um(j) - 0.99d0 * u0) / (-Um(j-1) + Um(j))
        bltre =   y(j) - (-y(j-1) + y(j)) * dudy
        exit
      endif
    enddo

    if (bltre > blt) then
      flag_re = flag_re + 1
    endif

    if (flag_re >= 1 .and. step >= start_rescale) then
      ufin(:,:)  = 0.d0
      vfin(:,:)  = 0.d0
      wfin(:,:)  = 0.d0
      Tfin(:,:)  = 0.d0
      pfin(:,:)  = 0.d0
      ufout(:,:) = 0.d0
      vfout(:,:) = 0.d0
      wfout(:,:) = 0.d0
      Tfout(:,:) = 0.d0
      pfout(:,:) = 0.d0
      do k = 1, nz
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset = 5 * (j-1)
          rhore = Qre(k_offset+j_offset+1)
          over_rhore = 1.d0 / rhore
          ure   = Qre(k_offset+j_offset+2) * over_rhore
          vre   = Qre(k_offset+j_offset+3) * over_rhore
          wre   = Qre(k_offset+j_offset+4) * over_rhore
          pre   = gamma_1 * (Qre(k_offset+j_offset+5) - 0.5d0 * rhore * (ure**2 + vre**2 + wre**2)) 
          Tre   = pre / (rhore * R)
          ufre(j,k) = ure - Um(j)
          vfre(j,k) = vre - Vm(j)
          wfre(j,k) = wre - Wm(j)
          Tfre(j,k) = Tre - Tm(j)
          pfre(j,k) = pre - pm(j)
      enddo;enddo

      ! friction velocity
      mu    = mu0_T0_S / (Tm(1) + 111.d0) * (Tm(1) * over_T0)**1.5
      nu    = mu / rhom(1)
      taure = mu * abs(-Um(1) + Um(2)) / (-y(1) + y(2))
      utre  = sqrt(taure / rhom(1))
      beta  = (bltre / blt)**0.1
      utin  = beta * utre

      do j = 1, ny
        ypin(j) = y(j) * utin / nu
        ypre(j) = y(j) * utre / nu
        etin(j) = y(j) / blt
        etre(j) = y(j) / bltre
        weight(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / (0.6d0 * etin(j) + 0.2d0)) / tanh(4.d0)))
      enddo

      jj_y(:) = -1
      do j = 1, ny
        do jj = 2, ny
          if (ypre(jj) > ypin(j)) then
            ady = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            ! mean
            Umin(j) = beta * ((1.d0 - ady) * Um(jj-1) + ady * Um(jj))!(Um(jj-1) + ady * (-Um(jj-1) + Um(jj)))
            Vmin(j) =         (1.d0 - ady) * Vm(jj-1) + ady * Vm(jj) !Vm(jj-1) + ady * (-Vm(jj-1) + Vm(jj))
            Tmin(j) =         (1.d0 - ady) * Tm(jj-1) + ady * Tm(jj) !Tm(jj-1) + ady * (-Tm(jj-1) + Tm(jj))
            pmin(j) =         (1.d0 - ady) * pm(jj-1) + ady * pm(jj) !pm(jj-1) + ady * (-pm(jj-1) + pm(jj))
            jj_y(j) = jj
            exit
          endif
      enddo;enddo

      jj_e(:) = -1
      do j = 1, ny
        do jj = 2, ny
          if (etre(jj) > etin(j)) then
            ade = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            ! mean
            Umout(j) = beta * ((1.d0 - ade) * Um(jj-1) + ade * Um(jj)) + (1.d0 - beta) * u0
            Vmout(j) =         (1.d0 - ade) * Vm(jj-1) + ade * Vm(jj) !Vm(jj-1) + ade * (-Vm(jj-1) + Vm(jj))
            Tmout(j) =         (1.d0 - ade) * Tm(jj-1) + ade * Tm(jj) !Tm(jj-1) + ade * (-Tm(jj-1) + Tm(jj))
            pmout(j) =         (1.d0 - ade) * pm(jj-1) + ade * pm(jj) !pm(jj-1) + ade * (-pm(jj-1) + pm(jj))
            jj_e(j)  = jj
            exit
          endif
      enddo;enddo
      
      do k = 1, nz
        do j = 1, ny
          jj = jj_y(j)
          if (jj > 0) then
            ady = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
            ufin(j,k) = beta * ((1.d0 - ady) * ufre(jj-1,k) + ady * ufre(jj,k))!(ufre(jj-1,k) + ady * (-ufre(jj-1,k) + ufre(jj,k)))
            vfin(j,k) = beta * ((1.d0 - ady) * vfre(jj-1,k) + ady * vfre(jj,k))!(vfre(jj-1,k) + ady * (-vfre(jj-1,k) + vfre(jj,k)))
            wfin(j,k) = beta * ((1.d0 - ady) * wfre(jj-1,k) + ady * wfre(jj,k))!(wfre(jj-1,k) + ady * (-wfre(jj-1,k) + wfre(jj,k)))
            Tfin(j,k) =         (1.d0 - ady) * Tfre(jj-1,k) + ady * Tfre(jj,k) !Tfre(jj-1,k) + ady * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfin(j,k) =         (1.d0 - ady) * pfre(jj-1,k) + ady * pfre(jj,k) !pfre(jj-1,k) + ady * (-pfre(jj-1,k) + pfre(jj,k))
          endif
          jj = jj_e(j)
          if (jj > 0) then
            ade = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
            ufout(j,k) = beta * ((1.d0 - ade) * ufre(jj-1,k) + ade * ufre(jj,k))!(ufre(jj-1,k) + ade * (-ufre(jj-1,k) + ufre(jj,k)))
            vfout(j,k) = beta * ((1.d0 - ade) * vfre(jj-1,k) + ade * vfre(jj,k))!(vfre(jj-1,k) + ade * (-vfre(jj-1,k) + vfre(jj,k)))
            wfout(j,k) = beta * ((1.d0 - ade) * wfre(jj-1,k) + ade * wfre(jj,k))!(wfre(jj-1,k) + ade * (-wfre(jj-1,k) + wfre(jj,k)))
            Tfout(j,k) =         (1.d0 - ade) * Tfre(jj-1,k) + ade * Tfre(jj,k) !Tfre(jj-1,k) + ade * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfout(j,k) =         (1.d0 - ade) * pfre(jj-1,k) + ade * pfre(jj,k) !pfre(jj-1,k) + ade * (-pfre(jj-1,k) + pfre(jj,k))
          endif
      enddo;enddo
  
      ! re-introducing
      do k = 1, nz
        kh = mod(k+nz/2,nz) + 1
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset = 5 * (j-1)
          weight_tmp   = weight(j)
          Jacobian_tmp = 1.d0 / Jacobian(nre2,j)
          uin = (Umin(j) + ufin(j,kh)) * (1.d0 - weight_tmp) + (Umout(j) + ufout(j,kh)) * weight_tmp
          vin = (Vmin(j) + vfin(j,kh)) * (1.d0 - weight_tmp) + (Vmout(j) + vfout(j,kh)) * weight_tmp
          win =            wfin(j,kh)  * (1.d0 - weight_tmp) +             wfout(j,kh)  * weight_tmp
          Tin = (Tmin(j) + Tfin(j,kh)) * (1.d0 - weight_tmp) + (Tmout(j) + Tfout(j,kh)) * weight_tmp
          pin = (pmin(j) + pfin(j,kh)) * (1.d0 - weight_tmp) + (pmout(j) + pfout(j,kh)) * weight_tmp
          rhoin = pin / (R * Tin)
          Qre(k_offset+j_offset+1) = rhoin * Jacobian_tmp
          Qre(k_offset+j_offset+2) = rhoin * uin * Jacobian_tmp
          Qre(k_offset+j_offset+3) = rhoin * vin * Jacobian_tmp
          Qre(k_offset+j_offset+4) = rhoin * win * Jacobian_tmp
          Qre(k_offset+j_offset+5) = (pin * over_gamma_1 + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)) * Jacobian_tmp
      enddo;enddo
    else
      ! cyclic boundary condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = 1, nz
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset = 5 * (j-1)
          do l = 1, 5
            Qre(k_offset+j_offset+l) = Qre(k_offset+j_offset+l) / Jacobian(nre2,j)
      enddo;enddo;enddo
    endif
  end subroutine set_rescale
end module calc_rescale

