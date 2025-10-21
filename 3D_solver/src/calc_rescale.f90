module calc_rescale
  use cudafor
  use mpi
  use mod_globals, only : nre1, nre2, rerank, nt, dt, gamma , R, Pr, u0, rho0, p0, M0, blt, start_rescale
  use mod_constant, only : Cp, gamma_1, over_gamma_1, mu0_T0_S, over_T0
contains
  subroutine calc_mean(step, flag_re, nx, ny, nz, Jacobian, QJ, Qm)
    integer, intent(in)            :: step, flag_re, nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(inout), device :: Qm(5,ny)
    real(8) Q1, Q2, Q3, Q4, Q5, rhoinv, Jacobian_tmp, volinv, step1, step2
    integer i, k
    volinv = 1.d0 / dble((nre2 - nre1 + 1) * (nz - 6))
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
      Qm(1,j) = Q1 * volinv
      Qm(2,j) = Q2 * volinv
      Qm(3,j) = Q3 * volinv
      Qm(4,j) = Q4 * volinv
      Qm(5,j) = Q5 * volinv
    enddo
  end subroutine calc_mean


  subroutine flatten(nx, ny, nz, Jacobian, QJ, Qre)
    integer, intent(in)          :: nx, ny, nz
    real(8), intent(in), device  :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(out), device :: Qre(ny*(nz-6)*5)
    integer j, k, l, j_offset, k_offset
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz-6
      k_offset = ny * 5 * (k-1)
      do j = 1, ny
        j_offset = 5 * (j-1)
        do l = 1, 5
          Qre(k_offset+j_offset+l) = QJ(l,nre2,j,k+3) * Jacobian(nre2,j)
    enddo;enddo;enddo
  end subroutine flatten


  subroutine Q_over_J(nx, ny, nz, Jacobian, Qre)
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in), device    :: Jacobian(nx,ny)
    real(8), intent(inout), device :: Qre(ny*(nz-6)*5)
    integer i, j, k, l, j_offset, k_offset
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz-6
      k_offset = ny * 5 * (k-1)
      do j = 1, ny
        j_offset = 5 * (j-1)
        do l = 1, 5
          i = k_offset + j_offset + l
          Qre(i) = Qre(i) * Jacobian(nre2,j)
    enddo;enddo;enddo
  end subroutine Q_over_J


  subroutine step_rescale(num, myrank, step, nx, ny, nz, flag_re, ireq, ireq2, y, Jacobian, QJ, Qm, Qre)
    integer, intent(in)            :: num, myrank, step, nx, ny, nz
    integer, intent(inout)         :: flag_re, ireq, ireq2(2)
    real(8), intent(in)            :: y(ny)
    real(8), intent(in), device    :: Jacobian(nx,ny), QJ(5,nx,ny,nz)
    real(8), intent(inout), device :: Qm(5,ny), Qre(ny*(nz-6)*5)
    real(8) Qm_cpu(5,ny), dudy, bltre
    integer ierr, j
    if (myrank == rerank) then
      if (num == 1) then
        call calc_mean(step, flag_re, nx, ny, nz, Jacobian, QJ, Qm)
        call flatten(nx, ny, nz, Jacobian, QJ, Qre)
        Qm_cpu = Qm
        bltre  = 0.d0
        do j = 2, ny
          if (Qm_cpu(2,j) >= 0.99d0 * u0) then
            dudy  = (Qm_cpu(2,j) - 0.99d0 * u0) / (-Qm_cpu(2,j-1) + Qm_cpu(2,j))
            bltre = y(j) - (-y(j-1) + y(j)) * dudy
            exit
          endif
        enddo
        print *, "Boundary layer thickness = ", bltre
      endif
      if (bltre >= blt) then
        block
          real(8) beta
          real(8), dimension(ny), device :: ady, ade, ypre, ypin, etre, etin, weight
          integer, dimension(ny), device :: jj_y, jj_e
          call set_rescale_cpu(nx, ny, nz, y, Qm_cpu, beta, ady, ade, ypre, ypin, etre, etin, weight, jj_y, jj_e)
          !call set_rescale_gpu(nx, ny, nz, beta, ady, ade, ypre, ypin, etre, etin, weight, Qm, jj_y, jj_e, Qre)
        end block
      endif
      call Q_over_J(nx, ny, nz, Jacobian, Qre)
    endif
    if (rerank /= 0 .and. myrank == rerank) then
      call MPI_ISEND(Qre, 5*ny*(nz-6), MPI_REAL8, 0,      0, MPI_COMM_WORLD, ireq, ierr)
    elseif (rerank /= 0 .and. myrank == 0) then
      call MPI_IRECV(Qre, 5*ny*(nz-6), MPI_REAL8, rerank, 0, MPI_COMM_WORLD, ireq, ierr)
    endif
  end subroutine step_rescale


  subroutine set_rescale_cpu(nx, ny, nz, y, Qm_cpu, beta, ady_gpu, ade_gpu, ypre_gpu, ypin_gpu, etre_gpu, etin_gpu, &
                             weight_gpu, jj_y_gpu, jj_e_gpu)
    integer, intent(in), value                  :: nx, ny, nz
    real(8), intent(in)                         :: y(ny), Qm_cpu(5,ny)
    real(8), intent(out)                        :: beta
    real(8), intent(out), dimension(ny), device :: ady_gpu, ade_gpu, ypre_gpu, ypin_gpu, etre_gpu, etin_gpu, weight_gpu
    integer, intent(out), dimension(ny), device :: jj_y_gpu, jj_e_gpu
    real(8) rhow, Tw, mu, nu, taure, utre, utin 
    real(8) over_tanh4, over_nu, over_bltre
    real(8), parameter     :: over_blt = 1.d0 / blt
    real(8), dimension(ny) :: ypin, ypre, etin, etre, weight, ady_cpu, ade_cpu, weight_cpu
    integer, dimension(ny) :: jj_y_cpu, jj_e_cpu
    ! friction velocity
    rhow  = Qm_cpu(1,1)
    Tw    = Qm_cpu(5,1) / (R * rhow)
    mu    = mu0_T0_S / (Tw + 111.d0) * (Tw * over_T0)**1.5
    taure = mu * abs(-Qm_cpu(2,1) + Qm_cpu(2,2)) / (-y(1) + y(2))
    utre  = sqrt(taure / rhow)
    beta  = (bltre * over_blt)**0.1
    utin  = beta * utre
    ! calc reciprocal
    over_nu    = rhow / mu
    over_bltre = 1.d0 / bltre
    over_tanh4 = 1.d0 / tanh(4.d0)
    ! calc inner and outer coordinate
    do j = 1, ny
      ypin(j) = y(j) * utin * over_nu
      ypre(j) = y(j) * utre * over_nu
      etin(j) = y(j) * over_blt
      etre(j) = y(j) * over_bltre
      weight_cpu(j) = min(1.d0, 0.5d0 * (1.d0 + tanh(4.d0 * (etin(j) - 0.2d0) / (0.6d0 * etin(j) + 0.2d0)) * over_tanh4))
    enddo
    jj_y_cpu(:) = -1
    do j = 1, ny
      do jj = 2, ny
        if (ypre(jj) > ypin(j)) then
          ady_cpu(j)  = (-ypre(jj-1) + ypin(j)) / (-ypre(jj-1) + ypre(jj))
          jj_y_cpu(j) = jj
          exit
        endif
    enddo;enddo
    jj_e_cpu(:) = -1
    do j = 1, ny
      do jj = 2, ny
        if (etre(jj) > etin(j)) then
          ade_cpu(j)  = (-etre(jj-1) + etin(j)) / (-etre(jj-1) + etre(jj))
          jj_e_cpu(j) = jj
          exit
        endif
    enddo;enddo
    ady_gpu    = ady_cpu
    ade_gpu    = ade_cpu
    ypre_gpu   = ypre
    ypin_gpu   = ypin
    etre_gpu   = etre
    etin_gpu   = etin
    jj_y_gpu   = jj_y_cpu
    jj_e_gpu   = jj_e_cpu
    weight_gpu = weight_cpu
  end subroutine set_rescale_cpu


  subroutine set_rescale_gpu(nx, ny, nz, beta, ady, ade, ypre, ypin, etre, etin, weight, Qm, jj_y, jj_e, Qre)
    integer, intent(in)                        :: nx, ny, nz ! nz-6
    real(8), intent(in)                        :: beta
    real(8), intent(in), dimension(ny), device :: ady, ade, ypre, ypin, etre, etin, weight
    real(8), intent(in), device                :: Qm(5,ny)
    integer, intent(in), dimension(ny), device :: jj_y, jj_e
    real(8), intent(inout), device             :: Qre(ny*nz*5) ! Q / J
    ! mean properties at rescaling plane
    real(8), dimension(ny), device    :: Tm, pm
    ! fluctuating properties at rescaling plane
    real(8), dimension(ny,nz), device :: ufre,  vfre,  wfre,  Tfre,  pfre
    ! fluctuating properties at both inner and outer region
    real(8), dimension(ny,nz), device :: ufin,  vfin,  wfin,  Tfin,  pfin
    real(8), dimension(ny,nz), device :: ufout, vfout, wfout, Tfout, pfout
    ! mean properties at both inner and outer region
    real(8), dimension(ny), device    :: Umin,  Vmin,  pmin,  Tmin,  Umout, Vmout, pmout, Tmout
    integer j, jj, k, kh, l, j_offset, k_offset
    block
      real(8) u_tmp, v_tmp, p_tmp, T_tmp
      !$cuf kernel do <<<*,*>>>
      do j = 1, ny
        u_tmp    = Qm(2,j)
        v_tmp    = Qm(3,j)
        p_tmp    = Qm(5,j)
        T_tmp    = p_tmp / (R * Qm(1,j))
        Umin(j)  = u_tmp
        Umout(j) = u_tmp
        Vmin(j)  = v_tmp
        Vmout(j) = v_tmp
        pmin(j)  = p_tmp
        pmout(j) = p_tmp
        Tm(j)    = T_tmp
        Tmin(j)  = T_tmp
        Tmout(j) = T_tmp
      enddo
    end block
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        ufin(j,k)  = 0.d0
        vfin(j,k)  = 0.d0
        wfin(j,k)  = 0.d0
        Tfin(j,k)  = 0.d0
        pfin(j,k)  = 0.d0
        ufout(j,k) = 0.d0
        vfout(j,k) = 0.d0
        wfout(j,k) = 0.d0
        Tfout(j,k) = 0.d0
        pfout(j,k) = 0.d0
    enddo;enddo
    block
      ! properties at rescaling plane
      real(8) ure, vre, wre, rhore, Tre, pre, over_rhore
      !$cuf kernel do <<<*,*>>>
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
          ufre(j,k) = ure - Qm(2,j)
          vfre(j,k) = vre - Qm(3,j)
          wfre(j,k) = wre - Qm(4,j)
          Tfre(j,k) = Tre -   Tm(j)
          pfre(j,k) = pre - Qm(5,j)
      enddo;enddo
    end block
    !$cuf kernel do <<<*,*>>>
    do j = 1, ny
      do jj = 2, ny
        if (ypre(jj) > ypin(j)) then
          ! mean
          Umin(j) = beta * ((1.d0 - ady(j)) * Qm(2,jj-1) + ady(j) * Qm(2,jj))!(Um(jj-1) + ady * (-Um(jj-1) + Um(jj)))
          Vmin(j) =         (1.d0 - ady(j)) * Qm(3,jj-1) + ady(j) * Qm(3,jj) !Vm(jj-1) + ady * (-Vm(jj-1) + Vm(jj))
          Tmin(j) =         (1.d0 - ady(j)) *   Tm(jj-1) + ady(j) *   Tm(jj) !Tm(jj-1) + ady * (-Tm(jj-1) + Tm(jj))
          pmin(j) =         (1.d0 - ady(j)) * Qm(5,jj-1) + ady(j) * Qm(5,jj) !pm(jj-1) + ady * (-pm(jj-1) + pm(jj))
          exit
        endif
    enddo;enddo
    !$cuf kernel do <<<*,*>>>
    do j = 1, ny
      do jj = 2, ny
        if (etre(jj) > etin(j)) then
          ! mean
          Umout(j) = beta * ((1.d0 - ade(j)) * Qm(2,jj-1) + ade(j) * Qm(2,jj)) + (1.d0 - beta) * u0
          Vmout(j) =         (1.d0 - ade(j)) * Qm(3,jj-1) + ade(j) * Qm(3,jj) !Vm(jj-1) + ade * (-Vm(jj-1) + Vm(jj))
          Tmout(j) =         (1.d0 - ade(j)) *   Tm(jj-1) + ade(j) *   Tm(jj) !Tm(jj-1) + ade * (-Tm(jj-1) + Tm(jj))
          pmout(j) =         (1.d0 - ade(j)) * Qm(5,jj-1) + ade(j) * Qm(5,jj) !pm(jj-1) + ade * (-pm(jj-1) + pm(jj))
          exit
        endif
    enddo;enddo
    !$cuf kernel do <<<*,*>>>
    do k = 1, nz
      do j = 1, ny
        jj = jj_y(j)
        if (jj > 0) then
          ufin(j,k) = beta * ((1.d0 - ady(j)) * ufre(jj-1,k) + ady(j) * ufre(jj,k))!(ufre(jj-1,k) + ady * (-ufre(jj-1,k) + ufre(jj,k)))
          vfin(j,k) = beta * ((1.d0 - ady(j)) * vfre(jj-1,k) + ady(j) * vfre(jj,k))!(vfre(jj-1,k) + ady * (-vfre(jj-1,k) + vfre(jj,k)))
          wfin(j,k) = beta * ((1.d0 - ady(j)) * wfre(jj-1,k) + ady(j) * wfre(jj,k))!(wfre(jj-1,k) + ady * (-wfre(jj-1,k) + wfre(jj,k)))
          Tfin(j,k) =         (1.d0 - ady(j)) * Tfre(jj-1,k) + ady(j) * Tfre(jj,k) !Tfre(jj-1,k) + ady * (-Tfre(jj-1,k) + Tfre(jj,k))
          pfin(j,k) =         (1.d0 - ady(j)) * pfre(jj-1,k) + ady(j) * pfre(jj,k) !pfre(jj-1,k) + ady * (-pfre(jj-1,k) + pfre(jj,k))
        endif
        jj = jj_e(j)
        if (jj > 0) then
            ufout(j,k) = beta * ((1.d0 - ade(j)) * ufre(jj-1,k) + ade(j) * ufre(jj,k))!(ufre(jj-1,k) + ade * (-ufre(jj-1,k) + ufre(jj,k)))
            vfout(j,k) = beta * ((1.d0 - ade(j)) * vfre(jj-1,k) + ade(j) * vfre(jj,k))!(vfre(jj-1,k) + ade * (-vfre(jj-1,k) + vfre(jj,k)))
            wfout(j,k) = beta * ((1.d0 - ade(j)) * wfre(jj-1,k) + ade(j) * wfre(jj,k))!(wfre(jj-1,k) + ade * (-wfre(jj-1,k) + wfre(jj,k)))
            Tfout(j,k) =         (1.d0 - ade(j)) * Tfre(jj-1,k) + ade(j) * Tfre(jj,k) !Tfre(jj-1,k) + ade * (-Tfre(jj-1,k) + Tfre(jj,k))
            pfout(j,k) =         (1.d0 - ade(j)) * pfre(jj-1,k) + ade(j) * pfre(jj,k) !pfre(jj-1,k) + ade * (-pfre(jj-1,k) + pfre(jj,k))
        endif
    enddo;enddo
    block
      ! rescaled properties at inlet
      real(8) uin, vin, win, rhoin, Tin, pin, weight_tmp
      !$cuf kernel do <<<*,*>>>
      do k = 1, nz
        kh = mod(k+nz/2,nz) + 1
        k_offset = ny * 5 * (k-1)
        do j = 1, ny
          j_offset = 5 * (j-1)
          weight_tmp   = weight(j)
          uin = (Umin(j) + ufin(j,kh)) * (1.d0 - weight_tmp) + (Umout(j) + ufout(j,kh)) * weight_tmp
          vin = (Vmin(j) + vfin(j,kh)) * (1.d0 - weight_tmp) + (Vmout(j) + vfout(j,kh)) * weight_tmp
          win =            wfin(j,kh)  * (1.d0 - weight_tmp) +             wfout(j,kh)  * weight_tmp
          Tin = (Tmin(j) + Tfin(j,kh)) * (1.d0 - weight_tmp) + (Tmout(j) + Tfout(j,kh)) * weight_tmp
          pin = (pmin(j) + pfin(j,kh)) * (1.d0 - weight_tmp) + (pmout(j) + pfout(j,kh)) * weight_tmp
          rhoin = pin / (R * Tin)
          Qre(k_offset+j_offset+1) = rhoin
          Qre(k_offset+j_offset+2) = rhoin * uin
          Qre(k_offset+j_offset+3) = rhoin * vin
          Qre(k_offset+j_offset+4) = rhoin * win
          Qre(k_offset+j_offset+5) = pin * over_gamma_1 + 0.5d0 * rhoin * (uin**2 + vin**2 + win**2)
      enddo;enddo
    end block
  end subroutine set_rescale_gpu
end module calc_rescale

