module set_init_dhit
  use cufft
  use mod_globals, only : gamma, RHO0, p0, up0, pi, A0, k0, mu0, c0, &
                          Re_lambda_target, Mat_target
  use mod_constant, only : id_accuracy
  implicit none
contains
  subroutine init_spectral_velocity(nx, ny, nz, Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(out) :: Q(nx,ny,nz,5)
    !===========================================================
    ! parameters
    !===========================================================
    integer :: offset
    integer :: Nf
    integer :: i,j,k
    integer :: icomp

    integer :: kxi,kyi,kzi
    integer :: kxc,kyc,kzc

    integer(8) :: ikx,iky,ikz

    integer :: seed_size
    integer, allocatable :: seed_arr(:)

    integer :: plan
    integer :: istat

    real(8) :: dk, kmag2, kmag, kx, ky, kz, k12

    real(8) :: Ek
    real(8) :: amp

    real(8) :: du1dx1_sq, dx, u1p, u1m
    real(8) :: lambda_meas, Re_lambda_meas, Mat_meas
    integer :: ip, im

    real(8) :: phi1
    real(8) :: phi2
    real(8) :: phi3

    real(8) :: r(3)

    real(8) :: norm

    real(8) :: urms_sq
    real(8) :: uscale

    real(8) :: u,v,w

    complex(8) :: a,b

    complex(8), allocatable :: uk(:,:,:,:)
    real(8), allocatable :: vel_r(:,:,:,:)

    complex(8), device, allocatable :: uk_d(:,:,:)
    real(8), device, allocatable :: vk_d(:,:,:)

    !===========================================================
    ! setup
    !===========================================================
    if     (kind(id_accuracy) == 2) then
      offset = 1
    elseif (kind(id_accuracy) == 4) then
      offset = 2
    elseif (kind(id_accuracy) == 8) then
      offset = 3
    endif
    Nf = nx - 2 * offset
    
    allocate(uk(Nf/2+1,Nf,Nf,3))
    allocate(vel_r(Nf,Nf,Nf,3))
    allocate(uk_d(Nf/2+1,Nf,Nf))
    allocate(vk_d(Nf,Nf,Nf))
    
    uk = cmplx(0.d0,0.d0,8)
    
    !===========================================================
    ! random seed
    !===========================================================
    call random_seed(size=seed_size)
    allocate(seed_arr(seed_size))
    seed_arr = 42
    call random_seed(put=seed_arr)
    !===========================================================
    ! generate Fourier modes (Rogallo 1981 / Johnsen et al. 2010)
    !===========================================================
    do kzi = 1, Nf
      ikz = kzi - 1
      if (ikz > Nf/2) ikz = ikz - Nf
      do kyi = 1, Nf
        iky = kyi - 1
        if (iky > Nf/2) iky = iky - Nf
        do kxi = 1, Nf/2 + 1
          ikx = kxi - 1
          
          kmag2 = dble(ikx*ikx + iky*iky + ikz*ikz)
          if (kmag2 < 0.5d0) cycle
          kmag = sqrt(kmag2)
          
          kx = dble(ikx)/kmag
          ky = dble(iky)/kmag
          kz = dble(ikz)/kmag
          
          k12 = sqrt(dble(ikx*ikx + iky*iky))

          ! ref.tex "Compressible homogeneous turbulence": E(k) = A0*k^4*exp(-2*k^2/k0^2)
          Ek  = A0 * kmag**4 * exp(-2.d0*kmag**2/k0**2)

          amp = sqrt(Ek / (2.d0 * pi * kmag**2))
          amp = amp * dble(Nf)**3
          
          call random_number(r)
          phi1 = 2.d0*pi*r(1)
          phi2 = 2.d0*pi*r(2)
          phi3 = 2.d0*pi*r(3)
          a = amp * exp(cmplx(0.d0,phi1,8)) * cos(phi3)
          b = amp * exp(cmplx(0.d0,phi2,8)) * sin(phi3)
          
          if (k12 < 0.5d0) then
            uk(kxi,kyi,kzi,1) = a
            uk(kxi,kyi,kzi,2) = cmplx(0.d0,0.d0,8)
            uk(kxi,kyi,kzi,3) = cmplx(0.d0,0.d0,8)
          else
            uk(kxi,kyi,kzi,1) = (dble(iky)/k12)*a + (dble(ikx)/k12)*kz*b
            uk(kxi,kyi,kzi,2) = -(dble(ikx)/k12)*a + (dble(iky)/k12)*kz*b
            uk(kxi,kyi,kzi,3) = -(k12/kmag)*b
          endif
    enddo;enddo;enddo
    !===========================================================
    ! Hermitian symmetry WITHIN kx=0 and kx=Nyquist planes ONLY
    !===========================================================
    do kxi = 1, Nf/2+1, Nf/2
      do kzi = 1, Nf
        kzc = mod(Nf-kzi+1, Nf) + 1
        do kyi = 2, Nf/2
          kyc = mod(Nf-kyi+1, Nf) + 1
          do icomp = 1,3
            uk(kxi,kyc,kzc,icomp) = conjg(uk(kxi,kyi,kzi,icomp))
      enddo;enddo;enddo
      do kyi = 1, Nf/2+1, Nf/2
        do kzi = 2, Nf/2
          kzc = mod(Nf-kzi+1, Nf) + 1
          do icomp = 1,3
            uk(kxi,kyi,kzc,icomp) = conjg(uk(kxi,kyi,kzi,icomp))
      enddo;enddo;enddo
    enddo
    !===========================================================
    ! self-conjugate modes (DC and Nyquist points)
    !===========================================================
    do icomp = 1,3
      do kzi = 1, Nf, Nf/2
        do kyi = 1, Nf, Nf/2
          do kxi = 1, Nf/2+1, Nf/2
            uk(kxi,kyi,kzi,icomp) = dcmplx(real(uk(kxi,kyi,kzi,icomp)),0.d0)
    enddo;enddo;enddo;enddo
    !===========================================================
    ! inverse FFT (Complex to Real : Z2D)
    !===========================================================
    istat = cufftPlan3d(plan, Nf, Nf, Nf, CUFFT_Z2D)
    do icomp = 1,3
      uk_d = uk(:,:,:,icomp)
      istat = cufftExecZ2D(plan, uk_d, vk_d)
      vel_r(:,:,:,icomp) = vk_d
    enddo
    istat = cufftDestroy(plan)
    !===========================================================
    ! normalization
    !===========================================================
    norm = 1.d0 / dble(Nf)**3
    urms_sq = 0.d0
    do icomp = 1,3
      do k = 1, Nf
        do j = 1, Nf
          do i = 1, Nf
            urms_sq = urms_sq + (vel_r(i,j,k,icomp)*norm)**2
    enddo;enddo;enddo;enddo
    urms_sq = urms_sq / dble(3*Nf**3)
    uscale = up0 / sqrt(urms_sq)
    !===========================================================
    ! self-consistency check: measured Re_lambda, Ma_t from the
    ! actual synthesized field (independent of the Python
    ! post-processing check run on the first VTK snapshot)
    !===========================================================
    dx = 2.d0*pi / dble(Nf)   ! solver's actual periodic grid spacing (Lx/Nf, Nf = nx - 6)
    du1dx1_sq = 0.d0
    do k = 1, Nf
      do j = 1, Nf
        do i = 1, Nf
          ip = mod(i, Nf) + 1
          im = mod(i-2+Nf, Nf) + 1
          u1p = vel_r(ip,j,k,1) * norm * uscale
          u1m = vel_r(im,j,k,1) * norm * uscale
          du1dx1_sq = du1dx1_sq + ((u1p - u1m) / (2.d0*dx))**2
    enddo;enddo;enddo
    du1dx1_sq = du1dx1_sq / dble(Nf**3)

    lambda_meas    = up0 / sqrt(du1dx1_sq)
    Re_lambda_meas = RHO0 * up0 * lambda_meas / mu0
    Mat_meas       = sqrt(3.d0) * up0 / c0
    !===========================================================
    ! diagnostics
    !===========================================================
    print *, '[DHIT init]'
    print *, 'Nf                 = ', Nf
    print *, 'u_rms(raw)         = ', sqrt(urms_sq)
    print *, 'u_rms(target)      = ', up0
    print *, 'u_scale            = ', uscale
    print *, 'Re_lambda(measured)= ', Re_lambda_meas
    print *, 'Re_lambda(target)  = ', Re_lambda_target
    print *, 'Ma_t(measured)     = ', Mat_meas
    print *, 'Ma_t(target)       = ', Mat_target
    !===========================================================
    ! pack conservative variables
    !===========================================================
    do k = 1, Nf
      do j = 1, Nf
        do i = 1, Nf
          u = vel_r(i,j,k,1)*norm*uscale
          v = vel_r(i,j,k,2)*norm*uscale
          w = vel_r(i,j,k,3)*norm*uscale
          
          Q(i+offset,j+offset,k+offset,1) = RHO0
          Q(i+offset,j+offset,k+offset,2) = RHO0 * u
          Q(i+offset,j+offset,k+offset,3) = RHO0 * v
          Q(i+offset,j+offset,k+offset,4) = RHO0 * w
          Q(i+offset,j+offset,k+offset,5) = p0 / (gamma-1.d0) + 0.5d0 * RHO0 * (u*u + v*v + w*w)
    enddo;enddo;enddo
    !===========================================================
    ! cleanup
    !===========================================================
    deallocate(uk)
    deallocate(vel_r)
    deallocate(uk_d)
    deallocate(vk_d)
    deallocate(seed_arr)
  end subroutine init_spectral_velocity
end module set_init_dhit
