module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, gamma, threshold, threadsE, threadsF, threadsG
  use calc_keep
  use calc_slau
  use calc_roe
  use calc_hybrid
  use calc_muscl
  implicit none
  interface flux6
    module procedure flux_KEEP6, flux_Roe6, flux_SLAU6, flux_Weighted6, flux_Threshold6
  end interface flux6

  interface flux4
    module procedure flux_KEEP4, flux_Roe4, flux_SLAU4, flux_Weighted4, flux_Threshold4
  end interface flux4

  interface flux2
    module procedure flux_KEEP2, flux_Roe2, flux_SLAU2, flux_Weighted2, flux_Threshold2
  end interface flux2

  interface calc_E
    module procedure calc_E2, calc_E4, calc_E6
  end interface calc_E

  interface calc_F
    module procedure calc_F2, calc_F4, calc_F6
  end interface calc_F

  interface calc_G
    module procedure calc_G2, calc_G4, calc_G6
  end interface calc_G
contains
  attributes(device) subroutine flux_KEEP6(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(6), u(6), v(6), w(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(5)
    F = KEEP6(rho, u, v, w, uu, p, T, Normal)
  end subroutine flux_KEEP6


  attributes(device) subroutine flux_KEEP4(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(4), u(4), v(4), w(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(5)
    F = KEEP4(rho, u, v, w, uu, p, T, Normal)
  end subroutine flux_KEEP4


  attributes(device) subroutine flux_KEEP2(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(2), u(2), v(2), w(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(5)
    F = KEEP2(rho, u, v, w, uu, p, T, Normal)
  end subroutine flux_KEEP2


  attributes(device) subroutine flux_SLAU6(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value  :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), w(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(8) wiggle, rho2(2), p2(2), V2(2,3)
    wiggle = wiggle_detector(p(2:5))
    call calc_6points(sensor, rho, u, v, w, p, rho2, p2, V2)
    F = SLAU(id_slau, id, rho2, p2, V2, Normal, wiggle)
  end subroutine flux_SLAU6


  attributes(device) subroutine flux_SLAU4(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value  :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), w(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(8) wiggle, rho2(2), p2(2), V2(2,3)
    wiggle = wiggle_detector(p)
    call calc_4points(sensor, rho, u, v, w, p, rho2, p2, V2)
    F = SLAU(id_slau, id, rho2, p2, V2, Normal, wiggle)
  end subroutine flux_SLAU4


  attributes(device) subroutine flux_SLAU2(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value  :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), w(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(8) V2(2,3)
    V2(:,1) = u
    V2(:,2) = v
    V2(:,3) = w
    F = SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
  end subroutine flux_SLAU2


  attributes(device) subroutine flux_Roe6(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(6), u(6), v(6), w(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(5)
    real(8) rho2(2), p2(2), V2(2,3), Vtmp(2), Ftmp
    call calc_6points(sensor, rho, u, v, w, p, rho2, p2, V2)
    select case(id)
    case(1)
      continue
    case(2)
      Vtmp(:) = V2(:,1)
      V2(:,1) = V2(:,2) ! v
      V2(:,2) = V2(:,3) ! w
      V2(:,3) = Vtmp    ! u
    case(3)
      Vtmp(:) = V2(:,1) 
      V2(:,1) = V2(:,3) ! w
      V2(:,3) = V2(:,2) ! v
      V2(:,2) = Vtmp    ! u
    end select
    F = Roe(rho2, V2, p2, Normal_x)
    select case(id)
    case(1)
      continue
    case(2)
      Ftmp = F(2)
      F(2) = F(4) ! w
      F(4) = F(3) ! v
      F(3) = Ftmp ! u
    case(3)
      Ftmp = F(2)
      F(2) = F(3)! v
      F(3) = F(4)! w
      F(4) = Ftmp! u
    end select
  end subroutine flux_Roe6


  attributes(device) subroutine flux_Roe4(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(4), u(4), v(4), w(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(5)
    real(8) rho2(2), p2(2), V2(2,3), Vtmp(2), Ftmp
    call calc_4points(sensor, rho, u, v, w, p, rho2, p2, V2)
    select case(id)
    case(1)
      continue
    case(2)
      Vtmp(:) = V2(:,1)
      V2(:,1) = V2(:,2) ! v
      V2(:,2) = V2(:,3) ! w
      V2(:,3) = Vtmp    ! u
    case(3)
      Vtmp(:) = V2(:,1) 
      V2(:,1) = V2(:,3) ! w
      V2(:,3) = V2(:,2) ! v
      V2(:,2) = Vtmp    ! u
    end select
    F = Roe(rho2, V2, p2, Normal_x)
    select case(id)
    case(1)
      continue
    case(2)
      Ftmp = F(2)
      F(2) = F(4) ! w
      F(4) = F(3) ! v
      F(3) = Ftmp ! u
    case(3)
      Ftmp = F(2)
      F(2) = F(3)! v
      F(3) = F(4)! w
      F(4) = Ftmp! u
    end select
  end subroutine flux_Roe4


  attributes(device) subroutine flux_Roe2(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(2), u(2), v(2), w(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8), intent(out), contiguous   :: F(5)
    real(8) V2(2,3), Ftmp
    select case(id)
    case(1)
      V2(:,1) = u
      V2(:,2) = v
      V2(:,3) = w
    case(2)
      V2(:,1) = v
      V2(:,2) = w
      V2(:,2) = u
    case(3)
      V2(:,1) = w
      V2(:,2) = v
      V2(:,2) = u
    end select
    F = Roe(rho, V2, p, Normal_x)
    select case(id)
    case(1)
      continue
    case(2)
      Ftmp = F(2)
      F(2) = F(4) ! w
      F(4) = F(3) ! v
      F(3) = Ftmp ! u
    case(3)
      Ftmp = F(2)
      F(2) = F(3)! v
      F(3) = F(4)! w
      F(4) = Ftmp! u
    end select
  end subroutine flux_Roe2


  attributes(device) subroutine flux_Weighted6(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    real(4), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), w(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(2) slau
    call flux_SLAU6(slau, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    F = sensor * F + (1.d0 - sensor) * KEEP6(rho, u, v, w, uu, p, T, Normal) 
  end subroutine flux_Weighted6


  attributes(device) subroutine flux_Weighted4(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    real(4), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), w(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(2) slau
    call flux_SLAU4(slau, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    F = sensor * F! + (1.d0 - sensor) * KEEP4(rho, u, v, w, uu, p, T, Normal)
  end subroutine flux_Weighted4


  attributes(device) subroutine flux_Weighted2(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(4), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), w(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(8) V2(2,3)
    V2(:,1) = u
    V2(:,2) = v
    V2(:,3) = w
    F = (1.d0 - sensor) * KEEP2(rho, u, v, w, uu, p, T, Normal)
    F = F + sensor * SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
  end subroutine flux_Weighted2


  attributes(device) subroutine flux_Threshold6(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    real(8), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(6), u(6), v(6), w(6), uu(6), p(6), T(6)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP6(rho, u, v, w, uu, p, T, Normal)
    else
      call flux_SLAU6(slau, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    endif
  end subroutine flux_Threshold6


  attributes(device) subroutine flux_Threshold4(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    real(8), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(4), u(4), v(4), w(4), uu(4), p(4), T(4)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP4(rho, u, v, w, uu, p, T, Normal)
    else
      call flux_SLAU4(slau, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    endif
  end subroutine flux_Threshold4


  attributes(device) subroutine flux_Threshold2(id_scheme, id, rho, u, v, w, uu, p, T, Normal, sensor, F)
    use mod_globals, only : id_slau
    real(8), intent(in), value       :: id_scheme
    integer, intent(in), value       :: id
    real(8), intent(in), contiguous  :: rho(2), u(2), v(2), w(2), uu(2), p(2), T(2)
    real(8), intent(in), contiguous  :: Normal(5)
    real(8), intent(in), value       :: sensor
    real(8), intent(out), contiguous :: F(5)
    real(8) V2(2,3)
    V2(:,1) = u
    V2(:,2) = v
    V2(:,3) = w
    if (sensor < threshold) then
      F = KEEP2(rho, u, v, w, uu, p, T, Normal)
    else
      F = SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
    endif
  end subroutine flux_Threshold2


  attributes(global) subroutine calc_E6(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
    real(8) fdx, tmp(6)
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(1,i,j,k)
          u(ii,jt,kt) = Q(2,i,j,k)
          v(ii,jt,kt) = Q(3,i,j,k)
          w(ii,jt,kt) = Q(4,i,j,k)
          p(ii,jt,kt) = Q(5,i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    associate(uu => u)
    if (3 <= i .and. i <= nx-3 .and. 8 <= kind(id_accuracy)) then
      tmp(:) = T(i-2:i+3,j,k)
      call flux6(id_scheme,1,rho(it-2:it+3,jt,kt),u(it-2:it+3,jt,kt),v(it-2:it+3,jt,kt),w(it-2:it+3,jt,kt),&
                 uu(it-2:it+3,jt,kt),p(it-2:it+3,jt,kt),tmp,Normal_x,fdx,E(:,i,j-1,k-1))
    elseif (2 <= i .and. i <= nx-2) then
      tmp(2:5) = T(i-1:i+2,j,k)
      call flux4(id_scheme,1,rho(it-1:it+2,jt,kt),u(it-1:it+2,jt,kt),v(it-1:it+2,jt,kt),w(it-1:it+2,jt,kt),&
                 uu(it-1:it+2,jt,kt),p(it-1:it+2,jt,kt),tmp(2:5),Normal_x,fdx,E(:,i,j-1,k-1))
    else
      tmp(3:4) = T(i:i+1,j,k)
      call flux2(id_scheme,1,rho(it:it+1,jt,kt),u(it:it+1,jt,kt),v(it:it+1,jt,kt),w(it:it+1,jt,kt),&
                 uu(it:it+1,jt,kt),p(it:it+1,jt,kt),tmp(3:4),Normal_x,fdx,E(:,i,j-1,k-1))
    endif
    end associate
  end subroutine calc_E6


  attributes(global) subroutine calc_F6(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
    real(8) fdy, tmp(6)
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(1,i,j,k)
          u(jj,it,kt) = Q(2,i,j,k)
          v(jj,it,kt) = Q(3,i,j,k)
          w(jj,it,kt) = Q(4,i,j,k)
          p(jj,it,kt) = Q(5,i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    associate(vv => v)
    if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
      tmp(:) = T(i,j-2:j+3,k)
      call flux6(id_scheme,2,rho(jt-2:jt+3,it,kt),u(jt-2:jt+3,it,kt),v(jt-2:jt+3,it,kt),w(jt-2:jt+3,it,kt),&
                 vv(jt-2:jt+3,it,kt),p(jt-2:jt+3,it,kt),tmp,Normal_y,fdy,F(:,i-1,j,k-1))
    elseif (2 <= j .and. j <= ny-2) then
      tmp(2:5) = T(i,j-1:j+2,k)
      call flux4(id_scheme,2,rho(jt-1:jt+2,it,kt),u(jt-1:jt+2,it,kt),v(jt-1:jt+2,it,kt),w(jt-1:jt+2,it,kt),&
                 vv(jt-1:jt+2,it,kt),p(jt-1:jt+2,it,kt),tmp(2:5),Normal_y,fdy,F(:,i-1,j,k-1))
    else
      tmp(3:4) = T(i,j:j+1,k)
      call flux2(id_scheme,2,rho(jt:jt+1,it,kt),u(jt:jt+1,it,kt),v(jt:jt+1,it,kt),w(jt:jt+1,it,kt),&
                 vv(jt:jt+1,it,kt),p(jt:jt+1,it,kt),tmp(3:4),Normal_y,fdy,F(:,i-1,j,k-1))
    endif
    end associate
  end subroutine calc_F6


  attributes(global) subroutine calc_G6(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
    real(8) :: fdz, tmp(6)
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-2, threadsG%z+3, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(1,i,j,k)
          u(kk,jt,it) = Q(2,i,j,k)
          v(kk,jt,it) = Q(3,i,j,k)
          w(kk,jt,it) = Q(4,i,j,k)
          p(kk,jt,it) = Q(5,i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    associate(ww => w)
    if (3 <= k .and. k <= nz-3 .and. 8 <= kind(id_accuracy)) then
      tmp(:) = T(i,j,k-2:k+3)
      call flux6(id_scheme,3,rho(kt-2:kt+3,jt,it),u(kt-2:kt+3,jt,it),v(kt-2:kt+3,jt,it),w(kt-2:kt+3,jt,it),&
                 ww(kt-2:kt+3,jt,it),p(kt-2:kt+3,jt,it),tmp,Normal_z,fdz,G(:,i-1,j-1,k))
    elseif (2 <= k .and. k <= nz-2) then
      tmp(2:5) = T(i,j,k-1:k+2)
      call flux4(id_scheme,3,rho(kt-1:kt+2,jt,it),u(kt-1:kt+2,jt,it),v(kt-1:kt+2,jt,it),w(kt-1:kt+2,jt,it),&
                 ww(kt-1:kt+2,jt,it),p(kt-1:kt+2,jt,it),tmp(2:5),Normal_z,fdz,G(:,i-1,j-1,k))
    else
      tmp(3:4) = T(i,j,k:k+1)
      call flux2(id_scheme,3,rho(kt:kt+1,jt,it),u(kt:kt+1,jt,it),v(kt:kt+1,jt,it),w(kt:kt+1,jt,it),&
                 ww(kt:kt+1,jt,it),p(kt:kt+1,jt,it),tmp(3:4),Normal_z,fdz,G(:,i-1,j-1,k))
    endif
    end associate
  end subroutine calc_G6


  attributes(global) subroutine calc_E4(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
    real(8) fdx, tmp(4)
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(ii,jt,kt) = Q(1,i,j,k)
          u(ii,jt,kt) = Q(2,i,j,k)
          v(ii,jt,kt) = Q(3,i,j,k)
          w(ii,jt,kt) = Q(4,i,j,k)
          p(ii,jt,kt) = Q(5,i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    associate(uu => u)
    if (2 <= i .and. i <= nx-2) then
      tmp(:) = T(i-1:i+2,j,k)
      call flux4(id_scheme,1,rho(it-1:it+2,jt,kt),u(it-1:it+2,jt,kt),v(it-1:it+2,jt,kt),w(it-1:it+2,jt,kt),&
                 uu(it-1:it+2,jt,kt),p(it-1:it+2,jt,kt),tmp,Normal_x,fdx,E(:,i,j-1,k-1))
    else
      tmp(2:3) = T(i:i+1,j,k)
      call flux2(id_scheme,1,rho(it:it+1,jt,kt),u(it:it+1,jt,kt),v(it:it+1,jt,kt),w(it:it+1,jt,kt),&
                 uu(it:it+1,jt,kt),p(it:it+1,jt,kt),tmp(2:3),Normal_x,fdx,E(:,i,j-1,k-1))
    endif
    end associate
  end subroutine calc_E4


  attributes(global) subroutine calc_F4(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
    real(8) fdy, tmp(4)
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(jj,it,kt) = Q(1,i,j,k)
          u(jj,it,kt) = Q(2,i,j,k)
          v(jj,it,kt) = Q(3,i,j,k)
          w(jj,it,kt) = Q(4,i,j,k)
          p(jj,it,kt) = Q(5,i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    associate(vv => v)
    if (2 <= j .and. j <= ny-2) then
      tmp(:) = T(i,j-1:j+2,k)
      call flux4(id_scheme,2,rho(jt-1:jt+2,it,kt),u(jt-1:jt+2,it,kt),v(jt-1:jt+2,it,kt),w(jt-1:jt+2,it,kt),&
                 vv(jt-1:jt+2,it,kt),p(jt-1:jt+2,it,kt),tmp,Normal_y,fdy,F(:,i-1,j,k-1))
    else
      tmp(2:3) = T(i,j:j+1,k)
      call flux2(id_scheme,2,rho(jt:jt+1,it,kt),u(jt:jt+1,it,kt),v(jt:jt+1,it,kt),w(jt:jt+1,it,kt),&
                 vv(jt:jt+1,it,kt),p(jt:jt+1,it,kt),tmp(2:3),Normal_y,fdy,F(:,i-1,j,k-1))
    endif
    end associate
  end subroutine calc_F4


  attributes(global) subroutine calc_G4(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
    real(8) :: fdz, tmp(4)
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-1, threadsG%z+2, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        rho(kk,jt,it) = Q(1,i,j,k)
          u(kk,jt,it) = Q(2,i,j,k)
          v(kk,jt,it) = Q(3,i,j,k)
          w(kk,jt,it) = Q(4,i,j,k)
          p(kk,jt,it) = Q(5,i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    associate(ww => w)
    if (2 <= k .and. k <= nz-2) then
      tmp(:) = T(i,j,k-1:k+2)
      call flux4(id_scheme,3,rho(kt-1:kt+2,jt,it),u(kt-1:kt+2,jt,it),v(kt-1:kt+2,jt,it),w(kt-1:kt+2,jt,it),&
                 ww(kt-1:kt+2,jt,it),p(kt-1:kt+2,jt,it),tmp,Normal_z,fdz,G(:,i-1,j-1,k))
    else
      tmp(2:3) = T(i,j,k:k+1)
      call flux2(id_scheme,3,rho(kt:kt+1,jt,it),u(kt:kt+1,jt,it),v(kt:kt+1,jt,it),w(kt:kt+1,jt,it),&
                 ww(kt:kt+1,jt,it),p(kt:kt+1,jt,it),tmp(2:3),Normal_z,fdz,G(:,i-1,j-1,k))
    endif
    end associate
  end subroutine calc_G4


  attributes(global) subroutine calc_E2(id_accuracy, nx, ny, nz, Q, T, sensor, E)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    real(8) fdx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    rho = Q(1,i:i+1,j,k)
    u   = Q(2,i:i+1,j,k)
    v   = Q(3,i:i+1,j,k)
    w   = Q(4,i:i+1,j,k)
    p   = Q(5,i:i+1,j,k)
    tmp = T(i:i+1,j,k)
    call flux2(id_scheme, 1, rho, u, v, w, u, p, tmp, Normal_x, fdx, E(:,i,j-1,k-1))
  end subroutine calc_E2


  attributes(global) subroutine calc_F2(id_accuracy, nx, ny, nz, Q, T, sensor, F)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    real(8) fdy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    rho = Q(1,i,j:j+1,k)
    u   = Q(2,i,j:j+1,k)
    v   = Q(3,i,j:j+1,k)
    w   = Q(4,i,j:j+1,k)
    p   = Q(5,i,j:j+1,k)
    tmp = T(i,j:j+1,k)
    call flux2(id_scheme, 2, rho, u, v, w, v, p, tmp, Normal_y, fdy, F(:,i-1,j,k-1))
  end subroutine calc_F2


  attributes(global) subroutine calc_G2(id_accuracy, nx, ny, nz, Q, T, sensor, G)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: T, sensor
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    real(8) :: fdz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    rho = Q(1,i,j,k:k+1)
    u   = Q(2,i,j,k:k+1)
    v   = Q(3,i,j,k:k+1)
    w   = Q(4,i,j,k:k+1)
    p   = Q(5,i,j,k:k+1)
    tmp = T(i,j,k:k+1)
    call flux2(id_scheme, 3, rho, u, v, w, w, p, tmp, Normal_z, fdz, G(:,i-1,j-1,k))
  end subroutine calc_G2
end module calc_flux

