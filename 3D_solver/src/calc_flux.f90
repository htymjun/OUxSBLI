module calc_flux
  use mod_globals, only : id_scheme, id_sensor, id_muscl, gamma, threshold, threadsE, threadsF, threadsG
  use calc_keep
  use calc_slau
  use calc_hybrid
  use calc_muscl
  implicit none
  interface flux6
    module procedure flux_KEEP6, flux_SLAU6, flux_Weighted6, flux_Threshold6
  end interface flux6

  interface flux4
    module procedure flux_KEEP4, flux_SLAU4, flux_Weighted4, flux_Threshold4
  end interface flux4

  interface flux2
    module procedure flux_KEEP2, flux_SLAU2, flux_Weighted2, flux_Threshold2
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
  attributes(device) function flux_KEEP6(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(6), u(6), v(6), w(6), uu(6), p(6)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8) F(5)
    F = KEEP6(rho, u, v, w, uu, p, Normal)
  end function flux_KEEP6

  attributes(device) function flux_KEEP4(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(4), u(4), v(4), w(4), uu(4), p(4)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8) F(5)
    F = KEEP4(rho, u, v, w, uu, p, Normal)
  end function flux_KEEP4
  
  attributes(device) function flux_KEEP2(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    integer(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value         :: id
    real(8), intent(in), contiguous    :: rho(2), u(2), v(2), w(2), uu(2), p(2)
    real(8), intent(in), contiguous    :: Normal(5)
    real(8), intent(in), value         :: sensor
    real(8) F(5)
    F = KEEP2(rho, u, v, w, uu, p, Normal)
  end function flux_KEEP2

  attributes(device) function flux_SLAU6(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(6), u(6), v(6), w(6), uu(6), p(6)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) wiggle, rho2(2), p2(2), V2(2,3), F(5)
    wiggle = wiggle_detector(p(2:5))
    call calc_6points(sensor, rho, u, v, w, p, rho2, p2, V2)
    F = SLAU(id_slau, id, rho2, p2, V2, Normal, wiggle)
  end function flux_SLAU6

  attributes(device) function flux_SLAU4(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(4), u(4), v(4), w(4), uu(4), p(4)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) wiggle, rho2(2), p2(2), V2(2,3), F(5)
    wiggle = wiggle_detector(p)
    call calc_4points(sensor, rho, u, v, w, p, rho2, p2, V2)
    F = SLAU(id_slau, id, rho2, p2, V2, Normal, wiggle)
  end function flux_SLAU4

  attributes(device) function flux_SLAU2(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    use mod_globals, only : id_slau
    real(kind=2), intent(in), value :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(2), u(2), v(2), w(2), uu(2), p(2)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) V2(2,3), F(5)
    V2(:,1) = u
    V2(:,2) = v
    V2(:,3) = w
    F = SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
  end function flux_SLAU2

  attributes(device) function flux_Weighted6(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    real(4), intent(in), value      :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(6), u(6), v(6), w(6), uu(6), p(6)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) F(5)
    real(2) slau
    F = (1.d0 - sensor) * KEEP6(rho, u, v, w, uu, p, Normal) &
        + sensor * flux_SLAU6(slau, id, rho, u, v, w, uu, p, Normal, sensor)
  end function flux_Weighted6

  attributes(device) function flux_Weighted4(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    real(4), intent(in), value      :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(4), u(4), v(4), w(4), uu(4), p(4)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) F(5)
    real(2) slau
    F = (1.d0 - sensor) * KEEP4(rho, u, v, w, uu, p, Normal) &
        + sensor * flux_SLAU4(slau, id, rho, u, v, w, uu, p, Normal, sensor)
  end function flux_Weighted4

  attributes(device) function flux_Weighted2(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    use mod_globals, only : id_slau
    real(4), intent(in), value      :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(2), u(2), v(2), w(2), uu(2), p(2)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) V2(2,3), F(5)
    V2(:,1) = u
    V2(:,2) = v
    V2(:,3) = w
    F = (1.d0 - sensor) * KEEP2(rho, u, v, w, uu, p, Normal) &
        + sensor * SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
  end function flux_Weighted2

  attributes(device) function flux_Threshold6(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    real(8), intent(in), value      :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(6), u(6), v(6), w(6), uu(6), p(6)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) F(5)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP6(rho, u, v, w, uu, p, Normal)
    else
      F = flux_SLAU6(slau, id, rho, u, v, w, uu, p, Normal, sensor)
    endif
  end function flux_Threshold6

  attributes(device) function flux_Threshold4(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    real(8), intent(in), value      :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(4), u(4), v(4), w(4), uu(4), p(4)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) F(5)
    real(2) slau
    if (sensor < threshold) then
      F = KEEP4(rho, u, v, w, uu, p, Normal)
    else
      F = flux_SLAU4(slau, id, rho, u, v, w, uu, p, Normal, sensor)
    endif
  end function flux_Threshold4

  attributes(device) function flux_Threshold2(id_scheme, id, rho, u, v, w, uu, p, Normal, sensor) result(F)
    use mod_globals, only : id_slau
    real(8), intent(in), value      :: id_scheme
    integer, intent(in), value      :: id
    real(8), intent(in), contiguous :: rho(2), u(2), v(2), w(2), uu(2), p(2)
    real(8), intent(in), contiguous :: Normal(5)
    real(8), intent(in), value      :: sensor
    real(8) V2(2,3), F(5)
    V2(:,1) = u
    V2(:,2) = v
    V2(:,3) = w
    if (sensor < threshold) then
      F = KEEP2(rho, u, v, w, uu, p, Normal)
    else
      F = SLAU(id_slau, id, rho, p, V2, Normal, 1.d0)
    endif
  end function flux_Threshold2

  attributes(global) subroutine calc_E6(id_accuracy, nx, ny, nz, Q, sensor, E)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt, idx
    real(8), dimension(-1:threadsE%x+3,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
    real(8) fdx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(it,jt,kt) = Q(1,i,j,k)
      u(it,jt,kt) = Q(2,i,j,k)
      v(it,jt,kt) = Q(3,i,j,k)
      w(it,jt,kt) = Q(4,i,j,k)
      p(it,jt,kt) = Q(5,i,j,k)
    if (it == blockDim%x) then ! 2nd-order
      rho(it+1,jt,kt) = Q(1,i+1,j,k)
        u(it+1,jt,kt) = Q(2,i+1,j,k)
        v(it+1,jt,kt) = Q(3,i+1,j,k)
        w(it+1,jt,kt) = Q(4,i+1,j,k)
        p(it+1,jt,kt) = Q(5,i+1,j,k)
      if (i <= nx-3) then ! 6th-order
        do idx = 2, 3
          rho(it+idx,jt,kt) = Q(1,i+idx,j,k)
            u(it+idx,jt,kt) = Q(2,i+idx,j,k)
            v(it+idx,jt,kt) = Q(3,i+idx,j,k)
            w(it+idx,jt,kt) = Q(4,i+idx,j,k)
            p(it+idx,jt,kt) = Q(5,i+idx,j,k)
        enddo
      elseif (i <= nx-2) then ! 4th-order
        rho(it+2,jt,kt) = Q(1,i+2,j,k)
          u(it+2,jt,kt) = Q(2,i+2,j,k)
          v(it+2,jt,kt) = Q(3,i+2,j,k)
          w(it+2,jt,kt) = Q(4,i+2,j,k)
          p(it+2,jt,kt) = Q(5,i+2,j,k)
      endif
    elseif (it == 1) then
      if (3 <= i) then ! 6th-order
        do idx = -2, -1
          rho(it+idx,jt,kt) = Q(1,i+idx,j,k)
            u(it+idx,jt,kt) = Q(2,i+idx,j,k)
            v(it+idx,jt,kt) = Q(3,i+idx,j,k)
            w(it+idx,jt,kt) = Q(4,i+idx,j,k)
            p(it+idx,jt,kt) = Q(5,i+idx,j,k)
        enddo
      elseif (2 <= i) then ! 4th-order
        rho(it-1,jt,kt) = Q(1,i-1,j,k)
          u(it-1,jt,kt) = Q(2,i-1,j,k)
          v(it-1,jt,kt) = Q(3,i-1,j,k)
          w(it-1,jt,kt) = Q(4,i-1,j,k)
          p(it-1,jt,kt) = Q(5,i-1,j,k)
      endif
    endif
    call syncthreads()
    associate(uu => u)
    if (3 <= i .and. i <= nx-3 .and. 8 <= kind(id_accuracy)) then
      E(:,i,j-1,k-1) = flux6(id_scheme,1,rho(it-2:it+3,jt,kt),u(it-2:it+3,jt,kt),&
                          v(it-2:it+3,jt,kt),w(it-2:it+3,jt,kt),uu(it-2:it+3,jt,kt),p(it-2:it+3,jt,kt),Normal_x,fdx)
    elseif (2 <= i .and. i <= nx-2) then
      E(:,i,j-1,k-1) = flux4(id_scheme,1,rho(it-1:it+2,jt,kt),u(it-1:it+2,jt,kt),&
                          v(it-1:it+2,jt,kt),w(it-1:it+2,jt,kt),uu(it-1:it+2,jt,kt),p(it-1:it+2,jt,kt),Normal_x,fdx)
    else
      E(:,i,j-1,k-1) = flux2(id_scheme,1,rho(it:it+1,jt,kt),u(it:it+1,jt,kt),&
                          v(it:it+1,jt,kt),w(it:it+1,jt,kt),uu(it:it+1,jt,kt),p(it:it+1,jt,kt),Normal_x,fdx)
    endif
    end associate
  end subroutine calc_E6

  attributes(global) subroutine calc_F6(id_accuracy, nx, ny, nz, Q, sensor, F)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt, idy
    real(8), dimension(-1:threadsF%y+3,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
    real(8) fdy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i = (blockIdx%x-1)*blockDim%x + it + 1
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt + 1
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(jt,it,kt) = Q(1,i,j,k)
      u(jt,it,kt) = Q(2,i,j,k)
      v(jt,it,kt) = Q(3,i,j,k)
      w(jt,it,kt) = Q(4,i,j,k)
      p(jt,it,kt) = Q(5,i,j,k)
    if (jt == blockDim%y) then ! 2nd-order
      rho(jt+1,it,kt) = Q(1,i,j+1,k)
        u(jt+1,it,kt) = Q(2,i,j+1,k)
        v(jt+1,it,kt) = Q(3,i,j+1,k)
        w(jt+1,it,kt) = Q(4,i,j+1,k)
        p(jt+1,it,kt) = Q(5,i,j+1,k)
      if (j <= ny-3) then ! 6th-order
        do idy = 2, 3
          rho(jt+idy,it,kt) = Q(1,i,j+idy,k)
            u(jt+idy,it,kt) = Q(2,i,j+idy,k)
            v(jt+idy,it,kt) = Q(3,i,j+idy,k)
            w(jt+idy,it,kt) = Q(4,i,j+idy,k)
            p(jt+idy,it,kt) = Q(5,i,j+idy,k)
        enddo
      elseif (j <= ny-2) then ! 4th-order
        rho(jt+2,it,kt) = Q(1,i,j+2,k)
          u(jt+2,it,kt) = Q(2,i,j+2,k)
          v(jt+2,it,kt) = Q(3,i,j+2,k)
          w(jt+2,it,kt) = Q(4,i,j+2,k)
          p(jt+2,it,kt) = Q(5,i,j+2,k)
      endif
    elseif (3 <= j .and. jt == 1) then
      if (3 <= j) then ! 6th-order
        do idy = -2, -1
          rho(jt+idy,it,kt) = Q(1,i,j+idy,k)
            u(jt+idy,it,kt) = Q(2,i,j+idy,k)
            v(jt+idy,it,kt) = Q(3,i,j+idy,k)
            w(jt+idy,it,kt) = Q(4,i,j+idy,k)
            p(jt+idy,it,kt) = Q(5,i,j+idy,k)
        enddo
      elseif (2 <= j) then ! 4th-order
        rho(jt-1,it,kt) = Q(1,i,j-1,k)
          u(jt-1,it,kt) = Q(2,i,j-1,k)
          v(jt-1,it,kt) = Q(3,i,j-1,k)
          w(jt-1,it,kt) = Q(4,i,j-1,k)
          p(jt-1,it,kt) = Q(5,i,j-1,k)
      endif
    endif
    call syncthreads()
    associate(vv => v)
    if (3 <= j .and. j <= ny-3 .and. 8 <= kind(id_accuracy)) then
      F(:,i-1,j,k-1) = flux6(id_scheme,2,rho(jt-2:jt+3,it,kt),u(jt-2:jt+3,it,kt),&
                            v(jt-2:jt+3,it,kt),w(jt-2:jt+3,it,kt),vv(jt-2:jt+3,it,kt),p(jt-2:jt+3,it,kt),Normal_y,fdy)
    elseif (2 <= j .and. j <= ny-2) then
      F(:,i-1,j,k-1) = flux4(id_scheme,2,rho(jt-1:jt+2,it,kt),u(jt-1:jt+2,it,kt),&
                            v(jt-1:jt+2,it,kt),w(jt-1:jt+2,it,kt),vv(jt-1:jt+2,it,kt),p(jt-1:jt+2,it,kt),Normal_y,fdy)
    else
      F(:,i-1,j,k-1) = flux2(id_scheme,2,rho(jt:jt+1,it,kt),u(jt:jt+1,it,kt),&
                            v(jt:jt+1,it,kt),w(jt:jt+1,it,kt),vv(jt:jt+1,it,kt),p(jt:jt+1,it,kt),Normal_y,fdy)
    endif
    end associate
  end subroutine calc_F6

  attributes(global) subroutine calc_G6(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt, idz
    real(8), dimension(-1:threadsG%z+3,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
    real(8) :: fdz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i = (blockIdx%x-1)*blockDim%x + it + 1
    j = (blockIdx%y-1)*blockDim%y + jt + 1
    k = (blockIdx%z-1)*blockDim%z + kt
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(kt,jt,it) = Q(1,i,j,k)
      u(kt,jt,it) = Q(2,i,j,k)
      v(kt,jt,it) = Q(3,i,j,k)
      w(kt,jt,it) = Q(4,i,j,k)
      p(kt,jt,it) = Q(5,i,j,k)
    if (kt == blockDim%z) then ! 2nd-order
      rho(kt+1,jt,it) = Q(1,i,j,k+1)
        u(kt+1,jt,it) = Q(2,i,j,k+1)
        v(kt+1,jt,it) = Q(3,i,j,k+1)
        w(kt+1,jt,it) = Q(4,i,j,k+1)
        p(kt+1,jt,it) = Q(5,i,j,k+1)
      if (k <= nz-3) then ! 6th-order
        do idz = 2, 3
          rho(kt+idz,jt,it) = Q(1,i,j,k+idz)
            u(kt+idz,jt,it) = Q(2,i,j,k+idz)
            v(kt+idz,jt,it) = Q(3,i,j,k+idz)
            w(kt+idz,jt,it) = Q(4,i,j,k+idz)
            p(kt+idz,jt,it) = Q(5,i,j,k+idz)
        enddo
      elseif (k <= nz-2) then ! 4th-order
        rho(kt+2,jt,it) = Q(1,i,j,k+2)
          u(kt+2,jt,it) = Q(2,i,j,k+2)
          v(kt+2,jt,it) = Q(3,i,j,k+2)
          w(kt+2,jt,it) = Q(4,i,j,k+2)
          p(kt+2,jt,it) = Q(5,i,j,k+2)
      endif
    elseif (kt == 1) then
      if (3 <= k) then ! 6th-order
        do idz = -2, -1
          rho(kt+idz,jt,it) = Q(1,i,j,k+idz)
            u(kt+idz,jt,it) = Q(2,i,j,k+idz)
            v(kt+idz,jt,it) = Q(3,i,j,k+idz)
            w(kt+idz,jt,it) = Q(4,i,j,k+idz)
            p(kt+idz,jt,it) = Q(5,i,j,k+idz)
        enddo
      elseif (2 <= k) then ! 4th-order
        rho(kt-1,jt,it) = Q(1,i,j,k-1)
          u(kt-1,jt,it) = Q(2,i,j,k-1)
          v(kt-1,jt,it) = Q(3,i,j,k-1)
          w(kt-1,jt,it) = Q(4,i,j,k-1)
          p(kt-1,jt,it) = Q(5,i,j,k-1)
      endif
    endif
    call syncthreads()
    associate(ww => w)
    if (3 <= k .and. k <= nz-3 .and. 8 <= kind(id_accuracy)) then
      G(:,i-1,j-1,k) = flux6(id_scheme,3,rho(kt-2:kt+3,jt,it),u(kt-2:kt+3,jt,it), &
                        v(kt-2:kt+3,jt,it),w(kt-2:kt+3,jt,it),ww(kt-2:kt+3,jt,it),p(kt-2:kt+3,jt,it),Normal_z,fdz)
    elseif (2 <= k .and. k <= nz-2) then
      G(:,i-1,j-1,k) = flux4(id_scheme,3,rho(kt-1:kt+2,jt,it),u(kt-1:kt+2,jt,it), &
                        v(kt-1:kt+2,jt,it),w(kt-1:kt+2,jt,it),ww(kt-1:kt+2,jt,it),p(kt-1:kt+2,jt,it),Normal_z,fdz)
    else
      G(:,i-1,j-1,k) = flux2(id_scheme,3,rho(kt:kt+1,jt,it),u(kt:kt+1,jt,it), &
                        v(kt:kt+1,jt,it),w(kt:kt+1,jt,it),ww(kt:kt+1,jt,it),p(kt:kt+1,jt,it),Normal_z,fdz)
    endif
    end associate
  end subroutine calc_G6
  
  attributes(global) subroutine calc_E4(id_accuracy, nx, ny, nz, Q, sensor, E)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt
    real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
    real(8) fdx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(it,jt,kt) = Q(1,i,j,k)
      u(it,jt,kt) = Q(2,i,j,k)
      v(it,jt,kt) = Q(3,i,j,k)
      w(it,jt,kt) = Q(4,i,j,k)
      p(it,jt,kt) = Q(5,i,j,k)
    if (it == blockDim%x) then ! 2nd-order
      rho(it+1,jt,kt) = Q(1,i+1,j,k)
        u(it+1,jt,kt) = Q(2,i+1,j,k)
        v(it+1,jt,kt) = Q(3,i+1,j,k)
        w(it+1,jt,kt) = Q(4,i+1,j,k)
        p(it+1,jt,kt) = Q(5,i+1,j,k)
      if (i <= nx-2) then ! 4th-order
        rho(it+2,jt,kt) = Q(1,i+2,j,k)
          u(it+2,jt,kt) = Q(2,i+2,j,k)
          v(it+2,jt,kt) = Q(3,i+2,j,k)
          w(it+2,jt,kt) = Q(4,i+2,j,k)
          p(it+2,jt,kt) = Q(5,i+2,j,k)
      endif
    elseif (2 <= i .and. it == 1) then ! 4th-order
      rho(it-1,jt,kt) = Q(1,i-1,j,k)
        u(it-1,jt,kt) = Q(2,i-1,j,k)
        v(it-1,jt,kt) = Q(3,i-1,j,k)
        w(it-1,jt,kt) = Q(4,i-1,j,k)
        p(it-1,jt,kt) = Q(5,i-1,j,k)
    endif
    call syncthreads()
    associate(uu => u)
    if (2 <= i .and. i <= nx-2) then
      E(:,i,j-1,k-1) = flux4(id_scheme,1,rho(it-1:it+2,jt,kt),u(it-1:it+2,jt,kt),&
                          v(it-1:it+2,jt,kt),w(it-1:it+2,jt,kt),uu(it-1:it+2,jt,kt),p(it-1:it+2,jt,kt),Normal_x,fdx)
    else
      E(:,i,j-1,k-1) = flux2(id_scheme,1,rho(it:it+1,jt,kt),u(it:it+1,jt,kt),&
                          v(it:it+1,jt,kt),w(it:it+1,jt,kt),uu(it:it+1,jt,kt),p(it:it+1,jt,kt),Normal_x,fdx)
    endif
    end associate
  end subroutine calc_E4

  attributes(global) subroutine calc_F4(id_accuracy, nx, ny, nz, Q, sensor, F)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt
    real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
    real(8) fdy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i = (blockIdx%x-1)*blockDim%x + it + 1
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt + 1
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(jt,it,kt) = Q(1,i,j,k)
      u(jt,it,kt) = Q(2,i,j,k)
      v(jt,it,kt) = Q(3,i,j,k)
      w(jt,it,kt) = Q(4,i,j,k)
      p(jt,it,kt) = Q(5,i,j,k)
    if (jt == blockDim%y) then ! 2nd-order
      rho(jt+1,it,kt) = Q(1,i,j+1,k)
        u(jt+1,it,kt) = Q(2,i,j+1,k)
        v(jt+1,it,kt) = Q(3,i,j+1,k)
        w(jt+1,it,kt) = Q(4,i,j+1,k)
        p(jt+1,it,kt) = Q(5,i,j+1,k)
      if (j <= ny-2) then ! 4th-order
        rho(jt+2,it,kt) = Q(1,i,j+2,k)
          u(jt+2,it,kt) = Q(2,i,j+2,k)
          v(jt+2,it,kt) = Q(3,i,j+2,k)
          w(jt+2,it,kt) = Q(4,i,j+2,k)
          p(jt+2,it,kt) = Q(5,i,j+2,k)
      endif
    elseif (2 <= j .and. jt == 1) then ! 4th-order
      rho(jt-1,it,kt) = Q(1,i,j-1,k)
        u(jt-1,it,kt) = Q(2,i,j-1,k)
        v(jt-1,it,kt) = Q(3,i,j-1,k)
        w(jt-1,it,kt) = Q(4,i,j-1,k)
        p(jt-1,it,kt) = Q(5,i,j-1,k)
    endif
    call syncthreads()
    associate(vv => v)
    if (2 <= j .and. j <= ny-2) then
      F(:,i-1,j,k-1) = flux4(id_scheme,2,rho(jt-1:jt+2,it,kt),u(jt-1:jt+2,it,kt),&
                            v(jt-1:jt+2,it,kt),w(jt-1:jt+2,it,kt),vv(jt-1:jt+2,it,kt),p(jt-1:jt+2,it,kt),Normal_y,fdy)
    else
      F(:,i-1,j,k-1) = flux2(id_scheme,2,rho(jt:jt+1,it,kt),u(jt:jt+1,it,kt),&
                            v(jt:jt+1,it,kt),w(jt:jt+1,it,kt),vv(jt:jt+1,it,kt),p(jt:jt+1,it,kt),Normal_y,fdy)
    endif
    end associate
  end subroutine calc_F4

  attributes(global) subroutine calc_G4(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt
    real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
    real(8) :: fdz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i = (blockIdx%x-1)*blockDim%x + it + 1
    j = (blockIdx%y-1)*blockDim%y + jt + 1
    k = (blockIdx%z-1)*blockDim%z + kt
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(kt,jt,it) = Q(1,i,j,k)
      u(kt,jt,it) = Q(2,i,j,k)
      v(kt,jt,it) = Q(3,i,j,k)
      w(kt,jt,it) = Q(4,i,j,k)
      p(kt,jt,it) = Q(5,i,j,k)
    if (kt == blockDim%z) then ! 2nd-order
      rho(kt+1,jt,it) = Q(1,i,j,k+1)
        u(kt+1,jt,it) = Q(2,i,j,k+1)
        v(kt+1,jt,it) = Q(3,i,j,k+1)
        w(kt+1,jt,it) = Q(4,i,j,k+1)
        p(kt+1,jt,it) = Q(5,i,j,k+1)
      if (k <= nz-2) then ! 4th-order
        rho(kt+2,jt,it) = Q(1,i,j,k+2)
          u(kt+2,jt,it) = Q(2,i,j,k+2)
          v(kt+2,jt,it) = Q(3,i,j,k+2)
          w(kt+2,jt,it) = Q(4,i,j,k+2)
          p(kt+2,jt,it) = Q(5,i,j,k+2)
      endif
    elseif (2 <= k .and. kt == 1) then ! 4th-order
      rho(kt-1,jt,it) = Q(1,i,j,k-1)
        u(kt-1,jt,it) = Q(2,i,j,k-1)
        v(kt-1,jt,it) = Q(3,i,j,k-1)
        w(kt-1,jt,it) = Q(4,i,j,k-1)
        p(kt-1,jt,it) = Q(5,i,j,k-1)
    endif
    call syncthreads()
    associate(ww => w)
    if (2 <= k .and. k <= nz-2) then
      G(:,i-1,j-1,k) = flux4(id_scheme,3,rho(kt-1:kt+2,jt,it),u(kt-1:kt+2,jt,it), &
                        v(kt-1:kt+2,jt,it),w(kt-1:kt+2,jt,it),ww(kt-1:kt+2,jt,it),p(kt-1:kt+2,jt,it),Normal_z,fdz)
    else
      G(:,i-1,j-1,k) = flux2(id_scheme,3,rho(kt:kt+1,jt,it),u(kt:kt+1,jt,it), &
                        v(kt:kt+1,jt,it),w(kt:kt+1,jt,it),ww(kt:kt+1,jt,it),p(kt:kt+1,jt,it),Normal_z,fdz)
    endif
    end associate
  end subroutine calc_G4
  
  attributes(global) subroutine calc_E2(id_accuracy, nx, ny, nz, Q, sensor, E)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k, it, jt, kt
    real(8), dimension(0:threadsE%x+2,threadsE%y,threadsE%z), shared :: rho, u, v, w, p
    real(8) fdx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    fdx = 0.5d0 * (sensor(i,j,k) + sensor(i+1,j,k))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(it,jt,kt) = Q(1,i,j,k)
      u(it,jt,kt) = Q(2,i,j,k)
      v(it,jt,kt) = Q(3,i,j,k)
      w(it,jt,kt) = Q(4,i,j,k)
      p(it,jt,kt) = Q(5,i,j,k)
    if (it == blockDim%x) then
      rho(it+1,jt,kt) = Q(1,i+1,j,k)
        u(it+1,jt,kt) = Q(2,i+1,j,k)
        v(it+1,jt,kt) = Q(3,i+1,j,k)
        w(it+1,jt,kt) = Q(4,i+1,j,k)
        p(it+1,jt,kt) = Q(5,i+1,j,k)
    endif
    call syncthreads()
    associate(uu => u)
    E(:,i,j-1,k-1) = flux2(id_scheme,1,rho(it:it+1,jt,kt),u(it:it+1,jt,kt),&
                           v(it:it+1,jt,kt),w(it:it+1,jt,kt),uu(it:it+1,jt,kt),p(it:it+1,jt,kt),Normal_x,fdx)
    end associate
  end subroutine calc_E2

  attributes(global) subroutine calc_F2(id_accuracy, nx, ny, nz, Q, sensor, F)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k, it, jt, kt
    real(8), dimension(0:threadsF%y+2,threadsF%x,threadsF%z), shared :: rho, u, v, w, p
    real(8) fdy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i = (blockIdx%x-1)*blockDim%x + it + 1
    j = (blockIdx%y-1)*blockDim%y + jt
    k = (blockIdx%z-1)*blockDim%z + kt + 1
    fdy = 0.5d0 * (sensor(i,j,k) + sensor(i,j+1,k))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(jt,it,kt) = Q(1,i,j,k)
      u(jt,it,kt) = Q(2,i,j,k)
      v(jt,it,kt) = Q(3,i,j,k)
      w(jt,it,kt) = Q(4,i,j,k)
      p(jt,it,kt) = Q(5,i,j,k)
    if (jt == blockDim%y) then
      rho(jt+1,it,kt) = Q(1,i,j+1,k)
        u(jt+1,it,kt) = Q(2,i,j+1,k)
        v(jt+1,it,kt) = Q(3,i,j+1,k)
        w(jt+1,it,kt) = Q(4,i,j+1,k)
        p(jt+1,it,kt) = Q(5,i,j+1,k)
    endif
    call syncthreads()
    associate(vv => v)
    F(:,i-1,j,k-1) = flux2(id_scheme,2,rho(jt:jt+1,it,kt),u(jt:jt+1,it,kt),&
                           v(jt:jt+1,it,kt),w(jt:jt+1,it,kt),vv(jt:jt+1,it,kt),p(jt:jt+1,it,kt),Normal_y,fdy)
    end associate
  end subroutine calc_F2
  
  attributes(global) subroutine calc_G2(id_accuracy, nx, ny, nz, Q, sensor, G)
    use mod_globals, only  : id_scheme
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in), value                 :: id_accuracy
    integer, intent(in), value                         :: nx, ny, nz
    real(8), intent(in), dimension(5,nx,ny,nz), device :: Q
    real(8), intent(in), dimension(nx,ny,nz), device   :: sensor
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k, it, jt, kt
    real(8), dimension(0:threadsG%z+2,threadsG%y,threadsG%x), shared :: rho, u, v, w, p
    real(8) :: fdz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i = (blockIdx%x-1)*blockDim%x + it + 1
    j = (blockIdx%y-1)*blockDim%y + jt + 1
    k = (blockIdx%z-1)*blockDim%z + kt
    fdz = 0.5d0 * (sensor(i,j,k) + sensor(i,j,k+1))
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho(kt,jt,it) = Q(1,i,j,k)
      u(kt,jt,it) = Q(2,i,j,k)
      v(kt,jt,it) = Q(3,i,j,k)
      w(kt,jt,it) = Q(4,i,j,k)
      p(kt,jt,it) = Q(5,i,j,k)
    if (kt == blockDim%z) then
      rho(kt+1,jt,it) = Q(1,i,j,k+1)
        u(kt+1,jt,it) = Q(2,i,j,k+1)
        v(kt+1,jt,it) = Q(3,i,j,k+1)
        w(kt+1,jt,it) = Q(4,i,j,k+1)
        p(kt+1,jt,it) = Q(5,i,j,k+1)
    endif
    call syncthreads()
    associate(ww => w)
    G(:,i-1,j-1,k) = flux2(id_scheme,3,rho(kt:kt+1,jt,it),u(kt:kt+1,jt,it), &
                           v(kt:kt+1,jt,it),w(kt:kt+1,jt,it),ww(kt:kt+1,jt,it),p(kt:kt+1,jt,it),Normal_z,fdz)
    end associate
  end subroutine calc_G2
end module calc_flux

