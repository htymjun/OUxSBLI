module calc_visc2
  use mod_globals, only : id_visc, gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third
  use calc_visc_common
  use calc_rand
  implicit none
contains
  !$dir inline
  attributes(device) subroutine store_shared_x(nx, ny, nz, i, j, k, it, jt, kt, Q, u, v, w)
    integer, intent(in), value  :: nx, ny, nz, i, j, k, it, jt, kt
    real(8), intent(in), device :: Q(5,nx,ny,nz)
    real(8), intent(inout)      :: u(threadsEv%x+1,0:threadsEv%y+1,0:threadsEv%z+1)
    real(8), intent(inout)      :: v(threadsEv%x+1,0:threadsEv%y+1,threadsEv%z)
    real(8), intent(inout)      :: w(threadsEv%x+1,0:threadsEv%z+1,threadsEv%y)
    ! k-1 ##########################################
    if (kt == 1) then
      if (jt == 1) then
        u(it,jt-1,kt-1) = Q(2,i,j-1,k-1)
        if (it == blockDim%x) then
          u(it+1,jt-1,kt-1) = Q(2,i+1,j-1,k-1)
        endif
      endif
      u(it,jt,kt-1) = Q(2,i,j,k-1)
      w(it,kt-1,jt) = Q(4,i,j,k-1)
      if (it == blockDim%x) then
        u(it+1,jt,kt-1) = Q(2,i+1,j,k-1)
        w(it+1,kt-1,jt) = Q(4,i+1,j,k-1)
      endif
      if (jt == blockDim%y) then
        u(it,jt+1,kt-1) = Q(2,i,j+1,k-1)
        if (it == blockDim%x) then
          u(it+1,jt+1,kt-1) = Q(2,i+1,j+1,k-1)
        endif
      endif
    endif
    ! k ############################################
    if (jt == 1) then
      u(it,jt-1,kt) = Q(2,i,j-1,k)
      v(it,jt-1,kt) = Q(3,i,j-1,k)
      if (it == blockDim%x) then
        u(it+1,jt-1,kt) = Q(2,i+1,j-1,k)
        v(it+1,jt-1,kt) = Q(3,i+1,j-1,k)
      endif
    endif
    u(it,jt,kt) = Q(2,i,j,k)
    v(it,jt,kt) = Q(3,i,j,k)
    w(it,kt,jt) = Q(4,i,j,k)
    if (it == blockDim%x) then
      u(it+1,jt,kt) = Q(2,i+1,j,k)
      v(it+1,jt,kt) = Q(3,i+1,j,k)
      w(it+1,kt,jt) = Q(4,i+1,j,k)
    endif
    if (jt == blockDim%y) then
      u(it,jt+1,kt) = Q(2,i,j+1,k)
      v(it,jt+1,kt) = Q(3,i,j+1,k)
      if (it == blockDim%x) then
        u(it+1,jt+1,kt) = Q(2,i+1,j+1,k)
        v(it+1,jt+1,kt) = Q(3,i+1,j+1,k)
      endif
    endif
    ! k+1 ##########################################
    if (kt == blockDim%z) then
      if (jt == 1) then
        u(it,jt-1,kt+1) = Q(2,i,j-1,k+1)
        if (it == blockDim%x) then
          u(it+1,jt-1,kt+1) = Q(2,i+1,j-1,k+1)
        endif
      endif
      u(it,jt,kt+1) = Q(2,i,j,k+1)
      w(it,kt+1,jt) = Q(4,i,j,k+1)
      if (it == blockDim%x) then
        u(it+1,jt,kt+1) = Q(2,i+1,j,k+1)
        w(it+1,kt+1,jt) = Q(4,i+1,j,k+1)
      endif
      if (jt == blockDim%y) then
        u(it,jt+1,kt+1) = Q(2,i,j+1,k+1)
        if (it == blockDim%x) then
          u(it+1,jt+1,kt+1) = Q(2,i+1,j+1,k+1)
        endif
      endif
    endif
    call syncthreads()
  end subroutine store_shared_x


  !$dir inline
  attributes(device) subroutine store_shared_y(nx, ny, nz, i, j, k, it, jt, kt, Q, u, v, w)
    integer, intent(in), value  :: nx, ny, nz, i, j, k, it, jt, kt
    real(8), intent(in), device :: Q(5,nx,ny,nz)
    real(8), intent(inout)      :: u(threadsFv%y+1,0:threadsFv%x+1,threadsFv%z)
    real(8), intent(inout)      :: v(threadsFv%y+1,0:threadsFv%x+1,0:threadsFv%z+1)
    real(8), intent(inout)      :: w(threadsFv%y+1,0:threadsFv%z+1,threadsFv%x)
    ! k-1 ##########################################
    if (kt == 1) then
      if (it == 1) then
        v(jt,it-1,kt-1) = Q(3,i-1,j,k-1)
      endif
      v(jt,it,kt-1) = Q(3,i,j,k-1)
      w(jt,kt-1,it) = Q(4,i,j,k-1)
      if (it == blockDim%x) then
        v(jt,it+1,kt-1) = Q(3,i+1,j,k-1)
      endif
      if (jt == blockDim%y) then
        if (it == 1) then
          v(jt+1,it-1,kt-1) = Q(3,i-1,j+1,k-1)
        endif
        v(jt+1,it,kt-1) = Q(3,i,j+1,k-1)
        w(jt+1,kt-1,it) = Q(4,i,j+1,k-1)
          if (it == blockDim%x) then
        v(jt+1,it+1,kt-1) = Q(3,i+1,j+1,k-1)
        endif
      endif
    endif
    ! k ############################################
    if (it == 1) then
      u(jt,it-1,kt) = Q(2,i-1,j,k)
      v(jt,it-1,kt) = Q(3,i-1,j,k)
    endif
    u(jt,it,kt) = Q(2,i,j,k)
    v(jt,it,kt) = Q(3,i,j,k)
    w(jt,kt,it) = Q(4,i,j,k)
    if (it == blockDim%x) then
      u(jt,it+1,kt) = Q(2,i+1,j,k)
      v(jt,it+1,kt) = Q(3,i+1,j,k)
    endif
    if (jt == blockDim%y) then
      if (it == 1) then
        u(jt+1,it-1,kt) = Q(2,i-1,j+1,k)
        v(jt+1,it-1,kt) = Q(3,i-1,j+1,k)
      endif
      u(jt+1,it,kt) = Q(2,i,j+1,k)
      v(jt+1,it,kt) = Q(3,i,j+1,k)
      w(jt+1,kt,it) = Q(4,i,j+1,k)
      if (it == blockDim%x) then
        u(jt+1,it+1,kt) = Q(2,i+1,j+1,k)
        v(jt+1,it+1,kt) = Q(3,i+1,j+1,k)
      endif
    endif
    ! k+1 ##########################################
    if (kt == blockDim%z) then
      if (it == 1) then
        v(jt,it-1,kt+1) = Q(3,i-1,j,k+1)
      endif
      v(jt,it,kt+1) = Q(3,i,j,k+1)
      w(jt,kt+1,it) = Q(4,i,j,k+1)
      if (it == blockDim%x) then
        v(jt,it+1,kt+1) = Q(3,i+1,j,k+1)
      endif
      if (jt == blockDim%y) then
        if (it == 1) then
          v(jt+1,it-1,kt+1) = Q(3,i-1,j+1,k+1)
        endif
        v(jt+1,it,kt+1) = Q(3,i,j+1,k+1)
        w(jt+1,kt+1,it) = Q(4,i,j+1,k+1)
        if (it == blockDim%x) then
          v(jt+1,it+1,kt+1) = Q(3,i+1,j+1,k+1)
        endif
      endif
    endif
    call syncthreads()
  end subroutine store_shared_y


  !$dir inline
  attributes(device) subroutine store_shared_z(nx, ny, nz, i, j, k, it, jt, kt, Q, u, v, w)
    integer, intent(in), value  :: nx, ny, nz, i, j, k, it, jt, kt
    real(8), intent(in), device :: Q(5,nx,ny,nz)
    real(8), intent(inout)      :: u(threadsGv%z+1,0:threadsGv%x+1,threadsGv%y)
    real(8), intent(inout)      :: v(threadsGv%z+1,0:threadsGv%y+1,threadsGv%x)
    real(8), intent(inout)      :: w(threadsGv%z+1,0:threadsGv%x+1,0:threadsGv%y+1)
    ! k ############################################
    if (jt == 1) then
      if (it == 1) then
        w(kt,it-1,jt-1) = Q(4,i-1,j-1,k)
      endif
      v(kt,jt-1,it) = Q(3,i,j-1,k)
      w(kt,it,jt-1) = Q(4,i,j-1,k)
      if (it == blockDim%x) then
        w(kt,it+1,jt-1) = Q(4,i+1,j-1,k)
      endif
    endif
    if (it == 1) then
      u(kt,it-1,jt) = Q(2,i-1,j,k)
      w(kt,it-1,jt) = Q(4,i-1,j,k)
    endif
    u(kt,it,jt) = Q(2,i,j,k)
    v(kt,jt,it) = Q(3,i,j,k)
    w(kt,it,jt) = Q(4,i,j,k)
    if (it == blockDim%x) then
      u(kt,it+1,jt) = Q(2,i+1,j,k)
      w(kt,it+1,jt) = Q(4,i+1,j,k)
    endif
    if (jt == blockDim%y) then
      if (it == 1) then
        w(kt,it-1,jt+1) = Q(4,i-1,j+1,k)
      endif
      v(kt,jt+1,it) = Q(3,i,j+1,k)
      w(kt,it,jt+1) = Q(4,i,j+1,k)
      if (it == blockDim%x) then
        w(kt,it+1,jt+1) = Q(4,i+1,j+1,k)
      endif
    endif
    ! k+1 ##########################################
    if (kt == blockDim%z) then
      if (jt == 1) then
        if (it == 1) then
          w(kt+1,it-1,jt-1) = Q(4,i-1,j-1,k+1)
        endif
        v(kt+1,jt-1,it) = Q(3,i,j-1,k+1)
        w(kt+1,it,jt-1) = Q(4,i,j-1,k+1)
        if (it == blockDim%x) then
          w(kt+1,it+1,jt-1) = Q(4,i+1,j-1,k+1)
        endif
      endif
      if (it == 1) then
        u(kt+1,it-1,jt) = Q(2,i-1,j,k+1)
        w(kt+1,it-1,jt) = Q(4,i-1,j,k+1)
      endif
      u(kt+1,it,jt) = Q(2,i,j,k+1)
      v(kt+1,jt,it) = Q(3,i,j,k+1)
      w(kt+1,it,jt) = Q(4,i,j,k+1)
      if (it == blockDim%x) then
        u(kt+1,it+1,jt) = Q(2,i+1,j,k+1)
        w(kt+1,it+1,jt) = Q(4,i+1,j,k+1)
      endif
      if (jt == blockDim%y) then
        if (it == 1) then
          w(kt+1,it-1,jt+1) = Q(4,i-1,j+1,k+1)
        endif
        v(kt+1,jt+1,it) = Q(3,i,j+1,k+1)
        w(kt+1,it,jt+1) = Q(4,i,j+1,k+1)
        if (it == blockDim%x) then
          w(kt+1,it+1,jt+1) = Q(4,i+1,j+1,k+1)
        endif
      endif
    endif 
    call syncthreads()
  end subroutine store_shared_z


  attributes(global) subroutine calc_Ev2(nx, ny, nz, dx, dy, dz, Q, E, seed)
    use calc_sutherland, only : mu6, mu2, mu_23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), shared :: u(threadsEv%x+1,0:threadsEv%y+1,0:threadsEv%z+1)
    real(8), shared :: v(threadsEv%x+1,0:threadsEv%y+1,threadsEv%z)
    real(8), shared :: w(threadsEv%x+1,0:threadsEv%z+1,threadsEv%y)
    integer i, j, k, it, jt, kt
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    call store_shared_x(nx, ny, nz, i, j, k, it, jt, kt, Q, u, v, w)

    block
      real(8), device :: Tx(2)
      Tx(:) = Q(5,i:i+1,j,k) / (R * Q(1,i:i+1,j,k))
      mx    = mu2(Tx(:))
      kTx   = Cp_over_Pr * mx * (-Tx(1) + Tx(2)) * dx(i)
    end block
    block
      real(8) my1, my2
      real(8), device :: Ty(2,3)
      Ty(:,:) = Q(5,i:i+1,j-1:j+1,k) / (R * Q(1,i:i+1,j-1:j+1,k))
      call mu_23(Ty(:,:), my1, my2)
      !muy = 0.25d0 * (my1 * (-u(it,jt-1,kt) + u(it,jt,kt) - u(it+1,jt-1,kt) + u(it+1,jt,kt)) &
      !              + my2 * (-u(it,jt,kt) + u(it,jt+1,kt) - u(it+1,jt,kt) + u(it+1,jt+1,kt))) * dy(j)
      !mvy = 0.25d0 * (my1 * (-v(it,jt-1,kt) + v(it,jt,kt) - v(it+1,jt-1,kt) + v(it+1,jt,kt)) &
      !              + my2 * (-v(it,jt,kt) + v(it,jt+1,kt) - v(it+1,jt,kt) + v(it+1,jt+1,kt))) * dy(j)
      muy = 0.25d0 * (my1 * (-u(it,jt-1,kt) - u(it+1,jt-1,kt)) + (my1 - my2) * (u(it,jt,kt) + u(it+1,jt,kt)) &
                    + my2 * ( u(it,jt+1,kt) + u(it+1,jt+1,kt))) * dy(j)
      mvy = 0.25d0 * (my1 * (-v(it,jt-1,kt) - v(it+1,jt-1,kt)) + (my1 - my2) * (v(it,jt,kt) + v(it+1,jt,kt)) &
                    + my2 * ( v(it,jt+1,kt) + v(it+1,jt+1,kt))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      real(8), device :: Tz(2,3)
      Tz(:,:) = Q(5,i:i+1,j,k-1:k+1) / (R * Q(1,i:i+1,j,k-1:k+1))
      call mu_23(Tz(:,:), mz1, mz2)
      !muz = 0.25d0 * (mz1 * (-u(it,jt,kt-1) + u(it,jt,kt) - u(it+1,jt,kt-1) + u(it+1,jt,kt)) &
      !              + mz2 * (-u(it,jt,kt) + u(it,jt,kt+1) - u(it+1,jt,kt) + u(it+1,jt,kt+1))) * dz(k)
      !mwz = 0.25d0 * (mz1 * (-w(it,kt-1,jt) + w(it,kt,jt) - w(it+1,kt-1,jt) + w(it+1,kt,jt)) &
      !              + mz2 * (-w(it,kt,jt) + w(it,kt+1,jt) - w(it+1,kt,jt) + w(it+1,kt+1,jt))) * dz(k)
      muz = 0.25d0 * (mz1 * (-u(it,jt,kt-1) - u(it+1,jt,kt-1)) + (mz1 - mz2) * (u(it,jt,kt) + u(it+1,jt,kt)) &
                    + mz2 * ( u(it,jt,kt+1) + u(it+1,jt,kt+1))) * dz(k)
      mwz = 0.25d0 * (mz1 * (-w(it,kt-1,jt) - w(it+1,kt-1,jt)) + (mz1 - mz2) * (w(it,kt,jt) + w(it+1,kt,jt)) &
                    + mz2 * ( w(it,kt+1,jt) + w(it+1,kt+1,jt))) * dz(k)
    end block
    mux = mx * (-u(it,jt,kt) + u(it+1,jt,kt)) * dx(i)
    mvx = mx * (-v(it,jt,kt) + v(it+1,jt,kt)) * dx(i)
    mwx = mx * (-w(it,kt,jt) + w(it+1,kt,jt)) * dx(i)
    txx  = 2.d0 * (2.d0 * mux - mvy - mwz) * one_third
    txy  = muy + mvx
    txz  = mwx + muz
    if (present(seed)) then
      block
        real(8) std_t, std_q, over_V, Zq, T1, T2
        real(8), device    :: rand(4), Z(6), Zx(6)
        real(8), parameter :: kb_over_dt = 1.380649d-23 / dt
        T1  = Q(5,i,j,k)   / (R * Q(1,i,j,k))
        T2  = Q(5,i+1,j,k) / (R * Q(1,i+1,j,k))
        over_V = dx(i) * dy(j) * dz(k)
        std_t  = sqrt(kb_over_dt * over_V * mx * (T1 + T2))
        std_q  = sqrt(kb_over_dt * over_V * mx * Cp_over_Pr * (T1**2 + T2**2))
        Z   = Z_tilde(seed(i,j,k))
        Zx  = Z_tilde(seed(i+1,j,k))
        Z   = 0.5d0 * (Z + Zx)
        Zq  = Zq_x(seed(i,j,k))
        txx = txx + std_t * (2.d0 * Z(1) - Z(4) - Z(6)) * one_third
        txy = txy + std_t * Z(2)
        txz = txz + std_t * Z(3)
        kTx = kTx + std_q * Zq
      end block
    endif
    utxx = 0.5d0 * (u(it,jt,kt) + u(it+1,jt,kt)) * txx
    vtxy = 0.5d0 * (v(it,jt,kt) + v(it+1,jt,kt)) * txy
    wtxz = 0.5d0 * (w(it,kt,jt) + w(it+1,kt,jt)) * txz
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev2
 

  attributes(global) subroutine calc_Ev_LES2(nx, ny, nz, dx, dy, dz, Q, mut, qc2, E)
    use calc_sutherland, only : mu6, mu2, mu23
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx, Hsgs
    real(8), dimension(2), device :: my, mysgs, mz, mzsgs
    real(8) mx, mxsgs, mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs, muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    ! SGS
    mx    = 0.5d0 * (mut(i,j,k) + mut(i+1,j,k))
    my(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k)), &
              0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))/)
    mz(:) = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k)), &
              0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))/)
    block
      real(8), device :: Tx(2)
      Tx(:) = Q(5,i:i+1,j,k) / (R * Q(1,i:i+1,j,k))
      mx    = mu2(Tx(:))
      kTx   = Cp_over_Pr * mx * (-Tx(1) + Tx(2)) * dx(i)
    end block
    block
      real(8), device :: Ty(2,3)
      Ty(:,:) = Q(5,i:i+1,j-1:j+1,k) / (R * Q(1,i:i+1,j-1:j+1,k))
      my(:)   = mu23(Ty(:,:))
      muy    = 0.25d0 * (my(1)    * (-Q(2,i,j-1,k) + Q(2,i,j,k) - Q(2,i+1,j-1,k) + Q(2,i+1,j,k)) &
                       + my(2)    * (-Q(2,i,j,k) + Q(2,i,j+1,k) - Q(2,i+1,j,k) + Q(2,i+1,j+1,k))) * dy(j)
      muysgs = 0.25d0 * (mysgs(1) * (-Q(2,i,j-1,k) + Q(2,i,j,k) - Q(2,i+1,j-1,k) + Q(2,i+1,j,k)) &
                       + mysgs(2) * (-Q(2,i,j,k) + Q(2,i,j+1,k) - Q(2,i+1,j,k) + Q(2,i+1,j+1,k))) * dy(j)
      mvy    = 0.25d0 * (my(1)    * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i+1,j-1,k) + Q(3,i+1,j,k)) &
                       + my(2)    * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i+1,j,k) + Q(3,i+1,j+1,k))) * dy(j)
      mvysgs = 0.25d0 * (mysgs(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i+1,j-1,k) + Q(3,i+1,j,k)) &
                       + mysgs(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i+1,j,k) + Q(3,i+1,j+1,k))) * dy(j)
    end block
    block
      real(8), device :: Tz(2,3)
      Tz(:,:) = Q(5,i:i+1,j,k-1:k+1) / (R * Q(1,i:i+1,j,k-1:k+1))
      mz(:)   = mu23(Tz(:,:))
      muz    = 0.25d0 * (mz(1)    * (-Q(2,i,j,k-1) + Q(2,i,j,k) - Q(2,i+1,j,k-1) + Q(2,i+1,j,k)) &
                       + mz(2)    * (-Q(2,i,j,k) + Q(2,i,j,k+1) - Q(2,i+1,j,k) + Q(2,i+1,j,k+1))) * dz(k)
      muzsgs = 0.25d0 * (mzsgs(1) * (-Q(2,i,j,k-1) + Q(2,i,j,k) - Q(2,i+1,j,k-1) + Q(2,i+1,j,k)) &
                       + mzsgs(2) * (-Q(2,i,j,k) + Q(2,i,j,k+1) - Q(2,i+1,j,k) + Q(2,i+1,j,k+1))) * dz(k)
      mwz    = 0.25d0 * (mz(1)    * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i+1,j,k-1) + Q(4,i+1,j,k)) &
                       + mz(2)    * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i+1,j,k) + Q(4,i+1,j,k+1))) * dz(k)
      mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i+1,j,k-1) + Q(4,i+1,j,k)) &
                       + mzsgs(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i+1,j,k) + Q(4,i+1,j,k+1))) * dz(k)
    end block
    mux    = mx    * (-Q(2,i,j,k) + Q(2,i+1,j,k)) * dx(i)
    muxsgs = mxsgs * (-Q(2,i,j,k) + Q(2,i+1,j,k)) * dx(i)
    mvx    = mx    * (-Q(3,i,j,k) + Q(3,i+1,j,k)) * dx(i)
    mvxsgs = mxsgs * (-Q(3,i,j,k) + Q(3,i+1,j,k)) * dx(i)
    mwx    = mx    * (-Q(4,i,j,k) + Q(4,i+1,j,k)) * dx(i)
    mwxsgs = mxsgs * (-Q(4,i,j,k) + Q(4,i+1,j,k)) * dx(i)
    txx  = 2.d0 * (2.d0 * mux - mvy - mwz) * one_third
    txy  = muy + mvx
    txz  = mwx + muz
    utxx = 0.5d0 * (Q(2,i,j,k) + Q(2,i+1,j,k)) * txx
    vtxy = 0.5d0 * (Q(3,i,j,k) + Q(3,i+1,j,k)) * txy
    wtxz = 0.5d0 * (Q(4,i,j,k) + Q(4,i+1,j,k)) * txz
    txx = txx + 2.d0 * (2.d0 * muxsgs - mvysgs - mwzsgs) * one_third
    txy = txy + muysgs + mvxsgs
    txz = txz + mwxsgs + muzsgs
    block
      real(8), device :: H(2)
      H(:) = (gamma * Q(5,i:i+1,j,k) / (Q(1,i:i+1,j,k) * gamma_1)) &
             + 0.5d0 * (Q(2,i:i+1,j,k)**2 + Q(3,i:i+1,j,k)**2 + Q(4,i:i+1,j,k)**2) + qc2(i:i+1,j,k)
      Hsgs = -mx * (-H(1) + H(2)) * dx(i) / Prt
    end block
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx + Hsgs)
  end subroutine calc_Ev_LES2
 

  attributes(global) subroutine calc_Fv2(nx, ny, nz, dy, dx, dz, Q, F, seed)
    use calc_sutherland, only : mu6, mu2, mu_23, mu_32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), shared :: u(threadsFv%y+1,0:threadsFv%x+1,threadsFv%z)
    real(8), shared :: v(threadsFv%y+1,0:threadsFv%x+1,0:threadsFv%z+1)
    real(8), shared :: w(threadsFv%y+1,0:threadsFv%z+1,threadsFv%x)
    integer i, j, k, it, jt, kt
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    call store_shared_y(nx, ny, nz, i, j, k, it, jt, kt, Q, u, v, w)

    block
      real(8) mx1, mx2
      real(8), device :: Tx(3,2)
      Tx(:,:) = Q(5,i-1:i+1,j:j+1,k) / (R * Q(1,i-1:i+1,j:j+1,k))
      call mu_23(Tx(:,:), mx1, mx2)
      !mux = 0.25d0 * (mx1 * (-u(jt,it-1,kt) + u(jt,it,kt) - u(jt+1,it-1,kt) + u(jt+1,it,kt)) &
      !              + mx2 * (-u(jt,it,kt) + u(jt,it+1,kt) - u(jt+1,it,kt) + u(jt+1,it+1,kt))) * dx(i)
      !mvx = 0.25d0 * (mx1 * (-v(jt,it-1,kt) + v(jt,it,kt) - v(jt+1,it-1,kt) + v(jt+1,it,kt)) &
      !              + mx2 * (-v(jt,it,kt) + v(jt,it+1,kt) - v(jt+1,it,kt) + v(jt+1,it+1,kt))) * dx(i)
      mux = 0.25d0 * (mx1 * (-u(jt,it-1,kt) - u(jt+1,it-1,kt)) + (mx1 - mx2) * (u(jt,it,kt) + u(jt+1,it,kt)) &
                    + mx2 * ( u(jt,it+1,kt) + u(jt+1,it+1,kt))) * dx(i)
      mvx = 0.25d0 * (mx1 * (-v(jt,it-1,kt) - v(jt+1,it-1,kt)) + (mx1 - mx2) * (v(jt,it,kt) + v(jt+1,it,kt)) &
                    + mx2 * ( v(jt,it+1,kt) + v(jt+1,it+1,kt))) * dx(i)
    end block
    block
      real(8), device :: Ty(2)
      Ty(:) = Q(5,i,j:j+1,k) / (R * Q(1,i,j:j+1,k))
      my    = mu2(Ty(:))
      kTy   = Cp_over_Pr * my * (-Ty(1) + Ty(2)) * dy(j)
    end block
    block
      real(8) mz1, mz2
      real(8), device :: Tz(2,3)
      Tz(:,:) = Q(5,i,j:j+1,k-1:k+1) / (R * Q(1,i,j:j+1,k-1:k+1))
      call mu_32(Tz(:,:), mz1, mz2)
      !mvz = 0.25d0 * (mz1 * (-v(jt,it,kt-1) + v(jt,it,kt) - v(jt+1,it,kt-1) + v(jt+1,it,kt)) &
      !              + mz2 * (-v(jt,it,kt) + v(jt,it,kt+1) - v(jt+1,it,kt) + v(jt+1,it,kt+1))) * dz(k)
      !mwz = 0.25d0 * (mz1 * (-w(jt,kt-1,it) + w(jt,kt,it) - w(jt+1,kt-1,it) + w(jt+1,kt,it)) &
      !              + mz2 * (-w(jt,kt,it) + w(jt,kt+1,it) - w(jt+1,kt,it) + w(jt+1,kt+1,it))) * dz(k)
      mvz = 0.25d0 * (mz1 * (-v(jt,it,kt-1) - v(jt+1,it,kt-1)) + (mz1 - mz2) * (v(jt,it,kt) + v(jt+1,it,kt)) &
                    + mz2 * ( v(jt,it,kt+1) + v(jt+1,it,kt+1))) * dz(k)
      mwz = 0.25d0 * (mz1 * (-w(jt,kt-1,it) - w(jt+1,kt-1,it)) + (mz1 - mz2) * (w(jt,kt,it) + w(jt+1,kt,it)) &
                    + mz2 * ( w(jt,kt+1,it) + w(jt+1,kt+1,it))) * dz(k)
    end block
    muy = my * (-u(jt,it,kt) + u(jt+1,it,kt)) * dy(j)
    mvy = my * (-v(jt,it,kt) + v(jt+1,it,kt)) * dy(j)
    mwy = my * (-w(jt,kt,it) + w(jt+1,kt,it)) * dy(j)
    tyx  = muy + mvx
    tyy  = 2.d0 * (2.d0 * mvy - mwz - mux) * one_third
    tyz  = mvz + mwy
    if (present(seed)) then
      block
        real(8) std_t, std_q, over_V, Zq, T1, T2
        real(8), device    :: rand(4), Z(6), Zy(6)
        real(8), parameter :: kb_over_dt = 1.380649d-23 / dt
        T1  = Q(5,i,j,k) / (R * Q(1,i,j,k))
        T2  = Q(5,i,j+1,k) / (R * Q(1,i,j+1,k))
        over_V = dx(i) * dy(j) * dz(k)
        std_t  = sqrt(kb_over_dt * over_V * my * (T1 + T2))
        std_q  = sqrt(kb_over_dt * over_V * my * Cp_over_Pr * (T1**2 + T2**2))
        Z   = Z_tilde(seed(i,j,k))
        Zy  = Z_tilde(seed(i,j+1,k))
        Z   = 0.5d0 * (Z + Zy)
        Zq  = Zq_y(seed(i,j,k))
        tyx = tyx + std_t * Z(2)
        tyy = tyy + std_t * (2.d0 * Z(4) - Z(6) - Z(1)) * one_third
        tyz = tyz + std_t * Z(5)
        kTy = kTy + std_q * Zq
      end block
    endif
    utyx = 0.5d0 * (u(jt,it,kt) + u(jt+1,it,kt)) * tyx
    vtyy = 0.5d0 * (v(jt,it,kt) + v(jt+1,it,kt)) * tyy
    wtyz = 0.5d0 * (w(jt,kt,it) + w(jt+1,kt,it)) * tyz
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv2
 

  attributes(global) subroutine calc_Fv_LES2(nx, ny, nz, dy, dx, dz, Q, mut, qc2, F)
    use calc_sutherland, only : mu6, mu2, mu23, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy, Hsgs
    real(8), dimension(2), device :: u2, v2, w2, mz, mzsgs, mx, mxsgs
    real(8) my, mysgs, muy, muysgs, mvy, mvysgs, mwy, mwysgs, mvz, mvzsgs, mwz, mwzsgs, mux, muxsgs, mvx, mvxsgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    ! SGS
    my     = 0.5d0 * (mut(i,j,k) + mut(i,j+1,k))
    mz(:)  = (/0.25d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k)), &
               0.25d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))/)
    mx(:)  = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k)), &
               0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))/)
    block
      real(8), device :: Tx(3,2)
      Tx(:,:) = Q(5,i-1:i+1,j:j+1,k) / (R * Q(1,i-1:i+1,j:j+1,k))
      mx(:)   = mu23(Tx(:,:))
      mux    = 0.25d0 * (mx(1)    * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                       + mx(2)    * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
      muxsgs = 0.25d0 * (mxsgs(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                       + mxsgs(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
      mvx    = 0.25d0 * (mx(1)    * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                       + mx(2)    * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
      mvxsgs = 0.25d0 * (mxsgs(1) * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                       + mxsgs(2) * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
    end block
    block
      real(8), device :: Ty(2)
      Ty(:) = Q(5,i,j:j+1,k) / (R * Q(1,i,j:j+1,k))
      my    = mu2(Ty(:))
      kTy   = Cp_over_Pr * my * (-Ty(1) + Ty(2)) * dy(j)
    end block
    block
      real(8), device :: Tz(2,3)
      Tz(:,:) = Q(5,i,j:j+1,k-1:k+1) / (R * Q(1,i,j:j+1,k-1:k+1))
      mz(:)   = mu32(Tz(:,:))
      mvz    = 0.25d0 * (mz(1)    * (-Q(3,i,j,k-1) + Q(3,i,j,k) - Q(3,i,j+1,k-1) + Q(3,i,j+1,k)) &
                       + mz(2)    * (-Q(3,i,j,k) + Q(3,i,j,k+1) - Q(3,i,j+1,k) + Q(3,i,j+1,k+1))) * dz(k)
      mvzsgs = 0.25d0 * (mzsgs(1) * (-Q(3,i,j,k-1) + Q(3,i,j,k) - Q(3,i,j+1,k-1) + Q(3,i,j+1,k)) &
                       + mzsgs(2) * (-Q(3,i,j,k) + Q(3,i,j,k+1) - Q(3,i,j+1,k) + Q(3,i,j+1,k+1))) * dz(k)
      mwz    = 0.25d0 * (mz(1)    * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i,j+1,k-1) + Q(4,i,j+1,k)) &
                       + mz(2)    * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i,j+1,k) + Q(4,i,j+1,k+1))) * dz(k)
      mwzsgs = 0.25d0 * (mzsgs(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i,j+1,k-1) + Q(4,i,j+1,k)) &
                       + mzsgs(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i,j+1,k) + Q(4,i,j+1,k+1))) * dz(k)
    end block
    muy    = my    * (-Q(2,i,j,k) + Q(2,i,j+1,k)) * dy(j)
    muysgs = mysgs * (-Q(2,i,j,k) + Q(2,i,j+1,k)) * dy(j)
    mvy    = my    * (-Q(3,i,j,k) + Q(3,i,j+1,k)) * dy(j)
    mvysgs = mysgs * (-Q(3,i,j,k) + Q(3,i,j+1,k)) * dy(j)
    mwy    = my    * (-Q(4,i,j,k) + Q(4,i,j+1,k)) * dy(j)
    mwysgs = mysgs * (-Q(4,i,j,k) + Q(4,i,j+1,k)) * dy(j)
    tyx    = muy + mvx
    tyy    = 2.d0 * (2.d0 * mvy - mwz - mux) * one_third
    tyz    = mvz + mwy
    utyx   = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j+1,k)) * tyx
    vtyy   = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j+1,k)) * tyy
    wtyz   = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j+1,k)) * tyz
    tyx    = tyx + muysgs + mvxsgs
    tyy    = tyy + 2.d0 * (2.d0 * mvysgs - mwzsgs - muxsgs) * one_third
    tyz    = tyz + mvzsgs + mwysgs
    block
      real(8), device :: H(2)
      H(:) = (gamma * Q(5,i,j:j+1,k) / (Q(1,i,j:j+1,k) * gamma_1)) &
             + 0.5d0 * (u2(:)**2 + v2(:)**2 + w2(:)**2) + qc2(i,j:j+1,k)
      Hsgs = -my * (-H(1) + H(2)) * dy(j) / Prt
    end block
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy + Hsgs)
  end subroutine calc_Fv_LES2
 

  attributes(global) subroutine calc_Gv2(nx, ny, nz, dx, dy, dz, Q, G, seed)
    use calc_sutherland, only : mu6, mu2, mu_32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), shared :: u(threadsGv%z+1,0:threadsGv%x+1,threadsGv%y)
    real(8), shared :: v(threadsGv%z+1,0:threadsGv%y+1,threadsGv%x)
    real(8), shared :: w(threadsGv%z+1,0:threadsGv%x+1,0:threadsGv%y+1)
    integer i, j, k, it, jt, kt
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    call store_shared_z(nx, ny, nz, i, j, k, it, jt, kt, Q, u, v, w)

    block
      real(8) mx1, mx2
      real(8), device :: Tx(3,2)
      Tx(:,:) = Q(5,i-1:i+1,j,k:k+1) / (R * Q(1,i-1:i+1,j,k:k+1))
      call mu_32(Tx(:,:), mx1, mx2)
      !mux = 0.25d0 * (mx1 * (-u(kt,it-1,jt) + u(kt,it,jt) - u(kt+1,it-1,jt) + u(kt+1,it,jt)) &
      !              + mx2 * (-u(kt,it,jt) + u(kt,it+1,jt) - u(kt+1,it,jt) + u(kt+1,it+1,jt))) * dx(i)
      !mwx = 0.25d0 * (mx1 * (-w(kt,jt,it-1) + w(kt,jt,it) - w(kt+1,jt,it-1) + w(kt+1,jt,it)) &
      !              + mx2 * (-w(kt,jt,it) + w(kt,jt,it+1) - w(kt+1,jt,it) + w(kt+1,jt,it+1))) * dx(i)
      mux = 0.25d0 * (mx1 * (-u(kt,it-1,jt) - u(kt+1,it-1,jt)) + (mx1 - mx2) * (u(kt,it,jt) + u(kt+1,it,jt)) &
                    + mx2 * ( u(kt,it+1,jt) + u(kt+1,it+1,jt))) * dx(i)
      mwx = 0.25d0 * (mx1 * (-w(kt,it-1,jt) - w(kt+1,it-1,jt)) + (mx1 - mx2) * (w(kt,it,jt) + w(kt+1,it,jt)) &
                    + mx2 * ( w(kt,it+1,jt) + w(kt+1,it+1,jt))) * dx(i)
    end block
    block
      real(8) my1, my2
      real(8), device :: Ty(3,2)
      Ty(:,:) = Q(5,i,j-1:j+1,k:k+1) / (R * Q(1,i,j-1:j+1,k:k+1))
      call mu_32(Ty(:,:), my1, my2)
      !mvy = 0.25d0 * (my1 * (-v(kt,jt-1,it) + v(kt,jt,it) - v(kt+1,jt-1,it) + v(kt+1,jt,it)) &
      !              + my2 * (-v(kt,jt,it) + v(kt,jt+1,it) - v(kt+1,jt,it) + v(kt+1,jt+1,it))) * dy(j)
      !mwy    = 0.25d0 * (my1    * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
      !                 + my2    * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
      mvy = 0.25d0 * (my1 * (-v(kt,jt-1,it) - v(kt+1,jt-1,it)) + (my1 - my2) * (v(kt,jt,it) + v(kt+1,jt,it)) &
                    + my2 * ( v(kt,jt+1,it) + v(kt+1,jt+1,it))) * dy(j)
      mwy = 0.25d0 * (my1 * (-w(kt,it,jt-1) - w(kt+1,it,jt-1)) + (my1 - my2) * (w(kt,it,jt) + w(kt+1,it,jt)) &
                    + my2 * ( w(kt,it,jt+1) + w(kt+1,it,jt+1))) * dy(j)
    end block
    block
      real(8), device :: Tz(2)
      Tz(:)   = Q(5,i,j,k:k+1) / (R * Q(1,i,j,k:k+1))
      mz      = mu2(Tz(:))
      kTz     = Cp_over_Pr * mz * (-Tz(1) + Tz(2)) * dz(k)
    end block
    muz = mz * (-u(kt,it,jt) + u(kt+1,it,jt)) * dz(k)
    mvz = mz * (-v(kt,jt,it) + v(kt+1,jt,it)) * dz(k)
    mwz = mz * (-w(kt,it,jt) + w(kt+1,it,jt)) * dz(k)
    tzx  = mwx + muz
    tzy  = mvz + mwy
    tzz  = 2.d0 * (2.d0 * mwz - mux - mvy) * one_third
    if (present(seed)) then
      block
        real(8) std_t, std_q, over_V, Zq, T1, T2
        real(8), device    :: rand(4), Z(6), Zz(6)
        real(8), parameter :: kb_over_dt = 1.380649d-23 / dt
        T1  = Q(5,i,j,k) / (R * Q(1,i,j,k))
        T2  = Q(5,i,j,k+1) / (R * Q(1,i,j,k+1))
        over_V = dx(i) * dy(j) * dz(k)
        std_t  = sqrt(kb_over_dt * over_V * mz * (T1 + T2))
        std_q  = sqrt(kb_over_dt * over_V * mz * Cp_over_Pr * (T1**2 + T2**2))
        Z   = Z_tilde(seed(i,j,k))
        Zz  = Z_tilde(seed(i,j,k+1))
        Z   = 0.5d0 * (Z + Zz)
        Zq  = Zq_y(seed(i,j,k))
        tzx = tzx + std_t * Z(3)
        tzy = tzy + std_t * Z(5)
        tzz = tzz + std_t * (2.d0 * Z(6) - Z(1) - Z(4)) * one_third
        kTz = kTz + std_q * Zq
      end block
    endif
    utzx = 0.5d0 * (u(kt,it,jt) + u(kt+1,it,jt)) * tzx
    vtzy = 0.5d0 * (v(kt,jt,it) + v(kt+1,jt,it)) * tzy
    wtzz = 0.5d0 * (w(kt,it,jt) + w(kt+1,it,jt)) * tzz
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv2


  attributes(global) subroutine calc_Gv_LES2(nx, ny, nz, dx, dy, dz, Q, mut, qc2, G)
    use calc_sutherland, only : mu6, mu2, mu32
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz)
    real(8), intent(in), device    :: mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz, Hsgs
    real(8), dimension(2), device :: mx, mxsgs, my, mysgs
    real(8) mz, mzsgs, muz, muzsgs, mvz, mvzsgs, mwz, mwzsgs, mwx, mwxsgs, mux, muxsgs, mvy, mvysgs, mwy, mwysgs
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x + 1
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y + 1
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    ! SGS
    mz    = 0.5d0 * (mut(i,j,k) + mut(i,j,k+1))
    mx(:) = (/0.25d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1)), &
              0.25d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))/)
    my(:) = (/0.25d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1)), &
              0.25d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))/)
    block
      real(8), device :: Tx(3,2)
      Tx(:,:) = Q(5,i-1:i+1,j,k:k+1) / (R * Q(1,i-1:i+1,j,k:k+1))
      mx(:)   = mu32(Tx(:,:))
      mux    = 0.25d0 * (mx(1)    * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j,k+1) + Q(2,i,j,k+1)) &
                       + mx(2)    * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j,k+1) + Q(2,i+1,j,k+1))) * dx(i)
      muxsgs = 0.25d0 * (mxsgs(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j,k+1) + Q(2,i,j,k+1)) &
                       + mxsgs(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j,k+1) + Q(2,i+1,j,k+1))) * dx(i)
      mwx    = 0.25d0 * (mx(1)    * (-Q(4,i-1,j,k) + Q(4,i,j,k) - Q(4,i-1,j,k+1) + Q(4,i,j,k+1)) &
                       + mx(2)    * (-Q(4,i,j,k) + Q(4,i+1,j,k) - Q(4,i,j,k+1) + Q(4,i+1,j,k+1))) * dx(i)
      mwxsgs = 0.25d0 * (mxsgs(1) * (-Q(4,i-1,j,k) + Q(4,i,j,k) - Q(4,i-1,j,k+1) + Q(4,i,j,k+1)) &
                       + mxsgs(2) * (-Q(4,i,j,k) + Q(4,i+1,j,k) - Q(4,i,j,k+1) + Q(4,i+1,j,k+1))) * dx(i)
    end block
    block
      real(8), device :: Ty(3,2)
      Ty(:,:) = Q(5,i,j-1:j+1,k:k+1) / (R * Q(1,i,j-1:j+1,k:k+1))
      my(:)   = mu32(Ty(:,:))
      mvy    = 0.25d0 * (my(1)    * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                       + my(2)    * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
      mvysgs = 0.25d0 * (mysgs(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                       + mysgs(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
      mwy    = 0.25d0 * (my(1)    * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                       + my(2)    * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
      mwysgs = 0.25d0 * (mysgs(1) * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                       + mysgs(2) * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
    end block
    block
      real(8), device :: Tz(2)
      Tz(:)   = Q(5,i,j,k:k+1) / (R * Q(1,i,j,k:k+1))
      mz      = mu2(Tz(:))
      kTz     = Cp_over_Pr * mz * (-Tz(1) + Tz(2)) * dz(k)
    end block
    muz    = mz    * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
    muzsgs = mzsgs * (-Q(2,i,j,k) + Q(2,i,j,k+1)) * dz(k)
    mvz    = mz    * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
    mvzsgs = mzsgs * (-Q(3,i,j,k) + Q(3,i,j,k+1)) * dz(k)
    mwz    = mz    * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
    mwzsgs = mzsgs * (-Q(4,i,j,k) + Q(4,i,j,k+1)) * dz(k)
    tzx    = mwx + muz
    tzy    = mvz + mwy
    tzz    = 2.d0 * (2.d0 * mwz - mux - mvy) * one_third
    utzx = 0.5d0 * (Q(2,i,j,k) + Q(2,i,j,k+1)) * tzx
    vtzy = 0.5d0 * (Q(3,i,j,k) + Q(3,i,j,k+1)) * tzy
    wtzz = 0.5d0 * (Q(4,i,j,k) + Q(4,i,j,k+1)) * tzz
    tzx    = tzx + mwx + muz
    tzy    = tzy + mvz + mwy
    tzz    = tzz + 2.d0 * (2.d0 * mwz - mux - mvy) * one_third
    block
      real(8), device :: H(2)
      H(:) = (gamma * Q(5,i,j,k:k+1) / (Q(1,i,j,k:k+1) * gamma_1)) &
             + 0.5d0 * (Q(2,i,j,k:k+1)**2 + Q(3,i,j,k:k+1)**2 + Q(4,i,j,k:k+1)**2) + qc2(i,j,k:k+1)
      Hsgs = -mz * (-H(1) + H(2)) * dz(k) / Prt
    end block
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz + Hsgs)
  end subroutine calc_Gv_LES2
end module calc_visc2

