module calc_visc2_curv
  use mod_globals, only : gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third
  use load_smem_visc2_curv
  implicit none
  private
  public calc_Ev2_curv, calc_Ev_LES2_curv, calc_Fv2_curv, calc_Fv_LES2_curv, calc_Gv2_curv, calc_Gv_LES2_curv
contains
  attributes(global) subroutine calc_Ev2_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, n_xi_x, n_xi_y, Q_2, Q_3, Q_4, T, mu, E)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), value                 :: dz
    real(8), intent(in), device, contiguous    :: xi_x(nx,ny)
    real(8), intent(in), device, contiguous    :: xi_y(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_x(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_y(nx,ny)
    real(8), intent(in), device, contiguous    :: n_xi_x(nx-1,ny-2)
    real(8), intent(in), device, contiguous    :: n_xi_y(nx-1,ny-2)
    real(8), intent(in), device, contiguous    :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(nx-1,ny-2,nz-2,5)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, txx, txy, txz
    real(8) mux, mvx, mwx, muy, mvy, muz, mwz, mwy, mvz
    real(8) nxx, nxy, S, xi_x_f, xi_y_f, eta_x_f, eta_y_f, mu_f, dT_dxi
    real(8) du_dxi, dv_dxi, dw_dxi, my1, my2, mz1, mz2
    real(8) mu_eta_u, mu_eta_v, mu_eta_w, dT_deta, muz_v, dTdx, dTdy
    real(8) txx_p, txy_p, tyy_p, txz_p, tyz_p
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc2_curv_x(it, jt, kt, j, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_xi_x(i, j-1)
    nxy = n_xi_y(i, j-1)
    S   = sqrt(nxx*nxx + nxy*nxy)
    
    ! Block A: Metrics at ξ-face
    xi_x_f  = 0.5d0 * (xi_x(i,j) + xi_x(i+1,j))
    xi_y_f  = 0.5d0 * (xi_y(i,j) + xi_y(i+1,j))
    eta_x_f = 0.5d0 * (eta_x(i,j) + eta_x(i+1,j))
    eta_y_f = 0.5d0 * (eta_y(i,j) + eta_y(i+1,j))
    mu_f    = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
    du_dxi  = u(idx+1) - u(idx)
    dv_dxi  = v(idx+1) - v(idx)
    dw_dxi  = w(idx+1) - w(idx)
    dT_dxi  = T(i+1,j,k) - T(i,j,k)
    
    ! Block B: η-tangential stencil (4-corner mu, unit Δη)
    my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,k) + mu(i+1,j-1,k) + mu(i+1,j,k))
    my2 = 0.0625d0 * (mu(i,j,k) + mu(i,j+1,k) + mu(i+1,j,k) + mu(i+1,j+1,k))
    mu_eta_u = my1*(-Q_2(i,j-1,k)-Q_2(i+1,j-1,k)) + (my1-my2)*(u(idx)+u(idx+1)) &
             + my2*(Q_2(i,j+1,k)+Q_2(i+1,j+1,k))
    mu_eta_v = my1*(-Q_3(i,j-1,k)-Q_3(i+1,j-1,k)) + (my1-my2)*(v(idx)+v(idx+1)) &
             + my2*(Q_3(i,j+1,k)+Q_3(i+1,j+1,k))
    mu_eta_w = my1*(-Q_4(i,j-1,k)-Q_4(i+1,j-1,k)) + (my1-my2)*(w(idx)+w(idx+1)) &
             + my2*(Q_4(i,j+1,k)+Q_4(i+1,j+1,k))
    dT_deta = 0.25d0 * ((T(i,j+1,k)+T(i+1,j+1,k)) - (T(i,j-1,k)+T(i+1,j-1,k)))
    
    ! Block C: z-tangential stencil (4-corner mu × 1/dz)
    mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k) + mu(i+1,j,k-1) + mu(i+1,j,k))
    mz2 = 0.0625d0 * (mu(i,j,k) + mu(i,j,k+1) + mu(i+1,j,k) + mu(i+1,j,k+1))
    muz = (mz1*(-Q_2(i,j,k-1)-Q_2(i+1,j,k-1)) + (mz1-mz2)*(u(idx)+u(idx+1)) &
         + mz2*(Q_2(i,j,k+1)+Q_2(i+1,j,k+1))) / dz
    mwz = (mz1*(-Q_4(i,j,k-1)-Q_4(i+1,j,k-1)) + (mz1-mz2)*(w(idx)+w(idx+1)) &
         + mz2*(Q_4(i,j,k+1)+Q_4(i+1,j,k+1))) / dz
    muz_v = (mz1*(-Q_3(i,j,k-1)-Q_3(i+1,j,k-1)) + (mz1-mz2)*(v(idx)+v(idx+1)) &
           + mz2*(Q_3(i,j,k+1)+Q_3(i+1,j,k+1))) / dz
    
    ! Block D: Physical gradients via chain rule
    mux = mu_f * du_dxi * xi_x_f + mu_eta_u * eta_x_f
    muy = mu_f * du_dxi * xi_y_f + mu_eta_u * eta_y_f
    mvx = mu_f * dv_dxi * xi_x_f + mu_eta_v * eta_x_f
    mvy = mu_f * dv_dxi * xi_y_f + mu_eta_v * eta_y_f
    mwx = mu_f * dw_dxi * xi_x_f + mu_eta_w * eta_x_f
    mwy = mu_f * dw_dxi * xi_y_f + mu_eta_w * eta_y_f
    
    ! Block E: Stress tensor projection and output
    txx_p = two_third * (2.d0*mux - mvy - mwz)
    txy_p = muy + mvx
    tyy_p = two_third * (2.d0*mvy - mux - mwz)
    txz_p = mwx + muz
    tyz_p = mwy + muz_v
    
    txx = (txx_p*nxx + txy_p*nxy) / S
    txy = (txy_p*nxx + tyy_p*nxy) / S
    txz = (txz_p*nxx + tyz_p*nxy) / S
    
    dTdx = xi_x_f*dT_dxi + eta_x_f*dT_deta
    dTdy = xi_y_f*dT_dxi + eta_y_f*dT_deta
    viscous_work = Cp_over_Pr*mu_f*(dTdx*nxx+dTdy*nxy)/S &
                 + 0.5d0*((u(idx)+u(idx+1))*txx+(v(idx)+v(idx+1))*txy+(w(idx)+w(idx+1))*txz)
    
    E(i,j-1,k-1,2) = E(i,j-1,k-1,2) - txx * S
    E(i,j-1,k-1,3) = E(i,j-1,k-1,3) - txy * S
    E(i,j-1,k-1,4) = E(i,j-1,k-1,4) - txz * S
    E(i,j-1,k-1,5) = E(i,j-1,k-1,5) - viscous_work * S
  end subroutine calc_Ev2_curv


  attributes(global) subroutine calc_Fv2_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, n_eta_x, n_eta_y, Q_2, Q_3, Q_4, T, mu, F)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), value                 :: dz
    real(8), intent(in), device, contiguous    :: xi_x(nx,ny)
    real(8), intent(in), device, contiguous    :: xi_y(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_x(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_y(nx,ny)
    real(8), intent(in), device, contiguous    :: n_eta_x(nx-2,ny-1)
    real(8), intent(in), device, contiguous    :: n_eta_y(nx-2,ny-1)
    real(8), intent(in), device, contiguous    :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(nx-2,ny-1,nz-2,5)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, tyx, tyy, tyz
    real(8) mux, mvx, mwx, muy, mvy, muz, mvz, mwz, mwx_temp, mwy, muz_w
    real(8) nex, ney, S, xi_x_f, xi_y_f, eta_x_f, eta_y_f, mu_f, du_deta
    real(8) dv_deta, dw_deta, dT_deta, mx1, mx2, mz1, mz2
    real(8) mu_xi_u, mu_xi_v, mu_xi_w, dT_dxi, dTdx, dTdy
    real(8) txx_p, txy_p, tyy_p, txz_p, tyz_p
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc2_curv_y(it, jt, kt, i, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    
    nex = n_eta_x(i-1, j)
    ney = n_eta_y(i-1, j)
    S   = sqrt(nex*nex + ney*ney)
    
    ! Metrics at η-face
    xi_x_f  = 0.5d0 * (xi_x(i,j) + xi_x(i,j+1))
    xi_y_f  = 0.5d0 * (xi_y(i,j) + xi_y(i,j+1))
    eta_x_f = 0.5d0 * (eta_x(i,j) + eta_x(i,j+1))
    eta_y_f = 0.5d0 * (eta_y(i,j) + eta_y(i,j+1))
    mu_f    = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
    du_deta = u(idx+1) - u(idx)
    dv_deta = v(idx+1) - v(idx)
    dw_deta = w(idx+1) - w(idx)
    dT_deta = T(i,j+1,k) - T(i,j,k)
    
    ! ξ-cross stencil (mx1, mx2 at ξ-face corners, unit Δξ)
    mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,j,k) + mu(i-1,j+1,k) + mu(i,j+1,k))
    mx2 = 0.0625d0 * (mu(i,j,k) + mu(i+1,j,k) + mu(i,j+1,k) + mu(i+1,j+1,k))
    mu_xi_u = mx1*(-Q_2(i-1,j,k)-Q_2(i-1,j+1,k)) + (mx1-mx2)*(u(idx)+u(idx+1)) &
            + mx2*(Q_2(i+1,j,k)+Q_2(i+1,j+1,k))
    mu_xi_v = mx1*(-Q_3(i-1,j,k)-Q_3(i-1,j+1,k)) + (mx1-mx2)*(v(idx)+v(idx+1)) &
            + mx2*(Q_3(i+1,j,k)+Q_3(i+1,j+1,k))
    mu_xi_w = mx1*(-Q_4(i-1,j,k)-Q_4(i-1,j+1,k)) + (mx1-mx2)*(w(idx)+w(idx+1)) &
            + mx2*(Q_4(i+1,j,k)+Q_4(i+1,j+1,k))
    dT_dxi = 0.25d0 * ((T(i+1,j,k)+T(i+1,j+1,k)) - (T(i-1,j,k)+T(i-1,j+1,k)))
    
    ! z-tangential stencil
    mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k) + mu(i,j+1,k-1) + mu(i,j+1,k))
    mz2 = 0.0625d0 * (mu(i,j,k) + mu(i,j,k+1) + mu(i,j+1,k) + mu(i,j+1,k+1))
    mvz = (mz1*(-Q_3(i,j,k-1)-Q_3(i,j+1,k-1)) + (mz1-mz2)*(v(idx)+v(idx+1)) &
         + mz2*(Q_3(i,j,k+1)+Q_3(i,j+1,k+1))) / dz
    muz = (mz1*(-Q_2(i,j,k-1)-Q_2(i,j+1,k-1)) + (mz1-mz2)*(u(idx)+u(idx+1)) &
         + mz2*(Q_2(i,j,k+1)+Q_2(i,j+1,k+1))) / dz
    mwz = (mz1*(-Q_4(i,j,k-1)-Q_4(i,j+1,k-1)) + (mz1-mz2)*(w(idx)+w(idx+1)) &
         + mz2*(Q_4(i,j,k+1)+Q_4(i,j+1,k+1))) / dz
    
    ! Physical gradients
    mux = mu_xi_u * xi_x_f + mu_f * du_deta * eta_x_f
    muy = mu_xi_u * xi_y_f + mu_f * du_deta * eta_y_f
    mvx = mu_xi_v * xi_x_f + mu_f * dv_deta * eta_x_f
    mvy = mu_xi_v * xi_y_f + mu_f * dv_deta * eta_y_f
    mwx_temp = mu_xi_w * xi_x_f + mu_f * dw_deta * eta_x_f
    mwy = mu_xi_w * xi_y_f + mu_f * dw_deta * eta_y_f
    
    ! Stress tensor projection (at η-face)
    txx_p = two_third * (2.d0*mux - mvy - mwz)
    txy_p = muy + mvx
    tyy_p = two_third * (2.d0*mvy - mux - mwz)
    txz_p = mwx_temp + muz
    tyz_p = mwy + mvz
    
    tyx = (txx_p*nex + txy_p*ney) / S
    tyy = (txy_p*nex + tyy_p*ney) / S
    tyz = (txz_p*nex + tyz_p*ney) / S
    
    dTdx = xi_x_f*dT_dxi + eta_x_f*dT_deta
    dTdy = xi_y_f*dT_dxi + eta_y_f*dT_deta
    viscous_work = Cp_over_Pr*mu_f*(dTdx*nex+dTdy*ney)/S &
                 + 0.5d0*((u(idx)+u(idx+1))*tyx+(v(idx)+v(idx+1))*tyy+(w(idx)+w(idx+1))*tyz)
    
    F(i-1,j,k-1,2) = F(i-1,j,k-1,2) - tyx * S
    F(i-1,j,k-1,3) = F(i-1,j,k-1,3) - tyy * S
    F(i-1,j,k-1,4) = F(i-1,j,k-1,4) - tyz * S
    F(i-1,j,k-1,5) = F(i-1,j,k-1,5) - viscous_work * S
  end subroutine calc_Fv2_curv


  attributes(global) subroutine calc_Gv2_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, Q_2, Q_3, Q_4, T, mu, G)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), value                 :: dz
    real(8), intent(in), device, contiguous    :: xi_x(nx,ny)
    real(8), intent(in), device, contiguous    :: xi_y(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_x(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_y(nx,ny)
    real(8), intent(in), device, contiguous    :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(nx-2,ny-2,nz-1,5)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, tzx, tzy, tzz
    real(8) mux, muz, mvz, mwz, mwx, mvy, mwy, muz_u, muz_v
    real(8) xi_x_f, xi_y_f, eta_x_f, eta_y_f, mu_f, mu_u_z
    real(8) mu_v_z, mu_w_z, mx1, mx2, my1, my2
    real(8) mu_xi_u, mu_xi_v, mu_xi_w, mu_eta_u, mu_eta_v, mu_eta_w, dTdz
    real(8) txx_p, txy_p, tyy_p, txz_p, tyz_p
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc2_curv_z(it, jt, kt, i, j, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    
    ! Metrics at z-face (cell-center values since z is uniform)
    xi_x_f  = xi_x(i,j)
    xi_y_f  = xi_y(i,j)
    eta_x_f = eta_x(i,j)
    eta_y_f = eta_y(i,j)
    mu_f    = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
    mu_u_z = mu_f * (u(idx+1) - u(idx)) / dz
    mu_v_z = mu_f * (v(idx+1) - v(idx)) / dz
    mu_w_z = mu_f * (w(idx+1) - w(idx)) / dz
    
    ! ξ-cross (at z-face corners: mu at (i-1/2, j, k+1/2))
    mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,j,k) + mu(i-1,j,k+1) + mu(i,j,k+1))
    mx2 = 0.0625d0 * (mu(i,j,k) + mu(i+1,j,k) + mu(i,j,k+1) + mu(i+1,j,k+1))
    mu_xi_u = mx1*(-Q_2(i-1,j,k)-Q_2(i-1,j,k+1)) + (mx1-mx2)*(u(idx)+u(idx+1)) &
            + mx2*(Q_2(i+1,j,k)+Q_2(i+1,j,k+1))
    mu_xi_v = mx1*(-Q_3(i-1,j,k)-Q_3(i-1,j,k+1)) + (mx1-mx2)*(v(idx)+v(idx+1)) &
            + mx2*(Q_3(i+1,j,k)+Q_3(i+1,j,k+1))
    mu_xi_w = mx1*(-Q_4(i-1,j,k)-Q_4(i-1,j,k+1)) + (mx1-mx2)*(w(idx)+w(idx+1)) &
            + mx2*(Q_4(i+1,j,k)+Q_4(i+1,j,k+1))
    
    ! η-cross (at z-face corners: mu at (i, j-1/2, k+1/2))
    my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,k) + mu(i,j-1,k+1) + mu(i,j,k+1))
    my2 = 0.0625d0 * (mu(i,j,k) + mu(i,j+1,k) + mu(i,j,k+1) + mu(i,j+1,k+1))
    mu_eta_u = my1*(-Q_2(i,j-1,k)-Q_2(i,j-1,k+1)) + (my1-my2)*(u(idx)+u(idx+1)) &
             + my2*(Q_2(i,j+1,k)+Q_2(i,j+1,k+1))
    mu_eta_v = my1*(-Q_3(i,j-1,k)-Q_3(i,j-1,k+1)) + (my1-my2)*(v(idx)+v(idx+1)) &
             + my2*(Q_3(i,j+1,k)+Q_3(i,j+1,k+1))
    mu_eta_w = my1*(-Q_4(i,j-1,k)-Q_4(i,j-1,k+1)) + (my1-my2)*(w(idx)+w(idx+1)) &
             + my2*(Q_4(i,j+1,k)+Q_4(i,j+1,k+1))
    
    ! Physical gradients
    mux = mu_xi_u * xi_x_f + mu_eta_u * eta_x_f
    mvy = mu_eta_v * eta_y_f + mu_xi_v * xi_y_f
    mwx = mu_xi_w * xi_x_f + mu_eta_w * eta_x_f
    mwy = mu_xi_w * xi_y_f + mu_eta_w * eta_y_f
    
    ! Stress tensor (z-face, no projection since nz=(0,0,1))
    tzx = mwx + mu_u_z
    tzy = mwy + mu_v_z
    tzz = two_third * (2.d0*mu_w_z - mux - mvy)
    
    dTdz = T(i,j,k+1) - T(i,j,k)
    viscous_work = Cp_over_Pr * mu_f * dTdz / dz &
                 + 0.5d0 * ((u(idx)+u(idx+1))*tzx + (v(idx)+v(idx+1))*tzy + (w(idx)+w(idx+1))*tzz)
    
    G(i-1,j-1,k,2) = G(i-1,j-1,k,2) - tzx
    G(i-1,j-1,k,3) = G(i-1,j-1,k,3) - tzy
    G(i-1,j-1,k,4) = G(i-1,j-1,k,4) - tzz
    G(i-1,j-1,k,5) = G(i-1,j-1,k,5) - viscous_work
  end subroutine calc_Gv2_curv


  ! LES versions follow (calc_Ev_LES2_curv, calc_Fv_LES2_curv, calc_Gv_LES2_curv)
  ! For brevity, these follow the same pattern with SGS stress computed from mut

  attributes(global) subroutine calc_Ev_LES2_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, n_xi_x, n_xi_y, Q_2, Q_3, Q_4, T, mu, mut, qc2, E)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), value                 :: dz
    real(8), intent(in), device, contiguous    :: xi_x(nx,ny)
    real(8), intent(in), device, contiguous    :: xi_y(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_x(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_y(nx,ny)
    real(8), intent(in), device, contiguous    :: n_xi_x(nx-1,ny-2)
    real(8), intent(in), device, contiguous    :: n_xi_y(nx-1,ny-2)
    real(8), intent(in), device, contiguous    :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: E(nx-1,ny-2,nz-2,5)
    integer, parameter :: sx = threadsEv%x + 1
    integer, parameter :: sy = threadsEv%y
    integer, parameter :: sz = threadsEv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, txx, txy, txz
    real(8) mux, muxsgs, mvx, mvxsgs, mwx, mwxsgs
    real(8) muy, muysgs, mvy, mvysgs, muz, muzsgs, mwz, mwzsgs, mwy, mwysgs, mvz, mvzsgs
    real(8) mut_f_xi, muysgs_phys, mvysgs_phys, mwysgs_phys
    real(8) nxx, nxy, S, xi_x_f, xi_y_f, eta_x_f, eta_y_f, mu_f, dT_dxi
    real(8) du_dxi, dv_dxi, dw_dxi, my1, my2, my1sgs, my2sgs, mz1, mz2, mz1sgs, mz2sgs
    real(8) mu_eta_u, mu_eta_v, mu_eta_w, dT_deta, dTdx, dTdy
    real(8) txx_p, txy_p, tyy_p, txz_p, tyz_p
    real(8) H1, H2, muetsgsdy, mysgsdy_temp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (it-1) + (jt-1)*sx + (kt-1)*sx*sy
    call load_smem_visc2_curv_x(it, jt, kt, j, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    nxx = n_xi_x(i, j-1)
    nxy = n_xi_y(i, j-1)
    S   = sqrt(nxx*nxx + nxy*nxy)
    
    ! Metrics and derivatives
    xi_x_f  = 0.5d0 * (xi_x(i,j) + xi_x(i+1,j))
    xi_y_f  = 0.5d0 * (xi_y(i,j) + xi_y(i+1,j))
    eta_x_f = 0.5d0 * (eta_x(i,j) + eta_x(i+1,j))
    eta_y_f = 0.5d0 * (eta_y(i,j) + eta_y(i+1,j))
    mu_f    = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
    du_dxi  = u(idx+1) - u(idx)
    dv_dxi  = v(idx+1) - v(idx)
    dw_dxi  = w(idx+1) - w(idx)
    dT_dxi  = T(i+1,j,k) - T(i,j,k)
    
    ! SGS mu stencils
    my1sgs = 0.0625d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i+1,j-1,k) + mut(i+1,j,k))
    my2sgs = 0.0625d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i+1,j,k) + mut(i+1,j+1,k))
    mz1sgs = 0.0625d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i+1,j,k-1) + mut(i+1,j,k))
    mz2sgs = 0.0625d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i+1,j,k) + mut(i+1,j,k+1))
    
    ! Molecular viscous terms (same as NS)
    my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,k) + mu(i+1,j-1,k) + mu(i+1,j,k))
    my2 = 0.0625d0 * (mu(i,j,k) + mu(i,j+1,k) + mu(i+1,j,k) + mu(i+1,j+1,k))
    mu_eta_u = my1*(-Q_2(i,j-1,k)-Q_2(i+1,j-1,k)) + (my1-my2)*(u(idx)+u(idx+1)) &
             + my2*(Q_2(i,j+1,k)+Q_2(i+1,j+1,k))
    mu_eta_v = my1*(-Q_3(i,j-1,k)-Q_3(i+1,j-1,k)) + (my1-my2)*(v(idx)+v(idx+1)) &
             + my2*(Q_3(i,j+1,k)+Q_3(i+1,j+1,k))
    mu_eta_w = my1*(-Q_4(i,j-1,k)-Q_4(i+1,j-1,k)) + (my1-my2)*(w(idx)+w(idx+1)) &
             + my2*(Q_4(i,j+1,k)+Q_4(i+1,j+1,k))
    dT_deta = 0.25d0 * ((T(i,j+1,k)+T(i+1,j+1,k)) - (T(i,j-1,k)+T(i+1,j-1,k)))
    
    mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k) + mu(i+1,j,k-1) + mu(i+1,j,k))
    mz2 = 0.0625d0 * (mu(i,j,k) + mu(i,j,k+1) + mu(i+1,j,k) + mu(i+1,j,k+1))
    muz = (mz1*(-Q_2(i,j,k-1)-Q_2(i+1,j,k-1)) + (mz1-mz2)*(u(idx)+u(idx+1)) &
         + mz2*(Q_2(i,j,k+1)+Q_2(i+1,j,k+1))) / dz
    mwz = (mz1*(-Q_4(i,j,k-1)-Q_4(i+1,j,k-1)) + (mz1-mz2)*(w(idx)+w(idx+1)) &
         + mz2*(Q_4(i,j,k+1)+Q_4(i+1,j,k+1))) / dz
    mvz = (mz1*(-Q_3(i,j,k-1)-Q_3(i+1,j,k-1)) + (mz1-mz2)*(v(idx)+v(idx+1)) &
         + mz2*(Q_3(i,j,k+1)+Q_3(i+1,j,k+1))) / dz
    
    ! SGS viscous terms
    muysgs = my1sgs*(-Q_2(i,j-1,k)-Q_2(i+1,j-1,k)) + (my1sgs-my2sgs)*(u(idx)+u(idx+1)) &
           + my2sgs*(Q_2(i,j+1,k)+Q_2(i+1,j+1,k))
    mvysgs = my1sgs*(-Q_3(i,j-1,k)-Q_3(i+1,j-1,k)) + (my1sgs-my2sgs)*(v(idx)+v(idx+1)) &
           + my2sgs*(Q_3(i,j+1,k)+Q_3(i+1,j+1,k))
    muzsgs = (mz1sgs*(-Q_2(i,j,k-1)-Q_2(i+1,j,k-1)) + (mz1sgs-mz2sgs)*(u(idx)+u(idx+1)) &
            + mz2sgs*(Q_2(i,j,k+1)+Q_2(i+1,j,k+1))) / dz
    mvzsgs = (mz1sgs*(-Q_3(i,j,k-1)-Q_3(i+1,j,k-1)) + (mz1sgs-mz2sgs)*(v(idx)+v(idx+1)) &
            + mz2sgs*(Q_3(i,j,k+1)+Q_3(i+1,j,k+1))) / dz
    mwzsgs = (mz1sgs*(-Q_4(i,j,k-1)-Q_4(i+1,j,k-1)) + (mz1sgs-mz2sgs)*(w(idx)+w(idx+1)) &
            + mz2sgs*(Q_4(i,j,k+1)+Q_4(i+1,j,k+1))) / dz
    mwysgs = my1sgs*(-Q_4(i,j-1,k)-Q_4(i+1,j-1,k)) + (my1sgs-my2sgs)*(w(idx)+w(idx+1)) &
           + my2sgs*(Q_4(i,j+1,k)+Q_4(i+1,j+1,k))

    ! Heat flux SGS
    H1 = Cp*T(i,j,k) + 0.5d0*(u(idx)**2+v(idx)**2+w(idx)**2) + qc2(i,j,k)
    H2 = Cp*T(i+1,j,k) + 0.5d0*(u(idx+1)**2+v(idx+1)**2+w(idx+1)**2) + qc2(i+1,j,k)
    Hsgs = -0.5d0*(mut(i,j,k)+mut(i+1,j,k))*(-H1+H2)/Prt
    
    ! Physical gradients
    mux = mu_f*du_dxi*xi_x_f + mu_eta_u*eta_x_f
    muy = mu_f*du_dxi*xi_y_f + mu_eta_u*eta_y_f
    mvx = mu_f*dv_dxi*xi_x_f + mu_eta_v*eta_x_f
    mvy = mu_f*dv_dxi*xi_y_f + mu_eta_v*eta_y_f
    mwx = mu_f*dw_dxi*xi_x_f + mu_eta_w*eta_x_f
    mwy = mu_f*dw_dxi*xi_y_f + mu_eta_w*eta_y_f
    
    ! SGS gradients
    muxsgs = 0.5d0*(mut(i,j,k)+mut(i+1,j,k))*du_dxi*xi_x_f + muysgs*eta_x_f
    mvxsgs = 0.5d0*(mut(i,j,k)+mut(i+1,j,k))*dv_dxi*xi_x_f + mvysgs*eta_x_f
    mwxsgs = 0.5d0*(mut(i,j,k)+mut(i+1,j,k))*dw_dxi*xi_x_f + mwysgs*eta_x_f
    mut_f_xi    = 0.5d0*(mut(i,j,k)+mut(i+1,j,k))
    muysgs_phys = mut_f_xi*du_dxi*xi_y_f + muysgs*eta_y_f
    mvysgs_phys = mut_f_xi*dv_dxi*xi_y_f + mvysgs*eta_y_f
    mwysgs_phys = mut_f_xi*dw_dxi*xi_y_f + mwysgs*eta_y_f

    ! Stress tensor projection
    txx_p = two_third*(2.d0*mux-mvy-mwz)
    txy_p = muy + mvx
    tyy_p = two_third*(2.d0*mvy-mux-mwz)
    txz_p = mwx + muz
    tyz_p = mwy + mvz
    
    txx = (txx_p*nxx + txy_p*nxy) / S
    txy = (txy_p*nxx + tyy_p*nxy) / S
    txz = (txz_p*nxx + tyz_p*nxy) / S
    
    ! Add SGS contribution (full curvilinear projection, mirroring molecular block)
    block
      real(8) txx_sgs, txy_sgs, tyy_sgs, txz_sgs, tyz_sgs
      txx_sgs = two_third*(2.d0*muxsgs    - mvysgs_phys - mwzsgs)
      txy_sgs = muysgs_phys + mvxsgs
      tyy_sgs = two_third*(2.d0*mvysgs_phys - muxsgs    - mwzsgs)
      txz_sgs = mwxsgs + muzsgs
      tyz_sgs = mwysgs_phys + mvzsgs
      txx = txx + (txx_sgs*nxx + txy_sgs*nxy) / S
      txy = txy + (txy_sgs*nxx + tyy_sgs*nxy) / S
      txz = txz + (txz_sgs*nxx + tyz_sgs*nxy) / S
    end block
    
    dTdx = xi_x_f*dT_dxi + eta_x_f*dT_deta
    dTdy = xi_y_f*dT_dxi + eta_y_f*dT_deta
    viscous_work = Cp_over_Pr*mu_f*(dTdx*nxx+dTdy*nxy)/S &
                 + 0.5d0*((u(idx)+u(idx+1))*txx+(v(idx)+v(idx+1))*txy+(w(idx)+w(idx+1))*txz)
    
    E(i,j-1,k-1,2) = E(i,j-1,k-1,2) - txx * S
    E(i,j-1,k-1,3) = E(i,j-1,k-1,3) - txy * S
    E(i,j-1,k-1,4) = E(i,j-1,k-1,4) - txz * S
    E(i,j-1,k-1,5) = E(i,j-1,k-1,5) - (viscous_work + Hsgs) * S
  end subroutine calc_Ev_LES2_curv


  attributes(global) subroutine calc_Fv_LES2_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, n_eta_x, n_eta_y, Q_2, Q_3, Q_4, T, mu, mut, qc2, F)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), value                 :: dz
    real(8), intent(in), device, contiguous    :: xi_x(nx,ny)
    real(8), intent(in), device, contiguous    :: xi_y(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_x(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_y(nx,ny)
    real(8), intent(in), device, contiguous    :: n_eta_x(nx-2,ny-1)
    real(8), intent(in), device, contiguous    :: n_eta_y(nx-2,ny-1)
    real(8), intent(in), device, contiguous    :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: F(nx-2,ny-1,nz-2,5)
    integer, parameter :: sx = threadsFv%x
    integer, parameter :: sy = threadsFv%y + 1
    integer, parameter :: sz = threadsFv%z
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, tyx, tyy, tyz
    real(8) mux, mvx, mwx, muy, mvy, muz, mvz, mwz, mwx_temp, mwy
    real(8) muxsgs, mvxsgs, mwxsgs, muysgs, mvysgs, mwysgs, muzsgs, mvzsgs, mwzsgs
    real(8) nex, ney, S, xi_x_f, xi_y_f, eta_x_f, eta_y_f, mu_f, du_deta
    real(8) dv_deta, dw_deta, dT_deta, mx1, mx2, mx1sgs, mx2sgs, mz1, mz2, mz1sgs, mz2sgs
    real(8) mu_xi_u, mu_xi_v, mu_xi_w, dT_dxi, dTdx, dTdy
    real(8) txx_p, txy_p, tyy_p, txz_p, tyz_p
    real(8) H1, H2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    idx = (jt-1) + (it-1)*sy + (kt-1)*sy*sx
    call load_smem_visc2_curv_y(it, jt, kt, i, k, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    
    nex = n_eta_x(i-1, j)
    ney = n_eta_y(i-1, j)
    S   = sqrt(nex*nex + ney*ney)
    
    ! Metrics at η-face
    xi_x_f  = 0.5d0 * (xi_x(i,j) + xi_x(i,j+1))
    xi_y_f  = 0.5d0 * (xi_y(i,j) + xi_y(i,j+1))
    eta_x_f = 0.5d0 * (eta_x(i,j) + eta_x(i,j+1))
    eta_y_f = 0.5d0 * (eta_y(i,j) + eta_y(i,j+1))
    mu_f    = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
    du_deta = u(idx+1) - u(idx)
    dv_deta = v(idx+1) - v(idx)
    dw_deta = w(idx+1) - w(idx)
    dT_deta = T(i,j+1,k) - T(i,j,k)
    
    ! ξ-cross stencil (molecular)
    mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,j,k) + mu(i-1,j+1,k) + mu(i,j+1,k))
    mx2 = 0.0625d0 * (mu(i,j,k) + mu(i+1,j,k) + mu(i,j+1,k) + mu(i+1,j+1,k))
    mu_xi_u = mx1*(-Q_2(i-1,j,k)-Q_2(i-1,j+1,k)) + (mx1-mx2)*(u(idx)+u(idx+1)) &
            + mx2*(Q_2(i+1,j,k)+Q_2(i+1,j+1,k))
    mu_xi_v = mx1*(-Q_3(i-1,j,k)-Q_3(i-1,j+1,k)) + (mx1-mx2)*(v(idx)+v(idx+1)) &
            + mx2*(Q_3(i+1,j,k)+Q_3(i+1,j+1,k))
    mu_xi_w = mx1*(-Q_4(i-1,j,k)-Q_4(i-1,j+1,k)) + (mx1-mx2)*(w(idx)+w(idx+1)) &
            + mx2*(Q_4(i+1,j,k)+Q_4(i+1,j+1,k))
    dT_dxi = 0.25d0 * ((T(i+1,j,k)+T(i+1,j+1,k)) - (T(i-1,j,k)+T(i-1,j+1,k)))
    
    ! ξ-cross stencil (SGS)
    mx1sgs = 0.0625d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j+1,k) + mut(i,j+1,k))
    mx2sgs = 0.0625d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j+1,k) + mut(i+1,j+1,k))
    
    ! z-tangential stencil (molecular)
    mz1 = 0.0625d0 * (mu(i,j,k-1) + mu(i,j,k) + mu(i,j+1,k-1) + mu(i,j+1,k))
    mz2 = 0.0625d0 * (mu(i,j,k) + mu(i,j,k+1) + mu(i,j+1,k) + mu(i,j+1,k+1))
    mvz = (mz1*(-Q_3(i,j,k-1)-Q_3(i,j+1,k-1)) + (mz1-mz2)*(v(idx)+v(idx+1)) &
         + mz2*(Q_3(i,j,k+1)+Q_3(i,j+1,k+1))) / dz
    muz = (mz1*(-Q_2(i,j,k-1)-Q_2(i,j+1,k-1)) + (mz1-mz2)*(u(idx)+u(idx+1)) &
         + mz2*(Q_2(i,j,k+1)+Q_2(i,j+1,k+1))) / dz
    mwz = (mz1*(-Q_4(i,j,k-1)-Q_4(i,j+1,k-1)) + (mz1-mz2)*(w(idx)+w(idx+1)) &
         + mz2*(Q_4(i,j,k+1)+Q_4(i,j+1,k+1))) / dz
    
    ! z-tangential stencil (SGS)
    mz1sgs = 0.0625d0 * (mut(i,j,k-1) + mut(i,j,k) + mut(i,j+1,k-1) + mut(i,j+1,k))
    mz2sgs = 0.0625d0 * (mut(i,j,k) + mut(i,j,k+1) + mut(i,j+1,k) + mut(i,j+1,k+1))
    
    ! Physical gradients (molecular)
    mux = mu_xi_u * xi_x_f + mu_f * du_deta * eta_x_f
    muy = mu_xi_u * xi_y_f + mu_f * du_deta * eta_y_f
    mvx = mu_xi_v * xi_x_f + mu_f * dv_deta * eta_x_f
    mvy = mu_xi_v * xi_y_f + mu_f * dv_deta * eta_y_f
    mwx_temp = mu_xi_w * xi_x_f + mu_f * dw_deta * eta_x_f
    mwy = mu_xi_w * xi_y_f + mu_f * dw_deta * eta_y_f
    
    ! SGS gradients
    muysgs = (mx1sgs*(-Q_2(i-1,j,k)-Q_2(i-1,j+1,k)) + (mx1sgs-mx2sgs)*(u(idx)+u(idx+1)) &
            + mx2sgs*(Q_2(i+1,j,k)+Q_2(i+1,j+1,k)))
    mvysgs = (mx1sgs*(-Q_3(i-1,j,k)-Q_3(i-1,j+1,k)) + (mx1sgs-mx2sgs)*(v(idx)+v(idx+1)) &
            + mx2sgs*(Q_3(i+1,j,k)+Q_3(i+1,j+1,k)))
    mwysgs = (mx1sgs*(-Q_4(i-1,j,k)-Q_4(i-1,j+1,k)) + (mx1sgs-mx2sgs)*(w(idx)+w(idx+1)) &
            + mx2sgs*(Q_4(i+1,j,k)+Q_4(i+1,j+1,k)))
    muzsgs = (mz1sgs*(-Q_2(i,j,k-1)-Q_2(i,j+1,k-1)) + (mz1sgs-mz2sgs)*(u(idx)+u(idx+1)) &
            + mz2sgs*(Q_2(i,j,k+1)+Q_2(i,j+1,k+1))) / dz
    mvzsgs = (mz1sgs*(-Q_3(i,j,k-1)-Q_3(i,j+1,k-1)) + (mz1sgs-mz2sgs)*(v(idx)+v(idx+1)) &
            + mz2sgs*(Q_3(i,j,k+1)+Q_3(i,j+1,k+1))) / dz
    mwzsgs = (mz1sgs*(-Q_4(i,j,k-1)-Q_4(i,j+1,k-1)) + (mz1sgs-mz2sgs)*(w(idx)+w(idx+1)) &
            + mz2sgs*(Q_4(i,j,k+1)+Q_4(i,j+1,k+1))) / dz
    
    ! Heat flux SGS
    H1 = Cp*T(i,j,k) + 0.5d0*(u(idx)**2+v(idx)**2+w(idx)**2) + qc2(i,j,k)
    H2 = Cp*T(i,j+1,k) + 0.5d0*(u(idx+1)**2+v(idx+1)**2+w(idx+1)**2) + qc2(i,j+1,k)
    Hsgs = -0.5d0*(mut(i,j,k)+mut(i,j+1,k))*(-H1+H2)/Prt
    
    ! Stress tensor projection (at η-face, molecular)
    txx_p = two_third * (2.d0*mux - mvy - mwz)
    txy_p = muy + mvx
    tyy_p = two_third * (2.d0*mvy - mux - mwz)
    txz_p = mwx_temp + muz
    tyz_p = mwy + mvz
    
    tyx = (txx_p*nex + txy_p*ney) / S
    tyy = (txy_p*nex + tyy_p*ney) / S
    tyz = (txz_p*nex + tyz_p*ney) / S
    
    ! Add SGS contribution (full curvilinear projection, mirroring molecular block)
    block
      real(8) mut_f_eta, muxsgs_phys, mvxsgs_phys, mwxsgs_phys
      real(8) muysgs_phys, mvysgs_phys, mwysgs_phys
      real(8) txx_sgs, txy_sgs, tyy_sgs, txz_sgs, tyz_sgs
      mut_f_eta   = 0.5d0*(mut(i,j,k)+mut(i,j+1,k))
      muxsgs_phys = muysgs*xi_x_f + mut_f_eta*du_deta*eta_x_f
      mvxsgs_phys = mvysgs*xi_x_f + mut_f_eta*dv_deta*eta_x_f
      mwxsgs_phys = mwysgs*xi_x_f + mut_f_eta*dw_deta*eta_x_f
      muysgs_phys = muysgs*xi_y_f + mut_f_eta*du_deta*eta_y_f
      mvysgs_phys = mvysgs*xi_y_f + mut_f_eta*dv_deta*eta_y_f
      mwysgs_phys = mwysgs*xi_y_f + mut_f_eta*dw_deta*eta_y_f
      txx_sgs = two_third*(2.d0*muxsgs_phys - mvysgs_phys - mwzsgs)
      txy_sgs = muysgs_phys + mvxsgs_phys
      tyy_sgs = two_third*(2.d0*mvysgs_phys - muxsgs_phys - mwzsgs)
      txz_sgs = mwxsgs_phys + muzsgs
      tyz_sgs = mwysgs_phys + mvzsgs
      tyx = tyx + (txx_sgs*nex + txy_sgs*ney) / S
      tyy = tyy + (txy_sgs*nex + tyy_sgs*ney) / S
      tyz = tyz + (txz_sgs*nex + tyz_sgs*ney) / S
    end block
    
    dTdx = xi_x_f*dT_dxi + eta_x_f*dT_deta
    dTdy = xi_y_f*dT_dxi + eta_y_f*dT_deta
    viscous_work = Cp_over_Pr*mu_f*(dTdx*nex+dTdy*ney)/S &
                 + 0.5d0*((u(idx)+u(idx+1))*tyx+(v(idx)+v(idx+1))*tyy+(w(idx)+w(idx+1))*tyz)
    
    F(i-1,j,k-1,2) = F(i-1,j,k-1,2) - tyx * S
    F(i-1,j,k-1,3) = F(i-1,j,k-1,3) - tyy * S
    F(i-1,j,k-1,4) = F(i-1,j,k-1,4) - tyz * S
    F(i-1,j,k-1,5) = F(i-1,j,k-1,5) - (viscous_work + Hsgs) * S
  end subroutine calc_Fv_LES2_curv


  attributes(global) subroutine calc_Gv_LES2_curv(nx, ny, nz, dz, xi_x, xi_y, eta_x, eta_y, Q_2, Q_3, Q_4, T, mu, mut, qc2, G)
    integer, intent(in), value                 :: nx
    integer, intent(in), value                 :: ny
    integer, intent(in), value                 :: nz
    real(8), intent(in), value                 :: dz
    real(8), intent(in), device, contiguous    :: xi_x(nx,ny)
    real(8), intent(in), device, contiguous    :: xi_y(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_x(nx,ny)
    real(8), intent(in), device, contiguous    :: eta_y(nx,ny)
    real(8), intent(in), device, contiguous    :: Q_2(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_3(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: Q_4(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: T(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mu(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: mut(nx,ny,nz)
    real(8), intent(in), device, contiguous    :: qc2(nx,ny,nz)
    real(8), intent(inout), device, contiguous :: G(nx-2,ny-2,nz-1,5)
    integer, parameter :: sx = threadsGv%x
    integer, parameter :: sy = threadsGv%y
    integer, parameter :: sz = threadsGv%z + 1
    real(8), dimension(0:sx*sy*sz-1), shared :: u, v, w
    integer i, j, k, it, jt, kt, idx
    real(8) viscous_work, Hsgs, tzx, tzy, tzz
    real(8) mux, muz, mvz, mwz, mwx, mvy, mwy
    real(8) muxsgs, mu_eta_usgs, mvxsgs, mwxsgs, mvysgs, mwysgs
    real(8) xi_x_f, xi_y_f, eta_x_f, eta_y_f, mu_f, mu_u_z
    real(8) mu_v_z, mu_w_z, mx1, mx2, mx1sgs, mx2sgs, my1, my2, my1sgs, my2sgs
    real(8) mu_xi_u, mu_xi_v, mu_xi_w, mu_eta_u, mu_eta_v, mu_eta_w, dTdz
    real(8) txx_p, txy_p, tyy_p, txz_p, tyz_p
    real(8) H1, H2
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    idx = (kt-1) + (jt-1)*sz + (it-1)*sz*sy
    call load_smem_visc2_curv_z(it, jt, kt, i, j, idx, nx, ny, nz, Q_2, Q_3, Q_4, u, v, w)
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    
    ! Metrics at z-face (cell-center values since z is uniform)
    xi_x_f  = xi_x(i,j)
    xi_y_f  = xi_y(i,j)
    eta_x_f = eta_x(i,j)
    eta_y_f = eta_y(i,j)
    mu_f    = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
    mu_u_z = mu_f * (u(idx+1) - u(idx)) / dz
    mu_v_z = mu_f * (v(idx+1) - v(idx)) / dz
    mu_w_z = mu_f * (w(idx+1) - w(idx)) / dz
    
    ! ξ-cross (at z-face corners: mu at (i-1/2, j, k+1/2))
    mx1 = 0.0625d0 * (mu(i-1,j,k) + mu(i,j,k) + mu(i-1,j,k+1) + mu(i,j,k+1))
    mx2 = 0.0625d0 * (mu(i,j,k) + mu(i+1,j,k) + mu(i,j,k+1) + mu(i+1,j,k+1))
    mu_xi_u = mx1*(-Q_2(i-1,j,k)-Q_2(i-1,j,k+1)) + (mx1-mx2)*(u(idx)+u(idx+1)) &
            + mx2*(Q_2(i+1,j,k)+Q_2(i+1,j,k+1))
    mu_xi_v = mx1*(-Q_3(i-1,j,k)-Q_3(i-1,j,k+1)) + (mx1-mx2)*(v(idx)+v(idx+1)) &
            + mx2*(Q_3(i+1,j,k)+Q_3(i+1,j,k+1))
    mu_xi_w = mx1*(-Q_4(i-1,j,k)-Q_4(i-1,j,k+1)) + (mx1-mx2)*(w(idx)+w(idx+1)) &
            + mx2*(Q_4(i+1,j,k)+Q_4(i+1,j,k+1))
    
    ! ξ-cross (SGS)
    mx1sgs = 0.0625d0 * (mut(i-1,j,k) + mut(i,j,k) + mut(i-1,j,k+1) + mut(i,j,k+1))
    mx2sgs = 0.0625d0 * (mut(i,j,k) + mut(i+1,j,k) + mut(i,j,k+1) + mut(i+1,j,k+1))
    
    ! η-cross (at z-face corners: mu at (i, j-1/2, k+1/2))
    my1 = 0.0625d0 * (mu(i,j-1,k) + mu(i,j,k) + mu(i,j-1,k+1) + mu(i,j,k+1))
    my2 = 0.0625d0 * (mu(i,j,k) + mu(i,j+1,k) + mu(i,j,k+1) + mu(i,j+1,k+1))
    mu_eta_u = my1*(-Q_2(i,j-1,k)-Q_2(i,j-1,k+1)) + (my1-my2)*(u(idx)+u(idx+1)) &
             + my2*(Q_2(i,j+1,k)+Q_2(i,j+1,k+1))
    mu_eta_v = my1*(-Q_3(i,j-1,k)-Q_3(i,j-1,k+1)) + (my1-my2)*(v(idx)+v(idx+1)) &
             + my2*(Q_3(i,j+1,k)+Q_3(i,j+1,k+1))
    mu_eta_w = my1*(-Q_4(i,j-1,k)-Q_4(i,j-1,k+1)) + (my1-my2)*(w(idx)+w(idx+1)) &
             + my2*(Q_4(i,j+1,k)+Q_4(i,j+1,k+1))
    
    ! η-cross (SGS)
    my1sgs = 0.0625d0 * (mut(i,j-1,k) + mut(i,j,k) + mut(i,j-1,k+1) + mut(i,j,k+1))
    my2sgs = 0.0625d0 * (mut(i,j,k) + mut(i,j+1,k) + mut(i,j,k+1) + mut(i,j+1,k+1))
    
    ! Physical gradients (molecular)
    mux = mu_xi_u * xi_x_f + mu_eta_u * eta_x_f
    mvy = mu_eta_v * eta_y_f + mu_xi_v * xi_y_f
    mwx = mu_xi_w * xi_x_f + mu_eta_w * eta_x_f
    mwy = mu_xi_w * xi_y_f + mu_eta_w * eta_y_f
    
    ! SGS gradients (ξ and η cross-stencils)
    muxsgs = (mx1sgs*(-Q_2(i-1,j,k)-Q_2(i-1,j,k+1)) + (mx1sgs-mx2sgs)*(u(idx)+u(idx+1)) &
            + mx2sgs*(Q_2(i+1,j,k)+Q_2(i+1,j,k+1)))
    mu_eta_usgs = (my1sgs*(-Q_2(i,j-1,k)-Q_2(i,j-1,k+1)) + (my1sgs-my2sgs)*(u(idx)+u(idx+1)) &
                + my2sgs*(Q_2(i,j+1,k)+Q_2(i,j+1,k+1)))
    mvysgs = (my1sgs*(-Q_3(i,j-1,k)-Q_3(i,j-1,k+1)) + (my1sgs-my2sgs)*(v(idx)+v(idx+1)) &
            + my2sgs*(Q_3(i,j+1,k)+Q_3(i,j+1,k+1)))
    mvxsgs = (mx1sgs*(-Q_3(i-1,j,k)-Q_3(i-1,j,k+1)) + (mx1sgs-mx2sgs)*(v(idx)+v(idx+1)) &
            + mx2sgs*(Q_3(i+1,j,k)+Q_3(i+1,j,k+1)))
    mwxsgs = (mx1sgs*(-Q_4(i-1,j,k)-Q_4(i-1,j,k+1)) + (mx1sgs-mx2sgs)*(w(idx)+w(idx+1)) &
            + mx2sgs*(Q_4(i+1,j,k)+Q_4(i+1,j,k+1)))
    mwysgs = (my1sgs*(-Q_4(i,j-1,k)-Q_4(i,j-1,k+1)) + (my1sgs-my2sgs)*(w(idx)+w(idx+1)) &
            + my2sgs*(Q_4(i,j+1,k)+Q_4(i,j+1,k+1)))

    ! Heat flux SGS
    H1 = Cp*T(i,j,k) + 0.5d0*(u(idx)**2+v(idx)**2+w(idx)**2) + qc2(i,j,k)
    H2 = Cp*T(i,j,k+1) + 0.5d0*(u(idx+1)**2+v(idx+1)**2+w(idx+1)**2) + qc2(i,j,k+1)
    Hsgs = -0.5d0*(mut(i,j,k)+mut(i,j,k+1))*(H2-H1) / (dz * Prt)

    ! Stress tensor (z-face, no projection since nz=(0,0,1))
    tzx = mwx + mu_u_z
    tzy = mwy + mu_v_z
    tzz = two_third * (2.d0*mu_w_z - mux - mvy)

    ! Add SGS contribution (physical gradients via chain rule, mirroring molecular block)
    block
      real(8) mut_f_z, muzsgs, mvzsgs, mwzsgs_direct
      real(8) muxsgs_phys, mvysgs_phys, mwxsgs_phys, mwysgs_phys
      mut_f_z      = 0.5d0*(mut(i,j,k)+mut(i,j,k+1))
      muzsgs       = mut_f_z*(u(idx+1)-u(idx))/dz
      mvzsgs       = mut_f_z*(v(idx+1)-v(idx))/dz
      mwzsgs_direct= mut_f_z*(w(idx+1)-w(idx))/dz
      muxsgs_phys  = muxsgs*xi_x_f + mu_eta_usgs*eta_x_f
      mvysgs_phys  = mvxsgs*xi_y_f + mvysgs*eta_y_f
      mwxsgs_phys  = mwxsgs*xi_x_f + mwysgs*eta_x_f
      mwysgs_phys  = mwxsgs*xi_y_f + mwysgs*eta_y_f
      tzx = tzx + (mwxsgs_phys + muzsgs)
      tzy = tzy + (mwysgs_phys + mvzsgs)
      tzz = tzz + two_third*(2.d0*mwzsgs_direct - muxsgs_phys - mvysgs_phys)
    end block
    
    dTdz = T(i,j,k+1) - T(i,j,k)
    viscous_work = Cp_over_Pr * mu_f * dTdz / dz &
                 + 0.5d0 * ((u(idx)+u(idx+1))*tzx + (v(idx)+v(idx+1))*tzy + (w(idx)+w(idx+1))*tzz)
    
    G(i-1,j-1,k,2) = G(i-1,j-1,k,2) - tzx
    G(i-1,j-1,k,3) = G(i-1,j-1,k,3) - tzy
    G(i-1,j-1,k,4) = G(i-1,j-1,k,4) - tzz
    G(i-1,j-1,k,5) = G(i-1,j-1,k,5) - (viscous_work + Hsgs)
  end subroutine calc_Gv_LES2_curv

end module calc_visc2_curv
