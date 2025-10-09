module calc_flux_base_c
  use iso_c_binding
  use calc_flux_base
  implicit none
contains
  subroutine calc_EFG_Euler_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, Q_ptr, &
                              T_ptr, mu_ptr, mut_ptr, qc2_ptr, E_ptr, F_ptr, G_ptr) bind(c, name="calc_EFG_Euler_c")
    integer, intent(in), value                   :: nx, ny, nz
    real(8), intent(inout), target               :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target               :: QJ_ptr, Q_ptr, T_ptr, mu_ptr, mut_ptr, qc2_ptr, E_ptr, F_ptr, G_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:), device, pointer     :: Jacobian
    real(8), dimension(:,:,:), device, pointer   :: T, mu, mut, qc2
    real(8), dimension(:,:,:,:), device, pointer :: QJ, ruvwp, E, F, G
    integer(kind=2) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [nx,ny])
    call c_f_pointer(c_loc(QJ_ptr),   QJ, [5,nx,ny,nz])
    call c_f_pointet(c_loc(Q_ptr), ruvwp, [5,nx,ny,nz])
    call c_f_pointer(c_loc(T_ptr),     T, [1,1,1])
    call c_f_pointer(c_loc(mu_ptr),   mu, [1,1,1])
    call c_f_pointer(c_loc(mut_ptr), mut, [1,1,1])
    call c_f_pointer(c_loc(qc2_ptr), qc2, [1,1,1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
  end subroutine calc_EFG_Euler_c


  subroutine calc_EFG_NS_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, Q_ptr, &
                           T_ptr, mu_ptr, mut_ptr, qc2_ptr, E_ptr, F_ptr, G_ptr) bind(c, name="calc_EFG_NS_c")
    integer, intent(in), value                   :: nx, ny, nz
    real(8), intent(inout), target               :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target               :: QJ_ptr, Q_ptr, T_ptr, mu_ptr, mut_ptr, qc2_ptr, E_ptr, F_ptr, G_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:), device, pointer     :: Jacobian
    real(8), dimension(:,:,:), device, pointer   :: T, mu, mut, qc2
    real(8), dimension(:,:,:,:), device, pointer :: QJ, ruvwp, E, F, G
    integer(kind=4) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [nx,ny])
    call c_f_pointer(c_loc(QJ_ptr),   QJ, [5,nx,ny,nz])
    call c_f_pointet(c_loc(Q_ptr), ruvwp, [5,nx,ny,nz])
    call c_f_pointer(c_loc(T_ptr),     T, [nx,ny,nz])
    call c_f_pointer(c_loc(mu_ptr),   mu, [nx,ny,nz])
    call c_f_pointer(c_loc(mut_ptr), mut, [1,1,1])
    call c_f_pointer(c_loc(qc2_ptr), qc2, [1,1,1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
  end subroutine calc_EFG_NS_c


  subroutine calc_EFG_LES_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr, QJ_ptr, Q_ptr, &
                            T_ptr, mu_ptr, mut_ptr, qc2_ptr, E_ptr, F_ptr, G_ptr) bind(c, name="calc_EFG_LES_c")
    integer, intent(in), value                   :: nx, ny, nz
    real(8), intent(inout), target               :: dx_ptr, dy_ptr, dz_ptr, Jacobian_ptr
    real(8), intent(inout), target               :: QJ_ptr, Q_ptr, T_ptr, mu_ptr, mut_ptr, qc2_ptr, E_ptr, F_ptr, G_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:), device, pointer     :: Jacobian
    real(8), dimension(:,:,:), device, pointer   :: T, mu, mut, qc2
    real(8), dimension(:,:,:,:), device, pointer :: QJ, ruvwp, E, F, G
    integer(kind=8) id_visc
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(Jacobian_ptr), Jacobian, [nx,ny])
    call c_f_pointer(c_loc(QJ_ptr),   QJ, [5,nx,ny,nz])
    call c_f_pointet(c_loc(Q_ptr), ruvwp, [5,nx,ny,nz])
    call c_f_pointer(c_loc(T_ptr),     T, [nx,ny,nz])
    call c_f_pointer(c_loc(mu_ptr),   mu, [nx,ny,nz])
    call c_f_pointer(c_loc(mut_ptr), mut, [nx,ny,nz])
    call c_f_pointer(c_loc(qc2_ptr), qc2, [nx,ny,nz])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call calc_EFG(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
  end subroutine calc_EFG_LES_c
end module calc_flux_base_c

