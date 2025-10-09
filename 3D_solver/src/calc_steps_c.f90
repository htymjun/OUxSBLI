module calc_steps_c
  use iso_c_binding
  use cudafor
  use calc_steps
  implicit none
contains
  subroutine calc_R_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, E_ptr, F_ptr, G_ptr, R_ptr) bind(c, name="calc_R_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), target :: dx_ptr, dy_ptr, dz_ptr
    real(8), intent(inout), target :: E_ptr, F_ptr, G_ptr, R_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:,:,:), device, pointer :: E, F, G, R
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(R_ptr), R, [5,nx-2,ny-2,nz-2])
    call calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R)
  end subroutine calc_R_c

  
  subroutine calc_step1_c(nx, ny, nz, coef, dx_ptr, dy_ptr, dz_ptr, E_ptr, F_ptr, G_ptr, &
                         Q_ptr, Q2_ptr) bind(c, name="calc_step1_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), value     :: coef
    real(8), intent(inout), target :: dx_ptr, dy_ptr, dz_ptr
    real(8), intent(inout), target :: E_ptr, F_ptr, G_ptr, Q_ptr, Q2_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:,:,:), device, pointer :: E, F, G, Q, Q2
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(Q_ptr), Q, [5,nx,ny,nz])
    call c_f_pointer(c_loc(Q2_ptr), Q2, [5,nx,ny,nz])
    call calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)
  end subroutine calc_step1_c
  
   
  subroutine calc_step_c(nx, ny, nz, coef1, coef2, dx_ptr, dy_ptr, dz_ptr, E_ptr, F_ptr, G_ptr, &
                         Q_ptr, Q2_ptr, Rs_ptr) bind(c, name="calc_step_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), value     :: coef1, coef2
    real(8), intent(inout), target :: dx_ptr, dy_ptr, dz_ptr
    real(8), intent(inout), target :: E_ptr, F_ptr, G_ptr, Q_ptr, Q2_ptr, Rs_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:,:,:), device, pointer :: E, F, G, Q, Q2, Rs
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(Q_ptr), Q, [5,nx,ny,nz])
    call c_f_pointer(c_loc(Q2_ptr), Q2, [5,nx,ny,nz])
    call c_f_pointer(c_loc(Rs_ptr), Rs, [5,nx-2,ny-2,nz-2])
    call calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
  end subroutine calc_step_c

  
  subroutine calc_step2_3_c(nx, ny, nz, coef1, coef2, coef3, coef4, dx_ptr, dy_ptr, dz_ptr, &
                            E_ptr, F_ptr, G_ptr, Qin_ptr, Qinout_ptr) bind(c, name="calc_step2_3_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), value     :: coef1, coef2, coef3, coef4
    real(8), intent(inout), target :: dx_ptr, dy_ptr, dz_ptr
    real(8), intent(inout), target :: E_ptr, F_ptr, G_ptr, Qin_ptr, Qinout_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:,:,:), device, pointer :: E, F, G, Qin, Qinout
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(Qin_ptr), Qin, [5,nx,ny,nz])
    call c_f_pointer(c_loc(Qinout_ptr), Qinout, [5,nx,ny,nz])
    call calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qinout)
  end subroutine calc_step2_3_c
 
   
  subroutine calc_step4_c(nx, ny, nz, dx_ptr, dy_ptr, dz_ptr, E_ptr, F_ptr, G_ptr, Rs_ptr, Q_ptr) bind(c, name="calc_step4_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), target :: dx_ptr, dy_ptr, dz_ptr
    real(8), intent(inout), target :: E_ptr, F_ptr, G_ptr, Rs_ptr, Q_ptr
    real(8), dimension(:), device, pointer       :: dx, dy, dz
    real(8), dimension(:,:,:,:), device, pointer :: E, F, G, Rs, Q
    call c_f_pointer(c_loc(dx_ptr), dx, [nx-1])
    call c_f_pointer(c_loc(dy_ptr), dy, [ny-1])
    call c_f_pointer(c_loc(dz_ptr), dz, [nz-1])
    call c_f_pointer(c_loc(E_ptr), E, [5,nx-1,ny-2,nz-2])
    call c_f_pointer(c_loc(F_ptr), F, [5,nx-2,ny-1,nz-2])
    call c_f_pointer(c_loc(G_ptr), G, [5,nx-2,ny-2,nz-1])
    call c_f_pointer(c_loc(Rs_ptr), Rs, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(Q_ptr), Q, [5,nx,ny,nz])
    call calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
  end subroutine calc_step4_c

  
  subroutine calc_error_c(nx, ny, nz, R1_ptr, R2_ptr, R1_new_ptr, R2_new_ptr, err) bind(c, name="calc_error_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), target :: R1_ptr, R2_ptr, R1_new_ptr, R2_new_ptr
    real(8), intent(out)           :: err
    real(8), dimension(:,:,:,:), device, pointer :: R1, R2, R1_new, R2_new 
    call c_f_pointer(c_loc(R1_ptr), R1, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(R2_ptr), R2, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(R1_new_ptr), R1_new, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(R2_new_ptr), R2_new, [5,nx-2,ny-2,nz-2])
    call calc_error(nx, ny, nz, R1, R2, R1_new, R2_new, err)
  end subroutine calc_error_c


  subroutine calc_Gauss_step_c(nx, ny, nz, a1, a2, &
                               R1_ptr, R2_ptr, Q_ptr, Q2_ptr) bind(c, name="calc_Gauss_step_c")
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), value     :: a1, a2
    real(8), intent(inout), target :: R1_ptr, R2_ptr, Q_ptr, Q2_ptr
    real(8), dimension(:,:,:,:), device, pointer :: R1, R2, Q, Q2
    call c_f_pointer(c_loc(R1_ptr), R1, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(R2_ptr), R2, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(Q_ptr),   Q, [5,nx,ny,nz])
    call c_f_pointer(c_loc(Q2_ptr), Q2, [5,nx,ny,nz])
    call calc_Gauss_step(nx, ny, nz, a1, a2, R1, R2, Q, Q2)
  end subroutine calc_Gauss_step_c
 

  subroutine calc_Gauss_step_Q_c(nx, ny, nz, a1, a2, &
                                 R1_ptr, R2_ptr, Q_ptr) bind(c, name="calc_Gauss_step_Q_c")
    integer, intent(in), value                                                    :: nx, ny, nz
    real(8), intent(in), value                                                    :: a1, a2
    real(8), intent(inout), target :: R1_ptr, R2_ptr, Q_ptr
    real(8), dimension(:,:,:,:), device, pointer :: R1, R2, Q
    call c_f_pointer(c_loc(R1_ptr), R1, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(R2_ptr), R2, [5,nx-2,ny-2,nz-2])
    call c_f_pointer(c_loc(Q_ptr),   Q, [5,nx,ny,nz])
    call calc_Gauss_step_Q(nx, ny, nz, a1, a2, R1, R2, Q)
  end subroutine calc_Gauss_step_Q_c
end module calc_steps_c

