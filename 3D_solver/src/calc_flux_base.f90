module calc_flux_base
  use mod_globals, only : id_accuracy, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG, &
  & blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv
  use calc_physical_quantities
  use calc_hybrid
  use calc_flux
  use calc_visc2
  use calc_visc4
  use calc_les
  use set
  implicit none
contains
  subroutine calc_EFG_Euler(nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, E, F, G)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(8), dimension(nx,ny,nz), device :: sensor
    integer stat
    print *, "calc quantities"
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, ruvwp)
    print *, "calc Ducros"
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    print *, "calc conv"
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, G)
    stat = cudaDeviceSynchronize()
    print *, "finish calc EFG"
  end subroutine calc_EFG_Euler

  
  subroutine calc_EFG_NS(nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, E, F, G, seed)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), dimension(nx,ny,nz), device :: sensor
    integer stat
    print *, "calc quantities"
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T, mu)
    print *, "calc Ducros"
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    print *, "calc conv"
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, G)
    stat = cudaDeviceSynchronize()
    print *, "calc visc"
    if (present(seed)) then
      if (id_visc == 2) then
        call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E, seed)
        call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F, seed)
        call calc_Gv4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G, seed)
      else
        call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E, seed)
        call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F, seed)
        call calc_Gv2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G, seed)
      endif
      call update_seed(nx, ny, nz, seed)
    else
      if (id_visc == 2) then
        call calc_Ev4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E)
        call calc_Fv4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F)
        call calc_Gv4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G)
      else
        call calc_Ev2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, E)
        call calc_Fv2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, F)
        call calc_Gv2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, G)
      endif
    endif
    print *, "finish calc EFG"
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_NS

  
  subroutine calc_EFG_LES(nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    real(8), dimension(nx,ny,nz), device :: sensor
    integer stat
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T, mu)
    call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, sensor, G)
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx, ny, nz, mut, qc2)
    if (id_visc == 2) then
      call calc_Ev_LES4<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, E)
      call calc_Fv_LES4<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, mut, qc2, F)
      call calc_Gv_LES4<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, G)
    else
      call calc_Ev_LES2<<<blocksEv,threadsEv,1>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, E)
      call calc_Fv_LES2<<<blocksFv,threadsFv,2>>>(nx, ny, nz, dy, dx, dz, ruvwp, T, mu, mut, qc2, F)
      call calc_Gv_LES2<<<blocksGv,threadsGv,3>>>(nx, ny, nz, dx, dy, dz, ruvwp, T, mu, mut, qc2, G)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_LES
end module calc_flux_base

