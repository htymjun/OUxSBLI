module calc_flux_base
  use mod_globals, only : blocksE,  blocksF,  blocksG,  threadsE,  threadsF,  threadsG, blocks, threads, &
                          blocksEv, blocksFv, blocksGv, threadsEv, threadsFv, threadsGv, id_accuracy
  use calc_physical_quantities
  use calc_hybrid
  use calc_flux
  use calc_visc2
  use calc_visc4
  use calc_visc
  use calc_les
  use set
  implicit none
  interface calc_EFG
    module procedure calc_EFG_Euler, calc_EFG_visc, calc_EFG_LES
  end interface calc_EFG
contains
  subroutine calc_EFG_Euler(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, over_Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, Rv, seed)
    integer(kind=2), intent(in), value   :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny), over_Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(1,1,1), mut(1,1,1), qc2(1,1,1)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), Rv(1,1,1,1)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), dimension(nx,ny,nz), device   :: sensor
    integer stat
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T)
    if (kind(id_tvd) == 8) then
      call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    else
      sensor = 0.d0
    endif
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_Euler

  
  subroutine calc_EFG_visc(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, over_Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, Rv, seed)
    integer(kind=4), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny), over_Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(nx,ny,nz), mut(1,1,1), qc2(1,1,1)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), Rv(4,nx-2,ny-2,nz-2)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), dimension(nx,ny,nz), device :: sensor
    integer stat
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T, mu)
    if (kind(id_tvd) == 8) then
      call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    else
      sensor = 0.d0
    endif
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, G)
    stat = cudaDeviceSynchronize()
    if (id_visc == 1) then
      call calc_visc_Laplacian2nd<<<blocks,threads,4>>>(nx, ny, nz, dx, dy, dz, over_Jacobian, ruvwp, T, mu, Rv)
    elseif (id_visc == 2) then
      call calc_visc_Laplacian4th<<<blocks,threads,4>>>(nx, ny, nz, dx, dy, dz, over_Jacobian, ruvwp, T, mu, Rv)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_visc

  
  subroutine calc_EFG_LES(id_visc, nx, ny, nz, dx, dy, dz, Jacobian, over_Jacobian, QJ, ruvwp, T, mu, mut, qc2, E, F, G, Rv, seed)
    integer(kind=8), intent(in), value :: id_visc
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny), over_Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz), mu(nx,ny,nz), mut(nx,ny,nz), qc2(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), Rv(4,nx-2,ny-2,nz-2)
    integer(8), intent(inout), device, optional :: seed(nx,ny,nz)
    real(8), dimension(nx,ny,nz), device :: sensor
    integer stat
    mut = 0.d0
    qc2 = 0.d0
    call calc_quantities_T_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T, mu)
    if (kind(id_tvd) == 8) then
      call calc_Ducros<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, sensor)
    else
      sensor = 0.d0
    endif
    call calc_E<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, E)
    call calc_F<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, F)
    call calc_G<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, sensor, G)
    call calc_mut<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, ruvwp, mut, qc2)
    stat = cudaDeviceSynchronize()
    call set_bc_mut(nx, ny, nz, mut, qc2)
    if (id_visc == 1) then
      call calc_visc_Laplacian2nd<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, over_Jacobian, ruvwp, T, mu, Rv)
    else
      call calc_visc_Laplacian4th<<<blocks,threads>>>(nx, ny, nz, dx, dy, dz, over_Jacobian, ruvwp, T, mu, Rv)
    endif
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG_LES
end module calc_flux_base

