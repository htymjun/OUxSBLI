module set_init_common
  use mod_globals, only : gamma, R
  use mod_constant, only : Cp, gamma_1, over_gamma_1
  implicit none
contains
  subroutine set_init_tbl(nx, ny, nz, xs, ys, zs, blt0, blt, rf, u0, p0, T0, M0, Q)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: xs(nx), ys(ny), zs(nz)
    real(8), intent(in)  :: blt0, blt, rf, u0, p0, T0, M0
    real(8), intent(out) :: Q(5,nx,ny,nz)
    integer i, j, k
    real(8) :: eta, rho, u, v, w, T, Tw, Taw, p_wall
    ! random
    real(8) :: std, ustd, vstd, wstd, Tstd
    real(8), allocatable :: randum(:,:,:,:)
    integer ir, jr, kr, nxr, nyr, nzr
    ! generate randum
    nxr = (nx+9) / 10
    nyr = (ny+1) / 2
    nzr = (nz+1) / 2
    allocate(randum(4,nxr,nyr,nzr))
    do k = 1, nzr
      do j = 1, nyr
        do i = 1, nxr
          call random_number(randum(1,i,j,k))
          call random_number(randum(2,i,j,k))
          call random_number(randum(3,i,j,k))
          call random_number(randum(4,i,j,k))
          randum(1,i,j,k) = 2.d0 * randum(1,i,j,k) - 1.d0
          randum(2,i,j,k) = 2.d0 * randum(2,i,j,k) - 1.d0
          randum(3,i,j,k) = 2.d0 * randum(3,i,j,k) - 1.d0
          randum(4,i,j,k) = 2.d0 * randum(4,i,j,k) - 1.d0
    enddo;enddo;enddo
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          eta = 5.d0 * ys(j) / blt0
          u   = min(u0, u0 * (0.0015d0 * eta**4 - 0.0181d0 * eta**3 + 0.029d0 * eta**2 + 0.3192 * eta + 0.0003d0))
          v   = 0.d0
          Taw = T0 * (1.d0 + rf * 0.5d0 * gamma_1 * M0**2)
          Tw  = Taw
          T   = Tw + (Taw - Tw) * u / u0 - rf * u**2 / (2.d0 * Cp)
          !call calc_Blasius(eta,d,u,v)
          if (10 < j .and. ys(j) <= blt) then
            ir = i / 10 + 1
            jr = j / 2  + 1
            kr = k / 2  + 1
            ustd = 0.2d0 * u0 * randum(1,ir,jr,kr)
            vstd = 0.1d0 * u0 * randum(2,ir,jr,kr)
            wstd = 0.1d0 * u0 * randum(3,ir,jr,kr)
            Tstd = T0 * gamma_1 * M0**2 * randum(4,ir,jr,kr) * 0.1d0
          else
            ustd = 0.d0
            vstd = 0.d0
            wstd = 0.d0
            Tstd = 0.d0
          endif
          u   = u + ustd
          v   = v + vstd
          w   = wstd
          T   = T + Tstd
          rho = p0 / (R * T)
          Q(1,i,j,k) = rho
          Q(2,i,j,k) = Q(1,i,j,k) * u
          Q(3,i,j,k) = Q(1,i,j,k) * v
          Q(4,i,j,k) = Q(1,i,j,k) * w
          Q(5,i,j,k) = p0 * over_gamma_1 + 0.5d0 * (Q(2,i,j,k)**2 + Q(3,i,j,k)**2 + Q(4,i,j,k)**2) / Q(1,i,j,k)
    enddo;enddo;enddo
    deallocate(randum)
    ! bottom
    Q(1,:,1,:) = Q(1,:,2,:)
    Q(2,:,1,:) = 0.d0
    Q(3,:,1,:) = 0.d0
    Q(4,:,1,:) = 0.d0
    p_wall = gamma_1 * (Q(5,2,2,2) - 0.5d0 * (Q(2,2,2,2)**2 + Q(3,2,2,2)**2 + Q(4,2,2,2)**2) / Q(1,2,2,2))
    Q(5,:,1,:) = p_wall * over_gamma_1
  end subroutine set_init_tbl
end module set_init_common

