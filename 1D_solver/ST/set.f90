module set
  use cudafor
  use mod_globals, only : dx, rho0, p0, rho1, p1
  use mod_constant, only : over_gamma_1
  implicit none
contains
  subroutine set_grid(myrank, nx, x)
    integer, intent(in)  :: myrank, nx
    real(8), intent(out) :: x(nx)
    integer i
    do i = 1, nx
      x(i) = dble(i-1) * dx
    enddo
  end subroutine set_grid


  subroutine set_init(myrank, nx, x, Q)
    integer, intent(in)  :: myrank, nx
    real(8), intent(in)  :: x(nx)
    real(8), intent(out) :: Q(nx,3)
    integer i
    do i = 1, nx/2
      Q(i,1) = rho0
      Q(i,2) = 0.d0
      Q(i,3) = p0 * over_gamma_1
    enddo
    do i = nx/2+1, nx
      Q(i,1) = rho1
      Q(i,2) = 0.d0
      Q(i,3) = p1 * over_gamma_1
    enddo
  end subroutine set_init


  !> Zero-gradient (flat extension) boundary, `ng` cells deep at each end.
  !>
  !> `ng` is max(ORDER, VISC_ORDER)/2, passed in from calc_time_dev rather than
  !> derived here so that this file stays out of the fypp pipeline. Filling the
  !> full stencil depth is what lets the flux kernels run one uniform interior
  !> stencil with no order-degrading branch near the edge -- see
  !> 1D_solver/CLAUDE.md. Cells 1..ng and nx-ng+1..nx are ghosts; the physical
  !> domain is ng+1..nx-ng.
  !>
  !> The sources (Q(ng+1), Q(nx-ng)) are physical cells, so they are never among
  !> the destinations and no staging buffer is needed.
  subroutine set_bc(myrank, nx, ng, Q)
    integer, intent(in), value     :: myrank, nx, ng
    real(8), intent(inout), device :: Q(nx,3)
    integer l, k
    !$cuf kernel do(2) <<<*,*>>>
    do l = 1, 3
      do k = 1, ng
        Q(k,l)      = Q(ng+1,l)
        Q(nx-k+1,l) = Q(nx-ng,l)
      enddo
    enddo
  end subroutine set_bc
end module set
