module set
  use mod_globals, only : nx, Lx, dx, gamma, R, rho0, rho1, p0, p1
  implicit none
contains
  subroutine set_grid(nx, x)
    integer, intent(in)  :: nx
    real(8), intent(out) :: x(nx)
    integer i
    x(1) = 0.d0
    do i = 1, nx-1
      x(i+1) = x(i) + dx
    enddo
  end subroutine set_grid
  
  subroutine set_init(nx,x,Q)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: x(nx)
    real(8), intent(out) :: Q(nx,3)
    integer i
    do i = 1, nx
      if (i < int(0.5 * nx)) then
        Q(i,1) = rho0
        Q(i,2) = 0.d0
        Q(i,3) = p0 / (gamma - 1.d0)
      else
        Q(i,1) = rho1
        Q(i,2) = 0.d0
        Q(i,3) = p1 / (gamma - 1.d0)
      endif
    enddo
  end subroutine set_init

  subroutine set_bc(nx,Q)
    integer, intent(in), value     :: nx
    real(8), intent(inout), device :: Q(nx,3)
    Q(1,1)  = rho0!Q(2,1)
    Q(1,2)  = 0.d0
    Q(1,3)  = p0 / (gamma - 1.d0)!Q(2,3)
    Q(nx,1) = rho1!Q(nx-1,1)
    Q(nx,2) = 0.d0
    Q(nx,3) = p1 / (gamma - 1.d0)!Q(nx-1,3)
  end subroutine set_bc
end module set

