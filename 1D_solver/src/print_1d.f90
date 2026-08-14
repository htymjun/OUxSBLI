!> ASCII output for the 1D solver -- writes Q.dat, read directly by plot.gnu.
module print_1d
  use mod_globals, only : gamma, R, rho0, Tlr, a, Lx
  implicit none
contains
  subroutine write_Q_dat(nx, x, Q)
    integer, intent(in) :: nx
    real(8), intent(in) :: x(nx), Q(nx,3)
    real(8) rho, u, p
    integer n
    open(10, file="Q.dat", status="replace", action="write")
    do n = 1, nx
      rho = Q(n,1)
      u   = Q(n,2) / rho
      p   = (gamma - 1.d0) * (Q(n,3) - 0.5d0 * rho * u**2)
      ! es17, not es16: a negative value fills es16.8e3 exactly, leaving no
      ! separator, which silently merges it with the preceding column.
      write(10,"(4es17.8e3)") x(n) / Lx, rho / rho0, u / a, p / (rho0 * R * Tlr)
    enddo
    close(10)
  end subroutine write_Q_dat
end module print_1d
