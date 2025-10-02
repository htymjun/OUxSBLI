module mod_constant
  use cudafor
  use mod_globals, only : R, gamma, Pr
  implicit none
  real(8), parameter :: gamma_1      = gamma - 1.d0
  real(8), parameter :: over_gamma_1 = 1.d0 / (gamma - 1.d0)
  real(8), parameter :: Cp           = gamma * R * over_gamma_1
  real(8), parameter :: Cp_over_Pr   = Cp / Pr
  real(8), parameter :: over_T0      = 1.d0 / 273.2d0
  real(8), parameter :: mu0_T0_S_over_T0_2_3 = 1.716d-5 * 384.2d0 * 273.2d0**(-1.5d0)
  real(8), parameter :: one_third    = 1.d0 / 3.d0
  real(8), parameter :: one_sixth    = 1.d0 / 6.d0
  real(8), parameter :: one_twelfth  = 1.d0 / 12.d0
  real(8), parameter :: one_sixty    = 1.d0 / 60.d0
  real(8), constant  :: Normal_x(5)  = (/0.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
  real(8), constant  :: Normal_y(5)  = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0/)
  real(8), constant  :: Normal_z(5)  = (/0.d0, 0.d0, 0.d0, 1.d0, 0.d0/)
end module mod_constant

