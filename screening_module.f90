
module screening_module
  use constants             ! so you can use kind=dbl if you want
  implicit none
  private
  public :: screening_potential

contains

  ! Linear screening potential between r=0.1 fm and r=6.0 fm
  real(kind=dbl) function screening_potential(r, V_screen1, V_screen2)
    implicit none
    real(kind=dbl), intent(in) :: r          ! radius in fm
    real(kind=dbl), intent(in) :: V_screen1  ! V at r = 0.1 fm [MeV]
    real(kind=dbl), intent(in) :: V_screen2  ! V at r = 6.0 fm [MeV]
    real(kind=dbl) :: r_min, r_max, m, b

    r_min = 0.1_dbl
    r_max = 6.0_dbl

    m = (V_screen2 - V_screen1) / (r_max - r_min)
    b = V_screen1 - m * r_min
    screening_potential = m * r + b
  end function screening_potential

end module screening_module
