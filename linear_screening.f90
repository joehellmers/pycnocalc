
program linear_screening

  implicit none
  integer :: i, n_points
  real :: r_min, r_max, delta_r, r
  real :: V_screen1, V_screen2
  real :: V
  integer :: outfile
  
  r_min = 0.1
  r_max = 6.0
  n_points = 21
  delta_r = (r_max - r_min) / real(n_points - 1)

  ! Open output file for writing
  outfile = 10
  open(unit=outfile, file='linear_screening.dat', status='replace')  

  ! Input parameters
  print *, 'Enter V_screening(0.1 fm) [MeV]:'
  read *, V_screen1
  print *, 'Enter V_screening(6 fm) [MeV]:'
  read *, V_screen2

  print *
  print *, '   r (fm)      V_screening (MeV)'
  print *, '----------------------------------'

  do i = 1, n_points
     r = r_min + (i - 1) * delta_r
     V = screening_potential(r, V_screen1, V_screen2)
     print '(F8.3, 5X, F8.3)', r, V
     write(outfile, '(F8.3, 2X, F8.3)') r, V
  end do

  close(outfile)
  print *, 'Data written to linear_screening.dat'

contains

  real function screening_potential(r, V_screen1, V_screen2)
    implicit none
    real, intent(in) :: r, V_screen1, V_screen2
    real :: r_min, r_max, m, b

    r_min = 0.1
    r_max = 6.0

    m = (V_screen2 - V_screen1) / (r_max - r_min)
    b = V_screen1 - m * r_min
    screening_potential = m * r + b

  end function screening_potential

end program linear_screening

