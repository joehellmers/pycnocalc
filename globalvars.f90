!
!	globalvars.f90
!	PycnoCalc
!
!	Joe Hellmers
!
module globalvars

use constants

implicit none

save

! Global Variables
integer :: A1 ! Atomic Mass Number 1
integer :: A2 ! Atomic Mass Number 2
integer :: Z1 ! Atomic Number 1
integer :: Z2 ! Atomic Number 2

real(kind=dbl) :: matter_density ! Density in cm/g^3 



real(kind=dbl) :: char_r_factor  ! The Characteristic radius factor in cm
real(kind=dbl) :: char_E_factor  ! The Characteristic energy factor in ergs
real(kind=dbl) :: inv_len_factor ! The dimensionless Inverse Length factor



end module globalvars

