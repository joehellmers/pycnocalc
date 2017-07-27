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

character(:), allocatable :: cfgfile ! The configuration file used

contains


!******************************************
!
! Subroutine to Initialize Global variables
!
!******************************************

subroutine initialize_globals

implicit none

! Setting Characteristic Length
char_r_factor = ((hbar*hbar)/((amu)*elementary_charge*elementary_charge))

!Setting Characteristic Energy Factor
char_E_factor = (amu*elementary_charge**4)/(hbar*hbar*boltzmann)

!Setting Dimensionless Inverse Length Factor
inv_len_factor = 1.0/((char_r_factor*char_r_factor*char_r_factor)/(2.0_dbl * amu))


end subroutine initialize_globals




end module globalvars

