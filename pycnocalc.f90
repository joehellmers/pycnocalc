!
!	PycnoCalc.f90
!	PycnoCalc
!
!	Joe Hellmers
!
!	References
!	SPVH1969: E. E. Salpeter and H. M. Van Horn, Astrophys. J. 155, 183 (1969)
!   SPVH1967: H. M. Van Horn, E. E. Salpeter, Phys. Rev, 157, 751
!   CPHRG1997: L. C. Chamon, B. V. Pereira, M. S. Hussein, M. A. Candido Ribeiro and D. Galetti, Phys. Rev. Lett. 79, 5218 (1997)
!	ST1983: S. L. Shapiro, S. A. Teukolsky "Black Holes, White Dwarfs, and Neutron Stars
!


program PycnoCalc

use constants
use vectors
use nucleon_interactions
use folding_potential

implicit none

! Function Headers
real(kind=dbl) :: mean_wt_nucleus
real(kind=dbl) :: mean_wt_electron
real(kind=dbl) :: nbr_density_nucleus
real(kind=dbl) :: nbr_density_electron
real(kind=dbl) :: char_r
real(kind=dbl) :: char_E
real(kind=dbl) :: inv_len_param
type(vect3Dcart) :: rvect
real(kind=dbl) :: NNM3Y ! Nucleon-Nucleon M3Y interaction
integer :: r_iterator
real(kind=dbl) :: this_r
real(kind=dbl) :: this_vfold


! Local Variables
real(4) :: t1
real(kind=dbl) :: work1, work2, work3, work4, work5, work6

call printintro

call initialize_globals


! Do some calculations and show them

rvect = initVect3Dcart(1.0_dbl,1.0_dbl,1.0_dbl)

NNM3Y = nucleonM3Y (rvect%mag,0.00_dbl)

print *, rvect%mag
print *, NNM3Y

print *


print *,'Zero Temp Reaction Rates'
print *,'------------------------'
print *,'Inside a White Dwarf'
print *,'C-C'
call react_rate_zero_temp(rho_wd, 12, 6, 4.0e16_dbl, work1, work2)
print *,'    Lower: ',work1
print *,'    Upper: ',work2
print *,'O-O'
call react_rate_zero_temp(rho_wd, 16, 8, 1.0e27_dbl, work1, work2)
print *,'    Lower: ',work1
print *,'    Upper: ',work2
print *,'Inside a Neutron Star'
print *,'C-C'
call react_rate_zero_temp(rho_ns, 12, 6, 4.0e16_dbl, work1, work2)
print *,'    Lower: ',work1
print *,'    Upper: ',work2
print *,'O-O'
call react_rate_zero_temp(rho_ns, 16, 8, 1.0e27_dbl, work1, work2)
print *,'    Lower: ',work1
print *,'    Upper: ',work2

print *,'Neutron Star'
print *,'------------'
print *,'C-C'
call react_rate_zero_temp(rho_ns_surface, 12, 6, 4.0e16_dbl, work1, work2)
print *,'Surface'
print *,'    Lower: ',work1
print *,'    Upper: ',work2
call react_rate_zero_temp(rho_ns_outer_crust, 12, 6, 4.0e16_dbl, work1, work2)
print *,'Outer Crust'
print *,'    Lower: ',work1
print *,'    Upper: ',work2
call react_rate_zero_temp(rho_ns_inner_crust, 12, 6, 4.0e16_dbl, work1, work2)
print *,'Inner Crust'
print *,'    Lower: ',work1
print *,'    Upper: ',work2
print *,'O-O'
call react_rate_zero_temp(rho_ns_surface, 16, 8, 1.0e27_dbl, work1, work2)
print *,'Surface'
print *,'    Lower: ',work1
print *,'    Upper: ',work2
call react_rate_zero_temp(rho_ns_outer_crust, 16, 8, 1.0e27_dbl, work1, work2)
print *,'Outer Crust'
print *,'    Lower: ',work1
print *,'    Upper: ',work2
call react_rate_zero_temp(rho_ns_inner_crust, 16, 8, 1.0e27_dbl, work1, work2)
print *,'Inner Crust'
print *,'    Lower: ',work1
print *,'    Upper: ',work2


t1 = secnds(0.0);

! Oxygen-16 <-> Nickel-58

do r_iterator=1,120
	this_r = 0.1_dbl*r_iterator
	this_vfold = -1.0_dbl*vfold_spherically_symmetric(0.1_dbl*r_iterator,58,16,4.230848400_dbl,2.460993151_dbl,15,0.50_dbl,0.50_dbl,0.1606844454_dbl,0.1820414223_dbl,6.0_dbl,6.0_dbl)
	print *,this_r,',',log10(this_vfold)
end do
print *,'Time elapsed=',secnds(t1)
t1 = secnds(0.0);

! Oxygen-16 <-> Carbor-12
!print *,vfold_uniform_sphere(1.0_dbl,16,12,3.15_dbl,2.86_dbl,15,0.50,0.50,0.1820414223,0.1859855319)
!print *,'Time elapsed=',secnds(t1)
!t1 = secnds(0.0);

! Oxygen-16 <-> Oxygen-15
!print *,vfold_uniform_sphere(1.0_dbl,16,16,3.15_dbl,3.15_dbl,15,0.50,0.50,0.1820414223,0.1820414223)
!print *,'Time elapsed=',secnds(t1)
!t1 = secnds(0.0);


print *, ' -------  End Run: PycnoCalc --------'



end program PycnoCalc


!******************************************
!
! Subroutine to Initialize Global variables
!
!******************************************


subroutine initialize_globals
use constants
use globalvars
implicit none

char_r_factor = ((hbar*hbar)/((amu)*elementary_charge*elementary_charge))
print *, 'Setting Characteristic Length Factor to', char_r_factor
print *
char_E_factor = (amu*elementary_charge*elementary_charge*elementary_charge*elementary_charge)/(hbar*hbar*boltzmann)
print *, 'Setting Characteristic Energy Factor to', char_E_factor
print *
inv_len_factor = 1.0/((char_r_factor*char_r_factor*char_r_factor)/(2.0_dbl * amu))
print *, 'Setting Dimensionless Inverse Length Factor to', inv_len_factor
print *


end subroutine initialize_globals



!******************************************
!
! Subroutine to print a banner
!
!******************************************

subroutine printintro
use constants

use logging

implicit none

call scr_and_log ('This is PycnoCalc - Pycnonulear Reaction Calculation System')


print *,'***********************************************************'
print *,'   A Collaboration of: Dr. Fridonlin Weber' 
print *,'                       Graduate Student Barbara Golf'
print *,'                       Undergrad Joe Hellmers'
print *,'***********************************************************'
print *

end subroutine printintro


!******************************************
!
! Function to calculate the mean weight per nuclus based on Atomic Weight and Number
!
!******************************************


real(kind=dbl) function mean_wt_nucleus (A, Z)

use constants

implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number


mean_wt_nucleus = A*(1.0 + ((Z * electron_mass)/(A*amu)))

end function mean_wt_nucleus


!******************************************
!
! Function to calculate the mean weight per electron based on Atomic Weight and Number
!
!******************************************

real(kind=dbl) function mean_wt_electron(A,Z)

use constants

implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number
real(kind=dbl) :: mean_wt_nucleus


mean_wt_electron = mean_wt_nucleus(A,Z)/Z

end function mean_wt_electron



!******************************************
!
! Function to calculate the number density of nuclei based on 
!	- matter density
!	- Atomic Weight
!   - Atomic Number
!
!******************************************

real(kind=dbl) function nbr_density_nucleus (rho, A, Z)
use constants
implicit none



real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number
real(kind=dbl) :: mean_wt_nucleus

nbr_density_nucleus = rho/(mean_wt_nucleus(A,Z)*amu)

end function nbr_density_nucleus

!******************************************
!
! Function to calculate the number density of electrons based on 
!	- matter density
!	- Atomic Weight
!   - Atomic Number
!
!******************************************
real(kind=dbl) function nbr_density_electron (rho, A, Z)
use constants
implicit none



real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number
real(kind=dbl) :: mean_wt_electron

nbr_density_electron = rho/(mean_wt_electron(A,Z)*amu)

end function nbr_density_electron


!******************************************
!
! Function to calculate the characteristic radius based on 
!	- Atomic Weight
!   - Atomic Number
! 
!	REF: SPVH1969 (2)
!
!******************************************

real(kind=dbl) function char_r (A, Z)
use constants
use globalvars
implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

char_r = char_r_factor/(A*Z*Z)

end function char_r


!******************************************
!
! Function to calculate the characteristic energy based on 
!	- Atomic Weight
!   - Atomic Number
! 
!	REF: SPVH1969 (2)
!
!******************************************

real(kind=dbl) function char_E (A, Z)
use constants
use globalvars
implicit none


integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

char_E = A*Z*Z*Z*Z*boltzmann*char_E_factor

end function char_E

!******************************************
!
! Function to calculate the inverse length parameter based on
!	- Density
!	- Atomic Weight
!   - Atomic Number
! 
!	REF: SPVH1969 (3)
!
!******************************************


real(kind=dbl) function inv_len_param (rho, A, Z)
use constants
use globalvars
implicit none



real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number
real(kind=dbl) :: mean_wt_nucleus
real(kind=dbl) :: nbr_density_nucleus
real(kind=dbl) :: char_r

!inv_len_param = char_r(A,Z)*(nbr_density_nucleus(rho,A,Z)/2.0)**(1.0/3.0)
inv_len_param = ((rho/(mean_wt_nucleus(A,Z)*inv_len_factor))**(1.0/3.0))/(A*Z*Z)


end function inv_len_param

!*****************************************************************
!
! Subroutine to calculate Pycnonuclear Reaction Rates at Zero Temp
!
!	INPUTS
!	- Density (g/cm^3)
!	- Atomic Weight
!   - Atomic Number
!   - SFactor (MeV*barns)f
! 
!	OUTPUTs 
!	- Lower Bound of the Rate (Static Approximation of Potential)
!     in (Reactions/sec/cm^3)
!   - Upper Bound of the Rate (Relaxed Approximation of Potential)
!     in (Reactions/sec/cm^3)
!
!	REF: SPVH1969 (39)
!
!*****************************************************************
subroutine react_rate_zero_temp(rho, A, Z, S, lower_rate, upper_rate)
use constants
use globalvars
implicit none


! Parameter declarations
real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number
real(kind=dbl), intent(in) :: S ! SFactor
real(kind=dbl), intent(out) :: lower_rate ! Lower Bound of Reaction Rate
real(kind=dbl), intent(out) :: upper_rate ! Upper Bound of Reaction Rate


! Function Headers
real(kind=dbl) :: inv_len_param
real(kind=dbl) :: mean_wt_nucleus


! local variables
real(kind=dbl) :: lambda
real(kind=dbl) :: common_part


! get the inverse length parameter

lambda = inv_len_param(rho, A, Z)

! Calculate the common part first to reduce load on system
common_part = (rho/mean_wt_nucleus(A,Z))*A*A*Z*Z*Z*Z*S*1.00e46_dbl*(lambda**(7.0/4.0))

lower_rate = common_part*3.90_dbl*exp(-2.638_dbl/sqrt(lambda))

upper_rate = common_part*4.76_dbl*exp(-2.516_dbl/sqrt(lambda))

end subroutine react_rate_zero_temp



