!
!	general_nuclear.f90
!	PycnoCalc
!
!	Created by hellmersjl on 05/05/2008.
!
!
module general_nuclear

use constants

implicit none

save


contains 



!******************************************
!
! Nuclear Radius Approximation
!
!******************************************


	real(kind=dbl) function nuclear_radius (A)
	use constants
	implicit none
	
	integer, intent(in) :: A
	
	nuclear_radius = 1.31_dbl * (A**(1.0_dbl/3.0_dbl)) - 0.84_dbl
	
	end function nuclear_radius


!******************************************
!
! Function to integegrate in order normalize to get central density
!
!******************************************
	
	real(kind=dbl) function normalize_rho(r,a,r0)

		implicit none

		real(kind=dbl), intent(in) :: r
		real(kind=dbl), intent(in) :: a
		real(kind=dbl), intent(in) :: r0

		normalize_rho = (r*r)/(1.0_dbl+exp((r-r0)/a))

	end function normalize_rho
	


!******************************************
!
! Standard Fermi density
!
!******************************************

	real(kind=dbl) function density_2pF  (rho0, r, r0, a)
	
		implicit none
		
		real(kind=dbl), intent(in) :: rho0, r, r0, a
		
		density_2pF = rho0/(1.0_dbl + exp((r-r0)/a))	
	
	end function density_2pF
	
end module general_nuclear

