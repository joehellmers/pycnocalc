!
!	vectors.f90
!	PycnoCalc
!
!	Created by hellmersjl on 9/23/07.
!
!
module vectors

use constants

implicit none

save

type :: vect3Dcart
	real(kind=dbl) :: x		! x-coordinate in cm
	real(kind=dbl) :: y		! y-coordinate in cm
    real(kind=dbl) :: z		! z-coordinate in cm
	real(kind=dbl) :: mag	! magnitude in cm
end type vect3Dcart

contains 
	type(vect3Dcart) function initVect3Dcart (x, y, z)
	use constants
	implicit none
	
	real(kind=dbl) x, y, z
	
		initVect3Dcart%x = x
		initVect3Dcart%y = y
		initVect3Dcart%z = z
		initVect3Dcart%mag = sqrt(x*x + y*y + z*z)
	
	end function initVect3Dcart
	
end module vectors


