!
!	mathdiff.f90
!	PycnoCalc
!
!	Created by hellmersjl on 10/19/18.
!
!   Taken mostly from Computational Physics Koonin and Meredith
!
module mathdiff

use constants

implicit none

save

contains
 
	real(kind=dbl) function ndiff5pt(f1, f2, f3, f4, h)
	
        real(kind=dbl), intent(in)   :: f1, f2, f3, f4, h
        
        ndiff5pt = (f1 - 8.0_dbl*f2 + 8.0_dbl*f3 - f4)/(h*12.0_dbl)

	end function ndiff5pt
	
end module mathdiff
