!
!	nucleon_interactions.f90
!	PycnoCalc
!
!	Created by hellmersjl on 9/23/07.
!
!	CPHRG1997 (6)
!
module nucleon_interactions

use constants

contains
	real(kind=dbl) function nucleonM3Y (r,delta_0)
	use vectors
	implicit none
	
	real(kind=dbl) :: r
	real(kind=dbl) :: delta_0
	
	nucleonM3Y = (7999.0*exp(-4.0*r)/(4*r) - 2134.0*exp(-2.5*r)/(2.5*r)) 	
	if (r < delta_0) then
		nucleonM3Y = nucleonM3Y - 262.0_dbl
	end if

	end function nucleonM3Y

end module nucleon_interactions
