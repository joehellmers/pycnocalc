!
!	integration.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/13/08.
!	Copyright 2008 __MyCompanyName__. All rights reserved.
!

module integration

use constants

implicit none

save


contains 
	
	real(kind=dbl) function trapezoid (func, x0, xN, N, param1, param2)
	
		use constants
		
		implicit none
	
		
		! Input Parameters
		real(kind=dbl), external :: func
		real(kind=dbl), intent(in) :: x0, xN
		integer, intent(in) :: N
		real(kind=dbl), intent(in), optional :: param1
		real(kind=dbl), intent(in), optional :: param2
		
		
		integer :: i
		real(kind=dbl) :: accumulator
		real(kind=dbl) :: dx
		
		accumulator = 0.00_dbl
		dx = (xN - x0)/N
				
		do i=0,N
			if ((i .eq. 0) .or. (i .eq. N)) then
				if (present(param1) .AND. present(param2)) then
					accumulator = accumulator + func(x0 + i*dx, param1, param2)
				else
					if (present(param1)) then
						accumulator = accumulator + func(x0 + i*dx, param1)
					else
						if (present(param2)) then
							accumulator = accumulator + func(x0 + i*dx, param2)
						else
							accumulator = accumulator + func(x0 + i*dx)
						end if
					end if
				end if
			else
				if (present(param1) .AND. present(param2)) then
					accumulator = accumulator + 2.0_dbl*func(x0 + i*dx, param1, param2)
				else
					if (present(param1)) then
						accumulator = accumulator + 2.0_dbl*func(x0 + i*dx, param1)
					else
						if (present(param2)) then
							accumulator = accumulator + 2.0_dbl*func(x0 + i*dx, param2)
						else
							accumulator = accumulator + 2.0_dbl*func(x0 + i*dx)
						end if
					end if
				end if
			end if
		end do

		trapezoid = (1.0_dbl/2.0_dbl)*accumulator*dx
		
	end function trapezoid
	
end module integration
