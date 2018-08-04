!
!	mathintegration.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/13/08.
!

module mathintegration

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

! -----------------------------------------------------------------------------------------------------------------------
    real(kind=dbl) function trapezoidArray(Np, N, array, h)
!		DOCUMENTATION NOTE:  THIS SUBROUTINE IS - IN ITS ENTIRETY - FROM DR. CALVIN JOHNSON'S COMPUTATIONAL 
! 		COURSE - FALL 2006 AT SDSU  --- CWJ - SDSU - August 2006
!		ALL CREDIT AND ERROR FROM THIS SECTION ARE TO EITHER HIS CREDIT OR HIS FAULT (and mine, I suppose, for not recognizing it...)
!  subroutine to compute integral using trapezoidal rule
! INPUT:
!   Np: declared dimension of array
!   N : used dimension of array
!   array(0:Np): array of points of function to be integrated
!   h :  dx
! OUTPUT:
!   ans : integral
        implicit none
!..........INPUT...............
        integer, intent(in)         :: Np, N        ! dimensions
        real(kind=dbl), intent(in)  :: array(0:Np)     ! this is legal, since being fed from outside
        real(kind=dbl), intent(in)  :: h               ! = dx

!.........INTERMEDIATE..........
        integer i
        real(kind=dbl)  :: ans

        ans = 0.0_dbl            ! initialize output
        ans = 0.5_dbl*(array(0) + array(N))
        do i = 1, N-1
            ans = ans+array(i)
        end do
        ans = ans*h
        trapezoidArray = ans

    end function trapezoidArray
	
end module mathintegration
