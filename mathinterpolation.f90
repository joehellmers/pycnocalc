!
!	mathinterpolation.f90
!	PycnoCalc
!
!	Created by hellmersjl on 8/30/18.
!

module mathinterpolation

use constants

implicit none

save

contains 

    subroutine cubicSpline1DArray(N, x, y, ypp, a, b, c, d)

        implicit none
        
        integer, intent(in)                         :: N    ! Number of data points
        real(kind=dbl), dimension(N), intent(in)    :: x    ! abcissa values
        real(kind=dbl), dimension(N), intent(in)    :: y    ! ordinate values
        real(kind=dbl), dimension(N), intent(out)   :: ypp  ! approximation for 2nd derivative of ordinate
        
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
	
end module mathinterpolation
