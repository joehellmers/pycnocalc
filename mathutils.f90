!
!	mathutils.f90
!	PycnoCalc
!
!	Created by hellmersjl on 8/15/18.
!
!   Taken mostly from Numerical Recipes 3rd Edition: Press, Teukolsky, et. al.
!
module mathutils

use constants

implicit none

save

contains 


	function arth_d(first,increment,n)
	
	    real(kind=dbl), intent(in)      :: first,increment
	    integer, intent(in)             :: n

        integer, parameter              :: NPAR_ARTH = 16
        integer, parameter              :: NPAR2_ARTH = 8

	    real(kind=dbl), dimension(n)    :: arth_d
	    integer                         :: k,k2
	    real(kind=dbl)                  :: temp
	
	    if (n .gt. 0) then
	        arth_d(1) = first
	    else
	        if (n .le. NPAR_ARTH) then
		        do k = 2,n
			        arth_d(k) = arth_d(k-1) + increment
		        end do
	        else
		        do k = 2,NPAR2_ARTH
			        arth_d(k) = arth_d(k-1) + increment
		        end do
		        temp = increment*NPAR2_ARTH
		        k = NPAR2_ARTH
		        do
			        if (k .ge. n) then
			            exit
			        end if
			        k2 = k + k
			        arth_d(k+1:min(k2,n)) = temp + arth_d(1:min(k,n-k))
			        temp=temp + temp
			        k = k2
		        end do
	        end if
	    end if
	end function arth_d

	integer function locate(xx,x)

        implicit none

        real(kind=dbl), dimension(:), intent(in) :: xx
        real(kind=dbl), intent(in) :: x

        integer :: n,jl,jm,ju
        logical :: ascnd
    
        n=size(xx)
    
        ascnd = (xx(n) >= xx(1))
        jl=0
        ju=n+1
        do
            if (ju-jl <= 1) exit
            jm=(ju+jl)/2
            if (ascnd .eqv. (x >= xx(jm))) then
                jl=jm
            else
                ju=jm
            end if
        end do
        if (x == xx(1)) then
            locate=1
        else if (x == xx(n)) then
            locate=n-1
        else
            locate=jl
        end if
	end function locate

	
end module mathutils
