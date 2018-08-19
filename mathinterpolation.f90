!
!	mathinterpolation.f90
!	PycnoCalc
!
!	Created by hellmersjl on 8/3/18.
!
!   Taken mostly from Numerical Recipes 3rd Edition: Press, Teukolsky, et. al.
!
module mathinterpolation

use constants

implicit none

save

contains 

	subroutine spline(x,y,yp1,ypn,y2)
	    
	    use mathlinearalg
	    
	    implicit none

    	real(kind=dbl), dimension(:), intent(in)    :: x,y
	    real(kind=dbl), intent(in)                  :: yp1,ypn
	    real(kind=dbl), dimension(:), intent(out)   :: y2
	    
	    integer                                     :: n
	    real(kind=dbl), dimension(size(x))          :: a,b,c,r
	    n = size(x)
	    ! TODO: verify size of input arrays
	
	    c(1:n-1) = x(2:n) - x(1:n-1)
	    r(1:n-1) = 6.0_dbl*((y(2:n) - y(1:n-1))/c(1:n-1))
	    r(2:n-1) = r(2:n-1)-r(1:n-2)
	    a(2:n-1) = c(1:n-2)
	    b(2:n-1) = 2.0_dbl*(c(2:n-1) + a(2:n-1))
	    b(1)=1.0_dbl
	    b(n)=1.0_dbl
	    if (yp1 > 0.99d30) then
		    r(1)=0.0_dbl
		    c(1)=0.0_dbl
	    else
		    r(1) = (3.0_dbl/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
		    c(1) = 0.5_dbl
	    end if
	    if (ypn .gt. 0.99d30) then
		    r(n) = 0.0_dbl
		    a(n) = 0.0_dbl
	    else
		    r(n) = (-3.0_dbl/(x(n) - x(n-1)))*((y(n) - y(n-1))/(x(n) - x(n-1)) - ypn)
		    a(n) = 0.5_dbl
	    end if
	    call tridag_par(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	end subroutine spline
    
    real(kind=dbl) function splint(xa,ya,y2a,x)
	
	    use mathutils
	    
	    implicit none
		
		real(kind=dbl), dimension(:), intent(in)    :: xa,ya,y2a
	    real(kind=dbl), intent(in)                  :: x

        integer         :: khi,klo,n
        real(kind=dbl)  :: a,b,h

        n = size(xa)
	    ! TODO: verify size of input arrays
        
        klo=max(min(locate(xa,x),n-1),1)
        
        khi=klo+1
	    h=xa(khi)-xa(klo)
	
	    if (h == 0.0) then
	        print *,'bad xa input in splint'
	    end if
	    a=(xa(khi)-x)/h
	    b=(x-xa(klo))/h
	    
	    splint=a*ya(klo)+b*ya(khi)+((a**3.0_dbl-a)*y2a(klo)+(b**3.0_dbl-b)*y2a(khi))*(h**2.0_dbl)/6.0_dbl
    
    end function        
        
end module mathinterpolation
