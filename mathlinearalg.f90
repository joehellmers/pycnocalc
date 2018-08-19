!
!	mathlinearalg.f90
!	PycnoCalc
!
!	Created by hellmersjl on 8/12/18.
!   
!   Taken mostly from Numerical Recipes 3rd Edition: Press, Teukolsky, et. al.
!

module mathlinearalg

use constants

implicit none

save

contains 
    
    subroutine tridag_ser(a,b,c,r,u)

        implicit none

	    real(kind=dbl), dimension(:), intent(in)    :: a,b,c,r
	    real(kind=dbl), dimension(:), intent(out)   :: u
        real(kind=dbl), dimension(size(b))          :: gam
        integer                                     :: n,j
        real(kind=dbl)                              :: bet

        n = size(b)
        ! TODO: verify all arrays are of the correct size
        
        ! Need to verify the pivot isn't zero
        bet=b(1)
        if (bet .eq. 0.0_dbl) then
            print *, "tridag_ser: Pivot point is zero!"
            return
        end if

	    u(1) = r(1)/bet
	    do j = 2, n
		    gam(j) = c(j-1)/bet
		    bet = b(j) - a(j-1)*gam(j)
		    if (bet == 0.0) then
			    print *, "tridag_ser: Error at code stage 2!"
			    return
			end if
		    u(j) = (r(j) - a(j-1)*u(j-1))/bet
	    end do
	
	    do j = n-1, 1, -1
		    u(j) = u(j) - gam(j+1)*u(j+1)
	    end do
       
    end subroutine tridag_ser
    
    recursive subroutine tridag_par(a,b,c,r,u)
    
        implicit none

	    real(kind=dbl), dimension(:), intent(in)    :: a,b,c,r
	    real(kind=dbl), dimension(:), intent(out)   :: u

        integer, parameter                          :: NPAR_TRIDAG=4
	    integer                                     :: n,n2,nm,nx
	    real(kind=dbl), dimension(size(b)/2)        :: y,q,piva
	    real(kind=dbl), dimension(size(b)/2-1)      :: x,z
	    real(kind=dbl), dimension(size(a)/2)        :: pivc

	    n = size(b)
	
	    if (n .lt. NPAR_TRIDAG) then
		    call tridag_ser(a,b,c,r,u)
	    else
		    if (maxval(abs(b(1:n))) == 0.0) then
			    print *,"tridag_par: possible singular matrix"
			end if
		    n2 = size(y)
		    nm = size(pivc)
		    nx = size(x)
		    piva = a(1:n-1:2)/b(1:n-1:2)
		    pivc = c(2:n-1:2)/b(3:n:2)
		    y(1:nm) = b(2:n-1:2) - piva(1:nm)*c(1:n-2:2) - pivc*a(2:n-1:2)
		    q(1:nm) = r(2:n-1:2) - piva(1:nm)*r(1:n-2:2) - pivc*r(3:n:2)
		    if (nm < n2) then
			    y(n2) = b(n)-piva(n2)*c(n-1)
			    q(n2) = r(n)-piva(n2)*r(n-1)
		    end if
		    x = -piva(2:n2)*a(2:n-2:2)
		    z = -pivc(1:nx)*c(3:n-1:2)
		    call tridag_par(x,y,z,q,u(2:n:2))
		    u(1) = (r(1) - c(1)*u(2))/b(1)
		    u(3:n-1:2) = (r(3:n-1:2) - a(2:n-2:2)*u(2:n-2:2) - c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
		    if (nm == n2) then
		        u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
		    end if
	    end if
	    
    end subroutine tridag_par
	
end module mathlinearalg
