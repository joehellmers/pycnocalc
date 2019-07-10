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
	real(kind=dbl) function nucleonM3Y (r, delta_0, inExcludeCore)
	use vectors
	implicit none
	
	real(kind=dbl), intent(in)      :: r
	real(kind=dbl), intent(in)      :: delta_0
	logical, intent(in), optional   :: inExcludeCore

	! Where did I get this number?
	
	real(kind=dbl), parameter   :: core_cutoff = 0.56754316782759184345508174374117515981197357177734375_dbl
   logical                     :: excludeCore = .false. 
	
	
	if (present(inExcludeCore)) then
        excludeCore = inExcludeCore
    end if	     
	
	if (excludeCore .and. (r .lt. core_cutoff)) then
	    nucleonM3Y = 0.0_dbl
	else
        if (r < delta_0) then
            nucleonM3Y = (7999.0_dbl*exp(-4.0_dbl*delta_0)/(4.0_dbl*delta_0) - 2134.0_dbl*exp(-2.5_dbl*delta_0)/(2.5_dbl*delta_0)) - 262.0_dbl
            !nucleonM3Y = - 262.0_dbl
        else    
            nucleonM3Y = (7999.0_dbl*exp(-4.0_dbl*r)/(4.0_dbl*r) - 2134.0_dbl*exp(-2.5_dbl*r)/(2.5_dbl*r)) 	
        end if
    end if
    
	end function nucleonM3Y

    real(kind=dbl) function nucleonSaoPaulo(r,E0,mu)

        implicit none

        real(kind=dbl), intent(in) :: r     ! distance
        real(kind=dbl), intent(in) :: mu    ! reduced mass
        real(kind=dbl), intent(in) :: E0    ! first particle zero-point vibrational energy

        real(kind=dbl) :: a_sub_m = 0.3_dbl ! average matter diffuseness  parameter IN NUCLEON in fm

        ! USING THIS (4July):  using nucleon - nucleon interaction EQUATION 46, multiple times exp (-v2) term to get velocity
        ! dependent nucleon-nucleon interaction
        nucleonSaoPaulo=-456.0_dbl*exp(-8.0_dbl*E0/mu)*exp(-r/a_sub_m)*(1.0_dbl+(r/a_sub_m)+0.33_dbl*((r/a_sub_m)**2))/(64.0_dbl*pi*(a_sub_m**3))
        ! TODO: Need to get this reference
        ! from eq. 46, Phys Rev C66, p11  AND from eq. 35, NL Desc Nuc Int. pdf, the units here are MeV*fm3

    end function nucleonSaoPaulo

end module nucleon_interactions
