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

    real(kind=dbl) function nucleonSaoPaulo(r,E0,mu)

        implicit none

        real(kind=dbl), intent(in) :: r     ! distance
        real(kind=dbl), intent(in) :: mu    ! reduced mass
        real(kind=dbl), intent(in) :: E0    ! first particle zero-point vibrational energy

        real(kind=dbl) :: a_sub_m = 0.3_dbl ! average matter diffuseness  parameter IN NUCLEON in fm

        ! USING THIS (4July):  using nucleon - nucleon interaction EQUATION 46, multiple times exp (-v2) term to get velocity
        ! dependent nucleon-nucleon interaction
        nucleonSaoPaulo=-456.0_dbl*exp(-8.0_dbl*E0/mu)*exp(-r/a_sub_m)*(1.0_dbl+(r/a_sub_m)+0.33_dbl*((r/a_sub_m)**2))/(64.0_dbl*pi*(a_sub_m**3))
        ! from eq. 46, Phys Rev C66, p11  AND from eq. 35, NL Desc Nuc Int. pdf, the units here are MeV*fm3

    end function nucleonSaoPaulo

end module nucleon_interactions
