!
!	general_nuclear.f90
!	PycnoCalc
!
!	Created by hellmersjl on 05/05/2008.
!
!
module general_nuclear

use constants

implicit none

save


contains 



!******************************************
!
! Nuclear Radius Approximation
!
!******************************************


	real(kind=dbl) function nuclear_radius (A)
	use constants
	implicit none
	
	integer, intent(in) :: A
	
	nuclear_radius = 1.31_dbl * (A**(1.0_dbl/3.0_dbl)) - 0.84_dbl
	
	end function nuclear_radius


!******************************************
!
! Function to integrate in order normalize to get central density
!
!******************************************
	
	real(kind=dbl) function normalize_rho(r,a,r0)

		implicit none

		real(kind=dbl), intent(in) :: r
		real(kind=dbl), intent(in) :: a
		real(kind=dbl), intent(in) :: r0

		normalize_rho = (r*r)/(1.0_dbl+exp((r-r0)/a))

	end function normalize_rho
	


!******************************************
!
! Standard Fermi density
!
!******************************************

	real(kind=dbl) function density_2pF  (rho0, r, r0, a)
	
		implicit none
		
		real(kind=dbl), intent(in) :: rho0, r, r0, a
		
		density_2pF = rho0/(1.0_dbl + exp((r-r0)/a))	
	
	end function density_2pF
	

!******************************************
!
! Central nuclear density for 2 parameter Fermi
!
!******************************************

	real(kind=dbl) function rho0_2pF(A, radius, diffuseness)

	use integration

	implicit none
	
	integer, intent(in) :: A 		! number of nucleons
	real(kind=dbl), intent(in) :: radius		! radius of nucleus
	real(kind=dbl), intent(in) :: diffuseness	! diffusivity

	rho0_2pF = A/((4.0_dbl*pi)*trapezoid(normalize_rho,0.0_dbl, 100.0_dbl, 1000, diffuseness, radius))

	end function rho0_2pF

!******************************************
!
! Spacing between nucleons at a particular density
!
! Function calculates the spacing between target particles - the exterior cube around a central body in the lattice of the crust
!
!******************************************

	real(kind=dbl) function bcc_spacing(rho, A, Z)   !output in units of fermi  

	use constants

	implicit none

	integer, intent(in) :: A, Z		!nucleon and proton number, respectively
	real(kind=dbl), intent(in) :: rho	!density of crust in GRAMS/CM^3
	real(kind=dbl) :: V			!volume around each nucleus 

      	V=log10(Z*m_proton_mev+(A-Z)*m_neutron_mev)-log10(rho)+12.25119	!in this line V is actually the log(V)
      	V=10.0**V                                               !in this line, log(V) becomes V
      	bcc_spacing = V**(1.0/3.0)  !output is in fermi 
	
	end function bcc_spacing

!******************************************
!
! Calculate reduced mass
!
!******************************************


    real(kind=dbl) function reduced_mass(A1,Z1,A2,Z2)   !output in units of MeV/c^2

    implicit none

    integer, intent(in) ::  A1, Z1, A2, Z2  !nucleon and proton number, respectively, for target (1) and projectile (2) nuclei, respectively

        reduced_mass=((Z1*m_proton_mev+(A1-Z1)*m_neutron_mev)*(Z2*m_proton_mev+(A2-Z2)*m_neutron_mev))/((Z1+Z2)*m_proton_mev+(A1+A2-Z1-Z2)*m_neutron_mev)

    end function reduced_mass


end module general_nuclear
