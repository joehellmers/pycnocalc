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


	real(kind=dbl) function nuclear_radius (A, inSQMFlag)
	
        use constants

        implicit none
	
	    real(kind=dbl), intent(in)      :: A
	    logical, intent(in), optional   :: inSQMFlag

        logical :: 	SQMFlag = .FALSE.

        if (present(inSQMFlag)) then
            SQMFlag = inSQMFlag
        end if
        
        if (SQMFlag) then
            nuclear_radius = A**(1.0_dbl/3.0_dbl)  ! p.5, Physical Props of Strangelets [Madsen]; first approximation... 
        else
            nuclear_radius = 1.31_dbl * (A**(1.0_dbl/3.0_dbl)) - 0.84_dbl
        end if

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

	real(kind=dbl) function density_2pF  (rho0, r, r0, A, inSQMFlag)
	
		implicit none
		
		real(kind=dbl), intent(in)      :: rho0         ! Central density of the nuclei/nugget
        real(kind=dbl), intent(in)      :: r            ! The current radius for which we want to calculate the density
        real(kind=dbl), intent(in)      :: r0           ! The total radius
        real(kind=dbl), intent(in)      :: A            ! The atomic mass in MeV
        logical, intent(in), optional   :: inSQMFlag    ! If the nuclei/nugget is SQM
	
        logical :: SQMFlag = .FALSE.

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        if (SQMFlag) then
            density_2pF = (3.0_dbl/4.0_dbl)*A/(pi*r0**3.0_dbl)
        else	
            density_2pF = rho0/(1.0_dbl + exp((r-r0)/a))	
	    end if

	end function density_2pF	

!******************************************
!
! Central nuclear density for 2 parameter Fermi
!
!******************************************

	real(kind=dbl) function rho0_2pF(A, radius, diffuseness)

	use mathintegration

	implicit none
	
	real(kind=dbl), intent(in)   :: A            ! number of nucleons
	real(kind=dbl), intent(in)  :: radius		! radius of nucleus
	real(kind=dbl), intent(in)  :: diffuseness	! diffusivity

    
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

    integer, intent(in)         ::  A1, Z1  !nucleon and proton number for target (1)
    real(kind=dbl), intent(in)  ::  A2, Z2  !nucleon and proton number for target (2)

        reduced_mass=((Z1*m_proton_mev+(A1-Z1)*m_neutron_mev)*(Z2*m_proton_mev+(A2-Z2)*m_neutron_mev))/((Z1+Z2)*m_proton_mev+(A1+A2-Z1-Z2)*m_neutron_mev)

    end function reduced_mass

!******************************************
!
! Adjust A2 and Z2 for SQM
!
!******************************************

    subroutine adjust_for_SQM(CFL,SQM_A2,SQM_mass, A2,Z2)   !output in units of MeV/c^2

        implicit none

        integer, intent(in)         :: CFL      ! The Color Flavor Locking indicator
        real(kind=dbl), intent(in)  :: SQM_A2   ! The Baryon number of the SQM "nugget"
        real(kind=dbl), intent(in)  :: SQM_mass ! The mass of the SQM "nugget"
        real(kind=dbl), intent(out) :: A2, Z2   ! The adjusted A2 and Z2 values

        A2 = SQM_A2
        if (CFL .eq. 1) then
            Z2 = 0.3_dbl*(SQM_mass/150.0_dbl)*((A2)**(2.0_dbl/3.0_dbl))                     ! charge of CFL state for strangelet...
        else
            if ((CFL .eq. 0) .and. (A2 .le. 1000)) then
                Z2 = 0.1_dbl*(A2)*((SQM_mass/150.0_dbl)**(2.0_dbl))                         ! charge of non-CFL stranglet state, low-baryon number
            else
                Z2 = 8.0_dbl*((SQM_mass/150.0_dbl)**(2.0_dbl))*((A2)**(1.0_dbl/3.0_dbl))    ! charge of non-CFL stranglet state, high-baryon number
            end if
        end if

    end subroutine adjust_for_SQM

    real(8) function number_density(A,radius)    !CURRENTLY IN USE ONLY FOR FIRST APPROX. FOR SQM
!	THIS DENSITY CALCULATION IS too simplistic!  NEED RHO(R)!	USE FINITE NUCLEI CALCULATION FOR NORMAL NUCLEAR MATTER AND USE
!	BAG MODEL TO CREATE SQM DENSITY IN NS CRUST!  

        implicit none

        real(kind=dbl), intent(in)  ::  A       !baryon number
        real(kind=dbl), intent(in)  ::  radius  !radius in fm

        number_density = (3.0_dbl/4.0_dbl)*A/(pi*radius**3.0_dbl)	  

    end function

end module general_nuclear
