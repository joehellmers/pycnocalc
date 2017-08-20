!
!	astrophysics.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

module astrophysics

use constants

implicit none

save


contains 
!******************************************
!
! Function to calculate the mean weight per nuclus based on Atomic Weight and Number
!
!******************************************


real(kind=dbl) function mean_wt_nucleus (A, Z)

use constants

implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number


mean_wt_nucleus = A*(1.0 + ((Z * electron_mass)/(A*amu)))

end function mean_wt_nucleus


!******************************************
!
! Function to calculate the mean weight per electron based on Atomic Weight and Number
!
!******************************************

real(kind=dbl) function mean_wt_electron(A,Z)

use constants

implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number


mean_wt_electron = mean_wt_nucleus(A,Z)/Z

end function mean_wt_electron

!*******************************d***********
!
! Function to calculate the number density of nuclei based on 
!	- matter density
!	- Atomic Weight
!   - Atomic Number
!
!******************************************

real(kind=dbl) function nbr_density_nucleus (rho, A, Z)
use constants
implicit none



real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

nbr_density_nucleus = rho/(mean_wt_nucleus(A,Z)*amu)

end function nbr_density_nucleus

!******************************************
!
! Function to calculate the number density of electrons based on 
!	- matter density
!	- Atomic Weight
!   - Atomic Number
!
!******************************************
real(kind=dbl) function nbr_density_electron (rho, A, Z)
use constants
implicit none



real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

nbr_density_electron = rho/(mean_wt_electron(A,Z)*amu)

end function nbr_density_electron

!******************************************
!
! Function to calculate the characteristic radius based on 
!	- Atomic Weight
!   - Atomic Number
! 
!	REF: SPVH1969 (2)
!
!******************************************

real(kind=dbl) function char_r (A, Z)
use constants
use globalvars
implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

char_r = char_r_factor/(A*Z*Z)

end function char_r

!******************************************
!
! Function to calculate the characteristic energy based on 
!	- Atomic Weight
!   - Atomic Number
! 
!	REF: SPVH1969 (2)
!
!******************************************

real(kind=dbl) function char_E (A, Z)
use constants
use globalvars
implicit none


integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

char_E = A*Z*Z*Z*Z*boltzmann*char_E_factor

end function char_E

!******************************************
!
! Function to calculate the inverse length parameter based on
!	- Density
!	- Atomic Weight
!   - Atomic Number
! 
!	REF: SPVH1969 (3)
!
!******************************************


real(kind=dbl) function inv_len_param (rho, A, Z)
use constants
use globalvars
implicit none



real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

inv_len_param = ((rho/(mean_wt_nucleus(A,Z)*inv_len_factor))**(1.0/3.0))/(A*Z*Z)


end function inv_len_param

!*****************************************************************
!
! Subroutine to calculate Pycnonuclear Reaction Rates at Zero Temp
!
!	INPUTS
!	- Density (g/cm^3)
!	- Atomic Weight
!   - Atomic Number
!   - SFactor (MeV*barns)f
! 
!	OUTPUTs 
!	- Lower Bound of the Rate (Static Approximation of Potential)
!     in (Reactions/sec/cm^3)
!   - Upper Bound of the Rate (Relaxed Approximation of Potential)
!     in (Reactions/sec/cm^3)
!
!	REF: SPVH1969 (39)
!
!*****************************************************************
subroutine react_rate_zero_temp(rho, A, Z, S, lower_rate, upper_rate)
use constants
use globalvars
implicit none


! Parameter declarations
real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number
real(kind=dbl), intent(in) :: S ! SFactor
real(kind=dbl), intent(out) :: lower_rate ! Lower Bound of Reaction Rate
real(kind=dbl), intent(out) :: upper_rate ! Upper Bound of Reaction Rate


! Function Headers
!real(kind=dbl) :: inv_len_param
!real(kind=dbl) :: mean_wt_nucleus


! local variables
real(kind=dbl) :: lambda
real(kind=dbl) :: common_part


! get the inverse length parameter

lambda = inv_len_param(rho, A, Z)

! Calculate the common part first to reduce load on system
common_part = (rho/mean_wt_nucleus(A,Z))*A*A*(Z**4)*S*1.00e46_dbl*(lambda**(7.0/4.0))

lower_rate = common_part*3.90_dbl*exp(-2.638_dbl/sqrt(lambda))

upper_rate = common_part*4.76_dbl*exp(-2.516_dbl/sqrt(lambda))

end subroutine react_rate_zero_temp

!*****************************************************************
!
! Coulomb Potential function - output in MEV
!
!*****************************************************************

    real(kind=dbl) function Vcoulomb(Z1_int,Z2,r,radius1,radius2)

        implicit none

        real(kind=dbl), intent(in)  :: r             ! distance (starts in fermi!)
        real(kind=dbl), intent(in)  :: radius1       ! radius of ion 1
        real(kind=dbl), intent(in)  :: radius2	    ! radius of ion 2
        integer, intent(in)         :: Z1_int                   ! proton count for nucleus 1
        real(kind=dbl), intent(in)  :: Z2  			        ! proton count for nucleus 2

        real(kind=dbl) :: Z1

        Z1 = dfloat(Z1_int)

        if (r .lt. (radius1+radius2)) then
            Vcoulomb=((3.0_dbl*(radius1+radius2)**2.0_dbl)-(r**2.0_dbl))*(Z1*Z2*1.4397_dbl)/(2.0_dbl*(radius1+radius2)**3.0_dbl)	! Vcoulomb is in MeV
        else
            Vcoulomb = (Z1*Z2*1.4397_dbl)/r	! Vcoulomb is in MeV
        end if

    end function Vcoulomb
	
!*****************************************************************
!
! Lattice + Vibrational Energy
! For Anisotropic Harmonic Oscillator
! Close to -0.5*hbar*w
! Output in unit of Rydberg Energy
!
!*****************************************************************
    real(kind=dbl) function E0_Energy(Z1_int,Z2, A1_int, A2, rho)

        implicit none

        integer, intent(in)         :: Z1_int           ! proton number for target (1)
        real(kind=dbl), intent(in)  :: Z2               ! proton number for target (2)
        integer, intent(in)         :: A1_int           ! nucleon number for target (1)
        real(kind=dbl), intent(in)  :: A2               ! nucleon number for target (2)
        real(kind=dbl), intent(in)  :: rho          ! density

        real(kind=dbl)  :: ln_rstar     ! natural log of Salpeter and VanHorn Bohr radius term
        real(kind=dbl)  :: ln_Estar     ! natural log of energy function from Salpeter and Vanhorn
        real(kind=dbl)  :: ln_lambda	! log of dimensionless length factor of Salpter and VanHorn
        real(kind=dbl)  :: ln_Ebcc		! natural log of the Salpeter bcc lattice energy - electrostatic interaction energy just from sitting in the lattice amidst other ions
        real(kind=dbl)  :: ln_Evib		! natural log of vibrational energy (salpeter)
        real(kind=dbl)  :: mn_mass, mn_chrg
        real(kind=dbl)  :: A1, Z1

        A1 = real(A1_int,dbl)
        Z1 = real(Z1_int,dbl)


        ln_rstar=log(real(A1+A2,dbl))-log(real(A1*A2*Z1*Z2,dbl))-31.87_dbl

        mn_mass=real((2*A1*A2)/(A1+A2),dbl)
        mn_chrg=real((Z1*A2+Z2*A1)/(A1+A2),dbl)

        ! I'm wondering if we can some how use the inverse length function in this module
        ln_lambda=-log(Z1**(1.0_dbl/3.0_dbl)+Z2**(1.0_dbl/3.0_dbl))+log((A1+A2)/(A1*A2*Z1*Z2))+(1.0_dbl/3.0_dbl)*log(rho)+(1.0_dbl/3.0_dbl)*log(mn_chrg)-(1.0_dbl/3.0_dbl)*log(mn_mass)-8.5455_dbl

        ! energy parameter, called Rydberg energy, based on ground state energy from Bohr model of the hydrogen atom
        ln_Estar=log(real(Z1*Z2,dbl))-ln_rstar-34.174_dbl

        ln_Ebcc = 0.5986_dbl + ln_lambda + ln_Estar
        ln_Evib = 0.6162_dbl + 0.5_dbl * ln_lambda + ln_Ebcc

        E0_Energy = exp(ln_Ebcc) + exp(ln_Evib)

    end function

    real(8) function lattice(rho,A,Z)   !output in units of fermi
!	Replace with bcc_spacing PYCNO
! 	function calculates the spacing between target particles - the exterior cube around a central body in the lattice of the crust

        implicit none

            real(kind=dbl), intent(in)    :: A, Z     !nucleon and proton number, respectively
            real(kind=dbl), intent(in)    :: rho      !density of crust in GRAMS/CM^3

            real(kind=dbl)    :: V                    !volume around each nucleus

            V=log10(Z*m_proton_mev+(A-Z)*m_neutron_mev) - log10(rho) + 12.25119_dbl	!in this line V is actually the log(V)
            V=10.0**V													!in this line, log(V) become V
            lattice = V**(1.0_dbl/3.0_dbl)  !output is in fermi

    end function

end module astrophysics

