!
!	astrophysics.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

module astrophysics

use constants
use screening_module

implicit none

save


contains 
!******************************************
!
! Function to calculate the mean weight per nuclus based on Atomic Weight and Number
! For 1 component matter
!******************************************


real(kind=dbl) function mean_wt_nucleus1comp (A, Z)

use constants

implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number


mean_wt_nucleus1comp = A*(1.0 + ((Z * electron_mass)/(A*amu)))

end function mean_wt_nucleus1comp

!******************************************
!
! Function to calculate the mean weight per nucleus 
! based on Atomic Weight and Number for two component matter
!
! SPVH1969 (15)
!
!******************************************

real(kind=dbl) function mean_wt_nucleus2comp(A1,Z1,X1,A2,Z2,X2)

use constants

implicit none

integer, intent(in) :: A1 ! Atomic Weight species 1
integer, intent(in) :: Z1 ! Atomic Number species 1
real(kind=dbl), intent(in) :: X1 ! Fraction of species 1
integer, intent(in) :: A2 ! Atomic Weight species 2
integer, intent(in) :: Z2 ! Atomic Number species 2
real(kind=dbl), intent(in) :: X2 ! Fraction of species 2

real(kind=dbl) :: part1
real(kind=dbl) :: part2
real(kind=dbl) :: mu_e

part1 = X1/A1
part1 = part1/(1.0_dbl + (Z1*electron_mass)/(A1*amu))

part2 = X2/A2
part2 = part2/(1.0_dbl + (Z2*electron_mass)/(A2*amu))

mu_e = 1.0_dbl/(part1+part2)

mean_wt_nucleus2comp = mu_e

end function mean_wt_nucleus2comp


!******************************************
!
! Function to calculate the mean weight per electron 
! based on Atomic Weight and Number for one component matter
!
! SPVH1969 (1)
!
!******************************************

real(kind=dbl) function mean_wt_electron1comp(A,Z)

use constants

implicit none

integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number


mean_wt_electron1comp = mean_wt_nucleus1comp(A,Z)/Z

end function mean_wt_electron1comp

!******************************************
!
! Function to calculate the mean weight per electron 
! based on Atomic Weight and Number for two component matter
!
! SPVH1969 (15)
!
!******************************************

real(kind=dbl) function mean_wt_electron2comp(A1,Z1,X1,A2,Z2,X2)

use constants

implicit none

integer, intent(in) :: A1 ! Atomic Weight species 1
integer, intent(in) :: Z1 ! Atomic Number species 1
real(kind=dbl), intent(in) :: X1 ! Fraction of species 1
integer, intent(in) :: A2 ! Atomic Weight species 2
integer, intent(in) :: Z2 ! Atomic Number species 2
real(kind=dbl), intent(in) :: X2 ! Fraction of species 2

real(kind=dbl) :: part1
real(kind=dbl) :: part2
real(kind=dbl) :: mu_e

part1 = (X1*Z1)/A1
part1 = part1/(1.0_dbl + (Z1*electron_mass)/(A1*amu))

part2 = (X2*Z2)/A2
part2 = part2/(1.0_dbl + (Z2*electron_mass)/(A2*amu))

mu_e = 1.0_dbl/(part1+part2)

mean_wt_electron2comp = mu_e

end function mean_wt_electron2comp

!******************************************
!
! Coulomb Excitation parameter 
!
! SPVH1969 (23)
!
!******************************************

real(kind=dbl) function beta_excitation2comp(A1_int,Z1_int,X1,A2_int,Z2_int,X2,rho,Temp)

use constants

implicit none

integer, intent(in)         :: A1_int   ! Atomic Weight species 1
integer, intent(in)         :: Z1_int   ! Atomic Number species 1
real(kind=dbl), intent(in)  :: X1       ! Fraction of species1
integer, intent(in)         :: A2_int   ! Atomic Weight species 2
integer, intent(in)         :: Z2_int   ! Atomic Number species 2
real(kind=dbl), intent(in)  :: X2       ! Fraction of species2
real(kind=dbl), intent(in)  :: rho      ! mass density in g/cm^3
real(kind=dbl), intent(in)  :: Temp     ! Temperature in Kelvin

real(kind=dbl)   :: A1
real(kind=dbl)   :: Z1
real(kind=dbl)   :: A2
real(kind=dbl)   :: Z2
real(kind=dbl)  :: part1
real(kind=dbl)  :: part2
real(kind=dbl)  :: part3
real(kind=dbl)  :: mu_e
real(kind=dbl)  :: fivethirds, onethird

A1 = real(A1_int,dbl)
Z1 = real(Z1_int,dbl)
A2 = real(A2_int,dbl)
Z2 = real(Z2_int,dbl)

fivethirds = (5.0_dbl/3.0_dbl)
onethird = (1.0_dbl/3.0_dbl)

mu_e = mean_wt_electron2comp(A1_int,Z1_int,X1,A2_int,Z2_int,X2)

part1 = (Z1 + Z2)**(fivethirds) - Z1**(fivethirds) - Z2**(fivethirds)

part2 = (Z1**2.0_dbl)*(Z2**2.0_dbl)*A1*A2
part2 = (part2/(A1+A2))**onethird

part1 = part1/part2


part3 = (42579000_dbl/Temp)*(((rho/mu_e)/16023000000_dbl)**onethird)

part1 = part1*part3

beta_excitation2comp = part1

end function beta_excitation2comp


!******************************************
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

nbr_density_nucleus = rho/(mean_wt_nucleus1comp(A,Z)*amu)

end function nbr_density_nucleus

!******************************************
!
! Function to calculate the number density of electrons for 2 component matters
!
!	- matter density
!	- Atomic Weight 1
!   - Atomic Number
!
! For 1 component matter
! 
! REF: SPVH1969 (1)
!
!******************************************
real(kind=dbl) function nbr_density_electron1comp (rho, A, Z)

use constants
implicit none

real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

nbr_density_electron1comp = rho/(mean_wt_electron1comp(A,Z)*amu)

end function nbr_density_electron1comp

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
real(kind=dbl) function inv_len_param1comp (rho, A, Z)
use constants
use globalvars
implicit none

real(kind=dbl), intent(in) :: rho ! mass density 
integer, intent(in) :: A ! Atomic Weight
integer, intent(in) :: Z ! Atomic Number

inv_len_param1comp = ((rho/(mean_wt_nucleus1comp(A,Z)*inv_len_factor))**(1.0/3.0))/(A*Z*Z)

end function inv_len_param1comp

!******************************************
!
! Function to calculate the inverse length parameter for two component material based on
!	- Density
!	- Atomic Weight1
!   - Atomic Number1
!   - Fraction of 1
!	- Atomic Weight2
!   - Atomic Number2
!   - Fraction of 2
! 
!	REF: SPVH1969 (16)
!
!******************************************
real(kind=dbl) function inv_len_param2comp (rho, A1, Z1, X1, A2, Z2, X2)
use constants
implicit none

real(kind=dbl), intent(in)  :: rho  ! mass density 
integer, intent(in)         :: A1   ! Atomic Weight
integer, intent(in)         :: Z1   ! Atomic Number
real(kind=dbl), intent(in)  :: X1   ! Fraction of species 1
integer, intent(in)         :: A2   ! Atomic Weight
integer, intent(in)         :: Z2   ! Atomic Number
real(kind=dbl), intent(in)  :: X2   ! Fraction of species 2

real(kind=dbl) :: mu_e

mu_e = mean_wt_electron2comp(A1, Z1, X1, A2, Z2, X2)

inv_len_param2comp = ((A1+A2)/(2.0_dbl*A1*A2*Z1*Z2))*((rho/((Z1*mu_e)*1.3574d11))**(1.0_dbl/3.0_dbl))


end function inv_len_param2comp


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

lambda = inv_len_param1comp(rho, A, Z)

! Calculate the common part first to reduce load on system
common_part = (rho/mean_wt_nucleus1comp(A,Z))*A*A*(Z**4)*S*1.00e46_dbl*(lambda**(7.0/4.0))

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
  ! Screened Coulomb Potential
  ! V_eff(r) = V_coulomb(r) - V_screen(r)
  !
  !*****************************************************************
  real(kind=dbl) function Vcoulomb_screened(Z1_int, Z2, r, radius1, radius2, &
                                            V_screen1, V_screen2)
    implicit none

    integer,      intent(in) :: Z1_int
    real(kind=dbl), intent(in) :: Z2
    real(kind=dbl), intent(in) :: r, radius1, radius2
    real(kind=dbl), intent(in) :: V_screen1, V_screen2

    real(kind=dbl) :: Vc, Vs

    ! Bare Coulomb
    Vc = Vcoulomb(Z1_int, Z2, r, radius1, radius2)

    ! Linear screening potential from screening_module
    Vs = screening_potential(r, V_screen1, V_screen2)

    ! Effective potential
    Vcoulomb_screened = Vc - Vs
  end function Vcoulomb_screened    
	
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

        ln_lambda=log(inv_len_param2comp (rho, A1_int, Z1_int, 0.5_dbl, int(A2), int(Z2), 0.5_dbl))
        
        ! energy parameter, called Rydberg energy, based on ground state energy from Bohr model of the hydrogen atom
        ln_Estar=log(real(Z1*Z2,dbl))-ln_rstar-34.174_dbl

        ! TODO: Where do these numbers come from?  They don't match SPVH
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
