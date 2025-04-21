!
!	Constants.f90
!	PycnoCalc
!
!	Joe Hellmers
!

module constants

implicit none

save

integer, parameter :: dbl=8
integer, parameter :: qdr=16


! Global Constants

real(kind=dbl), parameter :: lightspeed = 2.99792458e8_dbl ! Speed of light in meters/sec
real(kind=dbl), parameter :: amu = 1.660538782e-24_dbl ! Atomic Mass Unit in grams
real(kind=dbl), parameter :: electron_mass = 1.9093826e-28_dbl ! Electron rest mass in grams
real(kind=dbl), parameter :: elementary_charge = 1.602176487e-19_dbl * 3.0e9_dbl ! Elementary Charge -e in esu

real(kind=dbl), parameter :: rho_sun_mean = 1.4_dbl ! Central Density of the sun g/cm^3
real(kind=dbl), parameter :: rho_sun_central = 160.0_dbl ! Central Density of the sun g/cm^3
real(kind=dbl), parameter :: rho_wd = 1.0e7_dbl ! Mean density of a White Dwarf g/cm^3
real(kind=dbl), parameter :: rho_ns = 1.0e15_dbl ! Mean density of Neutron Star g/cm^3
real(kind=dbl), parameter :: rho_ns_surface = 1.0e6_dbl ! Neutron Star Surface g/cm^3
real(kind=dbl), parameter :: rho_ns_outer_crust = 4.3e11_dbl ! Mean density of Neutron Star g/cm^3
real(kind=dbl), parameter :: rho_ns_inner_crust = 2.4e14_dbl ! Mean density of Neutron Star g/cm^3

real(kind=dbl), parameter :: hbar = 1.054571596e-27_dbl ! Normalized Planck's constant in erg*sec
real(kind=dbl), parameter :: boltzmann = 1.3806503e-23_dbl * 1.0e7_dbl ! Boltzmann's constant in erg/degree Kelvin
real(kind=dbl), parameter :: pi = 3.1415926535897932_dbl ! PI
real(kind=dbl), parameter :: avogadro  = 6.02214179d23! Avogadro's constant (mole^-1)

real(kind=dbl), parameter :: m_proton_mev = 938.272046_dbl ! Proton pass in MeV
real(kind=dbl), parameter :: m_neutron_mev = 939.5654133_dbl ! Neutron pass in MeV

! Conversion Factors
real(kind=dbl), parameter :: ev2erg = 1.602176487d-12 ! eV to ergs

! Precomputed numerical values
real(kind=dbl), parameter :: one_third = 1d0/3d0
real(kind=dbl), parameter :: two_thirds = 2d0/3d0

real(kind=dbl), parameter :: new_screen_constant = 1d0

end module constants
