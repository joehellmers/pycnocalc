!
!	tests.f90
!	PycnoCalc
!
!	Created by hellmersjl on 08/11/17.
!
!   This is a module containing the various tests to be performed

module tests

use constants
use general_nuclear
use rate_calc

implicit none

save

contains 

!******************************************
!
! Test the turn_pt routine
!
!******************************************

subroutine turn_pt_001

    real(kind=dbl)  ::  mu
    real(kind=dbl)  ::  radius1
    real(kind=dbl)  ::  radius2
    real(kind=dbl)  ::  rho0_a1
    real(kind=dbl)  ::  rho0_a2
    integer         ::  turn1
    integer         ::  turn2
    real(kind=dbl)  ::  WKB

    integer :: A1
    integer :: A2
    integer :: Z1
    integer :: Z2

    call scr_and_log_str ('TEST: trun_pt_001:')

    A1 = 56
    A2 = 40
    Z1 = 26
    Z2 = 16

    mu = reduced_mass(A1, Z1, A2, Z2)
    radius1 = nuclear_radius(A1)
    radius2 = nuclear_radius(A2)
    rho0_A1 = rho0_2pF(A1,radius1,0.5_dbl)
    rho0_A2 = rho0_2pF(A2,radius2,0.5_dbl)

    ! turn_pt(R,      Rstep  ,Rmax   ,E0                   , mu,A1,A2,SQM_A2 , Z1, Z2, radius1, radius2, rho0_A1, rho0_A2, L, partition, turn1,turn2,WKB)
    call turn_pt(393.4_dbl,1.0_dbl,23,53.707226090017635_dbl,mu,A1,A2,0.0_dbl, Z1, Z2, radius1, radius2, rho0_A1, rho0_A2, 0, 15       , turn1,turn2,WKB)

    call scr_and_log(str='WKB =',nbr=WKB,fmt='(ES13.5)',lf=.TRUE.)    
    call scr_and_log(str='turn1 =',intval=turn1,fmt='(I5)',lf=.TRUE.)    
    call scr_and_log(str='turn2 =',intval=turn2,fmt='(I5)',lf=.TRUE.)    

end subroutine turn_pt_001

	
end module tests

