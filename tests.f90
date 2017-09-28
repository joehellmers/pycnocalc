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
! turn_pt subroutine: Test 001
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
    real(kind=dbl)  ::  E0
    real(kind=dbl)  ::  rho

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    integer         :: nucIntType = 1

    call scr_and_log_str ('TEST: trun_pt_001:')

    A1 = 56
    A2 = 40.0_dbl
    Z1 = 26
    Z2 = 16.0_dbl

    mu = reduced_mass(A1, Z1, A2, Z2)
    radius1 = nuclear_radius(real(A1,dbl),.FALSE.)
    radius2 = nuclear_radius(A2,.FALSE.)
    rho0_A1 = rho0_2pF(real(A1,dbl),radius1,0.5_dbl)
    rho0_A2 = rho0_2pF(A2,radius2,0.5_dbl)

    rho = 5500000000000.0029_dbl
    E0 = E0_Energy(Z1,Z2,A1,A2,rho)

    ! turn_pt(R,      Rstep  ,Rmax   ,E0,mu,A1,A2,SQM_A2 , Z1, Z2, radius1, radius2, rho0_A1, rho0_A2, L, partition, turn1,turn2,WKB)
    call turn_pt(393.4_dbl,1.0_dbl,23,E0,mu,A1,A2, Z1, Z2, radius1, radius2, rho0_A1, rho0_A2, 0, 15       , turn1,turn2,WKB,nucIntType, .FALSE.)

    call scr_and_log(str='WKB =',nbr=WKB,fmt='(ES13.5)',lf=.TRUE.)    
    call scr_and_log(str='turn1 =',intval=turn1,fmt='(I5)',lf=.TRUE.)    
    call scr_and_log(str='turn2 =',intval=turn2,fmt='(I5)',lf=.TRUE.)    

end subroutine turn_pt_001

!******************************************
!
! E0_Energy Function: Test 001
!
!******************************************

subroutine E0_Energy_001

    real(kind=dbl)  ::  E0
    real(kind=dbl)  ::  rho

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2

    call scr_and_log_str ('TEST: E0_Energy_001:')

    A1 = 56
    A2 = 40.0_dbl
    Z1 = 26
    Z2 = 26.0_dbl
    rho = 5500000000000.0029_dbl

    E0 = E0_Energy(Z1,Z2, A1, A2, rho)

    call scr_and_log(str='E0 = ',nbr=E0,fmt='(ES13.5)',lf=.TRUE.)

end subroutine E0_Energy_001
	

!******************************************
!
! Sfactor function: Test 001
!
!******************************************

subroutine Sfactor_001

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: Rstep = 1.0_dbl
    real(kind=dbl)  :: SQM_A2 = 0.0_dbl
    integer :: partition = 15
    integer :: nucIntType = 1

    real(kind=dbl)  ::  S

    call scr_and_log_str ('TEST: Sfactor_001:')

    A1 = 56
    A2 = 40_dbl
    Z1 = 26
    Z2 = 26_dbl

    rho = 5500000000000.0029_dbl

    S = Sfactor(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE.)

    call scr_and_log(str='S =',nbr=S,fmt='(ES13.5)',lf=.TRUE.)

end subroutine Sfactor_001

!******************************************
!
! pycnoRate function: Test 001
!
!******************************************
subroutine pycnoRate_001

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  ::  rho
    real(kind=dbl)  ::  Rstep = 1.0_dbl
    integer         :: partition = 15
    integer         :: nucIntType = 1

    real(kind=dbl)  ::  rate

    call scr_and_log_str ('TEST: pycnoRate_001:')

    A1 = 56
    A2 = 40.0_dbl
    Z1 = 26
    Z2 = 26.0_dbl

    rho = 550000000000000.0029_dbl

    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE.)

    call scr_and_log(str='Rate =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)

end subroutine pycnoRate_001

!******************************************
!
! pycnoRate function: Test 002
! 
! With SQM, CFL=1
!
!******************************************
subroutine pycnoRate_002

    integer         :: A1
    integer         :: Z1
    real(kind=dbl)  :: A2
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: Rstep = 1.0_dbl
    integer         :: partition = 15
    integer         :: CFL
    integer         :: nucIntType = 1
    real(kind=dbl)  :: SQM_A2
    real(kind=dbl)  :: SQM_mass

    
    real(kind=dbl)  ::  rate

    call scr_and_log_str ('TEST: pycnoRate_002:')

    A1 = 56
    Z1 = 26
    A2 = 40.0_dbl
    Z2 = 26.0_dbl
    CFL = 1    
    SQM_A2 = 40.0_dbl
    SQM_mass = 300.0_dbl

    rho = 5500000000000.0029_dbl

    call adjust_for_SQM(CFL, SQM_A2, SQM_mass, A2, Z2)

    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .TRUE.)

    call scr_and_log(str='Rate =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)

end subroutine pycnoRate_002

!******************************************
!
! pycnoRate function: Test 002
! 
! With SQM, CFL=1
!
!******************************************
subroutine pycnoRate_003

    integer         :: A1
    integer         :: Z1
    real(kind=dbl)  :: A2
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: Rstep = 1.0_dbl
    integer         :: partition = 15
    integer         :: nucIntType = 1
    integer         :: CFL
    real(kind=dbl)  :: SQM_A2
    real(kind=dbl)  :: SQM_mass

    
    real(kind=dbl)  ::  rate

    call scr_and_log_str ('TEST: pycnoRate_003:')

    A1 = 56
    Z1 = 26
    A2 = 40.0_dbl
    Z2 = 26.0_dbl
    CFL = 0
    SQM_A2 = 40.0_dbl
    SQM_mass = 300.0_dbl

    rho = 5500000000000.0029_dbl

    call adjust_for_SQM(CFL, SQM_A2, SQM_mass, A2, Z2)

    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .TRUE.)

    call scr_and_log(str='Rate =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)

end subroutine pycnoRate_003

!******************************************
!
! pycnoRate function: Test 004
!
!******************************************
subroutine pycnoRate_004

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  ::  rho
    real(kind=dbl)  ::  Rstep = 1.0_dbl
    integer         :: partition = 15
    integer         :: nucIntType = 2

    real(kind=dbl)  ::  rate

    call scr_and_log_str ('TEST: pycnoRate_004:')

    A1 = 56
    A2 = 40.0_dbl
    Z1 = 26
    Z2 = 26.0_dbl

    rho = 5500000000000.0029_dbl

    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE.)

    call scr_and_log(str='Rate =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)

end subroutine pycnoRate_004

!******************************************
!
!  subroutine adjust_for_SQM: Test 001
!
!******************************************
subroutine adjust_for_SQM_001

    real(kind=dbl)  :: A2
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: SQM_mass
    real(kind=dbl)  :: SQM_A2
    integer         :: CFL

    call scr_and_log_str ('TEST: adjust_for_SQM_001:')

    CFL = 1    
    SQM_A2 = 40
    SQM_mass = 300

    call adjust_for_SQM(CFL, SQM_A2, SQM_mass, A2, Z2)

    call scr_and_log(str='Adjusted A2 =',nbr=A2,fmt='(ES13.5)',lf=.TRUE.)
    call scr_and_log(str='Adjusted Z2 =',nbr=Z2,fmt='(ES13.5)',lf=.TRUE.)

end subroutine adjust_for_SQM_001

!******************************************
!
!  subroutine adjust_for_SQM: Test 002
!
!******************************************
subroutine adjust_for_SQM_002

    real(kind=dbl)  :: A2
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: SQM_mass
    real(kind=dbl)  :: SQM_A2
    integer         :: CFL

    call scr_and_log_str ('TEST: adjust_for_SQM_002:')

    CFL = 0
    SQM_A2 = 40_dbl
    SQM_mass = 300_dbl

    call adjust_for_SQM(CFL, SQM_A2, SQM_mass, A2, Z2)

    call scr_and_log_str("For Low Baryon Number")
    call scr_and_log(str='Adjusted A2 =',nbr=A2,fmt='(ES13.5)',lf=.TRUE.)
    call scr_and_log(str='Adjusted Z2 =',nbr=Z2,fmt='(ES22.14)',lf=.TRUE.)

    CFL = 0
    SQM_A2 = 1500_dbl
    SQM_mass = 300_dbl

    call adjust_for_SQM(CFL, SQM_A2, SQM_mass, A2, Z2)

    call scr_and_log_str("For High Baryon Number")
    call scr_and_log(str='Adjusted A2 =',nbr=A2,fmt='(ES13.5)',lf=.TRUE.)
    call scr_and_log(str='Adjusted Z2 =',nbr=Z2,fmt='(ES22.14)',lf=.TRUE.)

end subroutine adjust_for_SQM_002

subroutine vfold_cubes_001

    real(kind=dbl)  :: vfold
    real(kind=dbl)  :: E0

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: rho
    integer(kind=dbl) :: N1 = 10000000000_dbl

    call scr_and_log_str ('TEST: vfold_cubes_001:')

    A1 = 56
    A2 = 40.0_dbl
    Z1 = 26
    Z2 = 26.0_dbl

    rho = 5500000000000.0029_dbl
    E0 = E0_Energy(Z1,Z2,A1,A2,rho)

    vfold = vfold_cubes(5.0_dbl, A1, A2, Z1, Z2, N1, 1_dbl, 2.0_dbl, 0.1_dbl, E0, 1, .FALSE.)
    call scr_and_log(str='vfold =',nbr=vfold,fmt='(ES13.5)',lf=.TRUE.)


end subroutine vfold_cubes_001

subroutine mean_wt_001

    integer  :: A1
    integer  :: A2
    integer  :: Z1
    integer  :: Z2
    real(kind=dbl)  :: X1
    real(kind=dbl)  :: X2

    real(kind=dbl)  :: mean_wt

    A1 = 16
    Z1 = 8
    A2 = 16
    Z2 = 8
    X1 = 0.5_dbl
    X2 = 0.5_dbl

    mean_wt = mean_wt_electron2comp(A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='O-O 50-50: mu_e =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)

    mean_wt = mean_wt_electron1comp(A1,Z1)
    call scr_and_log(str='O-O 50-50 for 1 component: mu_e =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)

    mean_wt = mean_wt_nucleus2comp(A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='O-O 50-50: mu_A =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)

    mean_wt = mean_wt_nucleus1comp(A1,Z1)
    call scr_and_log(str='O-O 50-50 for 1 component: mu_A =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)


    A1 = 16
    Z1 = 8
    A2 = 12
    Z2 = 6
    X1 = 0.5_dbl
    X2 = 0.5_dbl

    mean_wt = mean_wt_electron2comp(A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='O-C 50-50: mu_e =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)

    mean_wt = mean_wt_nucleus2comp(A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='O-C 50-50: mu_A =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)

    A1 = 18
    Z1 = 8
    A2 = 12
    Z2 = 6
    X1 = 0.5_dbl
    X2 = 0.5_dbl

    mean_wt = mean_wt_electron2comp(A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='O18-C 50-50: mu_e =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)

    mean_wt = mean_wt_nucleus2comp(A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='O18-C 50-50: mu_A =',nbr=mean_wt,fmt='(ES13.5)',lf=.TRUE.)


end subroutine mean_wt_001

subroutine beta_excite_001

    integer  :: A1
    integer  :: A2
    integer  :: Z1
    integer  :: Z2
    real(kind=dbl)  :: X1
    real(kind=dbl)  :: X2
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: Temp

    real(kind=dbl)  :: beta

    A1 = 1
    Z1 = 1
    A2 = 1
    Z2 = 1
    X1 = 0.5_dbl
    X2 = 0.5_dbl
    Temp = 4.2579d7
    rho = 1.6203d10

    beta = beta_excitation2comp(A1,Z1,X1,A2,Z2,X2,rho,Temp)
    call scr_and_log(str='H-H 50-50 Unitary Temp Unitary rho: beta =',nbr=beta,fmt='(ES13.5)',lf=.TRUE.)

    A1 = 12
    Z1 = 6
    A2 = 12
    Z2 = 6
    X1 = 0.5_dbl
    X2 = 0.5_dbl
    Temp = 4.2579d7
    rho = 1.6203d10
    beta = beta_excitation2comp(A1,Z1,X1,A2,Z2,X2,rho,Temp)

    call scr_and_log(str='O-O 50-50 Unitary Temp Unitary rho: beta =',nbr=beta,fmt='(ES13.5)',lf=.TRUE.)


end subroutine beta_excite_001

subroutine inv_len_001

    integer         :: A1
    integer         :: Z1
    real(kind=dbl)  :: X1
    integer         :: A2
    integer         :: Z2
    real(kind=dbl)  :: X2
    real(kind=dbl)  :: rho

    real(kind=dbl)  :: lambda
    real(kind=dbl)  :: mu_e

    A1 = 1
    Z1 = 1
    A2 = 1
    Z2 = 1
    X1 = 0.5_dbl
    X2 = 0.5_dbl
    rho = 1.3574d11


    lambda =    inv_len_param2comp(rho,A1,Z1,X1,A2,Z2,X2)
    call scr_and_log(str='2comp: H-H 50-50 Unitary rho: lambda =',nbr=lambda,fmt='(ES13.5)',lf=.TRUE.)

    A1 = 12
    Z1 = 6
    A2 = 12
    Z2 = 6
    X1 = 0.5_dbl
    X2 = 0.5_dbl
    rho = 1.3574d11
    
    mu_e = mean_wt_electron2comp(A1,Z1,X1,A2,Z2,X2)
    
    lambda = inv_len_param2comp(rho,A1,Z1,X1,A2,Z2,X2)

    call scr_and_log(str='2comp:O-O 50-50 Unitary rho: lambda =',nbr=lambda,fmt='(ES13.5)',lf=.TRUE.)
    call scr_and_log(str='    mu_e =',nbr=mu_e,fmt='(ES13.5)',lf=.TRUE.)


end subroutine inv_len_001

end module tests


