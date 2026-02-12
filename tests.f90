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
use screening_module

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

    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl

    mu = reduced_mass(A1, Z1, A2, Z2)
    radius1 = nuclear_radius(real(A1,dbl),.FALSE.)
    radius2 = nuclear_radius(A2,.FALSE.)
    rho0_A1 = rho0_2pF(real(A1,dbl),radius1,0.5_dbl)
    rho0_A2 = rho0_2pF(A2,radius2,0.5_dbl)

    rho = 1000000000.0_dbl
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

    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl
    rho = 1000000000.0_dbl

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
    integer :: nucIntType = 3

    real(kind=dbl)  ::  S

    call scr_and_log_str ('TEST: Sfactor_001:')

    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl

    rho = 1000000000.0_dbl

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
    integer         :: nucIntType = 2

    real(kind=dbl)  ::  rate

    call scr_and_log_str ('TEST: pycnoRate_001:')

    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl

    rho = 1000000000.0_dbl

    
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 0)    
    call scr_and_log(str='Rate Original bcc =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 1)    
    call scr_and_log(str='Rate Schramm/Koonin bcc static =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 2)    
    call scr_and_log(str='Rate Schramm/Koonin bcc wigner-sietz =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 3)    
    call scr_and_log(str='Rate Schramm/Koonin bcc relaxed =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 4)    
    call scr_and_log(str='Rate Schramm/Koonin fcc static =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 5)    
    call scr_and_log(str='Rate Schramm/Koonin fcc wigner-sietz =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    rate = pycnoRate(A1, A2, Z1, Z2, rho, Rstep, partition, nucIntType, .FALSE., 6)    
    call scr_and_log(str='Rate Schramm/Koonin fcc relaxed =',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)

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

    call scr_and_log_str ('TEST: mean_wt_001:')

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

    call scr_and_log_str ('TEST: beta_excite_001:')

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

    call scr_and_log_str ('TEST: inv_len_001:')

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

subroutine rate_temp_adjust_002

    integer         :: A1
    integer         :: Z1
    real(kind=dbl)  :: X1
    integer         :: A2
    integer         :: Z2
    real(kind=dbl)  :: X2
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: temp
    real(kind=dbl)  :: adjust
    real(kind=dbl)  :: work1, work2
    integer         :: i
    
    call scr_and_log_str ('TEST: rate_temp_adjust_002:')

    A1 = 12
    Z1 = 6
    A2 = 12
    Z2 = 6
    X1 = 0.5_dbl
    X2 = 0.5_dbl
    rho = 33337662.176436048_dbl
    temp = 39317201.231798239_dbl

    adjust = tempRateAdjust(A1, Z1, X1, A2, Z2, X2, rho, temp)

    call scr_and_log(str='C-C 50-50 rho    = ',nbr=rho,fmt='(ES13.5)',lf=.TRUE.)
    call scr_and_log(str='          temp   = ',nbr=temp,fmt='(ES13.5)',lf=.TRUE.)
    call scr_and_log(str='          adjust = ',nbr=adjust,fmt='(ES13.5)',lf=.TRUE.)


end subroutine rate_temp_adjust_002

subroutine rate_temp_adjust_003

    integer         :: A1
    integer         :: Z1
    real(kind=dbl)  :: X1
    integer         :: A2
    integer         :: Z2
    real(kind=dbl)  :: X2
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: temp, temp_min, temp_max, delta_temp
    real(kind=dbl)  :: adjust
    real(kind=dbl)  :: work1, work2
    integer         :: i,j
    integer         :: n = 50
    
    call scr_and_log_str ('TEST: rate_temp_adjust_003')

    A1 = 12
    Z1 = 6
    A2 = 12
    Z2 = 6
    X1 = 0.5_dbl
    X2 = 0.5_dbl

    rho = 1d9
    temp_min = 2.6d7
    temp_max = 0.7d9
    delta_temp = (temp_max-temp_min)/n
    print *,"temp,adjust"
    do i = 0,n
        temp = temp_min + delta_temp*i
        adjust = tempRateAdjust(A1, Z1, X1, A2, Z2, X2, rho, temp)
        print *,temp,",",adjust    
    end do

    rho = 1d10
    temp_min = 6.65d7
    !temp_max = 0.7d9
    !delta_temp = (temp_max-temp_min)/n
    print *,"temp,adjust"
    do i = 0,n
        temp = temp_min + delta_temp*i
        adjust = tempRateAdjust(A1, Z1, X1, A2, Z2, X2, rho, temp)
        print *,temp,",",adjust   
    end do

    rho = 1d11
    temp_min = 14.25d7
    !temp_max = 0.7d9
    !delta_temp = (temp_max-temp_min)/n
    print *,"temp,adjust"
    do i = 0,n
        temp = temp_min + delta_temp*i
        adjust = tempRateAdjust(A1, Z1, X1, A2, Z2, X2, rho, temp)
        print *,temp,",",adjust   
    end do
    
    print *,""
    print *,"========================"
    rho = 34292563.783826306
    temp = 39909012.627657734
    adjust = tempRateAdjust(A1, Z1, X1, A2, Z2, X2, rho, temp, .true.)
    print *,adjust
    adjust = tempRateAdjust(A1, Z1, X1, A2, Z2, X2, rho, temp, .false.)
    print *,adjust
    print *,"========================"
    print *,""


end subroutine rate_temp_adjust_003

subroutine tridag_001

    use mathlinearalg
    
	integer, parameter :: NP=20
	integer                         :: k,n
	real(kind=dbl), dimension(NP)   :: diag,superd,subd,rhs,u
	character                       :: txt*3

    call scr_and_log_str ('TEST: tridag_001:')

	open(7,file='testdata/matrx2.dat',status='old')
	do
		read(7,'(a3)') txt
		if (txt == 'END') exit
		read(7,*)
		read(7,*) n
		read(7,*)
		read(7,*) (diag(k), k=1,n)
		read(7,*)
		read(7,*) (superd(k), k=1,n-1)
		read(7,*)
		read(7,*) (subd(k), k=2,n)
		read(7,*)
		read(7,*) (rhs(k), k=1,n)
!	carry out solution
!	use tridag_par since it also tests tridag_ser!
		call tridag_par(subd(2:n),diag(1:n), &
			superd(1:n-1),rhs(1:n),u(1:n))
		write(*,*) 'The solution vector is:'
		write(*,'(1x,6f12.6)') (u(k), k=1,n)
!	test solution
		write(*,*) '(matrix)*(sol''n vector) should be:'
		write(*,'(1x,6f12.6)') (rhs(k), k=1,n)
		write(*,*) 'Actual result is:'
		rhs(1)=diag(1)*u(1) + superd(1)*u(2)
		rhs(2:n-1)=subd(2:n-1)*u(1:n-2) + diag(2:n-1)*u(2:n-1) +&
					superd(2:n-1)*u(3:n)
		rhs(n)=subd(n)*u(n-1) + diag(n)*u(n)
		write(*,'(1x,6f12.6)') (rhs(k), k=1,n)
		write(*,*) '***********************************'
		write(*,*) 'Next problem...'
	end do
	close(7)

end subroutine tridag_001


subroutine spline_001
    
    use mathutils
    use mathinterpolation
    
    implicit none
    
	integer, parameter              :: N=40
	integer                         :: i
	real(kind=dbl)                  :: yp1,ypn
	real(kind=dbl), dimension(N)    :: x,y,y2

    call scr_and_log_str ('TEST: spline_001:')
	
	write(*,*) 'Second-derivatives for sin(x) from 0 to PI'
!	generate array for interpolation
	! x(1:N) = arth_d(1.0_dbl,1.0_dbl,N)*PI/N
	do i = 1,N
	    x(i) = (i-1)*2.0_dbl*PI/(N-1)
	end do
	y(:) = sin(x(:))
!	calculate 2nd derivative with SPLINE
	yp1 = cos(x(1))
	ypn = cos(x(N))
	call spline(x,y,yp1,ypn,y2)
!	test result
	write(*,'(t19,a,t35,a)') 'spline','actual'
	write(*,'(t6,a,t17,a,t33,a)') 'angle','2nd deriv','2nd deriv'
	do i = 1,N
		write(*,'(1x,f8.2,2f16.6)') x(i),y2(i),-sin(x(i))
	end do

end subroutine spline_001

subroutine spline_002

    use mathinterpolation

    implicit none

	integer, parameter :: NP=40
	integer :: i,nfunc
	real(kind=dbl)                  :: f,x,y,yp1,ypn
	real(kind=dbl), dimension(NP)   :: xa,ya,y2
	
	do nfunc=1,2
		if (nfunc == 1) then
			write(*,*) 'Sine function from 0 to 2*PI'
			!xa(1:NP)=arth(1,1,NP)*PI/NP
            do i = 1,NP
                xa(i) = (i-1)*2.0_dbl*PI/(NP-1) 
            end do
			ya(:)=sin(xa(:))
			yp1=cos(xa(1))
			ypn=cos(xa(NP))
		else if (nfunc == 2) then
			write(*,*) 'Exponential function from 0 to 1'
			!xa(1:NP)=arth(1.0_sp,1.0_sp,NP)/NP
            do i = 1,NP
                xa(i) = (i-1)*1.0/(NP-1) 
            end do
			ya(:)=exp(xa(:))
			yp1=exp(xa(1))
			ypn=exp(xa(NP))
		else
			stop
		end if
!	call SPLINE to get second derivatives
		call spline(xa,ya,yp1,ypn,y2)
!	call SPLINT for interpolations
		write(*,'(1x,t10,a1,t20,a4,t28,a13)') 'x','f(x)','interpolation'
		do i=1,10
			if (nfunc == 1) then
				x=(-0.05_dbl+i/10.0_dbl)*2.0*PI
				f=sin(x)
			else if (nfunc == 2) then
				x=-0.05_dbl+i/10.0_dbl
				f=exp(x)
			end if
			y=splint(xa,ya,y2,x)
			write(*,'(1x,3f18.12)') x,f,y
		end do
		write(*,*) '***********************************'
	end do

end subroutine spline_002

subroutine pycnoRateG05_CC_001

    real(kind=dbl)  :: T
    real(kind=dbl)  :: rho
    real(kind=dbl)  :: x12
    real(kind=dbl)  :: eps
    real(kind=dbl)  :: deps_dT
    real(kind=dbl)  :: deps_dRho
    real(kind=dbl)  :: rpyc
    
    
    T   = 9.2569312493209541D+07
    rho = 4.3765680728498936D+09
    x12 = 3.5666147092063111D-01
    rpyc = 0.0_dbl
    call pycnoRateG05_CC(T, rho, x12, eps, deps_dT, deps_drho, rpyc)
    write (*,*) 'T=',T,',rho=',rho,',Rpyc=',rpyc

    T   = 1.1209737001273090D+06
    rho = 1.4409212426081640D+10
    x12 = 3.5189091992105626D-01
    rpyc = 0.0_dbl
    call pycnoRateG05_CC(T, rho, x12, eps, deps_dT, deps_drho, rpyc)
    write (*,*) 'T=',T,',rho=',rho,',Rpyc=',rpyc

    T   = 1.0_dbl
    rho = 1.4409212426081640D+10
    x12 = 3.5189091992105626D-01
    rpyc = 0.0_dbl
    call pycnoRateG05_CC(T, rho, x12, eps, deps_dT, deps_drho, rpyc)
    write (*,*) 'T=',T,',rho=',rho,',Rpyc=',rpyc

    T   = 1.0D+9
    rho = 1.4409212426081640D+10
    x12 = 3.5189091992105626D-01
    rpyc = 0.0_dbl
    call pycnoRateG05_CC(T, rho, x12, eps, deps_dT, deps_drho, rpyc)
    write (*,*) 'T=',T,',rho=',rho,',Rpyc=',rpyc

end subroutine pycnoRateG05_CC_001

subroutine ndiff_001

    use mathdiff
    
    real(kind=dbl)  :: f1, f2, f3, f4, h
    real(kind=dbl)  :: x, ndiff
     
    integer ::  i

    call scr_and_log_str ('TEST: ndiff_001:')

    x = pi/6.0_dbl
    
    
    do i=1,10
        h = 1.0_dbl/(10.0_dbl**i)
        print *, "h=", h
        !print *, "x=", x
        !print *, "f(x)=", sin(x)
        !print *, "dfdx(x)=", cos(x)
        f1 = sin(x-2.0_dbl*h)
        f2 = sin(x-h)
        f3 = sin(x+h)
        f4 = sin(x+2.0_dbl*h)
        ndiff = ndiff5pt(f1, f2, f3, f4, h)
        print *, "ndiff(x)=", ndiff
        print *, "err=", ndiff-cos(x)
    end do
    
end subroutine ndiff_001

subroutine screening_potential_001

    real(kind=dbl) :: r
    real(kind=dbl) :: V_screen1
    real(kind=dbl) :: V_screen2
    real(kind=dbl) :: V

    call scr_and_log_str('TEST: screening_potential_001:')

    V_screen1 = -2.0_dbl    
    V_screen2 = -0.5_dbl    

    ! Test 1: at r_min = 0.1 fm
    r = 0.1_dbl
    V = screening_potential(r, V_screen1, V_screen2)
    call scr_and_log(str='r = 0.1 fm, V =', nbr=V, fmt='(ES13.5)', lf=.TRUE.)

    ! Test 2: at r_max = 6.0 fm
    r = 6.0_dbl
    V = screening_potential(r, V_screen1, V_screen2)
    call scr_and_log(str='r = 6.0 fm, V =', nbr=V, fmt='(ES13.5)', lf=.TRUE.)

    ! Test 3: mid-point to check for linearity
    r = 3.05_dbl
    V = screening_potential(r, V_screen1, V_screen2)
    call scr_and_log(str='r = 3.05 fm, V =', nbr=V, fmt='(ES13.5)', lf=.TRUE.)

end subroutine screening_potential_001

end module tests


