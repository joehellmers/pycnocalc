!
!	system_functions.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

module system_functions

use constants
use configuration

implicit none

save


contains 




!******************************************
!
! Subroutine to print a banner
!
!******************************************

subroutine printintro

use logging

implicit none

call scr_and_log (str='This is PycnoCalc - Pycnonulear Reaction Calculation System',fmt='(a)')

call scr_and_log(str=' ',fmt='(a)')
call scr_and_log(str='*****************************************************************',fmt='(a)')
call scr_and_log(str='  A Collaboration of: Dr. Fridolin Weber (fweber@mail.sdsu.edu)',fmt='(a)')
call scr_and_log(str='                      Joe Hellmers (hellmersjl@icloud.com)',fmt='(a)')
call scr_and_log(str='                      Barbara Golf (bsouci@gmail.com)',fmt='(a)')
call scr_and_log(str='                      Whitney Ryan (whitney.currier@gmail.com)',fmt='(a)')
call scr_and_log(str='*****************************************************************',fmt='(a)')
call scr_and_log(str=' ',fmt='(a)')

end subroutine printintro

!******************************************
!
! Do Some Sample Zero-Temp Calculations
!
!******************************************

subroutine sample_zero_temp_calcs

use astrophysics
use logging


real(kind=dbl) :: work1, work2


call scr_and_log_println


call scr_and_log_str('Zero Temp Reaction Rates')
call scr_and_log_str('------------------------')
call scr_and_log_str('Inside a White Dwarf')
call scr_and_log_str('C-C')
call react_rate_zero_temp(rho_wd, 12, 6, 4.0e16_dbl, work1, work2)
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('O-O')
call react_rate_zero_temp(rho_wd, 16, 8, 1.0e27_dbl, work1, work2)
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('Inside a Neutron Star')
call scr_and_log_str('C-C')
call react_rate_zero_temp(rho_ns, 12, 6, 4.0e16_dbl, work1, work2)
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('O-O')
call react_rate_zero_temp(rho_ns, 16, 8, 1.0e27_dbl, work1, work2)
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)

call scr_and_log_str('Neutron Star')
call scr_and_log_str('------------')
call scr_and_log_str('C-C')
call react_rate_zero_temp(rho_ns_surface, 12, 6, 4.0e16_dbl, work1, work2)
call scr_and_log_str('Surface')
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call react_rate_zero_temp(rho_ns_outer_crust, 12, 6, 4.0e16_dbl, work1, work2)
call scr_and_log_str('Outer Crust')
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call react_rate_zero_temp(rho_ns_inner_crust, 12, 6, 4.0e16_dbl, work1, work2)
call scr_and_log_str('Inner Crust')
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('O-O')
call react_rate_zero_temp(rho_ns_surface, 16, 8, 1.0e27_dbl, work1, work2)
call scr_and_log_str('Surface')
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call react_rate_zero_temp(rho_ns_outer_crust, 16, 8, 1.0e27_dbl, work1, work2)
call scr_and_log_str('Outer Crust')
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)
call react_rate_zero_temp(rho_ns_inner_crust, 16, 8, 1.0e27_dbl, work1, work2)
call scr_and_log_str('Inner Crust')
call scr_and_log_str('    Lower: ',lf=.FALSE.)
call scr_and_log(nbr=work1,fmt='(ES13.5)',lf=.TRUE.)
call scr_and_log_str('    Upper: ',lf=.FALSE.)
call scr_and_log(nbr=work2,fmt='(ES13.5)',lf=.TRUE.)

end subroutine sample_zero_temp_calcs 

!******************************************
!
! Do Some Sample Folding potential calculations
!
!******************************************


subroutine sample_folding_potential_calcs

use constants
use folding_potential
use utilities
use logging
use general_nuclear
use mathintegration

implicit none

! Local Variables
real(4) :: t1
integer :: r_iterator
real(kind=dbl) :: this_r
real(kind=dbl) :: this_vfold
real(kind=dbl) :: int_result

real(kind=dbl)		:: min_r, max_r, delta_r
integer			    :: A1, A2, Z1, Z2, n
real(kind=dbl)		:: diffuse1, diffuse2
real(kind=dbl)		:: rho0_1, rho0_2
real(kind=dbl)		:: tot_radius1, tot_radius2
integer			    :: unit,ierror
character(len=1)	:: delimiter
character(len=80)	:: outputfile
real(kind=dbl) 		:: start, finish, tot_time
real(kind=dbl)      :: total_walltime

integer count_0, count_1, count_rate, count_max
integer nn_int_type

t1 = secnds(0.0);

unit = 150

min_r		= ConvertStrToReal(getParamValue('FoldingSimple','r1'))
max_r		= ConvertStrToReal(getParamValue('FoldingSimple','r2'))
delta_r		= ConvertStrToReal(getParamValue('FoldingSimple','dr'))
A1		= ConvertStrToInt(getParamValue('FoldingSimple','A1')) 
A2		= ConvertStrToInt(getParamValue('FoldingSimple','A2')) 
Z1		= ConvertStrToInt(getParamValue('FoldingSimple','Z1')) 
Z2		= ConvertStrToInt(getParamValue('FoldingSimple','Z2')) 
n		= ConvertStrToInt(getParamValue('FoldingSimple','n')) 
diffuse1	= ConvertStrToReal(getParamValue('FoldingSimple','diffuse1')) 
diffuse2	= ConvertStrToReal(getParamValue('FoldingSimple','diffuse2')) 

! If the total radii are given, use them, otherwise calculate
 
if (checkParam('FoldingSimple','tot_radius1')) then
	tot_radius1	= ConvertStrToReal(getParamValue('FoldingSimple','tot_radius1')) + diffuse1
else
	tot_radius1 = nuclear_radius(real(A1,dbl))
end if

if (checkParam('FoldingSimple','tot_radius2')) then
	tot_radius2	= ConvertStrToReal(getParamValue('FoldingSimple','tot_radius2')) + diffuse2
else
	tot_radius2 = nuclear_radius(real(A2,dbl))
end if

! if the central densities are given use them, otherwise calculate

if (checkParam('FoldingSimple','rho0_1')) then
	rho0_1	= ConvertStrToReal(getParamValue('FoldingSimple','rho0_1'))
else
	rho0_1 = rho0_2pF(real(A1,dbl), tot_radius1, diffuse1)
end if

if (checkParam('FoldingSimple','rho0_2')) then
	rho0_2	= ConvertStrToReal(getParamValue('FoldingSimple','rho0_2'))
else
	rho0_2 = rho0_2pF(real(A2,dbl), tot_radius2, diffuse2)
end if

if (checkParam('FoldingSimple','nuc_interaction_type')) then
	nn_int_type	= ConvertStrToInt(getParamValue('FoldingSimple','nuc_interaction_type'))
else
	nn_int_type = 1
end if

print *, "Nuclear Interaction type: ", nn_int_type

!print *, 'Total Radius 1 = ', tot_radius1
!print *, 'Total Radius 2 = ', tot_radius2
!print *, 'Central Density 1 = ', rho0_1
!print *, 'Central Density 2 = ', rho0_2
!stop


delimiter	= getParamValue('FoldingSimple','delimiter')
outputfile	= getParamValue('FoldingSimple','outputfile')




open(unit, file='./' // results_dir // '/' // trim(outputfile), status='REPLACE', ACTION='WRITE', iostat=ierror)
if (ierror .NE. 0) then
	print *, 'Error opening data file for output'
	print *, ierror
	stop
end if



this_r = min_r

call system_clock(count_0, count_rate, count_max)
call cpu_time(start)

do 
	if (this_r > max_r) then
		exit
	end if
	this_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,nn_int_type,.FALSE.)
	call scr_and_log(nbr = this_r, fmt='(F10.5)')
    write (unit,*) this_r,delimiter,this_vfold,delimiter, log10(-1.0_dbl*this_vfold)
	this_r = this_r + delta_r
	
end do
call cpu_time(finish)
call system_clock(count_1, count_rate, count_max)
tot_time = finish - start
total_walltime = (count_1 * 1.0 / count_rate)/1000
call scr_and_log (str='Wall Time: ',nbr=total_walltime,fmt='(ES13.5)',lf=.true.)
call scr_and_log (str='Total Time: ',nbr=tot_time,fmt='(ES13.5)',lf=.true.)

close(unit)

end subroutine sample_folding_potential_calcs

!******************************************
!
! Generate C-C rates for a range of densities
!
!******************************************


subroutine genCCRates

use rate_calc

    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: initial_rho
    real(kind=dbl)  :: final_rho
    real(kind=dbl)  :: current_rho 
    real(kind=dbl)  :: delta_rho
    real(kind=dbl)  :: Rstep = 0.1 _dbl
    integer         :: partition = 15
    integer         :: nucIntType = 1
    character(len=1) :: nucIntTypeStr
    integer         :: i
    integer         :: rateCalcTypeIndex
    real(kind=dbl)  :: rate
    integer         :: N = 100 ! Number of intervals

    character(len=1)	:: delimiter
    character(len=120)	:: outputfile
    integer			    :: unit,ierror
    character(len=20)   :: nucIntTypeDesc
    character(len=20)   :: rxnDesc
    
    call scr_and_log_str ('CALCULATION: genCCRates:')

    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl

    initial_rho = 1000000000.0_dbl
    final_rho = 100000000000.0_dbl
    delta_rho = (final_rho-initial_rho)/N
    
    nucIntTypeStr = getParamValue('genCCRates','nucIntType')
    read(nucIntTypeStr,*) nucIntType
    delimiter	= getParamValue('genCCRates','delimiter')
    outputfile	= getParamValue('genCCRates','outputfile')

    select case (nucIntType)
        case (1)
            nucIntTypeDesc = 'SAOPAULO'
        case (2)
            nucIntTypeDesc = 'M3Y'
        case (3)
            nucIntTypeDesc = 'RMF'
        case default
            nucIntTypeDesc = 'M3Y'
    end select
    
    open(unit, file='./' // results_dir // '/' // trim(outputfile), status='REPLACE', ACTION='WRITE', iostat=ierror)
    if (ierror .NE. 0) then
        print *, 'Error opening data file for output'
        print *, ierror
        stop
    end if
    write(unit,*) 'nuclei,nninteraction,density,pycno_rate'
    call scr_and_log(str='nuclei,nninteraction,density,pycno_rate',fmt='(a)',lf=.TRUE.)
    ! TEMP: only do fcc rates
    do rateCalcTypeIndex = 4, 6
        do i = 1, N+1
            current_rho = initial_rho + delta_rho*(i-1)
            call scr_and_log(str='C-C,',fmt='(a)',lf=.FALSE.)
            if (rateCalcTypeIndex .eq. 0) then
                rxnDesc = trim(nucIntTypeDesc) // '-SPVH'
            end if
            if (rateCalcTypeIndex .eq. 1) then
                rxnDesc = trim(nucIntTypeDesc) // '-bcc-static'
            end if
            if (rateCalcTypeIndex .eq. 2) then
                rxnDesc = trim(nucIntTypeDesc) // '-bcc-ws'
            end if
            if (rateCalcTypeIndex .eq. 3) then
                rxnDesc = trim(nucIntTypeDesc) // '-bcc-relax'
            end if
            if (rateCalcTypeIndex .eq. 4) then
                rxnDesc = trim(nucIntTypeDesc) // '-fcc-static'
            end if
            if (rateCalcTypeIndex .eq. 5) then
                rxnDesc = trim(nucIntTypeDesc) // '-fcc-ws'
            end if
            if (rateCalcTypeIndex .eq. 6) then
                rxnDesc = trim(nucIntTypeDesc) // '-fcc-relax'
            end if
            rxnDesc = trim(rxnDesc)             
            if (nucIntType .eq. 1) then        
                call scr_and_log(str=rxnDesc,fmt='(a)',lf=.FALSE.)
            else
                call scr_and_log(str=rxnDesc,fmt='(a)',lf=.FALSE.)
            end if
            call scr_and_log(str='',nbr=current_rho,fmt='(ES13.5)',lf=.FALSE.)
            rate = pycnoRate(A1, A2, Z1, Z2, current_rho, Rstep, partition, nucIntType, .FALSE., rateCalcTypeIndex)
            write (unit,*) "C-C",delimiter,rxnDesc,delimiter,current_rho,delimiter,rate
            call scr_and_log(str=',',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
        end do
    end do
    
    close(unit)
    
                
end subroutine genCCRates	

subroutine graphM3Y

    use nucleon_interactions
    
    real(kind=dbl)  :: delta_r, r
    integer(kind=8) :: i
    integer(kind=8) :: N
    real(kind=dbl)  :: nn_potenergy
    real(kind=dbl)  :: r_bottom, r_top
    logical         :: change_pt_found = .false.
    logical         :: excludeCore = .true.
    real(kind=dbl)  :: r_start, r_end
    
    print *, "r,m3y_v_nn"
    N = 100
    r_start = 0.0_dbl
    r_end = 4.0_dbl
    delta_r = (r_end - r_start)/N
    do i = 0, N
        r = i*delta_r
        nn_potenergy = nucleonM3Y(r, 0.001_dbl, excludeCore)
        print *, r , ",", nn_potenergy
        !if (.not. (change_pt_found) .and. (nn_potenergy .lt. 0.0_dbl)) then
        !    change_pt_found = .true.
        !    r_top = r
        !    r_bottom = r - delta_r
        !end if
    end do
    
    ! Use this code to get the "root", where the interaction energy goes negative
    !if (change_pt_found) then
    !    change_pt_found = .false.
    !    N = 100
    !    r_bottom    = 0.56754316782759151_dbl
    !    r_top       = 0.56754316782759207_dbl
    !    delta_r = (r_top - r_bottom)/N        
    !    do i = 0, N
    !        r = r_bottom + i*delta_r
    !        nn_potenergy = nucleonM3Y(r,0.001_dbl)
    !        write (*,"(F100.98,2x,F100.98)") r, nn_potenergy
    !        !if (.not. (change_pt_found) .and. (nn_potenergy .le. 0.0_dbl)) then
            !    change_pt_found = .true.
            !    r_top = r
            !    r_bottom = r - delta_r
            !    !write (*,"(a,f100.98,1x,a,f100.98)") "r_bottom = ", r_bottom, "r_top = ", r_top
            !    !exit
            !end if
     !   end do
    !end if
    
    !print *, nucleonM3Y(0.56754316782759184345508174374117515981197357177734375_dbl,0.001_dbl)
    
end subroutine

subroutine graphSaoPaulo

    use nucleon_interactions
    use astrophysics
    use general_nuclear
    
    real(kind=dbl)  :: delta_r, r
    integer(kind=8) :: i
    integer(kind=8) :: N
    real(kind=dbl)  :: nn_potenergy_min, nn_potenergy_max
    real(kind=dbl)  :: r_bottom, r_top
    logical         :: change_pt_found = .false.
    real(kind=dbl)  :: r_start, r_end
    real(kind=dbl)  :: E0_min, E0_max
    real(kind=dbl)  :: rho_min = 1d7
    real(kind=dbl)  :: rho_max = 1d9
    integer         :: A1, Z1
    real(kind=dbl)  :: A2, Z2
    real(kind=dbl)  :: mu
    
    A1 = 12
    Z1 = 6
    A2 = 12.0_dbl
    Z2 = 12.0_dbl
    
    
    print *, "r,sao_paulo_v_nn_min,sao_paulo_v_nn_max"
    N = 100
    r_start = 0.0_dbl
    r_end = 4.0_dbl
    E0_min = E0_Energy(Z1,Z2, A1, A2, rho_min)
    E0_max = E0_Energy(Z1,Z2, A1, A2, rho_max)
    mu = reduced_mass(A1, Z1, A2, Z2)
    delta_r = (r_end - r_start)/N
    do i = 0, N
        r = i*delta_r
        nn_potenergy_min = nucleonSaoPaulo(r, E0_min, mu)
        nn_potenergy_max = nucleonSaoPaulo(r, E0_max, mu)
        print *, r , ",", nn_potenergy_min, nn_potenergy_max
        !if (.not. (change_pt_found) .and. (nn_potenergy .lt. 0.0_dbl)) then
        !    change_pt_found = .true.
        !    r_top = r
        !    r_bottom = r - delta_r
        !end if
    end do
    
    ! Use this code to get the "root", where the interaction energy goes negative
    !if (change_pt_found) then
    !    change_pt_found = .false.
    !    N = 100
    !    r_bottom    = 0.56754316782759151_dbl
    !    r_top       = 0.56754316782759207_dbl
    !    delta_r = (r_top - r_bottom)/N        
    !    do i = 0, N
    !        r = r_bottom + i*delta_r
    !        nn_potenergy = nucleonM3Y(r,0.001_dbl)
    !        write (*,"(F100.98,2x,F100.98)") r, nn_potenergy
    !        !if (.not. (change_pt_found) .and. (nn_potenergy .le. 0.0_dbl)) then
            !    change_pt_found = .true.
            !    r_top = r
            !    r_bottom = r - delta_r
            !    !write (*,"(a,f100.98,1x,a,f100.98)") "r_bottom = ", r_bottom, "r_top = ", r_top
            !    !exit
            !end if
     !   end do
    !end if
    
    !print *, nucleonM3Y(0.56754316782759184345508174374117515981197357177734375_dbl,0.001_dbl)
    
end subroutine

subroutine graphAllNN

    use nucleon_interactions
    use astrophysics
    use general_nuclear
    
    real(kind=dbl)  :: delta_r, r
    integer(kind=8) :: i
    integer(kind=8) :: N
    real(kind=dbl)  :: sp_potenergy_min, sp_potenergy_max
    real(kind=dbl)  :: m3y_potenergy
    real(kind=dbl)  :: rmf_potenergy
    real(kind=dbl)  :: r_bottom, r_top
    logical         :: change_pt_found = .false.
    real(kind=dbl)  :: r_start, r_end
    real(kind=dbl)  :: E0_min, E0_max
    real(kind=dbl)  :: rho_min = 1d7
    real(kind=dbl)  :: rho_max = 1d9
    integer         :: A1, Z1
    real(kind=dbl)  :: A2, Z2
    real(kind=dbl)  :: mu
    logical         :: excludeCore = .true.
    real(kind=dbl)  :: coreCutoff = 0.01_dbl

    A1 = 12
    Z1 = 6
    A2 = 12.0_dbl
    Z2 = 12.0_dbl
        
    print *, "Radius(fm),Vnn(SP_min),Vnn(SP_max),Vnn(M3Y),Vnn(RMF)"
    N = 100
    r_start = 0.3_dbl
    r_end = 2.5_dbl
    E0_min = E0_Energy(Z1,Z2, A1, A2, rho_min)
    E0_max = E0_Energy(Z1,Z2, A1, A2, rho_max)
    mu = reduced_mass(A1, Z1, A2, Z2)
    delta_r = (r_end - r_start)/N
    do i = 0, N
        r = r_start + i*delta_r
        sp_potenergy_min = nucleonSaoPaulo(r, E0_min, mu)
        sp_potenergy_max = nucleonSaoPaulo(r, E0_max, mu)
        m3y_potenergy = nucleonM3Y(r, 0.001_dbl, excludeCore, coreCutoff)
        rmf_potenergy = nucleonRMF(r,4) ! use L1 type of RMF nn-interaction
        print *, r , ",", sp_potenergy_min,",", sp_potenergy_max,",", m3y_potenergy,",",rmf_potenergy
    end do
        
end subroutine

subroutine genPycnoRateSplines
    
    ! Here we assume each rate type has 101 entries, and that they are in the csv file
    ! in the proper order
    
    use logging
    use mathinterpolation
    use rate_calc
    use csv_module

    type(csv_file) :: f
    
    integer, parameter              :: N = 101, M = 7      
    real(kind=dbl),  dimension(N)   :: densities, rates, rate_2nd_derivs
    integer			                :: unit,ierror
    real(kind=dbl)                  :: real_ierror
    character(:), allocatable       :: ratefile        
    logical                         :: firstLine = .true.
    real(kind=dbl)                  :: density, rate, rate_prime1, rate_primeN
    integer                         :: counter = 0
    integer                         :: i, j
    character(:), allocatable       :: species, nninteraction
    real(kind=dbl)                  :: thisrate, thisdensity
 
    character(len=60),dimension(:),allocatable :: rxnspecies
    character(len=60),dimension(:),allocatable :: rxntypes
    real(kind=dbl),dimension(:),allocatable    :: rxndensities
    real(kind=dbl),dimension(:),allocatable    :: rxnrates
     
    logical :: status_ok
    character(len=60),dimension(:),allocatable :: header 
    integer,dimension(:),allocatable :: itypes 
 
    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: Rstep = 0.1 _dbl
    integer         :: partition = 15
    integer         :: nucIntType
      ! 1 - M3Y
      ! 2 - Sao Paulo
    
    character(len=1) :: nucIntTypeStr
    character(len=120)	:: splineFile, densitiesFile, ratesFile, secondDerivsFile
    character(len=120)	:: outputfileStub
    character(len=120)	:: outputfileRaw
    character(len=120)	:: fileSuffix
    character(len=10)	:: fileSuffixNN
    integer             :: currIndex
    
    
    real(kind=dbl), dimension(10) :: A
    
    A = (/ &
    1D-10, 2.3D-14, 3.0D0, 0.6666D99, 573D14, &
    6.0D0, 7.0D0, 8.0D0, 9.0D0, 0.0D0 &
    /)

! Get the 
    nucIntTypeStr = getParamValue('genPycnoRateSplines','nucIntType')
    outputfileStub = getParamValue('genPycnoRateSplines','splinefile')
    read(nucIntTypeStr,*) nucIntType
    call scr_and_log_str('Spline Output File Stub: ' // trim(outputfileStub))        


! First we need to load the arrays

    ! Instantiate the CSV file object
    if (nucIntType .eq. 2) then
        call scr_and_log_str('Generating Splines for WD M3Y Rates')        
        call f%read('researchdata/pycno_rates_C-C_WD_M3Y.csv',header_row=1,status_ok=status_ok)
        fileSuffixNN = '.m3y'
    else
        call scr_and_log_str('Generating Splines for WD Sao Paulo Rates')        
        call f%read('researchdata/pycno_rates_C-C_WD_SaoPaulo.csv',header_row=1,status_ok=status_ok)
        fileSuffixNN = '.sp'
    end if
    
    ! get the header and type info
    call f%get_header(header,status_ok)
    call f%variable_types(itypes,status_ok)
    call f%get(1,rxnspecies,status_ok)
    call f%get(2,rxntypes,status_ok)
    call f%get(3,rxndensities,status_ok)
    call f%get(4,rxnrates,status_ok)
    ! Remove the CSV file object
    call f%destroy()
    
    
    ! Do a sanity check on the file contents
	do i = 1, M
        select case (i)
            case (1)
                fileSuffix = trim(fileSuffixNN) // '.spvh'    
            case (2)
                fileSuffix = trim(fileSuffixNN) // '.bcc-static'    
            case (3)
                fileSuffix = trim(fileSuffixNN) // '.bcc-ws'    
            case (4)
                fileSuffix = trim(fileSuffixNN) // '.bcc-relax'    
            case (5)
                fileSuffix = trim(fileSuffixNN) // '.fcc-static'    
            case (6)
                fileSuffix = trim(fileSuffixNN) // '.fcc-ws'    
            case (7)
                fileSuffix = trim(fileSuffixNN) // '.fcc-relax'    
        end select
        call scr_and_log_str('....file suffix =' // trim(fileSuffix))     
        splineFile = 'spline' // trim(fileSuffix)
        densitiesFile = 'densities' // trim(fileSuffix)
        ratesFile = 'rates' // trim(fileSuffix)
        secondDerivsFile = 'secondDerivs' // trim(fileSuffix)
        call scr_and_log_str('....spline file    = ' // trim(splineFile))     
        call scr_and_log_str('....density file    = ' // trim(densitiesFile))     
        call scr_and_log_str('....rates file      = ' // trim(ratesFile))     
        call scr_and_log_str('....2nd Derivs file = ' // trim(secondDerivsFile))     
        ! first load the work arrays
        do j = 1, N
            currIndex = (i-1) * N + j
            densities(j) = rxndensities(currIndex)
            rates(j) = rxnrates(currIndex)
        end do
        
        ! Now compute the spline
        rate_prime1 = (rates(2)-rates(1))/(densities(2)-densities(1))
        rate_primeN = (rates(N)-rates(N-1))/(densities(N)-densities(N-1))
        call spline(densities,rates,rate_prime1,rate_primeN,rate_2nd_derivs)
        
        ! For Sanity Check calculate the values at the "known" densities and rates
        !do j = 1, N
	    !    thisrate=splint(densities,rates,rate_2nd_derivs,densities(j))
	    !    print *, densities(j), rates(j), rate_2nd_derivs(j), thisrate
        !end do

        ! Ouput the spline
        open(unit, file='./' // results_dir // '/' // trim(splineFile), status='REPLACE', ACTION='WRITE', iostat=ierror)
        if (ierror .NE. 0) then
	        print *, 'Error opening spline file for output'
	        print *, ierror
	        stop
        end if

	    do j = 1, N
            write (unit,*) densities(j), rates(j), rate_2nd_derivs(j)
        end do
        close(unit)

        ! Ouput the densities
        open(unit, file='./' // results_dir // '/' // trim(densitiesFile), status='REPLACE', ACTION='WRITE', iostat=ierror)
        if (ierror .NE. 0) then
            print *, 'Error opening densities file for output'
            print *, ierror
            stop
        end if

        do j = 1, N
            write (unit,fmt="(ES13.5,A)",advance="no") densities(j),","
            if (mod(j,5) .eq. 0) then
                write (unit,fmt="(AA)") " &"
            end if
        end do
        close(unit)

        open(unit, file='./' // results_dir // '/' // trim(ratesFile), status='REPLACE', ACTION='WRITE', iostat=ierror)
        if (ierror .NE. 0) then
            print *, 'Error opening rates file for output'
            print *, ierror
            stop
        end if

        do j = 1, N
            write (unit,fmt="(ES13.5,A)",advance="no") rates(j),","
            if (mod(j,5) .eq. 0) then
                write (unit,fmt="(AA)") " &"
            end if
        end do
        close(unit)

        open(unit, file='./' // results_dir // '/' // trim(secondDerivsFile), status='REPLACE', ACTION='WRITE', iostat=ierror)
        if (ierror .NE. 0) then
            print *, 'Error opening 2nd Derivatives file for output'
            print *, ierror
            stop
        end if

        do j = 1, N
            write (unit,fmt="(ES13.5,A)",advance="no") rate_2nd_derivs(j),","
            if (mod(j,5) .eq. 0) then
                write (unit,fmt="(AA)") " &"
            end if
        end do
        close(unit)

    end do
    
    return

! Do some spot checking
    
    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl

    thisdensity = 74260000000.000000_dbl + ((75250000000.000000_dbl - 74260000000.000000_dbl)/3.0_dbl)
    thisrate = splint(densities,rates,rate_2nd_derivs,thisdensity)
    !print *, thisdensity, thisrate, pycnoRate(A1, A2, Z1, Z2, thisdensity, Rstep, partition, nucIntType, .FALSE.)
    
    print *, A
    
end subroutine genPycnoRateSplines

!*************************************************
!
! Graph M3Y, Sao Paulo, and RMF Folding potentials
!
!*************************************************

subroutine graphAllFolding

use constants
use folding_potential
use utilities
use logging
use general_nuclear
use mathintegration

implicit none

! Local Variables
real(4) :: t1
integer :: r_iterator
real(kind=dbl) :: this_r
real(kind=dbl) :: sp_vfold, m3y_vfold, rmf_vfold
real(kind=dbl) :: int_result

real(kind=dbl)		:: min_r, max_r, delta_r
integer			    :: A1, A2, Z1, Z2, n
real(kind=dbl)		:: diffuse1, diffuse2
real(kind=dbl)		:: rho0_1, rho0_2
real(kind=dbl)		:: tot_radius1, tot_radius2
integer			    :: unit,ierror
character(len=1)	:: delimiter
character(len=80)	:: outputfile
real(kind=dbl) 		:: start, finish, tot_time
real(kind=dbl)      :: total_walltime

integer count_0, count_1, count_rate, count_max
integer nn_int_type

t1 = secnds(0.0);

unit = 150

min_r		= ConvertStrToReal(getParamValue('graphAllFolding','r1'))
max_r		= ConvertStrToReal(getParamValue('graphAllFolding','r2'))
delta_r		= ConvertStrToReal(getParamValue('graphAllFolding','dr'))
A1		= ConvertStrToInt(getParamValue('graphAllFolding','A1')) 
A2		= ConvertStrToInt(getParamValue('graphAllFolding','A2')) 
Z1		= ConvertStrToInt(getParamValue('graphAllFolding','Z1')) 
Z2		= ConvertStrToInt(getParamValue('graphAllFolding','Z2')) 
n		= ConvertStrToInt(getParamValue('graphAllFolding','n')) 
diffuse1	= ConvertStrToReal(getParamValue('graphAllFolding','diffuse1')) 
diffuse2	= ConvertStrToReal(getParamValue('graphAllFolding','diffuse2')) 

! If the total radii are given, use them, otherwise calculate
 
if (checkParam('graphAllFolding','tot_radius1')) then
	tot_radius1	= ConvertStrToReal(getParamValue('graphAllFolding','tot_radius1')) + diffuse1
else
	tot_radius1 = nuclear_radius(real(A1,dbl))
end if

if (checkParam('graphAllFolding','tot_radius2')) then
	tot_radius2	= ConvertStrToReal(getParamValue('graphAllFolding','tot_radius2')) + diffuse2
else
	tot_radius2 = nuclear_radius(real(A2,dbl))
end if

! if the central densities are given use them, otherwise calculate

if (checkParam('graphAllFolding','rho0_1')) then
	rho0_1	= ConvertStrToReal(getParamValue('graphAllFolding','rho0_1'))
else
	rho0_1 = rho0_2pF(real(A1,dbl), tot_radius1, diffuse1)
end if

if (checkParam('graphAllFolding','rho0_2')) then
	rho0_2	= ConvertStrToReal(getParamValue('graphAllFolding','rho0_2'))
else
	rho0_2 = rho0_2pF(real(A2,dbl), tot_radius2, diffuse2)
end if

if (checkParam('graphAllFolding','nuc_interaction_type')) then
	nn_int_type	= ConvertStrToInt(getParamValue('graphAllFolding','nuc_interaction_type'))
else
	nn_int_type = 1
end if

delimiter	= getParamValue('graphAllFolding','delimiter')
outputfile	= getParamValue('graphAllFolding','outputfile')

open(unit, file='./' // results_dir // '/' // trim(outputfile), status='REPLACE', ACTION='WRITE', iostat=ierror)
if (ierror .NE. 0) then
	print *, 'Error opening data file for output'
	print *, ierror
	stop
end if

this_r = min_r

call system_clock(count_0, count_rate, count_max)
call cpu_time(start)

write (unit,*) "Radius(fm)",delimiter,"Sao Paulo",delimiter,"M3Y",delimiter,"RMF"
do 
	if (this_r > max_r) then
		exit
	end if
	sp_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,1,.FALSE.)
	m3y_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,2,.FALSE.)
	rmf_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,3,.FALSE.)
	call scr_and_log(nbr = this_r, fmt='(F10.5)')
    write (unit,*) this_r,delimiter,sp_vfold,delimiter,m3y_vfold,delimiter,rmf_vfold
	this_r = this_r + delta_r
	
end do
call cpu_time(finish)
call system_clock(count_1, count_rate, count_max)
tot_time = finish - start
total_walltime = (count_1 * 1.0 / count_rate)/1000
call scr_and_log (str='Wall Time: ',nbr=total_walltime,fmt='(ES13.5)',lf=.true.)
call scr_and_log (str='Total Time: ',nbr=tot_time,fmt='(ES13.5)',lf=.true.)

close(unit)

end subroutine graphAllFolding

!*******************************************************
!
! Graph M3Y, Sao Paulo, and RMF Folding total potentials
!
!*******************************************************

subroutine graphAllTotPotential

use constants
use folding_potential
use utilities
use logging
use general_nuclear
use mathintegration
use astrophysics

implicit none

! Local Variables
integer :: r_iterator
real(kind=dbl) :: this_r
real(kind=dbl) :: sp_vfold, m3y_vfold, rmf_vfold
real(kind=dbl) :: Vcoul
real(kind=dbl) :: int_result

real(kind=dbl)		:: min_r, max_r, delta_r
integer			    :: A1, A2, Z1, Z2, n
real(kind=dbl)		:: diffuse1, diffuse2
real(kind=dbl)		:: rho0_1, rho0_2
real(kind=dbl)		:: tot_radius1, tot_radius2
integer			    :: unit,ierror
character(len=1)	:: delimiter
character(len=120)	:: outputfile
real(kind=dbl) 		:: start, finish, tot_time
real(kind=dbl)      :: total_walltime

integer count_0, count_1, count_rate, count_max
integer nn_int_type

unit = 150

min_r		= ConvertStrToReal(getParamValue('graphAllTotPotential','r1'))
max_r		= ConvertStrToReal(getParamValue('graphAllTotPotential','r2'))
delta_r		= ConvertStrToReal(getParamValue('graphAllTotPotential','dr'))
A1		= ConvertStrToInt(getParamValue('graphAllTotPotential','A1')) 
A2		= ConvertStrToInt(getParamValue('graphAllTotPotential','A2')) 
Z1		= ConvertStrToInt(getParamValue('graphAllTotPotential','Z1')) 
Z2		= ConvertStrToInt(getParamValue('graphAllTotPotential','Z2')) 
n		= ConvertStrToInt(getParamValue('graphAllTotPotential','n')) 
diffuse1	= ConvertStrToReal(getParamValue('graphAllTotPotential','diffuse1')) 
diffuse2	= ConvertStrToReal(getParamValue('graphAllTotPotential','diffuse2')) 

! If the total radii are given, use them, otherwise calculate
 
if (checkParam('graphAllTotPotential','tot_radius1')) then
	tot_radius1	= ConvertStrToReal(getParamValue('graphAllTotPotential','tot_radius1')) + diffuse1
else
	tot_radius1 = nuclear_radius(real(A1,dbl))
end if

if (checkParam('graphAllTotPotential','tot_radius2')) then
	tot_radius2	= ConvertStrToReal(getParamValue('graphAllTotPotential','tot_radius2')) + diffuse2
else
	tot_radius2 = nuclear_radius(real(A2,dbl))
end if

! if the central densities are given use them, otherwise calculate

if (checkParam('graphAllTotPotential','rho0_1')) then
	rho0_1	= ConvertStrToReal(getParamValue('graphAllTotPotential','rho0_1'))
else
	rho0_1 = rho0_2pF(real(A1,dbl), tot_radius1, diffuse1)
end if

if (checkParam('graphAllTotPotential','rho0_2')) then
	rho0_2	= ConvertStrToReal(getParamValue('graphAllTotPotential','rho0_2'))
else
	rho0_2 = rho0_2pF(real(A2,dbl), tot_radius2, diffuse2)
end if

delimiter	= getParamValue('graphAllTotPotential','delimiter')
outputfile	= getParamValue('graphAllTotPotential','outputfile')

open(unit, file='./' // results_dir // '/' // trim(outputfile), status='REPLACE', ACTION='WRITE', iostat=ierror)
if (ierror .NE. 0) then
	print *, 'Error opening data file for output'
	print *, trim(outputfile)
	print *, ierror
	stop
end if

this_r = min_r

write (unit,*) "Radius(fm)",delimiter,"Vcoul",delimiter,"Sao Paulo",delimiter,"M3Y",delimiter,"RMF-L1"
do 
	if (this_r > max_r) then
		exit
	end if
	sp_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,1,.FALSE.)
	m3y_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,2,.FALSE.)
	rmf_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,3,.FALSE.)	
	Vcoul = Vcoulomb(Z1,real(Z2,dbl),this_r,tot_radius1,tot_radius2)
	call scr_and_log(nbr = this_r, fmt='(F10.5)')
    write (unit,*) this_r,delimiter,Vcoul,delimiter,sp_vfold+Vcoul,delimiter,m3y_vfold+Vcoul,delimiter,rmf_vfold+Vcoul
	this_r = this_r + delta_r
	
end do

close(unit)

end subroutine graphAllTotPotential


end module system_functions

