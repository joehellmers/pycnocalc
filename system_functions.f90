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
	this_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,2,.FALSE.)
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
    real(kind=dbl)  :: rate
    integer         :: N = 100 ! Number of intervals

    character(len=1)	:: delimiter
    character(len=120)	:: outputfile
    integer			    :: unit,ierror
    character(len=20)   :: nucIntTypeDesc
    
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

    if (nucIntType .eq. 1) then
        nucIntTypeDesc = 'SAOPAULO'
    else
        nucIntTypeDesc = 'M3Y'
    end if
    
    open(unit, file='./' // results_dir // '/' // trim(outputfile), status='REPLACE', ACTION='WRITE', iostat=ierror)
    if (ierror .NE. 0) then
	    print *, 'Error opening data file for output'
	    print *, ierror
	    stop
    end if


    write(unit,*) 'nuclei,nn-interation,density,pycno_rate'
    call scr_and_log(str='nuclei,nn-interation,density,pycno_rate',fmt='(a)',lf=.TRUE.)
    do i = 1, N+1
        current_rho = initial_rho + delta_rho*(i-1)
        call scr_and_log(str='C-C,',fmt='(a)',lf=.FALSE.)
        if (nucIntType .eq. 1) then        
            call scr_and_log(str='SAOPAULO,',fmt='(a)',lf=.FALSE.)
        else
            call scr_and_log(str='M3Y,',fmt='(a)',lf=.FALSE.)
        end if
        call scr_and_log(str='',nbr=current_rho,fmt='(ES13.5)',lf=.FALSE.)
        rate = pycnoRate(A1, A2, Z1, Z2, current_rho, Rstep, partition, nucIntType, .FALSE.)
        write (unit,*) "C-C",delimiter,nucIntTypeDesc,delimiter,current_rho,delimiter,rate
        call scr_and_log(str=',',nbr=rate,fmt='(ES13.5)',lf=.TRUE.)
    end do
            
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

subroutine genPycnoRateSplines
    
    use logging
    use mathinterpolation
    use rate_calc

    integer, parameter              :: N = 101      
    real(kind=dbl),  dimension(N)   :: densities, rates, rate_2nd_derivs
    integer			                :: unit,ierror
    real(kind=dbl)                  :: real_ierror
    character(:), allocatable       :: ratefile        
    logical                         :: firstLine = .true.
    real(kind=dbl)                  :: density, rate, rate_prime1, rate_primeN
    integer                         :: counter = 0
    integer                         :: i
    character(:), allocatable       :: species, nninteraction
    real(kind=dbl)                  :: thisrate, thisdensity
 
    integer         :: A1
    real(kind=dbl)  :: A2
    integer         :: Z1
    real(kind=dbl)  :: Z2
    real(kind=dbl)  :: Rstep = 0.1 _dbl
    integer         :: partition = 15
    integer         :: nucIntType = 2

! First we need to load the arrays

    allocate(character(13) :: ratefile)
    ratefile = 'researchdata/pycnorates_cc_m3y.csv'
    unit = 199
    open (unit,file=ratefile,status='OLD',action='READ', iostat=ierror)
	if (ierror .NE. 0) then
			call scr_and_log_str('Cannot open rate file')
	        call scr_and_log_str('File: ',lf=.FALSE.)	
			call scr_and_log_str(ratefile)
			call scr_and_log_str('Error: ',lf=.FALSE.)
			real_ierror =real(ierror)
			call scr_and_log(nbr=real_ierror,fmt='(f5.0)')
			stop
	end if

	readloop: Do
	
        if (firstLine) then
		    firstLine = .false.
		    read(unit,*,iostat=ierror) ! Ignore the headers in the first line
        end if
        read(unit,*,iostat=ierror) species, nninteraction, density, rate
	    
	    if (ierror .EQ. -1) then
			! End of file
			exit
		end if
		
		if (ierror .GT. 0) then
				print *, ierror
				call scr_and_log_str('Error reading rate file: ' // ratefile)
				stop
		end if
		counter = counter + 1
		densities(counter) = density
		rates(counter) = rate
	end do readloop

    close(unit)

! Now we need to calculate the spline (i.e. 2nd Derivates)

    rate_prime1 = (rates(2)-rates(1))/(densities(2)-densities(1))
    rate_primeN = (rates(N)-rates(N-1))/(densities(N)-densities(N-1))
	call spline(densities,rates,rate_prime1,rate_primeN,rate_2nd_derivs)

! For Sanity Check calculate the values at the "known" densities and rates

	do i = 1, N
	    thisrate=splint(densities,rates,rate_2nd_derivs,densities(i))
	    print *, densities(i), rates(i), thisrate
    end do

! Do some spot checking
    
    A1 = 12
    A2 = 12.0_dbl
    Z1 = 6
    Z2 = 6.0_dbl

    thisdensity = 74260000000.000000_dbl + ((75250000000.000000_dbl - 74260000000.000000_dbl)/3.0_dbl)
    thisrate = splint(densities,rates,rate_2nd_derivs,thisdensity)
    print *, thisdensity, thisrate, pycnoRate(A1, A2, Z1, Z2, thisdensity, Rstep, partition, nucIntType, .FALSE.)
    
end subroutine 

end module system_functions

