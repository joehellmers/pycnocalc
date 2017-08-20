!
!	system_functions.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

module system_functions

use constants

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
use configuration
use utilities
use logging
use general_nuclear
use integration

implicit none

! Local Variables
real(4) :: t1
integer :: r_iterator
real(kind=dbl) :: this_r
real(kind=dbl) :: this_vfold
real(kind=dbl) :: int_result

real(kind=dbl)		:: min_r, max_r, delta_r
integer			:: A1, A2, Z1, Z2, n
real(kind=dbl)		:: diffuse1, diffuse2
real(kind=dbl)		:: rho0_1, rho0_2
real(kind=dbl)		:: tot_radius1, tot_radius2
integer			:: unit,ierror
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
	this_vfold = vfold_spherically_symmetric(this_r,A1,real(A2,dbl),Z1,real(Z2,dbl),n,diffuse1,diffuse2,rho0_1,rho0_2,tot_radius1,tot_radius2,0.85_dbl,.FALSE.)
	!this_vfold =  -1.0_dbl*vfold_spherically_symmetric(0.1_dbl*r_iterator,58,16,4.230848400_dbl,2.460993151_dbl,15,0.50_dbl,0.50_dbl,0.1606844454_dbl,0.1820414223_dbl,6.0_dbl,6.0_dbl)
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

	
end module system_functions

