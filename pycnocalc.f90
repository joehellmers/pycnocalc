!
!	PycnoCalc.f90
!	PycnoCalc
!
!	Joe Hellmers
!
!	References
!	SPVH1969: E. E. Salpeter and H. M. Van Horn, Astrophys. J. 155, 183 (1969)
!   	SPVH1967: H. M. Van Horn, E. E. Salpeter, Phys. Rev, 157, 751
!   	CPHRG1997: L. C. Chamon, B. V. Pereira, M. S. Hussein, M. A. Candido Ribeiro and D. Galetti, Phys. Rev. Lett. 79, 5218 (1997)
!	ST1983: S. L. Shapiro, S. A. Teukolsky "Black Holes, White Dwarfs, and Neutron Stars
!


program PycnoCalc

use logging
use configuration
use system_functions
use globalvars
use astrophysics
use folding_potential
use utilities
use cmdline
use tests

implicit none

integer i

call process_cmdline
call initialize_logging
call load_config
call printintro

call initialize_globals
do i=1, function_cnt

	if (function_list(i)%function_name .EQ. 'LogConfig') then
		call LogConfig
	end if

	if (function_list(i)%function_name .EQ. 'SampPycno') then
		call sample_zero_temp_calcs
	end if
	
	if (function_list(i)%function_name .EQ. 'FoldingSimple') then
		call sample_folding_potential_calcs
	end if

	if (function_list(i)%function_name .EQ. 'turn_pt_001') then
		call turn_pt_001
	end if

end do

call printending 


end program PycnoCalc




