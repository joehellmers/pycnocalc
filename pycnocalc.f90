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

	if (function_list(i)%function_name .EQ. 'genCCRates') then
		call genCCRates
	end if

	if (function_list(i)%function_name .EQ. 'graphM3Y') then
		call graphM3Y
	end if

	if (function_list(i)%function_name .EQ. 'graphSaoPaulo') then
		call graphSaoPaulo
	end if

	if (function_list(i)%function_name .EQ. 'graphAllNN') then
		call graphAllNN
	end if

	if (function_list(i)%function_name .EQ. 'genPycnoRateSplines') then
		call genPycnoRateSplines
	end if

	if (function_list(i)%function_name .EQ. 'turn_pt_001') then
		call turn_pt_001
	end if

	if (function_list(i)%function_name .EQ. 'E0_Energy_001') then
		call E0_Energy_001
	end if

	if (function_list(i)%function_name .EQ. 'Sfactor_001') then
		call Sfactor_001
	end if

	if (function_list(i)%function_name .EQ. 'pycnoRate_001') then
		call pycnoRate_001
	end if

	if (function_list(i)%function_name .EQ. 'pycnoRate_002') then
		call pycnoRate_002
	end if

	if (function_list(i)%function_name .EQ. 'pycnoRate_003') then
		call pycnoRate_003
	end if

	if (function_list(i)%function_name .EQ. 'pycnoRate_004') then
		call pycnoRate_004
	end if

	if (function_list(i)%function_name .EQ. 'adjust_for_SQM_001') then
		call adjust_for_SQM_001
	end if

	if (function_list(i)%function_name .EQ. 'adjust_for_SQM_002') then
		call adjust_for_SQM_002
	end if

	if (function_list(i)%function_name .EQ. 'adjust_for_SQM_002') then
		call adjust_for_SQM_002
	end if

	if (function_list(i)%function_name .EQ. 'vfold_cubes_001') then
		call vfold_cubes_001
	end if

	if (function_list(i)%function_name .EQ. 'mean_wt_001') then
		call mean_wt_001
	end if

	if (function_list(i)%function_name .EQ. 'beta_excite_001') then
		call beta_excite_001
	end if

	if (function_list(i)%function_name .EQ. 'inv_len_001') then
		call inv_len_001
	end if

	if (function_list(i)%function_name .EQ. 'rate_temp_adjust_002') then
		call rate_temp_adjust_002
	end if

	if (function_list(i)%function_name .EQ. 'rate_temp_adjust_003') then
		call rate_temp_adjust_003
	end if

	if (function_list(i)%function_name .EQ. 'tridag_001') then
		call tridag_001
	end if

	if (function_list(i)%function_name .EQ. 'spline_001') then
		call spline_001
	end if

	if (function_list(i)%function_name .EQ. 'spline_002') then
		call spline_002
	end if

	if (function_list(i)%function_name .EQ. 'graphAllFolding') then
		call graphAllFolding
	end if

	if (function_list(i)%function_name .EQ. 'graphAllTotPotential') then
		call graphAllTotPotential
	end if

	if (function_list(i)%function_name .EQ. 'graphSFactor') then
		call graphSFactor
	end if

	if (function_list(i)%function_name .EQ. 'pycnoRateG05_CC_001') then
		call pycnoRateG05_CC_001
	end if

	if (function_list(i)%function_name .EQ. 'graphG05RateCC') then
		call graphG05RateCC
	end if

end do

call printending 


end program PycnoCalc




