!
!	cmdline.f90
!	PycnoCalc
!
!	Created by hellmersjl@icloud.com on 7/25/17.
!

module cmdline

use globalvars

implicit none

save


contains 	
	
	!************************************************
	!
	!  Get configuration file name from command line
	!
	!************************************************
	
	
	subroutine process_cmdline()
	implicit none
	
    integer arglen

	if (command_argument_count() .EQ. 1) then
        call get_command_argument(1, length = arglen)
        allocate(character(arglen) :: cfgfile)
        call get_command_argument(1, value = cfgfile)
    else
        allocate(character(13) :: cfgfile)
        cfgfile = 'pycnocalc.cfg'
    end if

	end subroutine process_cmdline


end module cmdline


