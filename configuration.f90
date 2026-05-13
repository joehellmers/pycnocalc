!
!	configuration.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

!
!	utilities.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

module configuration

use constants

implicit none

save

integer, parameter :: max_function_cnt = 25
integer, parameter :: max_param_cnt = 100
integer, parameter :: max_function_name_size = 64
integer, parameter :: max_param_name_size = 20
integer, parameter :: max_param_value_size = 30

type :: function_param_type
	character (len=max_param_name_size) :: param_name				!name of parameter
	character (len=max_param_value_size) :: param_value				!string value of parameter
end type function_param_type	

type :: function_list_type
	character (len=max_function_name_size) :: function_name			!name of function to perform
	type (function_param_type), dimension(max_param_cnt) :: params	!array of parameters
	integer :: param_cnt											! number of parameters
end type function_list_type


type (function_list_type), dimension(max_function_cnt) :: function_list
integer :: function_cnt

contains

	!*********************************
	!
	!  Load Configuration File
	!
	!*********************************

	subroutine load_config
	
	use constants
	use logging
	use utilities
    use globalvars
	
	implicit none

	integer :: unit, ierror
	real(kind=dbl) :: real_ierror
	character(len=80) :: cfg_line
	logical :: first_time 
	
	integer :: comment_cnt
	integer :: section_cnt 
	
	call system('cp ' // cfgfile // ' ' // results_dir)
	
	unit = 35
	
	first_time = .TRUE.
	function_cnt = 0
	
	open (unit,file=cfgfile,status='OLD',action='READ', iostat=ierror)
	if (ierror .NE. 0) then
			call scr_and_log_str('Cannot open configuration file')
	        call scr_and_log_str('File: ',lf=.FALSE.)	
			call scr_and_log_str(cfgfile)
			call scr_and_log_str('Error: ',lf=.FALSE.)
			real_ierror =real(ierror)
			call scr_and_log(nbr=real_ierror,fmt='(f5.0)')
			stop
	end if
	
	readloop: Do
	
		read(unit,'(a80)',iostat=ierror) cfg_line
		if (ierror .EQ. -1) then
			! End of file
			exit
		end if
		
		if (ierror .GT. 0) then
				print *, ierror
				call scr_and_log_str('Error Reading from Configuration file' // cfgfile)
				stop
		end if
		
		cfg_line = trim(adjustl(cfg_line))

		if ((.NOT. is_first_char(cfg_line,'#')) .AND. (len_trim(cfg_line) .NE. 0)) then
			if (first_time) then
				if (.NOT. is_first_char(cfg_line,'[')) then
					call scr_and_log_str('Invalid Configuration file format.  The first non-comment line must be a section.  ex. [section_name]')
					call scr_and_log_str(cfg_line)
					stop
				end if
				first_time = .FALSE.
			end if
			if (is_first_char(cfg_line,'[')) then
				function_cnt = function_cnt + 1
				function_list(function_cnt)%param_cnt = 0
				function_list(function_cnt)%function_name = parseFunctionName(cfg_line) 
			else
				if (index(cfg_line,'=') /= 0) then
					function_list(function_cnt)%param_cnt = function_list(function_cnt)%param_cnt + 1
					function_list(function_cnt)%params(function_list(function_cnt)%param_cnt)%param_name = parseParamName(cfg_line)
					function_list(function_cnt)%params(function_list(function_cnt)%param_cnt)%param_value = parseParamValue(cfg_line)
				end if
			end if
		end if

	end do readloop
		
	close(unit)

	end subroutine load_config


	!*********************************
	!
	!  Log Configuration File contents
	!
	!*********************************

		
	subroutine LogConfig
	
	use logging

	integer :: i, j
	
	
	i = 1
	j = 1
	
	call scr_and_log_str('Log of Configuration Contents')
	call scr_and_log_str('-----------------------------')
	call scr_and_log_str(' ')
	call scr_and_log_str('Function Count = ',lf=.FALSE.)
	call scr_and_log(nbr=real(function_cnt,kind=dbl),fmt='(f5.0)')
	call scr_and_log_str(' ')
	
	do i=1,function_cnt
		call scr_and_log_str('Function: ' // function_list(i)%function_name)
		call scr_and_log_str(' ')
		call scr_and_log_str('        Parameter Count = ',lf=.FALSE.)
		call scr_and_log(nbr=real(function_list(i)%param_cnt,kind=dbl),fmt='(f5.0)')
		call scr_and_log_str(' ')
		do j=1,function_list(i)%param_cnt
			call scr_and_log_str('     Parameter Name  : ' // function_list(i)%params(j)%param_name //' Value: ' // function_list(i)%params(j)%param_value)
		end do
	end do
	
	call scr_and_log_str(' ')
	call scr_and_log_str(' ')
	
	end subroutine LogConfig
	
	
	
	!*********************************
	!
	!  Parse out the function name
	!
	!*********************************
	
	character(len=max_function_name_size) function parseFunctionName(str)
	
	implicit none
	
	character(len=*), intent(in) :: str
	
	integer :: first_bracket, last_bracket
	
	first_bracket = index(str,'[')
	last_bracket = index(str,']',back=.TRUE.)
	
	
	
	parseFunctionName = str(first_bracket+1:last_bracket-1)
	
	end function parseFunctionName



	!*********************************
	!
	!  Parse out the parameter name
	!
	!*********************************
	
	character(len=max_param_name_size) function parseParamName(str)
	
	implicit none
	
	character(len=*), intent(in) :: str
	character(len=80) :: wrk_str
	
	wrk_str = trim(adjustl(str)) 
	
	
	parseParamName = str(1:index(str,'=')-1)
	
	end function parseParamName




	!*********************************
	!
	!  Parse out the parameter value
	!
	!*********************************
	
	character(len=max_param_value_size) function parseParamValue(str)
	
	implicit none
	
	character(len=*), intent(in) :: str
	
	parseParamValue = str(index(str,'=')+1:len(trim(str)))
	
	end function parseParamValue


	!******************************************
	!
	!  Get Parameter Value - Subroutine version
	!
	!*********************************
	
	!subroutine getParamValue(function_name, param_name, param_value, error_flg)
	
	!implicit none
	
	!character(len=*), intent(in)						:: function_name
	!character(len=*), intent(in)						:: param_name
	!character(len=max_param_value_size), intent(out)	:: param_value
	!logical, intent(out)								:: error_flg
	
	!integer :: i, j
	!logical :: fnd_flg


	!error_flg = .FALSE.
	!fnd_flg = .FALSE.
	
	!do i = 1, function_cnt		
	!	if (function_list(i)%function_name .EQ. function_name) then
	!		do j = 1, function_list(i)%param_cnt
	!			if (function_list(i)%params(j)%param_name .EQ. param_name) then
	!				fnd_flg = .TRUE.
	!				param_value = function_list(i)%params(j)%param_value
	!				exit
	!			end if
	!		end do
	!	end if
	!	if (fnd_flg) then
	!		exit
	!	end if
	!end do
	
	!if (fnd_flg) then
	!	error_flg = .FALSE.
	!else
	!	param_value = ''
	!	error_flg = .TRUE.
	!end if 
	
	
	!end subroutine getParamValue

	!******************************************
	!
	!  Get Parameter Value - Function version
	!
	!******************************************
	
	character(len=max_param_value_size) function getParamValue(function_name, param_name)
	
	implicit none
	
	character(len=*), intent(in)						:: function_name
	character(len=*), intent(in)						:: param_name
	
	integer :: i, j
	logical :: fnd_flg


	fnd_flg = .FALSE.
	
	do i = 1, function_cnt		
		if (function_list(i)%function_name .EQ. function_name) then
			do j = 1, function_list(i)%param_cnt
				if (function_list(i)%params(j)%param_name .EQ. param_name) then
					fnd_flg = .TRUE.
					getParamValue = function_list(i)%params(j)%param_value
					exit
				end if
			end do
		end if
		if (fnd_flg) then
			exit
		end if
	end do
	
	if (.NOT. fnd_flg) then
		getParamValue = ''
	end if 
	
	end function getParamValue

	!******************************************
	!
	!  Check for Parameter existence
	!
	!******************************************
	
	logical function checkParam(function_name, param_name)
	
	implicit none
	
	character(len=*), intent(in)	:: function_name
	character(len=*), intent(in)	:: param_name
	
	integer :: i, j
	logical :: fnd_flg


	fnd_flg = .FALSE.
	
	do i = 1, function_cnt		
		if (function_list(i)%function_name .EQ. function_name) then
			do j = 1, function_list(i)%param_cnt
				if (function_list(i)%params(j)%param_name .EQ. param_name) then
					fnd_flg = .TRUE.
					exit
				end if
			end do
		end if
		if (fnd_flg) then
			exit
		end if
	end do
	
	if (fnd_flg) then
		checkParam = .TRUE.
	else
		checkParam = .FALSE.
	end if 
	
	end function checkParam

	
end module configuration


