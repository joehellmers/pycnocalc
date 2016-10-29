module logging

use constants

implicit none
save

character(len=13), parameter :: log_file = 'PycnoCalc.log'

character(len=21)			:: results_dir  ! directory for results
character(len=2+21+1+13)	:: full_log_file	

contains



	! Initialize the log file
	
	subroutine initialize_logging
	
	integer :: unit, ierror
	character(len=8) :: currdate
	character(len=10) :: currtime
	
	call date_and_time(date=currdate,time=currtime)
	
	results_dir = 'results' // currdate // currtime(1:6)
	call system('mkdir ' // results_dir)
	unit = 25
	
	full_log_file = './' // results_dir // '/' // log_file
	
	open(unit, file=full_log_file, status='REPLACE', ACTION='WRITE', iostat=ierror)
	if (ierror .NE. 0) then
		print *, full_log_file
		print *, ierror
		stop
	end if

	write (unit,'(a)') '---- PycnoCalc Log ----'
	close(unit)
	
	end subroutine initialize_logging

    ! Print to log and screen

	subroutine scr_and_log (str, nbr, fmt, lf)
	implicit none
	
	character(len=*), intent(in), optional :: str
	real(kind=dbl), intent(in), optional :: nbr
	character(len=*), intent(in) :: fmt
	logical, intent(in), optional :: lf
	integer :: unit, ierror
	character(len=3) :: c_advance
	
	unit = 25
	
	if (present(lf)) then
		if (lf) then
			c_advance = 'YES'
		else
			c_advance = 'NO'
		end if
	else
		c_advance = 'YES'
	end if
	
	open(unit, file=full_log_file, status='OLD', ACTION='WRITE', iostat=ierror, access='APPEND')
	
	if (present(str) .and. present(nbr)) then
		write (unit,'(a)',advance='NO') str
		write (unit,fmt,advance=c_advance) nbr
		write(*, '(a)', advance='NO') str
		write(*, fmt, advance=c_advance) nbr
	else
		if (present(str)) then
			write (unit,fmt,advance=c_advance) str
			write(*, fmt, advance=c_advance) str
		else
			if (present(nbr)) then
				write (unit,fmt,advance=c_advance) nbr
				write(*, fmt,advance=c_advance) nbr
			else
				write (unit,'(a)',advance=c_advance) 'ERROR ----> No String or Number value given to Logging:scr_and_log'
				write(*, fmt, advance=c_advance) 'ERROR ----> No String or Number value given to Logging:scr_and_log'
			end if
		end if
	end if
	
	close(unit)

			
	end subroutine scr_and_log

	
	subroutine scr_and_log_println
	implicit none
		
	call scr_and_log(str=' ',fmt='(a)')
	
	
    end subroutine scr_and_log_println

	subroutine scr_and_log_str(in_str,lf)
	implicit none
	character(len=*), intent(in) :: in_str
	logical, intent(in), optional :: lf
	
	
	
	call scr_and_log(str=in_str,fmt='(a)',lf=lf)
	
	
    end subroutine scr_and_log_str
	

	subroutine printending
	
		call scr_and_log_str(' -------  End Run: PycnoCalc --------')
	
	end subroutine printending
	
end module logging
