!
!	utilities.f90
!	PycnoCalc
!
!	Created by hellmersjl on 5/24/08.
!

module utilities

use constants

implicit none

save


contains 
	
	
	!*********************************
	!
	!  Converts a string to an Integer
	!
	!*********************************
	
	
	integer (kind=8) function ConvertStrToInt(str, fmt)
	use constants
	implicit none
	
	character(len=*), intent(in) :: str
	character(len=*), intent(in), optional :: fmt
	
	if (.NOT. present(fmt)) then
		read(str,'(i19)') ConvertStrToInt
	else
		read(str,fmt) ConvertStrToInt
	end if
	
	end function ConvertStrToInt
		
	
	
	!************************************
	!
	!	Converts and integer to a string
	!
	!************************************
	
	
	character(len=19) function ConvertIntToStr(int, fmt)
	use constants
	implicit none
	
	integer(kind=8), intent(in) :: int
	character(len=*), intent(in), optional :: fmt
	
	if (.NOT. present(fmt)) then
		write(ConvertIntToStr,'(i19)') int
	else
		write(ConvertIntToStr,fmt) int
	end if
	
	end function ConvertIntToStr

	!*********************************
	!
	!	Converts a string to a real
	!
	!*********************************


	real (kind=dbl) function ConvertStrToReal(str, fmt)
	use constants
	implicit none
	
	character(len=*), intent(in) :: str
	character(len=*), intent(in), optional :: fmt
	
	if (.NOT. present(fmt)) then
		read(str,'(ES13.5)') ConvertStrToReal
	else
		read(str,fmt) ConvertStrToReal
	end if
	
	end function ConvertStrToReal

	!****************************************
	!
	!	Conversts a real to a string
	!
	!****************************************
	

	character(len=19) function ConvertRealtoStr(nbr, fmt)
	use constants
	implicit none
	
	integer(kind=dbl), intent(in) :: nbr
	character(len=*), intent(in), optional :: fmt
	
	if (.NOT. present(fmt)) then
		write(ConvertRealToStr,'(ES13.5)') nbr
	else
		write(ConvertRealToStr,fmt) nbr
	end if
	
	end function ConvertRealToStr


	
	!*********************************************************************************
	!
	!	Checks to see if the first non blank character of a string is equal to a value
	!
	!*********************************************************************************
	

	logical function is_first_char(str, chr)

	implicit none
	
	character(len=*), intent(in) :: str
	character(len=1), intent(in) :: chr
	logical :: save_result 
	integer :: i
	
	!print *,'is_first_char: str=', str
	!print *,'is_first_char: char=', chr
	!print *,'len(str)', len(str)
	
	do i=1, len(str)
		!print *,'i=', i
		!print *,'str(i:i)=', str(i:i)
		!print *,'achar(str(i:i))',aichar(str(i:i))
		if (iachar(str(i:i)) /= 32) then
			if (str(i:i) == chr) then
				save_result = .TRUE.
			else
				save_result = .FALSE.
			end if
			exit
		end if
	end do
	
	is_first_char = save_result
	
	end function is_first_char



end module utilities


