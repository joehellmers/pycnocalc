module logging
implicit none
save

contains

	subroutine scr_and_log (x)
	implicit none
	
	character(len=*) x

	print *, x
		
	end subroutine scr_and_log
	

end module logging
