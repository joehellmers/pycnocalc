!
!	rate_data.f90
!	PycnoCalc
!
!	Created by hellmersjl on 7/31/18.
!

module rate_data

use constants

implicit none

save

!enum :: rate_type
!    enumerator , pycno_saopaulo, pycno_m3y
!endenum
integer,parameter :: pycno_basic=1,pycno_saopaulo=2,pycno_m3y=3, rate_types(3) = [pycno_basic,pycno_saopaulo,pycno_m3y]

contains 

subroutine check_enum
    print *, "Checking enum: "
    !print *, "pycno_basic = ", pycno_basic
    !print *, "pycno_saopaulo = ", pycno_saopaulo
    !print *, "pycno_m3y = ", pycno_m3y
    
    !kind(pycno_m3y) rate_type = pycno_m3y
    !if (rate_type .eq. pycno_m3y) then
    !    print *, "Rate type is M3Y"
    !end if
    !if (rate_type .ne. pycno_basic) then
    !    print *, "Rate type is basic"
    !end if
    
        
end subroutine check_enum
	
end module rate_data
