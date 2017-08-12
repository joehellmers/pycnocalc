!
!	rate_calc.f90
!	PycnoCalc
!
!	Created by hellmersjl on 08/09/17.
!
module rate_calc

use constants
use folding_potential
use astrophysics
use integration
use logging

implicit none

save

contains


!*****************************************************************
!
! Subroutine for calculating turning points for S-factor integral
! finds the turning points by looking for the places where 
! Veff - E changes sign.  	
!
! REF: Golf Thesis
!
!*****************************************************************

    subroutine turn_pt(R,Rstep,Rmax,E0,mu,A1,A2,SQM_A2,Z1,Z2,radius1,radius2,rho0_A1,rho0_A2,L,partition,turn1,turn2,WKB)

        implicit none

        real(kind=dbl), intent(in)  :: R            ! Distance between nuclei in fermi (10e-15 m)
        real(kind=dbl), intent(in)  :: Rstep        ! lattice step size (user input)
        integer, intent(in)         :: Rmax         ! max size of array for R lattice
        real(kind=dbl), intent(in)  :: E0           ! zero-pt vibrational energy of incoming nucleon (1st)
        real(kind=dbl), intent(in)  :: mu           ! reduced mass value holder 
        integer, intent(in)  :: A1                  ! number of nucleons in the first nucleus
        integer, intent(in)  :: A2                  ! number of nucleons in the second nucleus
        real(kind=dbl), intent(in)  :: SQM_A2       ! baryon number for second nuclei if SQM involved
        integer, intent(in)  :: Z1                  ! proton count for nuclei 1
        integer, intent(in)  :: Z2                  ! proton count for nuclei 2
        real(kind=dbl), intent(in)  :: radius1      ! radius of the first nuclei in fermi
        real(kind=dbl), intent(in)  :: radius2      ! radius of the second nuclei in fermi
        real(kind=dbl), intent(in)  :: rho0_A1      ! central densities for first nucleus
        real(kind=dbl), intent(in)  :: rho0_A2      ! central densities for second nucleus
        integer, intent(in)         :: L            ! orbital angular momentum (not in use) 
        integer, intent(in)         :: partition    ! number of subdivisions to use for the calculation
        integer, intent(out)        :: turn1        ! position along R-axis of first turning point
        integer, intent(out)        :: turn2        ! position along R-axis of second turning point
        real(kind=dbl), intent(out) :: WKB			! value of WKB calculation to get through total barrier (not just coulomb barrier) between incoming and target nuclei
        
        real(kind=dbl) :: ln_WKB                    ! natural log of the value of the WKB calculation to get through the total barrier
        real(kind=dbl) :: V_2fold					! double folding potential calculated with distance of R between incoming and target nuclei
        real(kind=dbl) :: ndensity1, ndensity2		! number density of nuclei 1 and 2, calculated using baryon number
        real(kind=dbl) :: VEcheck(0:Rmax)	 	    ! array to hold difference between Veff and E of incoming particle at every point R (between incoming & target nuclei)
        real(kind=dbl) :: Vnucarray(0:Rmax)	    	! array to hold Veff between target & incoming particle at every point along R (between incoming & target nuclei)
        real(kind=dbl) :: Vcoulary(0:Rmax)			! array to hold Vcoulomb between traget and incoming particle at every point along R (between incoming & target nuclei)
        integer :: i, turn_counter                  ! ad hoc counter, counter for the number of turning points 
        real(kind=dbl) :: R_pos						! current position along R axis --this is the CURRENT separation between target and incoming particle	  
        real(kind=dbl) :: Energy					! function that calculates the Energy of incoming (or just second (projectile) nucleon)
        real(kind=dbl) :: Integrand(0:Rmax)  		! array to hold Integrand of S-factor calculation
        real(kind=dbl) :: S							! astrophysical S-factor - calculated using equation from PHYS REV C69 --rule of thumb model fitted to data
        real(kind=dbl) :: ln_S						! natural log of astrophysical S-factor
        integer        :: new_Rmax                  ! new Rmax based upon cutoff
         
        i = 0
        turn_counter = 0
        turn1 = -1
        turn2 = -1
      
!	Fill in VEcheck array first with values of Veff - E....  Veff is composed of Vnuc + Vcoul 
!	NOTE:  R_pos is the distance between the two nuclei centers!  You are essentially starting the incoming particle right next to the target particle and then backing up
!		to the starting separation distance of R ...   you are calculating VEcheck array IN REVERSE, starting where nuclei are touching (Rpos = 0) and then moving incoming
!		particle backward to Rpos = R...  So, for my visual sake, we are filling in VEcheck array from right to left (starting at Rmax and moving left to 0)

        do i = 0,Rmax
            R_pos = i * Rstep
            V_2fold = vfold_spherically_symmetric (R_pos, A1, A2, Z1, Z2, partition,0.5_dbl,0.5_dbl,rho0_A1, rho0_A2, radius1, radius2, E0)
            print *,"V_2fold at ", R_pos, " = ", V_2fold            
            Vnucarray(Rmax-i) = V_2fold
            Vcoulary(Rmax-i) = Vcoulomb(Z1,Z2,R_pos,radius1,radius2)
            VEcheck(Rmax-i) = V_2fold + Vcoulomb(Z1,Z2,R_pos,radius1,radius2) - E0
        end do
	  
!	now check the VEcheck array for turning point(s) - there may be more than one - especially if the total energy of the incoming particle starts out high (greater than the potential)
        do i = 1,Rmax
            if ((VEcheck(i).le.0).and.(VEcheck(i-1).ge.0)) then   		!if Veff - E(i) starts positive and goes negative....
                turn_counter = turn_counter + 1							!then you've found a turning point....	
                if (turn1.lt.0) then
                    turn1 = i - 1										  	!and the first turn just prior to this position
                else 													
                    turn2 = i - 1
                end if
            else 
                if ((VEcheck(i).ge.0).and.(VEcheck(i-1).le.0)) then 	!if Veff - E(i) starts negative and goes positive...
                    turn_counter = turn_counter + 1							!another turn....
                    if (turn1.lt.0) then
                        turn1 = i - 1										  	!and the turn is just prior to this position
                    else										
                        turn2 = i - 1
                    end if 
                end if            
            end if
        end do
	  
!		Equation from PHYS REV C69, 034603 (2005):   WKB = Integral from R1 to R2 [ sqrt( (8*mu)/hbar^2 *(Veff[R,E] - E)) ] dR
!			(Veff[R,E] - E) is the VEcheck array.
!			See notes for "Integrand - 24 April 2008"

        if (turn2.lt.0) then
        !	THERE MAY ONLY BE ONE TURNING POINT - STILL NEED TO INTEGRATE THROUGH BARRIER
            do i = 0,Rmax
                if (i.gt.turn1) then
                    Integrand(i) = 0.0
                else
                    Integrand(i) = .01432*sqrt(mu*VEcheck(i))    ! KEY DIFFERENCE:  hbar in SsubL def!
                end if
            end do		  	  
            WKB = trapezoidArray(Rmax,turn1,Integrand,Rstep)	
        else	  
        !	THERE ARE USUALLY TWO TURNING POINTS WHEN USING SALPETER AND VAN HORN... NEED TO INTEGRATE THROUGH BARRIER
            do i = 0,Rmax
                if ((i.le.turn1).or.(i.gt.turn2)) then
                    Integrand(i) = 0.0
                else
                    Integrand(i) = .01432*sqrt(mu*VEcheck(i))    ! KEY DIFFERENCE:  hbar in SsubL def!
                end if
            end do	  
            WKB = trapezoidArray(Rmax,turn2,Integrand,Rstep)
        end if

    end subroutine turn_pt
	
end module rate_calc
