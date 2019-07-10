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
use mathintegration
use logging
use general_nuclear

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

    subroutine turn_pt(R,Rstep,Rmax,E0,mu,A1,A2,Z1,Z2,radius1,radius2,rho0_A1,rho0_A2,L,partition,turn1,turn2,WKB,nucIntType,inSQMFlag)

        implicit none

        real(kind=dbl), intent(in)      :: R            ! Distance between nuclei in fermi (10e-15 m)
        real(kind=dbl), intent(in)      :: Rstep        ! lattice step size (user input)
        integer, intent(in)             :: Rmax         ! max size of array for R lattice
        real(kind=dbl), intent(in)      :: E0           ! zero-pt vibrational energy of incoming nucleon (1st)
        real(kind=dbl), intent(in)      :: mu           ! reduced mass value holder 
        integer, intent(in)             :: A1           ! number of nucleons in the first nucleus
        real(kind=dbl), intent(in)      :: A2           ! number of nucleons in the second nucleus
        integer, intent(in)             :: Z1           ! proton count for nuclei 1
        real(kind=dbl), intent(in)      :: Z2           ! proton count for nuclei 2
        real(kind=dbl), intent(in)      :: radius1      ! radius of the first nuclei in fermi
        real(kind=dbl), intent(in)      :: radius2      ! radius of the second nuclei in fermi
        real(kind=dbl), intent(in)      :: rho0_A1      ! central densities for first nucleus
        real(kind=dbl), intent(in)      :: rho0_A2      ! central densities for second nucleus
        integer, intent(in)             :: L            ! orbital angular momentum (not in use) 
        integer, intent(in)             :: partition    ! number of subdivisions to use for the calculation
        integer, intent(out)            :: turn1        ! position along R-axis of first turning point
        integer, intent(out)            :: turn2        ! position along R-axis of second turning point
        real(kind=dbl), intent(out)     :: WKB			! value of WKB calculation to get through total barrier (not just coulomb barrier) between incoming and target nuclei
        integer, intent(in)             :: nucIntType   ! Type of nuclear interaction to use
        logical, intent(in), optional   :: inSQMFlag    ! Indicate if we are using SQM for nuclei/nugget 2    

        real(kind=dbl)  :: ln_WKB                    ! natural log of the value of the WKB calculation to get through the total barrier
        real(kind=dbl)  :: V_2fold					! double folding potential calculated with distance of R between incoming and target nuclei
        real(kind=dbl)  :: ndensity1, ndensity2		! number density of nuclei 1 and 2, calculated using baryon number
        real(kind=dbl)  :: VEcheck(0:Rmax)	 	    ! array to hold difference between Veff and E of incoming particle at every point R (between incoming & target nuclei)
        real(kind=dbl)  :: Vnucarray(0:Rmax)	    	! array to hold Veff between target & incoming particle at every point along R (between incoming & target nuclei)
        real(kind=dbl)  :: Vcoulary(0:Rmax)			! array to hold Vcoulomb between traget and incoming particle at every point along R (between incoming & target nuclei)
        integer         :: i, turn_counter                  ! ad hoc counter, counter for the number of turning points 
        real(kind=dbl)  :: R_pos						! current position along R axis --this is the CURRENT separation between target and incoming particle	  
        real(kind=dbl)  :: Energy					! function that calculates the Energy of incoming (or just second (projectile) nucleon)
        real(kind=dbl)  :: Integrand(0:Rmax)  		! array to hold Integrand of S-factor calculation
        real(kind=dbl)  :: S							! astrophysical S-factor - calculated using equation from PHYS REV C69 --rule of thumb model fitted to data
        real(kind=dbl)  :: ln_S						! natural log of astrophysical S-factor
        integer         :: new_Rmax                  ! new Rmax based upon cutoff
        logical         :: SQMFlag = .FALSE.         ! Used internally

        logical         :: v_fold_threshold_flg
        logical         :: v_fold_first_time
        real(kind=dbl)  :: v_fold_max
        real(kind=dbl)  :: v_fold_min

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        v_fold_min = 1e-10
        v_fold_max = 0
        v_fold_threshold_flg = .false.
        v_fold_first_time = .true.

        i = 0
        turn_counter = 0
        turn1 = -1
        turn2 = -1
      
!	Fill in VEcheck array first with values of Veff - E....  Veff is composed of Vnuc + Vcoul 
!	NOTE:  R_pos is the distance between the two nuclei centers!  You are essentially starting the incoming particle right next to the target particle and then backing up
!		to the starting separation distance of R ...   you are calculating VEcheck array IN REVERSE, starting where nuclei are touching (Rpos = 0) and then moving incoming
!		particle backward to Rpos = R...  So, for my visual sake, we are filling in VEcheck array from right to left (starting at Rmax and moving left to 0)

        !print *,"r,TotEnergy"
        do i = 0,Rmax
            R_pos = i * Rstep
            if (.not. v_fold_threshold_flg) then 
                V_2fold = vfold_spherically_symmetric (R_pos, A1, A2, Z1, Z2, partition,0.5_dbl,0.5_dbl,rho0_A1, rho0_A2, radius1, radius2, E0, nucIntType, SQMFlag)
                !print *, V_2fold
                if (v_fold_max .LT. abs(v_2fold)) then
                    v_fold_max = abs(v_2fold)
                end if
                if (v_fold_first_time) then
                    v_fold_first_time = .false.
                else
                    if ((abs(v_2fold)/v_fold_max) .LT. v_fold_min) then
                        v_fold_threshold_flg = .true.			
                        v_2fold = 0.0_dbl
                    end if
                end if  
            else
                v_2fold = 0.0_dbl
            end if
            !print *,"V_2fold at ", R_pos, " = ", V_2fold            
            Vnucarray(Rmax-i) = V_2fold
            Vcoulary(Rmax-i) = Vcoulomb(Z1,Z2,R_pos,radius1,radius2)
            VEcheck(Rmax-i) = V_2fold + Vcoulary(Rmax-i) - E0
            !print *, Rmax-i, "Electrostatic ", Vcoulary(Rmax-i)
            !print *,R_pos, ",", VEcheck(Rmax-i)
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
                    Integrand(i) = 0.0_dbl
                else
                    if (VEcheck(i) .ge. 0.0_dbl) then
                        Integrand(i) = .01432_dbl*sqrt(mu*VEcheck(i))    ! KEY DIFFERENCE:  hbar in SsubL def!
                    else
                        Integrand(i) = -0.01432_dbl*sqrt(-1.0_dbl*mu*VEcheck(i))    ! KEY DIFFERENCE:  hbar in SsubL def!
                    end if
                end if
                !print *, Integrand(i),VEcheck(i)
            end do		  	  
            WKB = trapezoidArray(Rmax,turn1,Integrand,Rstep)	
        else	  
        !	THERE ARE USUALLY TWO TURNING POINTS WHEN USING SALPETER AND VAN HORN... NEED TO INTEGRATE THROUGH BARRIER
            do i = 0,Rmax
                if ((i.le.turn1).or.(i.gt.turn2)) then
                    Integrand(i) = 0.0_dbl
                else
                    if (VEcheck(i) .ge. 0.0_dbl) then
                        Integrand(i) = .01432_dbl*sqrt(mu*VEcheck(i))    ! KEY DIFFERENCE:  hbar in SsubL def!
                    else
                        Integrand(i) = -0.01432_dbl*sqrt(-1.0_dbl*mu*VEcheck(i))    ! KEY DIFFERENCE:  hbar in SsubL def!
                    end if
                end if
                !print *, Integrand(i),VEcheck(i)
            end do	  
            !print *,"mu=",mu
            WKB = trapezoidArray(Rmax,turn2,Integrand,Rstep)
        end if

        ! New, simpler way, just include positive energies
        !do i = 0, Rmax
        !    if (VEcheck(i) .gt. 0) then
        !        Integrand(i) = .01432_dbl*sqrt(mu*VEcheck(i))
        !    else
        !        Integrand(i) = 0.0_dbl
        !    end if
        !end do

        ! New, even simpler way, just integrate to first turning point
        !Integrand = 0.0_dbl
        !do i = 0, Rmax
        !    if (VEcheck(i) .gt. 0) then
        !        Integrand(i) = .01432_dbl*sqrt(mu*VEcheck(i))
        !    else
        !        exit
        !    end if
        !end do
        !write(*,*) Integrand
        WKB = trapezoidArray(Rmax,Rmax,Integrand,Rstep)
        
    end subroutine turn_pt


!*****************************************************************
!
! Calculate the S Factor
!
! REF: Golf Thesis
!
!*****************************************************************

    real(kind=dbl) function Sfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag)

        implicit none
        
        integer, intent(in)             :: A1_int, Z1_int
        real(kind=dbl), intent(in)      :: A2, Z2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Rstep
        integer, intent(in)             :: partition
        integer, intent(in)             :: nucIntType
        logical, intent(in), optional   :: inSQMFlag

        real(kind=dbl)  :: A1, Z1
        integer         :: L = 0
        real(kind=dbl)  :: E0
        real(kind=dbl)  :: ln_sigma
        real(kind=dbl)  :: R
        integer         :: Rmax
        real(kind=dbl)  :: mu
        real(kind=dbl)  :: radius1, radius2
        real(kind=dbl)  :: rho0_A1, rho0_A2
        integer         :: turn1, turn2
        real(kind=dbl)  :: WKB
        real(kind=dbl)  :: ln_S
        real(kind=dbl)  :: ln_Trans_total
        integer         :: i
        logical         :: SQMFlag = .FALSE.

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        A1 = real(A1_int,dbl)
        Z1 = real(Z1_int,dbl)

        R = sqrt(3.0_dbl)*0.5_dbl*lattice(rho,A1,Z1)
        Rmax = (R/Rstep)+1
        mu = reduced_mass(A1_int, Z1_int, A2, Z2)
        
        radius1 = nuclear_radius(A1,.FALSE.)
        radius2 = nuclear_radius(A2, SQMFlag)
        rho0_A1 = rho0_2pF(A1,radius1,0.5_dbl)
        rho0_A2 = rho0_2pF(A2,radius2,0.5_dbl)
        
        !print*, "radius1 = ", radius1
        !print*, "radius2 = ", radius2

!	Calculate E of "incoming" (ground state vibrating) particle coming toward lattice-bound target particle
        E0=E0_Energy(Z1_int,Z2, A1_int, A2, rho)
        !print*,'E0 =',E0

!	Calculate Veffective and WKB integration INSIDE turn_pts function - you should have all other input parameters at this point
!		AND you CAN'T carry the V_arrays back into the main program because Rmax is a DERIVED paramater - and used as the array dimension
        ln_sigma = 0.0
        do i = 0,L
            call turn_pt(R,Rstep,Rmax,E0,mu,A1_int,A2,Z1_int,Z2,radius1,radius2,rho0_A1,rho0_A2,L,partition,turn1,turn2,WKB,nucIntType,SQMFlag)
            ! ln of Total transmission Probability: ',ln_Trans_total
            ln_Trans_total = -WKB
            ln_sigma=ln_sigma+log(612.459_dbl)-log(mu*E0)+log(2.0_dbl*real(i,dbl)+1.0_dbl)+ln_Trans_total
        end do

        ln_S=ln_sigma+log(E0)+(Z1)*(Z2)*.0324_dbl*sqrt(mu/E0)
        !print*,'ln_S =',ln_S
        !print*,'log10_S =',(ln_S)*.4343
        Sfactor = exp(ln_S)
	end function Sfactor

!*****************************************************************
!
! Calculate Pycnonuclear reaction rates
! This version also calculates the Folding potential
!
! REF: Golf Thesis
!
!*****************************************************************

    real(kind=dbl) function pycnoRate(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag)

        implicit none

        integer, intent(in)             :: A1_int, Z1_int
        real(kind=dbl), intent(in)      :: A2, Z2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Rstep
        integer, intent(in)             :: partition
        integer, intent(in)             :: nucIntType
        logical, intent(in), optional   :: inSQMFlag

        real(kind=dbl)              :: ln_P0
        real(kind=dbl)              :: P0
        real(kind=dbl)              :: lambda
        real(kind=dbl)              :: ln_lambda
        real(kind=dbl)              :: ln_S
        real(kind=dbl)              :: S
        real(kind=dbl)              :: mn_mass, mn_chrg

        real(kind=dbl)              :: A1, Z1

        logical         :: SQMFlag = .FALSE.

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        A1 = real(A1_int,dbl)
        Z1 = real(Z1_int,dbl)

        S = Sfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, SQMFlag)
        ln_S = log(S)
        !print *,'ln_S = ', ln_S

        lambda = inv_len_param2comp (rho, A1_int, Z1_int, 0.5_dbl, int(A2), int(Z2), 0.5_dbl)
        ln_lambda = log(lambda)
        ln_P0=-2.638_dbl/(sqrt(lambda))+log(rho)+log(A1*A2)-log(A1+A2)+2.0_dbl*log(Z1*Z2)+ln_S+1.75_dbl*ln_lambda+109.36_dbl

        !print *, 'mn_mass = ', mn_mass
        !print *, 'mn_chrg = ', mn_chrg
        !print *, 'ln_lambda = ', ln_lambda
        !print *, 'lambda = ', lambda
        !print *, 'ln_P0 = ', ln_P0

        P0 = exp(ln_P0)

        pycnoRate = P0

    end function pycnoRate

!*****************************************************************
!
! Calculate the adjustment to the Pycno Reaction rate based on 
! Termparature
!
! REF: SPVH1969 (45)
!
!*****************************************************************

    real(kind=dbl) function tempRateAdjust(A1_int, Z1_int, X1, A2_int, Z2_int, X2, rho, Temp, inlowFlag)

        implicit none
        
        integer, intent(in)             :: A1_int, Z1_int, A2_int, Z2_int
        real(kind=dbl), intent(in)      :: X1, X2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Temp   
        logical, intent(in), optional   :: inlowFlag
                        
        real(kind=dbl)  :: beta_3halves_rt
        real(kind=dbl)  :: beta_factor
        real(kind=dbl)  :: inv_len_neg_sqrt
        logical         :: lowFlag
        real(kind=dbl)  :: param1, param2, param3, param4
        real(kind=dbl)  :: factor
                
        if (present(inlowFlag)) then
            lowFlag = inlowFlag
        else
            lowFlag = .true.
        end if

        if (lowFlag) then
            param1 = 0.0430_dbl
            param2 = 1.2624_dbl
            param3 = 1.2231_dbl
            param4 = 0.6310_dbl
        else
            param1 = 0.0485_dbl
            param2 = 2.9314_dbl
            param3 = 1.4331_dbl
            param4 = 1.4654_dbl
        end if

        beta_3halves_rt = beta_excitation2comp(A1_int,Z1_int,X1,A2_int,Z2_int,X2,rho,Temp)**(3.0_dbl/2.0_dbl)
        beta_factor = exp(-8.7833_dbl*beta_3halves_rt)
        inv_len_neg_sqrt = inv_len_param2comp (rho, A1_int, Z1_int, X1, A2_int, Z2_int, X2)**(-1.0_dbl/2.0_dbl)

        factor = 1.0_dbl - param4*beta_factor
        factor = inv_len_neg_sqrt*param3*beta_factor*factor 
        factor = exp(-7.272*beta_3halves_rt + factor)
        factor = ((1.0_dbl + param2*beta_factor)**(-1.0_dbl/2.0_dbl))*factor
        factor = 1.0_dbl + param1*inv_len_neg_sqrt*factor

        tempRateAdjust = factor
                    
    end function tempRateAdjust
    
end module rate_calc
