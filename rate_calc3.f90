!
!	rate_calc.f90
!	PycnoCalc
!
!	Created by hellmersjl on 08/09/17.
!
module rate_calc3

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

    subroutine turn_pt3(Rstep,Rmax,E0,mu,turn1,turn2,WKB,inSQMFlag)

        implicit none

        real(kind=dbl), intent(in)      :: Rstep        ! lattice step size (user input)
        integer, intent(in)             :: Rmax         ! max size of array for R lattice
        integer :: line
        real(kind=dbl), intent(in)      :: E0           ! zero-pt vibrational energy of incoming nucleon (1st)
        real(kind=dbl), intent(in)      :: mu           ! reduced mass value holder 
        integer, intent(out)            :: turn1        ! position along R-axis of first turning point
        integer, intent(out)            :: turn2        ! position along R-axis of second turning point
        real(kind=dbl), intent(out)     :: WKB			! value of WKB calculation to get through total barrier (not just coulomb barrier) between incoming and target nuclei
        logical, intent(in), optional   :: inSQMFlag    ! Indicate if we are using SQM for nuclei/nugget 2    

        real(kind=dbl)  :: ln_WKB                    ! natural log of the value of the WKB calculation to get through the total barrier
        real(kind=dbl)  :: V_2fold					! double folding potential calculated with distance of R between incoming and target nuclei
        real(kind=dbl)  :: ndensity1, ndensity2		! number density of nuclei 1 and 2, calculated using baryon number
        real(kind=dbl)  :: VEcheck(0:Rmax)	 	    ! array to hold difference between Veff and E of incoming particle at every point R (between incoming & target nuclei)
        real(kind=dbl)  :: Vnucarray(0:Rmax)	    	! array to hold Veff between target & incoming particle at every point along R (between incoming & target nuclei)
        real(kind=dbl)  :: Vcoulary(0:Rmax)			! array to hold Vcoulomb between traget and incoming particle at every point along R (between incoming & target nuclei)
        real(kind=dbl)  :: VEdata(0:Rmax)
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
        real(kind=dbl)  :: rcm, iter, UCfin, UNDfin, Uexfin, UNfin, Utotfin

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

            !print *,"V_2fold at ", R_pos, " = ", V_2fold 

            ! read data file containing all potentials, from DFMDEF2013 code
            open(unit=10, file='901.dat', status='old')
            line = 0
            do 
                read(10,*,end=100) rcm, iter, UCfin, UNDfin, Uexfin, UNfin, Utotfin
                if (line > Rmax) exit
                VEdata(line) = Utotfin
                line = line + 1

            end do
            100 continue
            close(10)

            do i=0, Rmax
                VEcheck(i) = VEdata(i) - E0
            end do

            ! Check the minimum and maximum value of Unuclear + Ucoulomb - E0
            print *, 'min VEcheck = ', minval(VEcheck)
            print *, 'max VEcheck = ', maxval(VEcheck)

            !print *, Rmax-i, "Electrostatic ", Vcoulary(Rmax-i)
            !print *,R_pos, ",", VEcheck(Rmax-i)

	  
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
        end if

        WKB = trapezoidArray(Rmax,Rmax,Integrand,Rstep)
        
    end subroutine turn_pt3
    
        real(kind=dbl) function lnSfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag)

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
        integer         :: Rmax
        real(kind=dbl)  :: mu
        integer         :: turn1, turn2
        real(kind=dbl)  :: WKB
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

        Rmax = 200
        mu = reduced_mass(A1_int, Z1_int, A2, Z2)
        E0 = E0_Energy(Z1_int, Z2, A1_int, A2, rho)

        ln_sigma = 0.0_dbl
        do i = 0,L
            call turn_pt3(Rstep, Rmax, E0, mu, turn1, turn2, WKB, SQMFlag)
            ln_Trans_total = -WKB
            ln_sigma = ln_sigma + log(612.459_dbl) - log(mu*E0) + &
                       log(2.0_dbl*real(i,dbl)+1.0_dbl) + ln_Trans_total
        end do

        lnSfactor = ln_sigma + log(E0) + (Z1)*(Z2)*0.0324_dbl*sqrt(mu/E0)

    end function lnSfactor

    real(kind=dbl) function Sfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag)

        implicit none

        integer, intent(in)             :: A1_int, Z1_int
        real(kind=dbl), intent(in)      :: A2, Z2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Rstep
        integer, intent(in)             :: partition
        integer, intent(in)             :: nucIntType
        logical, intent(in), optional   :: inSQMFlag

        real(kind=dbl) :: ln_S

        ln_S = lnSfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag)

        print *, 'ln_S    = ', ln_S
        print *, 'log10_S = ', ln_S / log(10.0d0)

        if (ln_S > log(huge(1.0d0))) then
            Sfactor = huge(1.0d0)
        else
            Sfactor = exp(ln_S)
        end if

    end function Sfactor

    real(kind=dbl) function pycnoRate(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag, inCalcType)

        implicit none

        integer, intent(in)             :: A1_int, Z1_int
        real(kind=dbl), intent(in)      :: A2, Z2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Rstep
        integer, intent(in)             :: partition
        integer, intent(in)             :: nucIntType
        logical, intent(in), optional   :: inSQMFlag
        integer, intent(in), optional   :: inCalcType

        real(kind=dbl)                  :: prepend, alpha1, alpha2, gamma
        real(kind=dbl)                  :: ln_P0
        real(kind=dbl)                  :: P0
        real(kind=dbl)                  :: lambda
        real(kind=dbl)                  :: ln_lambda
        real(kind=dbl)                  :: ln_S
        real(kind=dbl)                  :: A1, Z1
        logical                         :: SQMFlag = .FALSE.
        integer                         :: calcType = 0

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        if (present(inCalcType)) then
            calcType = inCalcType
        end if

        A1 = real(A1_int,dbl)
        Z1 = real(Z1_int,dbl)

        ln_S = lnSfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, SQMFlag)

        if (calcType .eq. 0) then

            lambda = inv_len_param2comp(rho, A1_int, Z1_int, 0.5_dbl, int(A2), int(Z2), 0.5_dbl)
            ln_lambda = log(lambda)

            ln_P0 = -2.638_dbl/(sqrt(lambda)) + log(rho) + log(A1*A2) - log(A1+A2) + &
                    2.0_dbl*log(Z1*Z2) + ln_S + (4.0_dbl/7.0_dbl)*ln_lambda + 109.36_dbl

            print *, 'ln_P0      = ', ln_P0
            print *, 'log10_rate = ', ln_P0 / log(10.0d0)

            if (ln_P0 > log(huge(1.0d0))) then
                P0 = huge(1.0d0)
            else if (ln_P0 < log(tiny(1.0d0))) then
                P0 = 0.0d0
            else
                P0 = exp(ln_P0)
            end if

            pycnoRate = P0

        else

            if (calcType .eq. 1 .or. calcType .eq. 2 .or. calcType .eq. 3) then
                prepend = 1.06_dbl
                gamma = 2.0_dbl
            else
                prepend = 2.69_dbl
                gamma = 4.0_dbl
            end if

            lambda = 0.0245_dbl*(A1**(-4.0_dbl/3.0_dbl))*(Z1**(-2.0_dbl))* &
                     (gamma**(-1.0_dbl/3.0_dbl))*((rho/(1.0d6))**(1.0_dbl/3.0_dbl))

            if (calcType .eq. 1) then
                alpha1 = 2.639_dbl
                alpha2 = -6.305_dbl
            end if

            if (calcType .eq. 2) then
                alpha1 = 2.516_dbl
                alpha2 = -6.793_dbl
            end if

            if (calcType .eq. 3) then
                alpha1 = 2.517_dbl
                alpha2 = -6.754_dbl
            end if

            if (calcType .eq. 4) then
                alpha1 = 2.401_dbl
                alpha2 = -6.315_dbl
            end if

            if (calcType .eq. 5) then
                alpha1 = 2.265_dbl
                alpha2 = -6.911_dbl
            end if

            if (calcType .eq. 6) then
                alpha1 = 2.260_dbl
                alpha2 = -6.923_dbl
            end if

            ln_P0 = log(prepend) + log(1.0d45) + log(rho) + log(A1*A2) + &
                    2.0_dbl*log(Z1) + 2.0_dbl*log(Z2) + ln_S + &
                    (4.0_dbl/7.0_dbl)*log(lambda) - alpha2 - alpha1*lambda**(-0.5_dbl)

            print *, 'ln_P0      = ', ln_P0
            print *, 'log10_rate = ', ln_P0 / log(10.0d0)

            if (ln_P0 > log(huge(1.0d0))) then
                P0 = huge(1.0d0)
            else if (ln_P0 < log(tiny(1.0d0))) then
                P0 = 0.0d0
            else
                P0 = exp(ln_P0)
            end if

            pycnoRate = P0

        end if

    end function pycnoRate

    end module rate_calc3
