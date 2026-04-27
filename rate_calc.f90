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

    subroutine turn_pt(R,Rstep,Rmax,E0,mu,A1,A2,Z1,Z2,radius1,radius2,rho0_A1,rho0_A2,L,partition,turn1,turn2,WKB,nucIntType,inSQMFlag,screening_model,screen_rho,screen_temp)

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
        real(kind=dbl), intent(out)     :: WKB          ! value of WKB calculation to get through total barrier (not just coulomb barrier) between incoming and target nuclei
        integer, intent(in)             :: nucIntType   ! Type of nuclear interaction to use
        logical, intent(in), optional   :: inSQMFlag    ! Indicate if we are using SQM for nuclei/nugget 2
        integer, intent(in), optional        :: screening_model ! screening selector: 0 = none, 2 = mean-field screened barrier
        real(kind=dbl), intent(in), optional :: screen_rho      ! mass density [g/cm^3] used by the screening model
        real(kind=dbl), intent(in), optional :: screen_temp     ! temperature [K] used by the screening model

        real(kind=dbl)  :: ln_WKB                    ! natural log of the value of the WKB calculation to get through the total barrier
        real(kind=dbl)  :: V_2fold                   ! double folding potential calculated with distance of R between incoming and target nuclei
        real(kind=dbl)  :: ndensity1, ndensity2      ! number density of nuclei 1 and 2, calculated using baryon number
        real(kind=dbl)  :: VEcheck(0:Rmax)           ! array to hold difference between Veff and E of incoming particle at every point R (between incoming & target nuclei)
        real(kind=dbl)  :: Vnucarray(0:Rmax)         ! array to hold Veff between target & incoming particle at every point along R (between incoming & target nuclei)
        real(kind=dbl)  :: Vcoulary(0:Rmax)          ! array to hold Vcoulomb between target and incoming particle at every point along R (between incoming & target nuclei)
        integer         :: i, turn_counter           ! ad hoc counter, counter for the number of turning points 
        real(kind=dbl)  :: R_pos                     ! current position along R axis --this is the CURRENT separation between target and incoming particle
        real(kind=dbl)  :: Energy                    ! function that calculates the Energy of incoming (or just second (projectile) nucleon)
        real(kind=dbl)  :: Integrand(0:Rmax)        ! array to hold Integrand of S-factor calculation
        real(kind=dbl)  :: S                         ! astrophysical S-factor - calculated using equation from PHYS REV C69 --rule of thumb model fitted to data
        real(kind=dbl)  :: ln_S                      ! natural log of astrophysical S-factor
        integer         :: new_Rmax                  ! new Rmax based upon cutoff
        logical         :: SQMFlag = .FALSE.         ! Used internally
        integer         :: local_screening_model     ! local copy of screening selector, defaults to 0 (unscreened)
        real(kind=dbl)  :: local_screen_rho          ! local density [g/cm^3] for screening calculations
        real(kind=dbl)  :: local_screen_temp         ! local temperature [K] for screening calculations

        logical         :: v_fold_threshold_flg
        logical         :: v_fold_first_time
        real(kind=dbl)  :: v_fold_max
        real(kind=dbl)  :: v_fold_min

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        ! Set defaults for optional screening inputs so older unscreened calls still work
        local_screening_model = 0
        if (present(screening_model)) local_screening_model = screening_model

        local_screen_rho = 0.0_dbl
        if (present(screen_rho)) local_screen_rho = screen_rho

        local_screen_temp = 0.0_dbl
        if (present(screen_temp)) local_screen_temp = screen_temp

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

        do i = 0,Rmax
            R_pos = i * Rstep
            if (.not. v_fold_threshold_flg) then 
                V_2fold = vfold_spherically_symmetric (R_pos, A1, A2, Z1, Z2, partition,0.5_dbl,0.5_dbl,rho0_A1, rho0_A2, radius1, radius2, E0, nucIntType, SQMFlag)
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

            Vnucarray(Rmax-i) = V_2fold

            if (local_screening_model == 2 .and. local_screen_rho > 0.0_dbl .and. local_screen_temp > 0.0_dbl) then
                Vcoulary(Rmax-i) = U_barrier_mean_field(Z1, Z2, R_pos, radius1, radius2, &
                                                        local_screen_rho, local_screen_temp, A1, Z1)
            else
                Vcoulary(Rmax-i) = Vcoulomb(Z1, Z2, R_pos, radius1, radius2)
            end if

            VEcheck(Rmax-i) = V_2fold + Vcoulary(Rmax-i) - E0
        end do
	  
!	now check the VEcheck array for turning point(s) - there may be more than one - especially if the total energy of the incoming particle starts out high (greater than the potential)
        do i = 1,Rmax
            if ((VEcheck(i).le.0).and.(VEcheck(i-1).ge.0)) then
                turn_counter = turn_counter + 1
                if (turn1.lt.0) then
                    turn1 = i - 1
                else
                    turn2 = i - 1
                end if
            else 
                if ((VEcheck(i).ge.0).and.(VEcheck(i-1).le.0)) then
                    turn_counter = turn_counter + 1
                    if (turn1.lt.0) then
                        turn1 = i - 1
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
            do i = 0,Rmax
                if (i.gt.turn1) then
                    Integrand(i) = 0.0_dbl
                else
                    if (VEcheck(i) .ge. 0.0_dbl) then
                        Integrand(i) = .01432_dbl*sqrt(mu*VEcheck(i))
                    else
                        Integrand(i) = -0.01432_dbl*sqrt(-1.0_dbl*mu*VEcheck(i))
                    end if
                end if
            end do
        else
            do i = 0,Rmax
                if ((i.le.turn1).or.(i.gt.turn2)) then
                    Integrand(i) = 0.0_dbl
                else
                    if (VEcheck(i) .ge. 0.0_dbl) then
                        Integrand(i) = .01432_dbl*sqrt(mu*VEcheck(i))
                    else
                        Integrand(i) = -0.01432_dbl*sqrt(-1.0_dbl*mu*VEcheck(i))
                    end if
                end if
            end do
        end if

        WKB = trapezoidArray(Rmax,Rmax,Integrand,Rstep)
        
    end subroutine turn_pt


!*****************************************************************
!
! Calculate the S Factor
!
! REF: Golf Thesis eqn 3.32
!
!*****************************************************************

    real(kind=dbl) function Sfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag, screening_model, screen_temp)

        implicit none
        
        integer, intent(in)             :: A1_int, Z1_int
        real(kind=dbl), intent(in)      :: A2, Z2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Rstep
        integer, intent(in)             :: partition
        integer, intent(in)             :: nucIntType
        logical, intent(in), optional   :: inSQMFlag
        integer, intent(in), optional        :: screening_model ! screening selector passed to turn_pt
        real(kind=dbl), intent(in), optional :: screen_temp     ! temperature [K] for screened barrier evaluation

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
        integer         :: local_screening_model          ! local copy of screening selector, defaults to 0
        real(kind=dbl)  :: local_screen_temp             ! local temperature [K] used in screened barrier calculation

        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        local_screening_model = 0
        if (present(screening_model)) local_screening_model = screening_model

        local_screen_temp = 0.0_dbl
        if (present(screen_temp)) local_screen_temp = screen_temp

        A1 = real(A1_int,dbl)
        Z1 = real(Z1_int,dbl)

        R = sqrt(3.0_dbl)*0.5_dbl*lattice(rho,A1,Z1)
        Rmax = (R/Rstep)+1
        mu = reduced_mass(A1_int, Z1_int, A2, Z2)
        
        radius1 = nuclear_radius(A1,.FALSE.)
        radius2 = nuclear_radius(A2, SQMFlag)
        rho0_A1 = rho0_2pF(A1,radius1,0.5_dbl)
        rho0_A2 = rho0_2pF(A2,radius2,0.5_dbl)

        E0 = E0_Energy(Z1_int,Z2, A1_int, A2, rho)

        ln_sigma = 0.0
        do i = 0,L
            call turn_pt(R,Rstep,Rmax,E0,mu,A1_int,A2,Z1_int,Z2,radius1,radius2, &
                         rho0_A1,rho0_A2,L,partition,turn1,turn2,WKB,nucIntType,SQMFlag, &
                         local_screening_model, rho, local_screen_temp)
            ln_Trans_total = -WKB
            ln_sigma=ln_sigma+log(612.459_dbl)-log(mu*E0)+log(2.0_dbl*real(i,dbl)+1.0_dbl)+ln_Trans_total
        end do

        ln_S=ln_sigma+log(E0)+(Z1)*(Z2)*.0324_dbl*sqrt(mu/E0)
        Sfactor = exp(ln_S)
    end function Sfactor

!*****************************************************************
!
! Calculate Pycnonuclear reaction rates
! This version also calculates the Folding potential
!
! REF: Original source Golf Thesis
!      Updated from Schramm and Koonin 1990, an equation near eq. 4, eq. 33, and table I
!
!*****************************************************************

    real(kind=dbl) function pycnoRate(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, inSQMFlag, inCalcType, screening_model, screen_temp)

        implicit none

        integer, intent(in)             :: A1_int, Z1_int
        real(kind=dbl), intent(in)      :: A2, Z2
        real(kind=dbl), intent(in)      :: rho
        real(kind=dbl), intent(in)      :: Rstep
        integer, intent(in)             :: partition
        integer, intent(in)             :: nucIntType
        logical, intent(in), optional   :: inSQMFlag
        integer, intent(in), optional   :: inCalcType 
            !  0 - Original Calculation from Golf
            !  1 - bcc, static
            !  2 - bcc, wigner-sietz cell
            !  3 - bcc, relaxed
            !  4 - fcc, static
            !  5 - fcc, wigner-sietz cell
            !  6 - fcc, relaxed 
        integer, intent(in), optional        :: screening_model ! screening selector passed down to Sfactor / turn_pt
        real(kind=dbl), intent(in), optional :: screen_temp     ! temperature [K] used by the mean-field screening model
        
        real(kind=dbl)                  :: prepend, alpha1, alpha2, gamma

        real(kind=dbl)              :: ln_P0
        real(kind=dbl)              :: P0
        real(kind=dbl)              :: lambda
        real(kind=dbl)              :: ln_lambda
        real(kind=dbl)              :: ln_S
        real(kind=dbl)              :: S
        real(kind=dbl)              :: mn_mass, mn_chrg

        real(kind=dbl)              :: A1, Z1

        logical         :: SQMFlag = .FALSE.
        integer         :: calcType = 0
        integer         :: local_screening_model          ! local copy of screening selector, defaults to 0
        real(kind=dbl)  :: local_screen_temp             ! local temperature [K] for screened barrier evaluation
        
        if (present(inSQMFlag)) then
            if (inSQMFlag) then
                SQMFlag = .TRUE.
            end if
        end if

        if (present(inCalcType)) then
            calcType = inCalcType
        end if

        local_screening_model = 0
        if (present(screening_model)) local_screening_model = screening_model

        local_screen_temp = 0.0_dbl
        if (present(screen_temp)) local_screen_temp = screen_temp

        A1 = real(A1_int,dbl)
        Z1 = real(Z1_int,dbl)

        S = Sfactor(A1_int, A2, Z1_int, Z2, rho, Rstep, partition, nucIntType, SQMFlag, &
                    local_screening_model, local_screen_temp)
        ln_S = log(S)
        
        if (calcType .eq. 0) then
            lambda = inv_len_param2comp (rho, A1_int, Z1_int, 0.5_dbl, int(A2), int(Z2), 0.5_dbl)
            ln_lambda = log(lambda)
            ln_P0=-2.638_dbl/(sqrt(lambda))+log(rho)+log(A1*A2)-log(A1+A2)+2.0_dbl*log(Z1*Z2)+ln_S+1.75_dbl*ln_lambda+109.36_dbl        
            P0 = exp(ln_P0)
            pycnoRate = P0
        else
        
            if (calcType .eq. 1 .or. calcType .eq. 2 .or. calcType .eq. 3) then
                prepend = 1.06_dbl
                gamma = 2.0_dbl
            else
                prepend = 2.69_dbl
                gamma = 4.0_dbl
            end if

            lambda = 0.0245*(A1**(-4.0_dbl/3.0_dbl))*(Z1**(-2.0_dbl))*(gamma**(-1.0_dbl/3.0_dbl))*((rho/(1.0d6))**(1.0_dbl/3.0_dbl))

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
                alpha2 = -6.841_dbl
            end if
        
            P0 = prepend * 1.0d45 * S * rho*(A1*A2)*(Z1*Z1)*(Z2*Z2)*(lambda**(7.0_dbl/4.0_dbl))*exp(-alpha2-alpha1*lambda**(-1.0_dbl/2.0_dbl))
        end if
        pycnoRate = P0 

    end function pycnoRate

!*****************************************************************
!
! Calculate the adjustment to the Pycno Reaction rate based on 
! Temperature
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

      real(dbl) function pow3(x) result(res) ! x^3
         real(dbl), intent(in) :: x
         res = x*x*x
      end function pow3

      real(dbl) function pow4(x) result(res) ! x^4
         real(dbl), intent(in) :: x
         res = x*x*x*x
      end function pow4

      real(kind=dbl) function log_cr(x) result(res)
         real(kind=dbl), intent(in) :: x
         res = log(x)
      end function log_cr

      real(dbl) function exp_cr(x) result(res) ! E^x
         real(dbl), intent(in) :: x
         res = exp(x)
      end function exp_cr

      real(kind=dbl) function pow_cr(x,y) result(res) ! x^y
         real(kind=dbl), intent(in) :: x,y
         integer :: iy, j
         if (x == 0d0) then
            res = 0d0
            return
         end if
         iy = floor(y)
         if (y == dble(iy) .and. abs(iy) < 100) then ! integer power of x
            res = 1d0
            do j=1,abs(iy)
               res = res*x
            end do
            if (iy < 0) res = 1d0/res
            return
         end if
         res = exp_cr(log_cr(x)*y)
      end function pow_cr

      subroutine pycnoRateG05_CC(T, Rho, X12, eps, deps_dT, deps_dRho, retRpyc)
         
         ! from Gasques, et al, Nuclear fusion in dense matter.  astro-ph/0506386.
         ! Phys Review C, 72, 025806 (2005)
         
         real(dbl), intent(in) :: T
         real(dbl), intent(in) :: Rho
         real(dbl), intent(in) :: X12 ! mass fraction of c12
         real(dbl), intent(out) :: eps ! rate in ergs/g/sec
         real(dbl), intent(out) :: deps_dT ! partial wrt temperature
         real(dbl), intent(out) :: deps_dRho ! partial wrt density
         real(dbl), intent(out) :: retRpyc ! Pycnonuclear Reaction rate
         
         real(dbl), parameter :: exp_cutoff = 200d0
         real(dbl), parameter :: exp_max_result = 1d200
         real(dbl), parameter :: exp_min_result = 1d-200

         real(dbl), parameter :: sqrt3 = 1.73205080756888d0 ! sqrt[3]
         
         real(dbl), parameter :: Z = 6 ! charge of C12
         real(dbl), parameter :: A = 12 ! atomic number for C12

         real(dbl), parameter :: b_ne20 = 160.64788d0
         real(dbl), parameter :: b_c12  =  92.1624d0
         real(dbl), parameter :: b_he4  =  28.2928d0
         
         ! coefficients from "optimal" model 1
         real(dbl), parameter :: Cexp = 2.638
         real(dbl), parameter :: Cpyc = 3.90
         real(dbl), parameter :: Cpl = 1.25
         real(dbl), parameter :: CT = 0.724
         real(dbl), parameter :: Csc = 1.0754
         real(dbl), parameter :: Lambda = 0.5  
         
         real(dbl) :: m, Zqe, Z2e2, Q1212, MeV_to_erg, barn_to_cm2, MeV_barn_to_erg_cm2, ergs_c12c12
         
         real(dbl) :: Tk, ni, a_ion, Gamma, omega_p, Tp
         real(dbl) :: Ea, tau, Epk, SEpk, lam, Rpyc, Ppyc, Fpyc
         real(dbl) :: Ttilda, Gammatilda, tautilda, phi, gam1, gam2, gam
         real(dbl) :: P, lnF, F, deltaR, R

         real(dbl) :: Epk1, Epk2x, SEpkx1, SEpkx2x, SEpkx2d, SEpkx2n, SEpkx2, phin, phid
         real(dbl) :: Epk2c, dEpk1_dT, dni_dRho, dTk_dT, da_ion_dRho, dGamma_dT, &
               dGamma_dRho, domega_p_dRho, dTp_dRho
         real(dbl) :: dtau_dT, dEpk1_dRho, dEpk2x_dT, dEpk2x_dRho, dEpk2c_dT, dEpk2c_dRho, dEpk_dT, dEpk_dRho
         real(dbl) :: dSEpkx1_dT, dSEpkx1_dRho, dSEpkx2x_dT, dSEpkx2x_dRho, dSEpkx2d_dT, dSEpkx2d_dRho
         real(dbl) :: dSEpkx2n_dT, dSEpkx2n_dRho, dSEpkx2_dT, dSEpkx2_dRho, dSEpk_dT, dSEpk_dRho
         real(dbl) :: dlam_dRho, dPpyc_dRho, dFpyc_dRho, dRpyc_dT, dRpyc_dRho, gam_n, gam_d, dgam_dT, dgam_dRho
         real(dbl) :: dphin_dGamma, dphin_dT, dphin_dRho, dphid_dGamma,  &
               dphid_dT, dphid_dRho, dphidT, dphidRho
         real(dbl) :: dTtilda_dT, dTtilda_dRho, dGammatilda_dT, dGammatilda_dRho,  &
               dtautilda_dTtilda, dtautilda_dT, dtautilda_dRho
         real(dbl) :: dgam_n_dT, dgam_n_dRho, dgam_d_dT, dgam_d_dRho, dP_dgam, dP_dTtilda, dP_dT, dP_dRho
         real(dbl) :: lnF1, lnF2x, lnF2, lnF3, dlnF1_dT, dlnF1_dRho, dlnF2x_dT, dlnF2x_dRho, dlnF2_dT, dlnF2_dRho
         real(dbl) :: dlnF3_dT, dlnF3_dRho, dlnF_dT, dlnF_dRho, dF_dT, dF_dRho, niterm, dniterm_dRho
         real(dbl) :: d_deltaR_dT, d_deltaR_dRho, dRdT, dRdRho, exponent
                   
         integer :: status, pwdlen
         integer :: unit = 45
         character(len=128) :: pwd, pycnoTypeFile 

 1       format(a55,99(1pd26.16))
                  
         m = A * amu
         Zqe = Z*elementary_charge
         Z2e2 = Zqe*Zqe
         Q1212 = b_ne20 + b_he4 - 2 * b_c12
         MeV_to_erg = 1d6 * ev2erg
         barn_to_cm2 = 1d-24
         MeV_barn_to_erg_cm2 = MeV_to_erg * barn_to_cm2
         ergs_c12c12 = MeV_to_erg * Q1212

         Tk = T * boltzmann
         dTk_dT = boltzmann

         ni = Rho * avogadro / A
         dni_dRho = avogadro / A
         
         a_ion = pow_cr(3 / (4 * Pi * ni), one_third)
         da_ion_dRho = -one_third * a_ion * dni_dRho / ni
         
         Gamma = Z2e2 / (a_ion * Tk)
         dGamma_dT = - Gamma * dTk_dT / Tk
         dGamma_dRho =  - Gamma * da_ion_dRho / a_ion
         
         omega_p = DSQRT(4 * Pi * Z2e2 * ni / m)
         domega_p_dRho = 0.5d0 * omega_p * dni_dRho / ni
         
         Tp = hbar * omega_p
         dTp_dRho = hbar * domega_p_dRho
         
         Ea = m * (Z2e2 / hbar)**2
         
         tau = 3 * pow_cr(Pi/2,two_thirds) * pow_cr(Ea/Tk, one_third)
         dtau_dT = -one_third * tau * dTk_dT / Tk
         
         Epk1 = hbar * omega_p
         dEpk1_dRho = hbar * domega_p_dRho
         dEpk1_dT = 0
         
         exponent = -Lambda*Tp/Tk
         if (exponent < -exp_cutoff) then
            Epk2x = exp_min_result
            dEpk2x_dT = 0
            dEpk2x_dRho = 0
         else if (exponent > exp_cutoff) then
            Epk2x = exp_max_result
            dEpk2x_dT = 0
            dEpk2x_dRho = 0
         else
            Epk2x = exp_cr(exponent)
            dEpk2x_dT = Epk2x * Lambda * Tp * dTk_dT / (Tk*Tk)
            dEpk2x_dRho = -Epk2x * Lambda * Tp * dTp_dRho / Tk
         end if
         
         Epk2c = Z2e2/a_ion + Tk*tau/3
         dEpk2c_dT = (dTk_dT*tau + Tk*dtau_dT)/3
         dEpk2c_dRho = -Z2e2*da_ion_dRho/(a_ion*a_ion)
         
         Epk = (Epk1 + Epk2c * Epk2x) / MeV_to_erg
         dEpk_dT = (dEpk1_dT + dEpk2c_dT * Epk2x + Epk2c * dEpk2x_dT) / MeV_to_erg
         dEpk_dRho = (dEpk1_dRho + dEpk2c_dRho * Epk2x + Epk2c * dEpk2x_dRho) / MeV_to_erg

         SEpkx1 = -0.428*Epk
         dSEpkx1_dT = -0.428*dEpk_dT
         dSEpkx1_dRho = -0.428*dEpk_dRho
         
         SEpkx2x = 0.613*(8-Epk)
         dSEpkx2x_dT = -0.613*dEpk_dT
         dSEpkx2x_dRho = -0.613*dEpk_dRho
         
         if (SEpkx2x < -exp_cutoff) then
            SEpkx2d = 1 + exp_min_result
            dSEpkx2d_dT = 0
            dSEpkx2d_dRho = 0
         else if (SEpkx2x > exp_cutoff) then
            SEpkx2d = 1 + exp_max_result
            dSEpkx2d_dT = 0
            dSEpkx2d_dRho = 0
         else
            SEpkx2d = 1 + exp_cr(SEpkx2x)
            dSEpkx2d_dT = (SEpkx2d-1)*dSEpkx2x_dT
            dSEpkx2d_dRho = (SEpkx2d-1)*dSEpkx2x_dRho
         end if
         
         SEpkx2n = 3 * pow_cr(Epk,0.308d0)
         dSEpkx2n_dT = 3 * 0.308 * dEpk_dT / Epk
         dSEpkx2n_dRho = 3 * 0.308 * dEpk_dRho / Epk
         
         SEpkx2 = SEpkx2n / SEpkx2d
         dSEpkx2_dT = (dSEpkx2n_dT - SEpkx2n * dSEpkx2d_dT / SEpkx2d) / SEpkx2d
         dSEpkx2_dRho = (dSEpkx2n_dRho - SEpkx2n * dSEpkx2d_dRho / SEpkx2d) / SEpkx2d
         
         exponent = SEpkx1+SEpkx2
         if (exponent < -exp_cutoff) then
            SEpk = 5.15d16 * exp_min_result
            dSEpk_dT = 0
            dSEpk_dRho = 0
         else if (exponent > exp_cutoff) then
            SEpk = 5.15d16 * exp_max_result
            dSEpk_dT = 0
            dSEpk_dRho = 0
         else
            SEpk = 5.15d16 * exp_cr(exponent)
            dSEpk_dT = SEpk * (dSEpkx1_dT + dSEpkx2_dT)
            dSEpk_dRho = SEpk * (dSEpkx1_dRho + dSEpkx2_dRho)
         end if

         lam = pow_cr(Rho/(1.3574d11*A),one_third) / (A*Z*Z)
         dlam_dRho = one_third * lam / Rho
         
         Ppyc = 8*Cpyc*11.515/pow_cr(lam,Cpl)
         dPpyc_dRho = - Ppyc * dlam_dRho / lam
         
         Fpyc = exp_cr(-Cexp/DSQRT(lam))
         dFpyc_dRho = Fpyc * Cexp * dlam_dRho / (2 * lam*sqrt(lam))
         
         Rpyc = Rho * A * pow4(Z) * (1 / (8 * 11.515)) * 1d46 * pow3(lam) * SEpk * Fpyc * Ppyc
         
         if (Rpyc < exp_min_result) then
            Rpyc = 0
            dRpyc_dT = 0
            dRpyc_dRho = 0
         else
            dRpyc_dT = Rpyc * (dSEpk_dT / SEpk)
            dRpyc_dRho = Rpyc * (3 * dlam_dRho / lam + dSEpk_dRho / SEpk + dFpyc_dRho / Fpyc + dPpyc_dRho / Ppyc)
         end if

         phin = DSQRT(Gamma)
         dphin_dGamma = 0.5 * phin / Gamma
         dphin_dT = dGamma_dT * dphin_dGamma
         dphin_dRho = dGamma_dRho * dphin_dGamma
         
         phid = pow_cr(Csc**4 / 9 + Gamma**2, 0.25d0)
         dphid_dGamma = 0.5d0 * Gamma / pow3(phid)
         dphid_dT = dGamma_dT * dphid_dGamma
         dphid_dRho = dGamma_dRho *dphid_dGamma
         
         phi = phin / phid
         dphidT = (dphin_dT - dphid_dT * phin / phid) / phid
         dphidRho = (dphin_dRho - dphid_dRho * phin / phid) / phid
         
         Ttilda = DSQRT(Tk*Tk + (CT*Tp)*(CT*Tp))
         dTtilda_dT = Tk * dTk_dT / Ttilda
         dTtilda_dRho = CT*CT * Tp * dTp_dRho / Ttilda
         
         Gammatilda = Z2e2 / (a_ion * Ttilda)
         dGammatilda_dT = -Gammatilda * dTtilda_dT / Ttilda
         dGammatilda_dRho = -Gammatilda * (dTtilda_dRho / Ttilda + da_ion_dRho / a_ion)
         
         tautilda = 3 * pow_cr(Pi/2,two_thirds) * pow_cr(Ea / Ttilda, one_third)
         dtautilda_dTtilda = - tautilda / (3 * Ttilda)
         dtautilda_dT = dtautilda_dTtilda * dTtilda_dT
         dtautilda_dRho = dtautilda_dTtilda * dTtilda_dRho
         
         gam1 = two_thirds
         gam2 = two_thirds * (Cpl + 0.5)
         gam_n = Tk*Tk * gam1 + Tp*Tp * gam2
         dgam_n_dT = 2 * Tk * dTk_dT * gam1
         dgam_n_dRho = 2 * Tp * dTp_dRho * gam2
         
         gam_d = Tk*Tk + Tp*Tp
         dgam_d_dT = 2 * Tk * dTk_dT
         dgam_d_dRho = 2 * Tp * dTp_dRho

         gam = gam_n / gam_d
         dgam_dT = (dgam_n_dT - gam_n * dgam_d_dT / gam_d)/gam_d
         dgam_dRho = (dgam_n_dRho - gam_n * dgam_d_dRho / gam_d)/gam_d
         
         P = (8 * pow_cr(Pi/2,one_third)/sqrt3) * pow_cr(Ea/Ttilda, gam)
         if (P < exp_min_result) then
            P = 0
            dP_dgam = 0
            dP_dTtilda = 0
            dP_dT = 0
            dP_dRho = 0
         else
            dP_dgam = P * log_cr(Ea / Ttilda)
            dP_dTtilda = -P * gam / Ttilda
            dP_dT = dP_dgam * dgam_dT + dP_dTtilda * dTtilda_dT
            dP_dRho = dP_dgam * dgam_dRho + dP_dTtilda * dTtilda_dRho
         end if
         
         lnF1 = -tautilda
         exponent = -Lambda * Tp / Tk
         if (exponent < -exp_cutoff) then
            lnF2x = exp_min_result
            dlnF2x_dT = 0
            dlnF2x_dRho = 0
         else if (exponent > exp_cutoff) then
            lnF2x = exp_max_result
            dlnF2x_dT = 0
            dlnF2x_dRho = 0
         else
            lnF2x = exp_cr(exponent)
            dlnF2x_dT = lnF2x * Lambda * Tp * dTk_dT / Tk**2
            dlnF2x_dRho = -lnF2x * Lambda * dTp_dRho / Tk
         end if
         
         lnF2 = Csc * Gammatilda * phi * lnF2x
         lnF3 = - Lambda * Tp / Tk
         lnF = lnF1 + lnF2 + lnF3
         if (lnF < -exp_cutoff) then
            F = exp_min_result
            dF_dT = 0
            dF_dRho = 0
         else if (lnF > exp_cutoff) then
            F = exp_max_result
            dF_dT = 0
            dF_dRho = 0
         else
            F = exp_cr(lnF)
            dlnF1_dT = -dtautilda_dT
            dlnF1_dRho = -dtautilda_dRho
            dlnF2_dT = lnF2 * (dGammatilda_dT / Gammatilda + dphidT / phi + dlnF2x_dT / lnF2x)
            dlnF2_dRho = lnF2 * (dGammatilda_dRho / Gammatilda + dphidRho / phi + dlnF2x_dRho / lnF2x)
            dlnF3_dT = -lnF3 * dTk_dT / Tk
            dlnF3_dRho = lnF3 * dTp_dRho / Tp
            dlnF_dT = dlnF1_dT + dlnF2_dT + dlnF3_dT
            dlnF_dRho = dlnF1_dRho + dlnF2_dRho + dlnF3_dRho
            dF_dT = dlnF_dT * F
            dF_dRho = dlnF_dRho * F
         end if

         niterm = ni*ni * (hbar / (2 * m * Z2e2))
         dniterm_dRho = 2 * niterm * dni_dRho / ni
         
         deltaR = niterm * SEpk * MeV_barn_to_erg_cm2 * P * F
         if (deltaR < 1d-100) then
            deltaR = 0
            d_deltaR_dT = 0
            d_deltaR_dRho = 0
         else
            d_deltaR_dT = deltaR * (dSEpk_dT / SEpk + dP_dT / P + dF_dT / F)
            d_deltaR_dRho = deltaR * (dniterm_dRho / niterm + dSEpk_dRho / SEpk + dP_dRho / P + dF_dRho / F)
         end if
         
         R = X12 * (Rpyc + deltaR)
         dRdT = X12 * (dRpyc_dT + d_deltaR_dT)
         dRdRho = X12 * (dRpyc_dRho + d_deltaR_dRho)

         eps = R * ergs_c12c12 / Rho
         deps_dT = dRdT * ergs_c12c12 / Rho
         deps_dRho = (dRdRho * ergs_c12c12 - eps) / Rho
         
         retRpyc = Rpyc
         
      end subroutine pycnoRateG05_CC
    
end module rate_calc
