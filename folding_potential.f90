!
!	folding_potential.f90
!	PycnoCalc
!
!	Created by hellmersjl on 10/27/07.
!
module folding_potential

use constants

implicit none

save

contains

	real(kind=dbl) function vfold_spherically_symmetric (R, A1, A2, Z1, Z2, partition,diffuse1,diffuse2,rho0_1, rho0_2,tot_radius1, tot_radius2, E0, nucIntType, inSQMFlag)
	
	! R  : Distance between nuclei in fermi (10e-15 m)
	! A1 : The number of nucleons in the first nucleus
	! A2 : The number of nucleons in the second nucleus
	! Z1 : The number of protons in the first nucleus
	! Z2 : The number of protons in the second nucleus
	! partition : The number of increments to divide the range of integrations
	! diffuse1 : the diffuseness parameter for the first nucleus
	! diffuse2 : the diffuseness paramter for the second nucleus
	! rho0_1 : the central density of the first nucleus
	! rho0_2 : the central density of the second nucleus
	! tot_radius1 : the effective radius of the first nucleus
	! tot_radius2 : the effective radius of the second nucleus
    ! E0: zero-point vibrational energy for first (incoming particle)
    ! inSQMFlag: An optional input parameter indicating if we are doing an SQM calculation
	
	! returns the Folding potential in MeV
	
	use constants
	use nucleon_interactions
	use general_nuclear
    use configuration
	
	implicit none
	
	real(kind=dbl), intent(in)      :: R
	integer, intent(in)             :: A1
	real(kind=dbl), intent(in)      :: A2
	integer, intent(in)             :: Z1
	real(kind=dbl), intent(in)      :: Z2
	integer, intent(in)             :: partition
	real(kind=dbl), intent(in)      :: diffuse1
	real(kind=dbl), intent(in)      :: diffuse2
	real(kind=dbl), intent(in)      :: rho0_1
	real(kind=dbl), intent(in)      :: rho0_2
	real(kind=dbl), intent(in)      :: tot_radius1
	real(kind=dbl), intent(in)      :: tot_radius2
	real(kind=dbl), intent(in)      :: E0
	integer, intent(in)             :: nucIntType
	logical, intent(in), optional   :: inSQMFlag
	
	integer :: r1
	integer :: r2
	integer :: theta1
	integer :: theta2
	integer :: phi1
	integer :: phi2
    character(len=1) :: foldingNuclearInteraction

	real(kind=dbl) :: delta_r1
	real(kind=dbl) :: delta_r2
	real(kind=dbl) :: delta_theta1
	real(kind=dbl) :: delta_phi1
	real(kind=dbl) :: delta_theta2
	real(kind=dbl) :: delta_phi2
	real(kind=dbl) :: dV1
	real(kind=dbl) :: dV2
	real(kind=dbl) :: d
	
	real(kind=dbl) :: accumulator1
	real(kind=dbl) :: accumulator2
	
	real(kind=dbl) :: radius1
	real(kind=dbl) :: radius2
	
	
	
	real(kind=dbl) :: ndensity1
	real(kind=dbl) :: ndensity2
	real(kind=dbl) :: this_r1
	real(kind=dbl) :: this_theta1
	real(kind=dbl) :: this_phi1
	real(kind=dbl) :: this_r2
	real(kind=dbl) :: this_theta2
	real(kind=dbl) :: this_phi2
	real(kind=dbl) :: dx
	real(kind=dbl) :: dy
	real(kind=dbl) :: dz
	real(kind=dbl) :: mu
	
	real(kind=dbl) :: my_cnt
	real(kind=dbl) :: pi_div_2
    logical :: SQMFlag = .FALSE.

    if (present(inSQMFlag)) then
        if (inSQMFlag) then
            SQMFlag = .TRUE.
        end if
    end if

    mu = reduced_mass(A1,Z1,A2,Z2)
	pi_div_2 = pi/2.0_dbl
	
    radius1 = nuclear_radius(real(A1,dbl), .FALSE.)
    radius2 = nuclear_radius(A2, SQMFlag)
	
	accumulator1 = 0
	accumulator2 = 0
		
	!print *,'"Exact" Volume1 = ',(1.0/3.0)*pi*tot_radius1**3
	
	delta_r1=tot_radius1/partition
	delta_r2=tot_radius2/partition
	delta_theta1 = pi/partition
	delta_phi1 = 2.0_dbl*pi/partition
    delta_theta2 = pi/partition
	delta_phi2 = 2.0_dbl*pi/partition
	
    ! If this is an SQM calculation we just need to do the 2nd nuclei/nugget density once
    if (SQMFlag) then
        ndensity2 = number_density(A2, radius2)
    end if

	my_cnt = 0
    ! foldingNuclearInteraction = getParamValue("FoldingSimple","nuc_interaction_type")

	do r1=1,partition
		this_r1 = r1*delta_r1-delta_r1/2.0
		ndensity1 = density_2pF(rho0_1,this_r1,radius1,diffuse1,.FALSE.)
		do theta1=1,partition
			this_theta1 = theta1*delta_theta1-delta_theta1/2.0
			do phi1=1,partition
				this_phi1 = phi1*delta_phi1 - delta_phi1/2.0 - pi_div_2
				accumulator2 = accumulator2 + this_r1*this_r1*sin(this_theta1)*delta_theta1*delta_phi1*delta_r1
				do r2=1,partition
					this_r2 = r2*delta_r2-delta_r2/2.0
                    if (.not. SQMFlag) then					
                        ndensity2 = density_2pF(rho0_2,this_r2,radius2,diffuse2,.FALSE.)
                    end if
					do theta2=1,partition
						this_theta2 = theta2*delta_theta2-delta_theta2/2.0
						do phi2=1, partition
							this_phi2 = phi2*delta_phi2 - delta_phi2/2.0
							dV1 = this_r1*this_r1*sin(this_theta1)*delta_theta1*delta_phi1*delta_r1
							dV2 = (this_r2)*(this_r2)*sin(this_theta2)*delta_theta2*delta_phi2*delta_r2
							dx = R + this_r2*sin(this_theta2)*cos(this_phi2) - this_r1*sin(this_theta1)*cos(this_phi1) 
							dy = this_r2*sin(this_theta2)*sin(this_phi2) - this_r1*sin(this_theta1)*sin(this_phi1)
							dz = this_r2*cos(this_theta2) - this_r1*cos(this_theta1)
							d = sqrt(dx*dx + dy*dy + dz*dz)
                            if (nucIntType .eq. 1) then
							    accumulator1 = accumulator1 + ndensity1*ndensity2*dV1*dV2*nucleonSaoPaulo(d,E0,mu)
                            else
							    accumulator1 = accumulator1 + ndensity1*ndensity2*dV1*dV2*nucleonM3Y(d,(delta_r1+delta_r2),.true.)
                            end if
							my_cnt = my_cnt + 1
						end do
					end do
				end do
			end do
		end do
	end do
	
	!print *,'Calculated Volume1 = ',accumulator2
	
    vfold_spherically_symmetric = accumulator1	
	
	end function vfold_spherically_symmetric
	

    real(kind=dbl) function vfold_cubes(R, A1, A2, Z1, Z2, N1, N2, side1, side2, E0, nucIntType, inSQMFlag)
    
	use nucleon_interactions
	use general_nuclear
 
    real(kind=dbl), intent(in)      :: R            ! Distance between centers of cubes
	integer, intent(in)             :: A1
	real(kind=dbl), intent(in)      :: A2
	integer, intent(in)             :: Z1
	real(kind=dbl), intent(in)      :: Z2
	integer(kind=dbl), intent(in)             :: N1, N2       ! Number of increments to split the cubes into 
	real(kind=dbl), intent(in)      :: side1        ! Length of side for cube 1
	real(kind=dbl), intent(in)      :: side2        ! Length of side for cube 2
	real(kind=dbl), intent(in)      :: E0
	integer, intent(in)             :: nucIntType
	logical, intent(in), optional   :: inSQMFlag

    real(kind=dbl)                  :: dx1, dy1, dz1
    real(kind=dbl)                  :: dx2, dy2, dz2
    real(kind=dbl)                  :: dV1, dV2
    real(kind=dbl)                  :: d
    integer                         :: i
    real(kind=dbl)                  :: accumulator

	real(kind=dbl) :: ndensity1
	real(kind=dbl) :: ndensity2
	real(kind=dbl) :: mu

    dx1 = side1/N1
    dy1 = side1/N1
    dz1 = side1/N1

    dx2 = side2/N2
    dy2 = side2/N2
    dz2 = side2/N2
    
    dV1 = dx1 * dy1 * dz1
	dV2 = dx2 * dy2 * dz2


    ndensity1 = A1/dV1
    ndensity2 = A2/dV2
   
    mu = reduced_mass(A1,Z1,A2,Z2)

    print *,"dx1=",dx1
    print *,"side1=",side1

!$OMP PARALLEL DO
    do i = 1, N1+1
        d = 0.0_dbl - side1/2.0_dbl + (i-1)*dx1
        ! print *,i,":d=",d
        if (nucIntType .eq. 1) then
	        accumulator = accumulator + ndensity1*ndensity2*dV1*dV2*nucleonSaoPaulo(d,E0,mu)
        else
	        accumulator = accumulator + ndensity1*ndensity2*dV1*dV2*nucleonM3Y(d,(dx1+dx2))
        end if
    end do
!$OMP END PARALLEL DO

    vfold_cubes = accumulator

    end function vfold_cubes

end module folding_potential
