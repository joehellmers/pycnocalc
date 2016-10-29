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

	real(kind=dbl) function vfold_spherically_symmetric (R, A1, A2, partition,diffuse1,diffuse2,rho0_1, rho0_2,tot_radius1, tot_radius2)
	
	! R  : Distance between nuclei in fermi (10e-15 m)
	! A1 : The number of nucleon's in the first nucleus
	! A2 : The number of nucleon's in the second nucleus
	! partition : The number of increments to divide the range of integrations
	! diffuse1 : the diffuseness parameter for the first nucleus
	! diffuse2 : the diffuseness paramter for the second nucleus
	! rho0_1 : the central density of the first nucleus
	! rho0_2 : the central density of the second nucleus
	! tot_radius1 : the effective radius of the first nucleus
	! tot_radius2 : the effective radius of the second nucleus
	
	! returns the Folding potential in MeV
	
	use constants
	use nucleon_interactions
	use general_nuclear
	
	implicit none
	
	real(kind=dbl) :: R
	integer :: A1
	integer :: A2
	integer :: partition
	real(kind=dbl) :: diffuse1
	real(kind=dbl) :: diffuse2
	real(kind=dbl) :: rho0_1
	real(kind=dbl) :: rho0_2
	real(kind=dbl) :: tot_radius1
	real(kind=dbl) :: tot_radius2
	
	integer :: r1
	integer :: r2
	integer :: theta1
	integer :: theta2
	integer :: phi1
	integer :: phi2
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
	
	
	real(kind=dbl) :: my_cnt
	real(kind=dbl) :: pi_div_2

	pi_div_2 = pi/2.0_dbl
	
	radius1 = nuclear_radius(A1)
	radius2 = nuclear_radius(A2)
	
	accumulator1 = 0
	accumulator2 = 0
	
	! Use these number densities for a uniform distribution
	!ndensity1 = (3.0/4.0)*A1/(pi*tot_radius1**3) ! the number density of nucleus 1 in #/fm^3
	!ndensity2 = (3.0/4.0)*A2/(pi*tot_radius2**3) ! the number density of nucleus 2 in #/fm^3
	
	!print *,'"Exact" Volume1 = ',(1.0/3.0)*pi*tot_radius1**3
	
	delta_r1=tot_radius1/partition
	delta_r2=tot_radius2/partition
	delta_theta1 = ((pi_div_2)/2.0_dbl)/partition
	delta_phi1 = pi/partition
    delta_theta2 = pi/partition
	delta_phi2 = 2*pi/partition
	
	my_cnt = 0

!!$omp parallel &
!!$omp shared ( delta_r1, delta_r2, delta_theta1, delta_phi1, delta_theta2, delta_phi2 ) &
!!$omp shared ( this_r1, ndensity1, this_theta1, this_phi1, this_r2, ndensity2, this_theta2) &
!!$omp private ( r1, theta1, phi1, r2, theta2, phi2 ) &
!!$omp reduction (+:accumulator1,accumulator2)

!!$omp do
	do r1=1,partition
		this_r1 = r1*delta_r1-delta_r1/2.0
		ndensity1 = density_2pF(rho0_1,this_r1,radius1,diffuse1)
		do theta1=1,partition
			this_theta1 = theta1*delta_theta1-delta_theta1/2.0
			do phi1=1,partition
				this_phi1 = phi1*delta_phi1 - delta_phi1/2.0 - pi_div_2
				accumulator2 = accumulator2 + this_r1*this_r1*sin(this_theta1)*delta_theta1*delta_phi1*delta_r1
				do r2=1,partition
					this_r2 = r2*delta_r2-delta_r2/2.0
					ndensity2 = density_2pF(rho0_2,this_r2,radius2,diffuse2)
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
							accumulator1 = accumulator1 + ndensity1*ndensity2*dV1*dV2*nucleonM3Y(d,(delta_r1+delta_r2))							
							my_cnt = my_cnt + 1
						end do
					end do
				end do
			end do
		end do
	end do
!!$omp end do

!!$omp end parallel	
	
	!print *,'Inner Loop Total Count=',my_cnt
	!print *,'Calculated Volume1 = ',accumulator2
	
	
	vfold_spherically_symmetric = -8.0_dbl*accumulator1
	
	
	end function vfold_spherically_symmetric
	
end module folding_potential
