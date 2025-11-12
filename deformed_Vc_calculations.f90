module calculations

    implicit none

    integer, parameter :: dbl = selected_real_kind(15, 307)
    real(kind=dbl), parameter :: pi = acos(-1.0_dbl)
    real(kind=dbl), parameter :: e2 = 1.439964_dbl
    
    contains

    ! Function that calculates the Deformed Coulomb potential
    real(kind=dbl) function Vc(R, B2, theta2, Zp, Zt, At)
    
        implicit none
    
        real(kind=dbl), intent(in) :: R, B2, theta2, Zp, Zt, At
        real(kind=dbl) :: rt, r0t, P2, c
    
        c = cos(theta2)
        P2 = 0.5_dbl*(3.0_dbl*c*c - 1.0_dbl)
        r0t = 1.28_dbl*(At**(1.0_dbl/3.0_dbl)) - 0.76_dbl + 0.8_dbl*(At**(-1.0_dbl/3.0_dbl))
        rt = r0t*(1.0_dbl + sqrt(5.0_dbl/(4.0_dbl*pi))*B2*P2)
    
        Vc = ((Zp*Zt*e2)/R)*(1.0_dbl + ((rt*rt*B2*P2)/(R*R))*(sqrt(9.0_dbl/(20.0_dbl*pi)) + (3.0_dbl*B2*P2)/(7.0_dbl*pi)))
    
    end function Vc

end module calculations

program deformedNuclei

    use calculations
    implicit none
  
    real(kind=dbl), parameter :: Rmin = 2.0_dbl, Rmax = 6.0_dbl, dR = 0.1_dbl         ! Loop parameters
    real(kind=dbl) :: R, V, B2, deg, angle0, angle547, angle90, Zp, Zt, At            ! Radius, Coulomb Potential, quadrupole 
                                                                                      ! deformation parameter, degrees, angles, nucleons,
                                                                                      ! protons, neutrons
    integer :: iounit, iounit2, iounit3, iounit4
    real(kind=dbl) :: Vc_usual                                                        ! Coulomb potential with no deformation


    ! Angles for each data file
    angle0 = 0.0
    deg = 54.7_dbl
    angle547 = deg*pi/180.0_dbl
    angle90 = pi/2.0_dbl


    ! Properties of 117 Lanthanum, from National Nuclear Data Center, https://www.nndc.bnl.gov/ensdf/DSfromRefServlet?kn=1976LI30&searchType=references
    B2     = 0.31_dbl                                  ! https://www.sciencedirect.com/science/article/pii/S1876610218311603?ref=pdf_download&fr=RR-2&rr=99628a2688de51f0      
    Zp     = 57.0_dbl     
    Zt     = 60.0_dbl        
    At     = 117.0_dbl       

    ! Open a data file for each angle
    open(newunit=iounit, file='Vc_curve_angle0.dat', status='replace', action='write')
    open(newunit=iounit2, file='Vc_curve_angle547.dat', status='replace', action='write')
    open(newunit=iounit3, file='Vc_curve_angle90.dat', status='replace', action='write')

  
    ! Writes the Deformed Coulomb Potential as R changes, into each data file
    R = Rmin
    do while (R <= Rmax + 0.5_dbl*dR)
        V = Vc(R, B2, angle0, Zp, Zt, At)
        write(iounit,'(F6.2,",",ES20.12)') R, V
        R = R + dR
    end do
    
    R = Rmin
    do while (R <= Rmax + 0.5_dbl*dR)
        V = Vc(R, B2, angle547, Zp, Zt, At)
        write(iounit2,'(F6.2,",",ES20.12)') R, V
        R = R + dR
    end do
    
    R = Rmin
    do while (R <= Rmax + 0.5_dbl*dR)
        V = Vc(R, B2, angle90, Zp, Zt, At)
        write(iounit3,'(F6.2,",",ES20.12)') R, V
        R = R + dR
    end do

    close(iounit)
    close(iounit2)
    close(iounit3)

    ! Open a data file for spherical coulomb potential
    open(newunit=iounit4, file='Vc_curve_spherical.dat', status='replace', action='write')
    
    R = Rmin
    
    ! Writes the Spherical (no deformation) Coulomb Potential as R changes and no angle dependence, into a data file
        do while (R <= Rmax + 0.5_dbl*dR)
        Vc_usual = (Zp * Zt * e2) / R
        write(iounit4,'(F6.2,",",ES20.12)') R, Vc_usual
        R = R + dR
    end do
    
    close(iounit4)
    
end program deformedNuclei
    
    
    