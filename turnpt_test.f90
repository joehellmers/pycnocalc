program run_turn_pt

        use rate_calc3
    implicit none


    real(kind=dbl) :: Rstep, E0, mu, WKB
    integer :: Rmax, turn1, turn2

    ! -------------------------
    ! SET INPUT VALUES 
    ! -------------------------
    Rstep = 0.1d0               ! radius step (needs to be consistent with data)
    Rmax  = 300                 ! (maximum radius - minimum radius) / Rstep
    E0    = 25.0d0              ! input energy
    mu    = 1.0d0               ! 

    ! -------------------------
    ! CALL SUBROUTINE that calculates turning points
    ! -------------------------
    call turn_pt3(Rstep, Rmax, E0, mu, turn1, turn2, WKB)

    ! -------------------------
    ! PRINT RESULTS (turning points & WKB, -1 turning point value means no turning point)
    ! -------------------------
    print *, 'turn1 = ', turn1
    print *, 'turn2 = ', turn2
    print *, 'WKB   = ', WKB

end program run_turn_pt
