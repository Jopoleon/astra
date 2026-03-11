module auxiliary

implicit none

contains

!---------------------------------------------------------------------
    double precision function dfj(f_old, f_new)

    double precision, intent(in) :: f_old, f_new
    double precision :: f_mid

    f_mid = 0.5*(abs(f_old) + abs(f_new))
    dfj = abs(f_old - f_new)
    if (f_mid >= 1.e-6) then
        dfj = dfj/f_mid
    endif
 
    return
    end function dfj

!---------------------------------------------------------------------
    integer function IFSTEP

! Input
!       LEQ(1)  LEQ(2)  LEQ(3)  LEQ(4)  LEQ(5)  LEQ(6-9)
!  NE  TE  TI  CU equil    free
!       LEQ(10) LEQ(11) LEQ(12), etc
!        F0  F1  F2
!       LEQ(j) - usage of j-th equation in the model
!        =-1 no equation (default)
!        = 0 type: AS
!        = 1 type: EQ
!        = 2 type: FU
!        = 3 type: heat conductivity flux + heat convection
! Returned value
!  1 - time step made, normal exit
!  0 - time step will be repeated
!---------------------------------------------------------------------

    use const_inc, only: TAUINC, DELVAR, TAU, TAUPRP, TAUMIN, TAUMAX, &
        DPOUT, NA, NB1, LEQ, NSTEPS, ROC, ROCO, FTN, FTO
    use status_inc, only: NE, NI, TE, TI, FP, VR, UPAR, &
        NEO, NIO, TEO, TIO, FPO, VRO, UPARO, &
        F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, &
        F0O, F1O, F2O, F3O, F4O, F5O, F6O, F7O, F8O, F9O
    use debugger, only: markloc
    use metrics, only: CUOFP

    integer :: j
    double precision :: CTAU, TAUO, YY, TAUN

    call markloc('IFSTEP')

    CTAU = 1./TAUINC

    do j=1, NA
        if (LEQ(1)  > 0) CTAU = MAX(CTAU, ABS(NEO(j)/NE(j) - 1.)/DELVAR)
        if (LEQ(2)  > 0) CTAU = MAX(CTAU, ABS(TEO(j)/TE(j) - 1.)/DELVAR)
        if (LEQ(3)  > 0) CTAU = MAX(CTAU, ABS(TIO(j)/TI(j) - 1.)/DELVAR)
        if (LEQ(10) > 0) CTAU = MAX(CTAU, dfj(F0O(j), F0(j))/DELVAR)
        if (LEQ(11) > 0) CTAU = MAX(CTAU, dfj(F1O(j), F1(j))/DELVAR)
        if (LEQ(12) > 0) CTAU = MAX(CTAU, dfj(F2O(j), F2(j))/DELVAR)
        if (LEQ(13) > 0) CTAU = MAX(CTAU, dfj(F3O(j), F3(j))/DELVAR)
        if (LEQ(14) > 0) CTAU = MAX(CTAU, dfj(F4O(j), F4(j))/DELVAR)
        if (LEQ(15) > 0) CTAU = MAX(CTAU, dfj(F5O(j), F5(j))/DELVAR)
        if (LEQ(16) > 0) CTAU = MAX(CTAU, dfj(F6O(j), F6(j))/DELVAR)
        if (LEQ(17) > 0) CTAU = MAX(CTAU, dfj(F7O(j), F7(j))/DELVAR)
        if (LEQ(18) > 0) CTAU = MAX(CTAU, dfj(F8O(j), F8(j))/DELVAR)
        if (LEQ(19) > 0) CTAU = MAX(CTAU, dfj(F9O(j), F9(j))/DELVAR)
    enddo

    TAUO   = TAU
    TAUPRP = TAUO
    TAUN   = TAU
    TAUN   = MIN(TAUMAX, TAUN/CTAU, DPOUT)
    TAU = MAX(TAUMIN, TAUN)   ! due to DELVAR & TAUINC

    if (TAU >= TAUO) then
        IFSTEP = 1
        NSTEPS = NSTEPS + 1
        return  ! -> proceed to the next time step
    endif

! Reset all functions F(t) to old values

    do j=1, NB1
        NE(j) = NEO(j)
        NI(j) = NIO(j)
        TE(j) = TEO(j)
        TI(j) = TIO(j)
        FP(j) = FPO(j)
        VR(j) = VRO(j)
        UPAR(j) = UPARO(j)
        F0(j) = F0O(j)
        F1(j) = F1O(j)
        F2(j) = F2O(j)
        F3(j) = F3O(j)
        F4(j) = F4O(j)
        F5(j) = F5O(j)
        F6(j) = F6O(j)
        F7(j) = F7O(j)
        F8(j) = F8O(j)
        F9(j) = F9O(j)
    enddo

! Reset ROC
    ROC = ROCO
    FTO = FTN

    call CUOFP  ! Restore CU & MU, to be checked after VR is done EFable
    IFSTEP = 0

    return
    end function IFSTEP

!---------------------------------------------------------------------
    integer function IFTREQ(YACC)
!---------------------------------------------------------------------
!    Control parameter:
! NTREQ == 1 iterations switched off (timing VR(t) downshifted)
!       IFTREQ == 2 is returned
! NTREQ >= 2 iterations switched on  (VR(t) properly included)
!
!    Input:
! YACC - requested accuracy
! 
! ITREQ = 0 is set in BLOCKDATA, reset in stepup, incremented in METRIC
! ITREQ == 1 save metric parameters (G11, G22, VR, ROC) and exit
!            IFTREQ == 2 is returned
!
! ITREQ >= 1 & ITREQ <= NITREQ - Mark time evolution section
!
!---------------------------------------------------------------------
!    Returned value:
! IFTREQ == 0 - No convergence
! IFTREQ == 1 - Convergence condition is achieved
! IFTREQ == 2 - No convergence and maximum iteration number is achieved
!---------------------------------------------------------------------

    use parameter_inc, only: NRD
    use status_inc, only: G11, G22, VR, FP
    use const_inc, only: IPART, ITREQ, ROC, NA1
    use debugger, only: markloc

    integer, parameter :: ITREQMIN=1, ITREQMAX=200
    double precision, intent(in) :: YACC

    integer :: j, NTREQ
    double precision :: Y1, Y2, YV, YI, YR, YER(ITREQMAX)
    double precision, dimension(NRD) :: YWA, YWB, YWC, YWD

    save YER, YR, YWA, YWB, YWC, YWD

    call markloc('IFTREQ')

    NTREQ = 1
! Disable NTREQ setting for the initial phase
    if (IPART == 1) NTREQ = ITREQMAX
    IFTREQ = 2
    if (NTREQ == 1) return ! Iterations are off, no check

    if (NTREQ > ITREQMAX) then  ! 
         write(*, '(2A,I5)') " >>> Warning >>> Max number of transport/equilibrium", &
             " iterations is reduced to ", ITREQMAX
         NTREQ = ITREQMAX
    endif
    IFTREQ = 0

    if (ITREQ /= 0) then ! From 2nd iteration

        if (ITREQ >= ITREQMAX) then
            write(*,'(2A,2I5)') " >>> Warning >>> Max number of initialization ", &
                 " iterations is overwhelmed, stopping ", ITREQ, ITREQMAX
            stop
        endif

! Convergence check:  enabled if NTREQ > 1
        YR = abs(YR/ROC - 1.d0)
        Y1 = 0.d0
        Y2 = 0.d0
        YV = 0.d0
        YI = 0.d0
        do j=1, NA1  ! Analize NA1new vs NA1old
            Y1 = max(Y1, abs(YWA(j)/G11(j) - 1.d0))
            Y2 = max(Y2, abs(YWB(j)/G22(j) - 1.d0))
            YV = max(YV, abs(YWC(j)/VR(j)  - 1.d0))
            YI = max(YI, abs(YWD(j)/FP(j)  - 1.d0))
        enddo
        if (Y1 > max(Y2, YV, YR, YI)) then
            YER(ITREQ) = Y1
            j = 1
        endif
        if (Y2 > max(Y1, YV, YR, YI)) then
            YER(ITREQ) = Y2
            j = 2
        endif
        if (YV > max(Y1, Y2, YR, YI)) then
            YER(ITREQ) = YV
            j = 3
        endif
        if (YR > max(Y1, Y2, YV, YI)) then
            YER(ITREQ) = YR
            j = 4
        endif
        if (YI > max(Y1, Y2, YV, YR)) then
            YER(ITREQ) = YI
            j = 5
        endif

! Check convergence:
! Note! the tolerance, YACC(alias ATREQ), should be >> than ACEQLB (1.d-6)
        IFTREQ = 0
        if (ITREQ == NTREQ    ) IFTREQ = 2 ! Go to next t-step anyway
        if (YACC >= YER(ITREQ)) IFTREQ = 1 ! Converged -> next t-step
        if (ITREQ < ITREQMIN  ) IFTREQ = 0 ! force another iteration if below minimum
        if (IFTREQ > 0) return ! converged

    endif

    ITREQ = ITREQ + 1 

    YR = ROC
    do j=1, NA1   ! Store some metric data
        YWA(j) = G11(j)
        YWB(j) = G22(j)
        YWC(j) = VR(j)
        YWD(j) = FP(j)
    enddo

    return
    end function IFTREQ

!---------------------------------------------------------------------
    subroutine OLDNEW

    use status_inc, only: NE, TE, TI, FP, VR, UPAR, UPS0, UPS1, UPS2, &
        NEO, TEO, TIO, FPO, VRO, UPARO, UPS0O, UPS1O, UPS2O, &
        F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, &
        F0O, F1O, F2O, F3O, F4O, F5O, F6O, F7O, F8O, F9O
    use const_inc, only: NB1, BTN, FTO, FTN, ROCO, ROC, BTOR
    use debugger, only: markloc

    integer :: j

    call markloc('OLDNEW')

    do j=1, NB1
        NEO(j) = NE(j)
        TEO(j) = TE(j)
        TIO(j) = TI(j)
        FPO(j) = FP(j)
        VRO(j) = VR(j)
        UPARO(j) = UPAR(j)
        UPS0O(j) = UPS0(j)
        UPS1O(j) = UPS1(j)
        UPS2O(j) = UPS2(j)
        F0O(j) = F0(j)
        F1O(j) = F1(j)
        F2O(j) = F2(j)
        F3O(j) = F3(j)
        F4O(j) = F4(j)
        F5O(j) = F5(j)
        F6O(j) = F6(j)
        F7O(j) = F7(j)
        F8O(j) = F8(j)
        F9O(j) = F9(j)
    enddo

    BTN  = BTOR
    FTN  = FTO
    ROCO = ROC

    return
    end subroutine OLDNEW

!---------------------------------------------------------------------
    double precision function LINEAV

! LINEAV [10#19/m#3]: Horizontal chord average density (r) [m]
! Integral {0, r} ( NE ) dl / a

    use const_inc, only: NA, ABC, NA1
    use status_inc, only: AMETR, NE

    integer :: j

    LINEAV = 2.*AMETR(1)*NE(1)
    do j=2, NA
        LINEAV = LINEAV + (AMETR(j) - AMETR(j-1))*(NE(j) + NE(j-1))
    enddo
    LINEAV = 0.5*(LINEAV + (ABC - AMETR(NA))*(NE(NA1) + NE(NA)))/ABC

    return
    end function lineav

end module auxiliary
