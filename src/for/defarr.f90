!---------------------------------------------------------------------
subroutine DEFARR
!---------------------------------------------------------------------
! 1) Check positiveness of Z_eff, n_e, n_i, T_e, T_i.  Stop if negative.
! 2) Extend definition of all standard arrays beyond NA1.
!Efable
! 3) Define upsilons for momentum transport equation
!---------------------------------------------------------------------

use outcmn_inc, only: exp_file, NSBR, DTNAME
use const_inc, only: GP2, RTOR, BTOR, HRO, ROC, ABC, NA1, NB1, NAB, &
    NSDELOUT, TIME, TAU, TSTART, WTE, WTI, WNE
use status_inc
use debugger, only: markloc, astra_stop, flightsim

implicit none

integer :: j, js
double precision :: YV, YF, YMU, YN, YNE, YNI, YTE, YTI, YZF
character(len=2) :: str
character(len=132) :: err_msg

call markloc('DEFARR')

YZF = 0.
YNE = 0.
YNI = 0.
YTE = 0.
YTI = 0.
YV = GP2*BTOR*HRO
do j=1, NA1
    SQEPS(j) = SQRT(AMETR(j)/(RTOR + SHIF(j)))
! VP = c<E_||/h>/B_p in Hinton & Hazeltine Rev.Mod.Phys.48(1976) p.297
    if (MU(j) <= 0.d0) then
        write(*,'(/A,F10.6,A)')'>>> ERROR >>>  Time =', TIME, ' sec'
        write(*,'(A/A,F10.6)') &
            '               The rotational transform is less or equal zero ', &
            '               at rho_N =', RHO(j)/ROC
        if (TIME <= TSTART + TAU/2.) write(*, '(2A,1H"/)') &
            '               Check if it is defined in the data file "', &
            TRIM(exp_file)
        if (flightsim >= 0) call err_catch_a
        call astra_stop
    endif
    VP(j) = ULON(j)/(YV*j*MU(j))
    YZF = max(YZF, ZEF(j))
    YNE = max(YNE, NE(j))
    YNI = max(YNI, NI(j))
    YTE = max(YTE, TE(j))
    YTI = max(YTI, TI(j))
enddo

str = '3M'
if (min(YZF, YNE, YNI, YTE, YTI) <= 0.) then

    if (YZF <= 0.) then
        str = "ZF"
    elseif (YNE <= 0.) then
        str = "NE"
    elseif (YNI <= 0.) then
        str = "NI"
    elseif (YTE <= 0.) then
        str = "TE"
    elseif (YTI <= 0.) then
        str = "TI"
    endif

    write(*, '(/A, F6.3, A)') '>>> ERROR >>>  Time =',TIME,'   DEFVAR'

    SELECT CASE(str)
    CASE('3M')
        err_msg = '               The 3M equilibrium solver does not converge'
    CASE('ZF')
        err_msg = '               The variable "ZEF" is less or equal zero'
    CASE DEFAULT
        err_msg = '               The variable "' // str // '" is less or equal zero'
    END SELECT
    call astra_stop(err_msg)

endif

if (NA1 /= NB1) then

! Determine arrays beyond ABC
    YF = GP2*BTOR
    YMU = (MU(NA1) - MV(NA1))*ROC**2
    js = 1
    do j=1, NSBR
        if (DTNAME(NSDELOUT+4*j) == 'NEUTAB') js = 3
    enddo

    do j=NA1+1, NB1
        FPO(j)   = FPO(NA1)
        MU(j)    = MU(NA1)
        ULON(j)  = ULON(NA1)
        TTRQ(j)  = TTRQ(NA1)
        CU(j)    = CU(NA1)
        VP(j)    = VP(NA1)
        GN(j)    = GN(NA1)
        SQEPS(j) = SQEPS(NA1)
        PELH(j)  = PELH(NA1)
        PETOT(j) = PETOT(NA1)
        PELON(j) = PELON(NA1)
        CULH(j)  = CULH(NA1)
        PITOT(j) = PITOT(NA1)
        SNTOT(j) = SNTOT(NA1)
        ZEF(j)   = ZEF(NA1)
        PEECR(j) = PEECR(NA1)
        PEFW(j)  = PEFW(NA1)
        PEICR(j) = PEICR(NA1)
        PIICR(j) = PIICR(NA1)
        PIFW(j)  = PIFW(NA1)
        PEPER(j) = PEPER(NA1)

        QE(j)    = QE(NA1)
        QI(j)    = QI(NA1)
        QN(j)    = QN(NA1)
        QF0(j)   = QF0(NA1)
        QF1(j)   = QF1(NA1)
        QF2(j)   = QF2(NA1)
        QF3(j)   = QF3(NA1)
        QF4(j)   = QF4(NA1)
        QF5(j)   = QF5(NA1)
        QF6(j)   = QF6(NA1)
        QF7(j)   = QF7(NA1)
        QF8(j)   = QF8(NA1)
        QF9(j)   = QF9(NA1)

        CUECR(j) = CUECR(NA1)
        CUFW(j)  = CUFW(NA1)
        CUICR(j) = CUICR(NA1)
        CUTOR(j) = CUTOR(NA1)
        UPL(j)   = UPL(NA1)
        ER(j)    = ER(NA1)
        VRS(j)   = VRS(NA1)
        NI(j)    = NI(NA1)
        NIO(j)   = NIO(NA1)
        ZMAIN(j) = ZMAIN(NA1)
        AMAIN(j) = AMAIN(NA1)

        VRO(j)   = VRO(NA1)
        VR(j)    = VR(NA1)
        SHIF(j)  = SHIF(NA1)
        SHIV(j)  = SHIV(NA1)
        ELON(j)  = ELON(NA1)
        TRIA(j)  = TRIA(NA1)
        G11(j)   = G11(NA1)
        G22(j)   = G22(NA1)
        G33(j)   = G33(NA1)

        G41(j)   = G41(NA1)
        G42(j)   = G42(NA1)
        G43(j)   = G43(NA1)
        G44(j)   = G44(NA1)
        G45(j)   = G45(NA1)

        VPFP(j)  = VPFP(NA1)
        IPOL(j)  = IPOL(NA1)
        CUBS(j)  = CUBS(NA1)
        VOLUM(j) = VOLUM(NA1)
        CV(j)    = CV(NA1)
        DRODA(j) = DRODA(NA1)
        PRES(j)  = PRES(NA1)
        PDE(j)   = PDE(NA1)
        PDI(j)   = PDI(NA1)

        SDN(j)   = SDN(NA1)
        SD0(j)   = SD0(NA1)
        SD1(j)   = SD1(NA1)
        SD2(j)   = SD2(NA1)
        SD3(j)   = SD3(NA1)
        SD4(j)   = SD4(NA1)
        SD5(j)   = SD5(NA1)
        SD6(j)   = SD6(NA1)
        SD7(j)   = SD7(NA1)
        SD8(j)   = SD8(NA1)
        SD9(j)   = SD9(NA1)

        UPS0(j)  = UPS0(NA1)
        UPS1(j)  = UPS1(NA1)
        UPS0O(j) = UPS0O(NA1)
        UPS1O(j) = UPS1O(NA1)

        if (js == 1) then
            NN(j) = NN(NA1)
            TN(j) = TN(NA1)
        endif
        NIZ1(j)  = NIZ1(NA1)
        NIZ2(j)  = NIZ2(NA1)
        NIZ3(j)  = NIZ3(NA1)
        NALF(j)  = NALF(NA1)
        PRAD(j)  = PRAD(NA1)
        VTOR(j)  = VTOR(NA1)
        NHYDR(j) = NHYDR(NA1)
        NDEUT(j) = NDEUT(NA1)
        NTRIT(j) = NTRIT(NA1)
        NHE3(j)  = NHE3(NA1)
        ZIM1(j)  = ZIM1(NA1)
        ZIM2(j)  = ZIM2(NA1)
        ZIM3(j)  = ZIM3(NA1)
        PBOL1(j) = PBOL1(NA1)
        PBOL2(j) = PBOL2(NA1)
        PBOL3(j) = PBOL3(NA1)
        PSXR1(j) = PSXR1(NA1)
        PSXR2(j) = PSXR2(NA1)
        PSXR3(j) = PSXR3(NA1)
        ZEF1(j)  = ZEF1(NA1)
        ZEF2(j)  = ZEF2(NA1)
        ZEF3(j)  = ZEF3(NA1)
        DIMP1(j) = DIMP1(NA1)
        DIMP2(j) = DIMP2(NA1)
        DIMP3(j) = DIMP3(NA1)
        VIMP1(j) = VIMP1(NA1)
        VIMP2(j) = VIMP2(NA1)
        VIMP3(j) = VIMP3(NA1)
        VPOL(j)  = VPOL(NA1)
        NMAIN(j) = NMAIN(NA1)

        if (js == 3 .and. j > NAB) then
            NN(j) = NN(NAB)
            TN(j) = TN(NAB)
        endif
        YN = exp((ABC - AMETR(j))/WNE)
        NE(j)    = NE(NA1)*YN
        NI(j)    = NI(NA1)*YN
        NALF(j)  = NALF(NA1)*YN
        NHE3(j)  = NHE3(NA1)*YN
        NHYDR(j) = NHYDR(NA1)*YN
        NDEUT(j) = NDEUT(NA1)*YN
        NTRIT(j) = NTRIT(NA1)*YN
        NMAIN(j) = NMAIN(NA1)*YN
        TE(j)    = TE(NA1)*exp((ABC - AMETR(j))/WTE)
        TI(j)    = TI(NA1)*exp((ABC - AMETR(j))/WTI)
! Efable: added UPAR decay in SOL with TI length scale
        UPAR(j) = UPAR(NA1)*exp((ABC-AMETR(j))/WTI)
! Efable: should be changed when a SOL equilibrium is implemented
        MV(j) = MV(NA1)
        FV(j) = FV(j-1) + YF*RHO(j)*MV(j)*(RHO(j) - RHO(j-1))
        MU(j) = MV(j) + YMU/RHO(j)**2
        FP(j) = FP(j-1) + YF*RHO(j)*MU(j)*(RHO(j) - RHO(j-1))
        CU(j) = 0.
        CUBS(j) = 0.
        CV(j) = 0.
        CD(j) = 0.
        CUTOR(j) = 0.
    enddo
endif

return
end subroutine DEFARR

!---------------------------------------------------------------------
double precision function GETXVAR(VAR, TIME1)
!---------------------------------------------------------------------
! Similar to INTVAR but returns a single quantity at arbitrary time
!---------------------------------------------------------------------
! The time evolution of the input data is taken from
! VARDAT(1,NTVAR) - time   |
! VARDAT(2,NTVAR) - value   |  for variables
! VARDAT(3,NTVAR) - error (not used) |
! NTVAR total number of the time slices for all variables in a data file
!---------------------------------------------------------------------
! Input:
! IVAR
! INDVAR
! IFDFVX
! VARDAT
! TIME1
!
! jj  - ordinal number of the record in VARDAT(1-2-3,jj)
! N=INDVAR(jj) - ordinal number of the quantity in DEVAR(N) & DEVARX(N)
! IFDFVX(N)  - type of variable
! IFDFVX = 0 - determined by the data file (independent on time),
!  = 1 - determined by the data file (dependent on time),
!  = 2 - determined by the MODEL,
!  = 3 - Keyboard
!  = 4 - any change of the variable is forbidden 
!    eg. (AB, RTOR, ELONM, TRICH or set interactively)
!---------------------------------------------------------------------

use parameter_inc, only: NTVAR
use outcmn_inc, only: PRNAME, IFDFVX
use expdat, only: IVAR, INDVAR, VARDAT
use debugger, only: markloc

implicit none

double precision, intent(in) :: TIME1
character(len=*), intent(in) :: VAR

integer :: jj, N1, N2
double precision :: YDT, YDTR, YDTL, XVAR

getxvar = 0.

call markloc('GETVARX')

N1 = 0
N2 = 0

do jj=1, IVAR
    N2 = N1
    N1 = INDVAR(jj)
    if (IFDFVX(N1) < 0) then
        write(*, *) "The quantity ", VAR, " is not defined in exp file"
        stop
    endif
    if (TRIM(PRNAME(N1)) == TRIM(VAR) ) EXIT
enddo

if (IFDFVX(N1) == 0 .or. N1 /= N2) XVAR = VARDAT(2, jj)
if (N1 == N2) then
    if (TIME1 - VARDAT(1, jj-1) >= 0.) then
        if (TIME1 - VARDAT(1, jj) < 0.) then
            YDT  =  VARDAT(1, jj) - VARDAT(1, jj-1)
            YDTR = (VARDAT(1, jj) - TIME1)/YDT
            YDTL = (TIME1 - VARDAT(1, jj-1))/YDT
            XVAR = VARDAT(2, jj)*YDTL + VARDAT(2, jj-1)*YDTR      
        else if (TIME1 - VARDAT(1, jj) == 0.) then
            XVAR = VARDAT(2, jj)
        endif
    endif
endif
if (IFDFVX(N1) <= 1) GETXVAR = XVAR

return
end function GETXVAR

!---------------------------------------------------------------------
double precision function ARRNA1(ARR, H)
!---------------------------------------------------------------------
! Compute edge value of the array ARR(NA1) assuming that d^3(ARR)/dr^3=0
! EFable shall be substituted by EXTRAP_EF with quadratic, or cubic, spline
!---------------------------------------------------------------------
! It should be
! ARRNA1 = ARR(2)+(ARR(2)-ARR(1))*H ! call ARRNA1(UPL(NA-1),HRO/HRO)
! ARRNA1 = ARR(3)+(ARR(3)-ARR(2))*H
! ARRNA1 = ARR(3)*(1+(1.5+0.5*H)*H)+ARR(1)*H*0.5*(1+H)-
!     -  ARR(2)*H*(2+H)
! ARRNA1 = ARR(3)*(1+H)-ARR(2)*H
! ARRNA1 = ARR(3)

double precision, intent(in) :: ARR, H

ARRNA1 = ARR !+H*(2.*ARR(1)+ARR(-2)-3.*ARR(-1))

return
end function ARRNA1

!---------------------------------------------------------------------
integer function IFSTEP(IFCONV, updwno)
! IFCONV dummy parameter
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

use const_inc, only: TIME, TAUINC, DELVAR, TAU, TAUPRP, TAUMIN, TAUMAX, &
    DTOUT, DPOUT, NA, NB1, LEQ, NSTEPS, UPDWN, ROC, ROCO, &
    BTN, IPLN, FTN, FTO
use status_inc, only: NEO, NIO, TEO, TIO, FJO, FPO, VRO, UPARO, &
    NE, NI, TE, TI, FJ, FP, VR, UPAR
use debugger, only: markloc, flightsim

implicit none

integer, intent(in) :: IFCONV
double precision, intent(in) :: updwno

integer :: j, jj
double precision :: CTAU, TAUO, YY, TAUN

call markloc('IFSTEP')

CTAU = 1./TAUINC

do j = 1, NA
    if (LEQ(1) > 0) CTAU = MAX(CTAU, ABS(NEO(j)/NE(j) - 1.)/DELVAR)
    if (LEQ(2) > 0) CTAU = MAX(CTAU, ABS(TEO(j)/TE(j) - 1.)/DELVAR)
    if (LEQ(3) > 0) CTAU = MAX(CTAU, ABS(TIO(j)/TI(j) - 1.)/DELVAR)
    if (flightsim >= 1 .and. j == 1) CTAU = MAX(CTAU, ABS(Updwno/UPDWN - 1.)/DELVAR)
    do jj=0, 9
        if (LEQ(jj+10) > 0) then
            YY = 0.5*(abs(FJO(j, jj)) + abs(FJ(j, jj)))
            if (YY < 1.E-6) then  ! Allow zero FJ
                YY = abs(FJO(j, jj) - FJ(j, jj))
            else
                YY = abs(FJO(j, jj) - FJ(j, jj))/YY
            endif
            CTAU = MAX(CTAU, YY/DELVAR)
        endif
    enddo
enddo

TAUO   = TAU
TAUPRP = TAUO
TAUN   = TAU
TAUN   = MIN(TAUMAX, TAUN/CTAU, DTOUT, DPOUT)

TAUN = MAX(TAUMIN, TAUN)   ! due to DELVAR & TAUINC
TAU = TAUN 
if (flightsim >= 1) tauo = 1.d-6*nint(tauo*1.d6)
if (flightsim >= 1) tau = 1.d-6*nint(tau*1.d6)

if (TAU >= TAUO) then
    IFSTEP = 1
    NSTEPS = NSTEPS+1
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
    do jj=0, 9
        FJ(j, jj) = FJO(j, jj)
    enddo
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
! NITREQ == 1 iterations switched off (timing VR(t) downshifted)
!       IFTREQ == 2 is returned
! NITREQ >= 2 iterations switched on  (VR(t) properly included)
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
use const_inc, only: NITREQ, IPART, ITREQ, ROC, PSIEXT, PSPLEX, NA1, NB1
use debugger, only: markloc

implicit none

integer, parameter :: ITREQMIN=1, ITREQMAX=200
double precision, intent(in) :: YACC

integer :: j, NTREQ

double precision :: Y1, Y2, YV, YI, YR, YER(ITREQMAX), YPSE, YPSP
double precision, dimension(NRD) :: YWA, YWB, YWC, YWD

save YER, YR, YWA, YWB, YWC, YWD
save YPSE, YPSP

call markloc('IFTREQ')

NTREQ = nint(NITREQ)
! Disable NITREQ setting for the initial phase
if (IPART == 1) NTREQ = ITREQMAX
IFTREQ = 2
if (NTREQ == 1) return ! Iterations are off, no check

if (NTREQ > ITREQMAX) then  ! 
     write(*,'(2A,I5)') " >>> Warning >>> Max number of transport/equilibrium", &
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
    YR   = abs(YR/ROC      - 1.d0)
    YPSE = abs(YPSE/PSIEXT - 1.d0)
    YPSP = abs(YPSP/PSPLEX - 1.d0)
    Y1 = 0.d0
    Y2 = 0.d0
    YV = 0.d0
    YI = 0.d0
    do j=1, NA1  ! Analize NA1new vs NA1old
        Y1 = max(Y1, abs(YWA(j)/G11(j) - 1.d0))
        Y2 = max(Y2, abs(YWB(j)/G22(j) - 1.d0))
        YV = max(YV, abs(YWC(j)/VR(j) - 1.d0))
        YI = max(YI, abs(YWD(j)/FP(j) - 1.d0))
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

YR   = ROC
YPSE = PSIEXT
YPSP = PSPLEX
do j=1, NB1   ! Store some metric data
    YWA(j) = G11(j)
    YWB(j) = G22(j)
    YWC(j) = VR(j)
    YWD(j) = FP(j)
enddo

return
end function IFTREQ

!---------------------------------------------------------------------
subroutine OLDNEW

use status_inc, only: NE, TE, TI, FP, VR, FJO, UPAR, &
    UPS0, UPS1, UPS2, NEO, TEO, TIO, FPO, VRO, FJ, &
    UPARO, UPS0O, UPS1O, UPS2O
use const_inc, only: NB1, BTN, TIME, TAU, IPLN, FTO, FTN, ROCO, ROC, BTOR, BTN
use debugger, only: markloc

implicit none

integer :: j, jj

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
    do jj=0, 9
        FJO(j, jj) = FJ(j, jj)
    enddo
enddo

BTN = BTOR
FTN = FTO
ROCO = ROC

return
end subroutine OLDNEW

!---------------------------------------------------------------------
subroutine SMOOTH(ALFA, NO, FO, XO, N, FN, XN)
!---------------------------------------------------------------------
!  Subroutine minimizes the value of functional
!  INTEGRAL(alfa*P(x)*(dU/dx)**2+(U-F)**2)*dx, 
!  where FO(NO) is a function, given on the grid XOld(NOld)
! P(x) is equal to unit now
! ALFA=alfa<<0.01*XO(NO)**2 is regularizator
! NO - number of old grid points
! N=<NRD - number of new grid points
! 0<=XO(NO) - old grid |     both grids are arbitrary
! 0<=XN(N)  - new grid |     but XO(NO)=XN(N)
! FO(NO) - origin function, given on the grid XO(NO)
! FN(N) - smoothed function on grid XN(N)
!  The result is function FN(XN), given on the new grid
! with additional conditions:
! dFN/dx(x=0)=0 - cylindrical case and
! FN(XN(N))=FO(XO(NO))
!---------------------------------------------------------------------

use parameter_inc, only: NRD
use debugger, only: astra_stop

implicit none

integer, intent(in) :: NO, N
double precision, intent(in) :: ALFA, XO(*), FO(*), XN(*)
double precision, intent(out) :: FN(*)

integer :: I, J
double precision :: YF, YX, YP, YQ, YD, FJ, P(NRD)
character(len=132) :: err_msg

if (N > NRD .or. NO <= 0) then
    write(err_msg, '(A, 1X, i, 1X, i)') ' >>> SMOOTH: array is out of limits', NRD, NO
    call err_catch_a
    call astra_stop(err_msg)
endif
if (NO == 1) then
    do j=1, N
        FN(j) = FO(1)
    enddo
    return
endif
if (NO == 2) then
    do j=1, N
        FN(j) = (FO(2)*(XN(j) - XO(1)) - FO(1)*(XN(j) - XO(2)))/(XO(2) - XO(1))
    enddo
    return
endif
if (N < 2) then
    call err_catch_a
    call astra_stop(' >>> SMOOTH: no output grid is provided')
endif
if (abs(XO(NO) - XN(N)) > XN(N)/N) then
    write(*, *) '>>> SMOOTH: grids are not aligned'
    write(*, '(1A23, I4, F8.4)') '     Old grid size/edge', NO, XO(NO)
    write(*, '(1A23, I4, F8.4)') '     New grid size/edge', N, XN(N)
    call err_catch_a
    call astra_stop
endif
do j=2, N
    YP = XN(j) - XN(j-1)
    if (YP <= 0.d0) then
        write(*, *)'>>> SMOOTH: new grid is not increasing monotonically'
        write(*, '(A, I4, A, F8.4)')'Node ', j-1, '   Value', XN(j-1)
        write(*, '(A, I4, A, F8.4)')'Node ', j  , '   Value', XN(j)
        call err_catch_a
        call astra_stop
    endif
    P(j) = ALFA/YP/XO(NO)**2
enddo
P(1)  = 0.
FN(1) = 0.
I = 1
YF = (FO(2) - FO(1))/(XO(2) - XO(1))
YX = 2./(XN(2) + XN(1))
YP = 0.
YQ = 0.
do j=1, N-1
    if (XO(I) <= XN(j)) then
        do
            I = I + 1
            if (I > NO) I = NO
            if (I == NO .or. XO(I) >= XN(j)) EXIT
        enddo
        YF = (FO(I) - FO(I-1))/(XO(I) - XO(I-1))
    endif
    FJ = FO(I) + YF*(XN(j) - XO(I))
    YD = 1. + YX*(YP + P(j+1))
    P(j) = YX*P(j+1)/YD
    FN(j) = (FJ + YX*YQ)/YD
    if (j /= N-1) then
        YX = 2./(XN(j+2) - XN(j))
        YP = (1. - P(j))*P(j+1)
        YQ = FN(j)*P(j+1)
    endif
enddo

FN(N) = FO(NO)
do j=N-1, 1, -1
    FN(j) = P(j)*FN(j+1) + FN(j)
enddo

return
end subroutine SMOOTH

!---------------------------------------------------------------------
subroutine SMAP(ALFA, NO, XO, N, XN, F)
! Similar to SMOOTH but the same array, F, is used for input and output

use parameter_inc, only: NRD

implicit none

integer, intent(in) :: NO, N
double precision, intent(in) :: ALFA, XO(*), XN(*)
double precision, intent(inout) :: F(*)

integer :: J
double precision :: P(NRD)

call SMOOTH(ALFA, NO, F, XO, N, P, XN)
do j=1, N
    F(j) = P(j)
enddo

return
end subroutine SMAP

!---------------------------------------------------------------------
subroutine TRANSF(NO, FO, XO, N, FN, XN)
!---------------------------------------------------------------------
!  The subroutine transmits a function FO(1:NO)
!  from an arbitrary grid XO(1:NO) to a functon F(1:N)
!  on an another arbitrary grid XN(1:N) by the method
!  of quadratic interpolation.
! Input: NO, FO(1:NO), XO(1:NO), N, XN(1:N)
! Output: FN(1:N)
!---------------------------------------------------------------------
!       Note: [XO(NO)-XO(1)]*[XN(N)-XN(1)] must be > 0
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: NO, N
double precision, intent(in) , dimension(*) :: XO, FO, XN
double precision, intent(out), dimension(*) :: FN

integer :: I, J
double precision :: YF1, YF2, YF3

I = 1
YF1 = FO(1)/((XO(1) - XO(2))*(XO(1) - XO(3)))
YF2 = FO(2)/((XO(2) - XO(1))*(XO(2) - XO(3)))
YF3 = FO(3)/((XO(3) - XO(1))*(XO(3) - XO(2)))
do j=1, N
    if(2.*XN(j) > XO(I+1) + XO(I+2)) then
        do
            I = I + 1
            I = min(I, NO - 2)
            if (2.*XN(j) <= XO(I+1) + XO(I+2) .or. I >= NO-2) EXIT
        enddo
        YF1 = FO(I  )/((XO(I  ) - XO(I+1))*(XO(I  ) - XO(I+2)))
        YF2 = FO(I+1)/((XO(I+1) - XO(I  ))*(XO(I+1) - XO(I+2)))
        YF3 = FO(I+2)/((XO(I+2) - XO(I  ))*(XO(I+2) - XO(I+1)))
    endif
    FN(j) = YF1*(XN(j) - XO(I+1))*(XN(j) - XO(I+2)) + &
            YF2*(XN(j) - XO(I  ))*(XN(j) - XO(I+2)) + &
            YF3*(XN(j) - XO(I  ))*(XN(j) - XO(I+1))
enddo

return
end subroutine TRANSF

!---------------------------------------------------------------------
subroutine QMAP(NO, XO, N, XN, F)
! Similar to TRANSF but uses the same array, F, for input and output.

use parameter_inc, only: NRD

implicit none

integer, intent(in) :: NO, N
double precision, intent(in), dimension(*) :: XO, XN
double precision, intent(inout), dimension(*) :: F

integer :: j
double precision :: FN(NRD)

call TRANSF(NO, F, XO, N, FN, XN)

do j=1, N
    F(j) = FN(j)
enddo

return
end subroutine QMAP

!---------------------------------------------------------------------
subroutine CHEBFT(NIN, FIN, XIN, NCHCF, CHEBCF)
!---------------------------------------------------------------------
! The subroutine returns NCHCF coefficients of a Chebyshev
! polynomial fit to the function FIN(1:NIN) given as a
! function of any "radial" variable on the grid XIN(1:NIN)
! Example:
! call CHEBFT(NA1,TE,FP,5,CHOUT)
! out = PFITN(FP,CHOUT,5)
!-----------------------------------------------------------------------

use debugger, only: markloc, astra_stop

implicit none

integer, parameter :: NMAX=10
double precision, parameter :: PI=3.141592654

integer, intent(in) :: NIN, NCHCF
double precision, intent(in) , dimension(NIN)   :: FIN, XIN
double precision, intent(out), dimension(NCHCF) :: CHEBCF

integer :: K, J
double precision :: YA, YB, FAC, YD, SUM, PIOVN, BMA, BPA
double precision, dimension(NMAX) :: YC, YF

call markloc('CHEBFT')

if (NCHCF > NMAX) call astra_stop('>>> Chebyshev fit error: too high power')

YA = max(XIN(1), XIN(NIN))
YB = min(XIN(1), XIN(NIN))
PIOVN = PI/NCHCF

! The following two lines require that at the boundary
!     the fit concides with the original function
YD = COS(0.5*PIOVN)
YA = (2.*YA - YB*(1. - YD))/(1. + YD)

BMA = 0.5*(YB - YA)
BPA = 0.5*(YB + YA)
do K=1,NCHCF
    YC(K) = BPA + BMA*COS(PIOVN*(K - 0.5))
enddo

call TRANSF(NIN, FIN, XIN, NCHCF, YF, YC)

FAC = 2./NCHCF
do J=1, NCHCF
    SUM = 0.
    do K=1, NCHCF
        SUM = SUM + YF(K)*COS(PIOVN*(K - 0.5)*(J - 1.))
    enddo
    YC(J) = FAC*SUM
enddo

call CHEBPC(YC, CHEBCF, YF, NCHCF)

FAC = 1./BMA
do J=2, NCHCF
    CHEBCF(J) = CHEBCF(J)*FAC
    FAC = FAC/BMA
enddo
do J=1, NCHCF-1
    do K=NCHCF-1, J, -1
        CHEBCF(K) = CHEBCF(K) - BPA*CHEBCF(K+1)
    enddo
enddo

return
end subroutine CHEBFT

!---------------------------------------------------------------------
subroutine CHEBPC(C, D, F, N)
!---------------------------------------------------------------------
! F(1:N) array for internal use
!
!   Input: N number of coefficients
!  C(1:N) 
!   Output: D(1:N) array of Chebyshev coefficients
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: N
double precision, intent(in) , dimension(N) :: C
double precision, intent(out), dimension(N) :: D, F

integer :: J, K
double precision :: SV

D = 0.
F = 0.

D(1) = C(N)
do J=N-1, 2, -1
    do K=N-J+1, 2, -1
        SV = D(K)
        D(K) = 2.*D(K-1) - F(K)
        F(K) = SV
    enddo
    SV = D(1)
    D(1) = -F(1) + C(J)
    F(1) = SV
enddo
do J=N, 2, -1
    D(J) = D(J-1) - F(J)
enddo
D(1) = -F(1) + 0.5*C(1)

return
end subroutine CHEBPC

!---------------------------------------------------------------------
double precision function PFITN(x, cp, N)
!---------------------------------------------------------------------
! Nth order Polynomial FIT to a function f(a)
! CP are the N polynomial coefficients given
! PFITN = f(a) = \Sum {CP_j * a^(j-1)} ;   1<=j<=N
!---------------------------------------------------------------------

implicit none

integer, intent(in) :: N
double precision, intent(in) :: x, cp(N)

integer :: j

PFITN = CP(N)
do j=N-1, 1, -1
    PFITN = PFITN*x + CP(j)
enddo

return
end function PFITN
