SUBROUTINE MIXINT(OPTION, RECOND)

!-----------------------------------------------------------------------
!     (G.Pereverzev 17-NOV-94)
!     (     updated 14-JUL-99)
!-----------------------------------------------------------------------
! ! Note !   The Astra compiler places this subroutine
!  AFTER all the transport equations solved
!-----------------------------------------------------------------------
! B.B.Kadomtsev, Disruptive instability in tokamaks,
!        Sov. J. Plasma Phys. Vol.1, No.5, Sept.-Oct. (1975) pp.389-391.
! V.V.Parail, G.V.Pereverzev, Internal disruption in a tokamak,
!        Sov. J. Plasma Phys. Vol.6, No.1, Jan-Feb. (1980) pp.14-17.
!-----------------------------------------------------------------------
!   OPTION = 0,10
!  Kadomtsev-type reconnection for a single q=1 resonance 
!  surface. Reconnection condition: dt > RECOND
!   RECOND[sec] determines a minimal time interval
!  between the successive disruptions.
!  For this type of reconnection it is additionally 
!  required that q = 1 surface exists in the plasma
!
!   OPTION = 1,11 
!  Kadomtsev-type reconnection for a single q=1 resonance 
!  surface. Reconnection condition: rho_s > RECOND*ROC
!   RECOND[d/l] minimal resonance radius
!
!  In both cases above, when multiple resonance encountered
!  the outermost resonance surface is taken into account
!  only. This can cause a contradiction when
!  (a) a resonance surface first appears long after
!      the prescribed time interval (OPTION = 0)
!  (b) a resonance surface appears far outside
!      rho = RECOND*ROC (OPTION = 1)
!  In those cases, message(*) can be printed (see below).
!
!   OPTION = 2,12
!  q = 1 double tearing mode central reconnection.
!   RECOND[s] is ignored when two resonance surfaces occur.
!  If only one resonance surface is found then
!  the OPTION = 0 type reconnection is done with
!  the reconnection condition saw_period >= RECOND.
!  This option is useful for the first call.
!  In this case, message(*) can be printed (see below).
!
! Note: The subroutine affects control times TAUMIN and TAUMAX 
!  so that TAUMAX = min(TAUMAX,PERIOD/5.)
!   TAUMIN = min(TAUMIN,TAUMAX/2.)
!  where PERIOD is the recent sawtooth duration
!
! To print energy conservation check enable lines marked with "Co"
!
!     (*) To enable message prints use 
!   OPTION = 10, 11, 12, respectively.
!-----------------------------------------------------------------------

use pi_const, only: GP, GP2
use scalars, only: BTOR, RTOR, HRO, ROC, NA1, NB1, &
    TIME, TSTART, TAU, TAUMIN, LEQ
use status, only: NRD, TE, TI, NE, NI, VR, FP, CU, MU, IPOL, G22, G33

implicit none

double precision, intent(in) :: OPTION, RECOND
double precision :: RHOS(3), PSIM(3), FPSTAR(NRD), WBPOLR, QRES
logical :: verbose
integer :: IS(3), IOPT, JNRES, j, jj, IMIX
double precision :: TMIX, YMURES, YMUC, YXS, YWP, YX, &
    HFP, YSIGN, YSIGNO, YFPC, YCM, YCM1, D1PS, D2PS, RHOMIX, YV0, DPS, &
    YV2, YV3, YV4, YV5, YV6, YNE, YNI, YTE, YTI, YX2, YCU, YM, YM1, YT, &
    DEMAX, DIMAX, TEMAX, TIMAX, DEMIN, DIMIN, TEMIN, TIMIN, DDE, DDI, &
    DECONT, DICONT, WECONT, WICONT, YTE0, YTI0, YNE0, YDTE, YDTI, delta_t
character(len=132) :: STRI

save TMIX
data TMIX/-99999./

! Store the ASTRA start time:
if (TMIX < -99998.) TMIX = TSTART

IOPT = nint(OPTION + 0.001)
if (iopt > 9) then
    verbose = .TRUE.
else
    verbose = .FALSE.
endif

! IOPT = 1, 2, 3 - reconnection option
IOPT = MOD(IOPT, 10)
delta_t = TIME - TMIX + 0.5*TAU

! Other value "n/m" can be used instead of 1 for the double tearing mode
!QRES = 1.5
QRES = 1.
YMURES = 1./QRES
YMUC = (4.*MU(1) - MU(2))/3.
YFPC = (9.*FP(1) - FP(2))/8.
YCM = GP2 * BTOR * HRO**2 * YMURES
YCM1 = 0.5*YCM
JNRES = 0
FPSTAR(1) = FP(1) - YFPC - 0.125*YCM
YSIGNO = sign(1.d0, YMUC - YMURES)
do j=2, NB1
! "Psi" for a homogeneous current can be calculated as
!    HFP(1) = 0.125*YMURES*YCM
!    HFP(j) = HFP(j-1)+YMURES*YCM*(j-1)
! or as
    HFP = YCM1*(j - 0.5)**2
    FPSTAR(j) = FP(j) - YFPC - HFP
    YSIGN = sign(1.d0, MU(j) - YMURES)
    if (YSIGN /= YSIGNO) then
        YSIGNO = YSIGN
        if (j > 2) then
            JNRES = JNRES + 1
            if (JNRES <= 3) then
                RHOS(JNRES) = (j - 1 + (YMURES - MU(j-1))/(MU(j) - MU(j-1)))*HRO
                IS(JNRES) = j
            endif
        endif
    endif
enddo

if (FPSTAR(NB1) > 0. .and. verbose)   then
    write(*, *) 'Inverse q-profile', MU(1:NA1)
endif

if (JNRES == 0) return  ! No resonance found

! One resonance surface
if (JNRES > 1) goto 20
! The two lines allow to switch the reconnection condition between
! the OPTION = 1 type with the condition rho_s > RECOND*ROC and
! the OPTION = 0 type with the condition saw_period > RECOND.
if (IOPT > 1 .and. delta_t > RECOND) then
    if (verbose) write(*, *) 'One resonance', RHOS(1)/ROC, &
        ' found instead of two expected. ', "Kadomtsev's reconnection done."
    goto 12
endif

11 continue

if ((IOPT == 1 .and. RHOS(JNRES) > RECOND*ROC) .or. (IOPT == 0 .and. delta_t > RECOND)) goto 12

return

12 continue

if (JNRES == 2 .and. verbose) write(*, *) '>>> MIXINT >>> Double tearing mode ignored'
j = IS(JNRES)
YXS = RHOS(JNRES)/HRO-j+0.5
D1PS = FPSTAR(j+1) -   FPSTAR(j-1)
D2PS = FPSTAR(j+1) -2.*FPSTAR(j) + FPSTAR(j-1)
PSIM(JNRES) = 0.5*YXS*(D1PS + YXS*D2PS) + FPSTAR(j)
do j=IS(JNRES)-1, NA1
    if (FPSTAR(j) > 0.) CYCLE
    RHOMIX = (j - 1.5 - FPSTAR(j-1)/(FPSTAR(j) - FPSTAR(j-1)))*HRO
    IMIX = j - 1
    goto 15
enddo

if (FPSTAR(NA1) > 0.) return

15 continue

YWP = WBPOLR(ROC)
YTE0 = TE(1)
YTI0 = TI(1)
YNE0 = NE(1)
YV0  = 0.
DECONT = 0.
DICONT = 0.
WECONT = 0.
WICONT = 0.

do j=1, IMIX
    FPSTAR(j) = PSIM(JNRES)*(1. - ((j - 0.5)*HRO/RHOMIX)**4)
    FP(j) = FPSTAR(j) + YFPC + YCM1*(j-0.5)**2
    YV0  = YV0 + VR(J)
    DECONT = DECONT + NE(J)*VR(J)
    DICONT = DICONT + NI(J)*VR(J)
    WECONT = WECONT + NE(J)*TE(J)*VR(J)
    WICONT = WICONT + NI(J)*TI(J)*VR(J)
enddo
do j=1, IMIX
    if (LEQ(1) > 0) NE(j) = DECONT/YV0
    if (LEQ(1) > 0) NI(j) = DICONT/YV0
    if (LEQ(2) > 0) TE(j) = WECONT/DECONT
    if (LEQ(3) > 0) TI(j) = WICONT/DICONT
enddo

goto 50
 
! Two resonance surfaces

20 continue

if (JNRES > 2) goto 30

do jj =1, JNRES
    j = IS(jj)
    YXS = RHOS(jj)/HRO - j + 0.5
    D1PS = FPSTAR(j+1) - FPSTAR(j-1)
    D2PS = FPSTAR(j+1) - 2.*FPSTAR(j) + FPSTAR(j-1)
    PSIM(jj) = 0.5*YXS*(D1PS + YXS*D2PS) + FPSTAR(j)
enddo

if (PSIM(2) < 0.) return

! Enforced switch to the OPTION = 0 or 1 type reconnection.
! The outermost resonance is taking into account only.
if (IOPT < 2) goto 11
 
! This "goto" allows OPTION = 0 or 1 type reconnection.
! The reconnection takes into account the outermost resonance only.

do j=IS(JNRES)-1, NA1
    if (FPSTAR(j) > PSIM(1)) CYCLE
    DPS = (PSIM(1) - FPSTAR(j-1))/(FPSTAR(j) - FPSTAR(j-1))
    RHOMIX = (j - 1.5 - DPS)*HRO
    IMIX = j-  1
    EXIT
enddo

do jj=1, JNRES
    j = IS(jj)
    YXS = RHOS(jj)/HRO - j + 0.5
    if (jj == 1) then
        DEMAX = NE(j-1) + YXS*(NE(j) - NE(j-1))
        DIMAX = NI(j-1) + YXS*(NI(j) - NI(j-1))
        TEMAX = TE(j-1) + YXS*(TE(j) - TE(j-1))
        TIMAX = TI(j-1) + YXS*(TI(j) - TI(j-1))
    endif
    if (jj == 2) then
        DEMIN = NE(j-1) + YXS*(NE(j) - NE(j-1))
        DIMIN = NI(j-1) + YXS*(NI(j) - NI(j-1))
        TEMIN = TE(j-1) + YXS*(TE(j) - TE(j-1))
        TIMIN = TI(j-1) + YXS*(TI(j) - TI(j-1))
    endif
enddo

YWP = WBPOLR(ROC)
YTE0 = TE(1)
YTI0 = TI(1)
YNE0 = NE(1)
YV0 = 0.
YV2 = 0.
YV3 = 0.
YV4 = 0.
YV5 = 0.
YV6 = 0.
DECONT = 0.
DICONT = 0.
WECONT = 0.
WICONT = 0.
do j=1, IMIX
    YX = (j-0.5)*HRO/RHOMIX
    FPSTAR(j) = (PSIM(1)-PSIM(2))*(1. - (1. - YX**4)**2) + PSIM(2)
    FP(j) = FPSTAR(j) + YFPC + YCM1*(j-0.5)**2
    YV0 = YV0 + VR(J)
    YV2 = YV2 + VR(J)*YX**2
    YV3 = YV3 + VR(J)*YX**3
    YV4 = YV4 + VR(J)*YX**4
    YV5 = YV5 + VR(J)*YX**5
    YV6 = YV6 + VR(J)*YX**6
    DECONT = DECONT + NE(J)*VR(J)
    DICONT = DICONT + NI(J)*VR(J)
    WECONT = WECONT + NE(J)*TE(J)*VR(J)
    WICONT = WICONT + NI(J)*TI(J)*VR(J)
enddo
DDE  = DEMAX - DEMIN
DDI  = DIMAX - DIMIN
YDTE = TEMAX - TEMIN
YDTI = TIMAX - TIMIN
YNE = (DECONT - DEMIN*YV0 - DDE*YV2)/(YV2 - YV3)
YNI = (DICONT - DIMIN*YV0 - DDI*YV2)/(YV2 - YV3)
YTE = (WECONT - DEMIN*TEMIN*YV0 - (DEMIN*YDTE + TEMIN*DDE)*YV2 &
      -YNE*(TEMIN*(YV2 - YV3) + YDTE*(YV4 - YV5)) - DDE*YDTE*YV4)/ &
           (DEMIN*(YV2 - YV3) + (DDE + YNE)*(YV4 - YV5) - YNE*(YV5 - YV6))
YTI = (WICONT - DIMIN*TIMIN*YV0 - (DIMIN*YDTI + TIMIN*DDI)*YV2 &
      -YNI*(TIMIN*(YV2 - YV3) + YDTI*(YV4 - YV5)) - DDI*YDTI*YV4)/ &
           (DIMIN*(YV2 - YV3) + (DDI + YNI)*(YV4 - YV5) - YNI*(YV5 - YV6))

do j=1, IMIX
    YX = (j-0.5)*HRO/RHOMIX
    YX2 = YX*YX
    if (LEQ(1) > 0) NE(j) = DEMIN + (DDE  + YNE*(1. - YX))*YX2
    if (LEQ(1) > 0) NI(j) = DIMIN + (DDI  + YNI*(1. - YX))*YX2
    if (LEQ(2) > 0) TE(j) = TEMIN + (YDTE + YTE*(1. - YX))*YX2
    if (LEQ(3) > 0) TI(j) = TIMIN + (YDTI + YTI*(1. - YX))*YX2
enddo

goto 50

! Three resonance surfaces
30 continue

if (JNRES > 3) goto 40

if (IOPT < 2) goto 11

return

40 continue

if (JNRES >= 4 .and. verbose) write(*, *) '>>> MIXER: too many resonant surfaces'

return

50 continue

YCU = 2.5*BTOR/(GP*RTOR*HRO)
YM = 0.

do J=1, IMIX
    MU(J) = min(YMURES, (FP(J+1) - FP(J))/(J*YCM))
    YM1 = YM
    YM = J*G22(J)*MU(J)
    CU(J) = YCU*G33(J)*IPOL(J)**3*(YM - YM1)/(J-0.5)
enddo

YDTE = max(0.d0, (YWP - WBPOLR(ROC))/(0.0024*DECONT*HRO))
do j=1, IMIX
    if (LEQ(2) > 0) TE(j) = TE(j) + YDTE
enddo

YT = TIME - TMIX
if (YT < 0.75) then
    write(STRI, '(1A9, 1F7.2, 1A4)')' Period =', 1000.*YT, ' ms;'
else
    write(STRI, '(1A9, 1F7.3, 1A4)')' Period =', YT, ' s; '
endif

YT = YTE0 - TE(1)
if (YT < 0.75) then
    write(STRI(21:), '(1A, 1F6.1, 1A, 1F5.2)') '  Te amplitude: ', &
        1000.*YT, ' eV;   r0/a =', RHOMIX/ROC
else
    write(STRI(21:), '(1A, 1F6.3, 1A, 1F5.2)') '  Te amplitude: ', &
        YT, ' keV;  r0/a =', RHOMIX/ROC
endif

if (verbose) write(*, '(A)') STRI
YT = max(1.d-4, (TIME - TMIX)/5.)

TAU = TAUMIN
TMIX = TIME

end subroutine mixint
