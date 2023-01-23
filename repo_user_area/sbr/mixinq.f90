subroutine MIXINQ(OPTION, RECOND)
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
!   OPTION = 0,10,20,30
!  Kadomtsev-type reconnection for a single q=1 resonance 
!  surface. Reconnection condition: dt > RECOND
!   RECOND[sec] determines a minimal time interval
!  between the successive disruptions.
!  For this type of reconnection it is additionally 
!  required that q = 1 surface exists in the plasma 
!
!   OPTION = 1,11,21,31 
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
!   OPTION = 2,12,22,32
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
!
! Additional information (involving CFs and CARs ) is available for 
!   OPTION = 20, 21, 22.
!  See the examples below:
!  CAR1 gives psi* quantity before reconnection
!  CAR2 gives psi* quantity after reconnection
!  CMHD1 - relative radius of the 1st resonance surface
!  CMHD2 - relative radius of the 2nd resonance surface
!  CF3 - relative radius of the 3rd resonance surface
!  CMHD3 - relative radius of the reconnection region
!  CMHD4 - psi* maximal value
!
! Options 30, 31, 32 include (10+20), (11+21), (12+22).
!-----------------------------------------------------------------------

use parameter_inc, only: NRD
use const_inc, only: BTOR, RTOR, HRO, ROC, GP, GP2, NA1, NB1, &
    TIME, TSTART, TAU, TAUMIN, CMHD1, CMHD2, CMHD3, CMHD4, LEQ
use status_inc, only: TE, TI, NE, NI, VR, FP, CU, MU, IPOL, G22, G33, CAR1, CAR2

implicit none

double precision RHOS(3),PSIM(3),FPSTAR(NRD),WBPOLR,OPTION,QRES
integer IS(3),IOPT,KOPT,JNRES,j,jj,IMIX
double precision TMIX,RECOND,NEAVR,YMURES,YMUC,YXS,YRS,YWP,YX,   &
   HFP,YSIGN,YSIGNO,YFPC,YCM,YCM1,D1PS,D2PS,RHOMIX,YV0,DPS, &
   YV2,YV3,YV4,YV5,YV6,YNE,YNI,YTE,YTI,YX2,YCU,YM,YM1,YT,   &
   DEMAX,DIMAX,TEMAX,TIMAX,DEMIN,DIMIN,TEMIN,TIMIN,DDE,DDI, &
   DECONT,DICONT,WECONT,WICONT,YTE0,YTI0,YNE0,YDTE,YDTI
character*132 STRI
save TMIX
data TMIX/-99999./

! Store the ASTRA start time:
if (TMIX .lt. -99998.) TMIX = TSTART
! or the time of the 1st call:
! if (TMIX .lt. -99998.) TMIX = TIME
IOPT = OPTION+0.001

! KOPT = 0,1,2,3 - output option
KOPT = IOPT/10

! IOPT = 1,2,3 - reconnection option
IOPT = IOPT-10*KOPT

! Other value "n/m" can be used instead of 1 for the double tearing mode
QRES = 1.5
QRES = 1.
YMURES = 1./QRES
YMUC = (4.*MU(1)-MU(2))/3.
YFPC = (9.*FP(1)-FP(2))/8.
YCM = GP2*BTOR*HRO*HRO*YMURES
YCM1 = 0.5*YCM
JNRES = 0
FPSTAR(1) = FP(1)-YFPC-0.125*YCM
if (KOPT .ge. 2) CAR1(1) = FPSTAR(1)
if (KOPT .ge. 2) CAR2(1) = FPSTAR(1)
YSIGNO = sign(1.d0,YMUC-YMURES)

do  j=2,NB1
! "Psi" for a homogeneous current can be calculated as
!    HFP(1) = 0.125*YMURES*YCM
!    HFP(j) = HFP(j-1)+YMURES*YCM*(j-1)
! or as
   HFP = YCM1*(j-0.5)**2
   FPSTAR(j) = FP(j)-YFPC-HFP
   if (KOPT .ge. 2) CAR1(j) = FPSTAR(j)
   if (KOPT .ge. 2) CAR2(j) = FPSTAR(j)
   YSIGN  = sign(1.d0,MU(j)-YMURES)
   if (YSIGN.eq.YSIGNO) CYCLE
   YSIGNO = YSIGN
   if (j .eq. 2) CYCLE
   JNRES = JNRES+1
   if (JNRES .gt. 3) CYCLE
   YRS = (j-1+(YMURES-MU(j-1))/(MU(j)-MU(j-1)))*HRO
   RHOS(JNRES) = YRS
   IS(JNRES) = j
enddo
if (FPSTAR(NB1).gt.0. .and. (KOPT.eq.1 .or. KOPT.eq.3))   then
   write(*,*)MU
   pause'Inverse q-profile'
   endif
if ( KOPT .ge. 2)   then
   CMHD1 = 0.
   CMHD2 = 0.
!    CF3 = 0.
   endif
if (JNRES .eq. 0)   return  ! No resonance found
if ( KOPT .ge. 2)   CMHD1 = RHOS(1)/ROC

! One resonance surface
if (JNRES .gt. 1)   goto   20
! The two lines allow to switch the reconnection condition between
! the OPTION = 1 type with the condition rho_s > RECOND*ROC and
! the OPTION = 0 type with the condition saw_period > RECOND.
! if ( IOPT .gt. 1 .and. RHOS(JNRES) .gt. RECOND*ROC)   then
if ( IOPT .gt. 1 .and. TIME-TMIX+.5*TAU .gt. RECOND)   then
   if (KOPT .eq. 1 .or. KOPT .eq. 3)  write(*,*) &
       'One resonance',RHOS(1)/ROC,' found instead of ', &
       'two expected. ',"Kadomtsev's reconnection done."
   goto  12
endif
 11 if ( IOPT .eq. 1 .and. RHOS(JNRES) .gt. RECOND*ROC .or. &
         IOPT .eq. 0 .and. TIME-TMIX+.5*TAU .gt. RECOND)   goto   12
return

 12 if (JNRES .eq. 2 .and. (KOPT.eq.1 .or. KOPT.eq.3)) &
     write(*,*)'>>> MIXINT >>> Double tearing mode ignored'
j = IS(JNRES)
YXS = RHOS(JNRES)/HRO-j+0.5
D1PS = FPSTAR(j+1)-FPSTAR(j-1)
D2PS = FPSTAR(j+1)-2.*FPSTAR(j)+FPSTAR(j-1)
PSIM(JNRES) = 0.5*YXS*(D1PS+YXS*D2PS)+FPSTAR(j)
do   14   j=IS(JNRES)-1,NA1
   if (FPSTAR(j) .gt. 0.)   goto   14
   RHOMIX=(j-1.5-FPSTAR(j-1)/(FPSTAR(j)-FPSTAR(j-1)))*HRO
   IMIX = j-1
   goto  15
 14 continue
if (FPSTAR(NA1) .gt. 0.)   return
 15 continue

if (KOPT .ge. 2)  CMHD4 = PSIM(JNRES)
if (KOPT .ge. 2)  CMHD3 = RHOMIX/ROC

! Conserving quantities:
!   YN  = NEAVR(ROC)
!   YWE = WER(ROC)
!   YWI = WIR(ROC)

YWP = WBPOLR(ROC)
YTE0 = TE(1)
YTI0 = TI(1)
YNE0 = NE(1)
YV0=0.
DECONT=0.
DICONT=0.
WECONT=0.
WICONT=0.
if (KOPT .ge. 1)write(*,*)
do j=1,IMIX
   FPSTAR(j) = PSIM(JNRES)*(1.-((j-0.5)*HRO/RHOMIX)**4)
   FP(j) = FPSTAR(j)+YFPC+YCM1*(j-0.5)**2
   if (KOPT .ge. 2) CAR2(j) = FPSTAR(j)
   YV0  = YV0+VR(J)
   DECONT = DECONT+NE(J)*VR(J)
   DICONT = DICONT+NI(J)*VR(J)
   WECONT = WECONT+NE(J)*TE(J)*VR(J)
   WICONT = WICONT+NI(J)*TI(J)*VR(J)
enddo

goto   50
     
! Two resonance surfaces
 20 if ( KOPT .ge. 2)   CMHD2 = RHOS(2)/ROC
if (JNRES .gt. 2)   goto   30
do jj =1,JNRES
   j = IS(jj)
   YXS = RHOS(jj)/HRO-j+0.5
   D1PS = FPSTAR(j+1)-FPSTAR(j-1)
   D2PS = FPSTAR(j+1)-2.*FPSTAR(j)+FPSTAR(j-1)
   PSIM(jj) = 0.5*YXS*(D1PS+YXS*D2PS)+FPSTAR(j)
enddo

if (PSIM(2) .lt. 0.)   return

! Enforced switch to the OPTION = 0 or 1 type reconnection.
! The outermost resonance is taking into account only.
if (IOPT .lt. 2)   goto   11
 
! This "goto" allows OPTION = 0 or 1 type reconnection.
! The reconnection takes into account the outermost resonance only.
! if (IOPT .eq. 2 .and. abs(PSIM(1)) .lt. .1*PSIM(2))  goto  12

do j=IS(JNRES)-1,NA1
    if (FPSTAR(j) .gt. PSIM(1)) CYCLE
    DPS = (PSIM(1)-FPSTAR(j-1))/(FPSTAR(j)-FPSTAR(j-1))
    RHOMIX=(j-1.5-DPS)*HRO
    IMIX = j-1
    EXIT
enddo

if (KOPT .ge. 2)  CMHD4 = PSIM(JNRES)
if (KOPT .ge. 2)  CMHD3 = RHOMIX/ROC

do jj =1,JNRES
   j = IS(jj)
   YXS = RHOS(jj)/HRO-j+0.5
   if (jj .eq. 1)   then
      DEMAX = NE(j-1)+YXS*(NE(j)-NE(j-1))
      DIMAX = NI(j-1)+YXS*(NI(j)-NI(j-1))
      TEMAX = TE(j-1)+YXS*(TE(j)-TE(j-1))
      TIMAX = TI(j-1)+YXS*(TI(j)-TI(j-1))
   endif
   if (jj .eq. 2)   then
      DEMIN = NE(j-1)+YXS*(NE(j)-NE(j-1))
      DIMIN = NI(j-1)+YXS*(NI(j)-NI(j-1))
      TEMIN = TE(j-1)+YXS*(TE(j)-TE(j-1))
      TIMIN = TI(j-1)+YXS*(TI(j)-TI(j-1))
   endif
enddo

YWP = WBPOLR(ROC)
YTE0 = TE(1)
YTI0 = TI(1)
YNE0 = NE(1)
YV0=0.
YV2=0.
YV3=0.
YV4=0.
YV5=0.
YV6=0.
DECONT=0.
DICONT=0.
WECONT=0.
WICONT=0.
do j=1,IMIX
   YX = (j-0.5)*HRO/RHOMIX
   FPSTAR(j) = (PSIM(1)-PSIM(2))*(1.-(1.-(YX)**4)**2)+PSIM(2)
   FP(j) = FPSTAR(j)+YFPC+YCM1*(j-0.5)**2
   if (KOPT .ge. 2) CAR2(j) = FPSTAR(j)
   YV0 = YV0+VR(J)
   YV2 = YV2+VR(J)*YX**2
   YV3 = YV3+VR(J)*YX**3
   YV4 = YV4+VR(J)*YX**4
   YV5 = YV5+VR(J)*YX**5
   YV6 = YV6+VR(J)*YX**6
   DECONT = DECONT+NE(J)*VR(J)
   DICONT = DICONT+NI(J)*VR(J)
   WECONT = WECONT+NE(J)*TE(J)*VR(J)
   WICONT = WICONT+NI(J)*TI(J)*VR(J)
enddo

DDE = DEMAX-DEMIN
DDI = DIMAX-DIMIN
YDTE = TEMAX-TEMIN
YDTI = TIMAX-TIMIN
YNE = (DECONT-DEMIN*YV0-DDE*YV2)/(YV2-YV3)
YNI = (DICONT-DIMIN*YV0-DDI*YV2)/(YV2-YV3)
YTE = (WECONT-DEMIN*TEMIN*YV0-(DEMIN*YDTE+TEMIN*DDE)*YV2 &
          -YNE*(TEMIN*(YV2-YV3)+YDTE*(YV4-YV5))-DDE*YDTE*YV4)/ &
          (DEMIN*(YV2-YV3)+(DDE+YNE)*(YV4-YV5)-YNE*(YV5-YV6))
YTI = (WICONT-DIMIN*TIMIN*YV0-(DIMIN*YDTI+TIMIN*DDI)*YV2 &
          -YNI*(TIMIN*(YV2-YV3)+YDTI*(YV4-YV5))-DDI*YDTI*YV4)/ &
          (DIMIN*(YV2-YV3)+(DDI+YNI)*(YV4-YV5)-YNI*(YV5-YV6))
do j=1,IMIX
   YX = (j-0.5)*HRO/RHOMIX
   YX2 = YX*YX
enddo

goto   50

! Three resonance surfaces
 30 continue

if (JNRES .gt. 3)   goto   40
if (IOPT .lt. 2) then
   goto 11
endif
return

 40 continue
if (JNRES .ge. 4 .and. (KOPT.eq.2 .or. KOPT.eq.3)) &
      write(*,*)'>>> MIXER: too many resonant surfaces'
return

 50 YCU=2.5*BTOR/(GP*RTOR*HRO)
YM=0.

do J=1,IMIX
   MU(J)=min(YMURES,(FP(J+1)-FP(J))/(J*YCM))
   YM1=YM
   YM=J*G22(J)*MU(J)
   CU(J)=YCU*G33(J)*IPOL(J)**3*(YM-YM1)/(J-0.5)
enddo

YDTE = max(0.d0,(YWP-WBPOLR(ROC))/(0.0024*DECONT*HRO))

YT = TIME-TMIX
if (YT .lt. 0.75) then
write(STRI,'(1A9,1F7.2,1A4)')' Period =',1000.*YT,' ms;'
else
write(STRI,'(1A9,1F7.3,1A4)')' Period =',YT,' s; '
endif

YT = YTE0-TE(1)
if (YT .lt. 0.75) then
   write(STRI(21:),'(1A,1F6.1,1A,1F5.2)') &
      '  Te amplitude: ',1000.*YT,' eV;   r0/a =',RHOMIX/ROC
else
   write(STRI(21:),'(1A,1F6.3,1A,1F5.2)') &
      '  Te amplitude: ',YT,' keV;  r0/a =',RHOMIX/ROC
endif

if (KOPT .ge. 1) write(*,'(A)') STRI
YT = max(1.d-4,(TIME-TMIX)/5.)

TAU = TAUMIN
TMIX = TIME

return
end subroutine mixinq
