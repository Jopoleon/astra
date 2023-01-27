!----------------------------------------------------------------------|
subroutine METRIC
!----------------------------------------------------------------------|

use outcmn_inc, only: CPT, CPTEQL
use status_inc, only: VRO, VR
use const_inc, only: IPART, FTO, FTN, ROC, GP, BTOR, ROCO, RTOR, SHIFT, &
   ABC, ELONG, TRIAN, NB1, NEQUIL, MEQUIL, LEQ, IPEQL, &
   TIME, TSTART, TIMEQL, DTEQL, BTN
use debugger, only: markloc, astra_stop

implicit none

integer :: jrho, jexit, NDTEQUILMY, equil_solver
double precision :: ROC3A
character(len=120) :: err_msg

call markloc('METRIC')

call ADDTIME(CPT)

if (IPART == 1) then ! do only at initiation

   FTN = FTO
   BTN = BTOR
   ROC = sqrt(FTO/GP/BTOR)
   ROCO = ROC
   do jrho=1, NB1
      VRO(jrho) = VR(jrho)
   enddo
endif

LEQ(5) = nint(IPEQL)

SELECT CASE(LEQ(5))

CASE(-2)  ! Cylindircal case, No equilibrium solver. No toroidicity
   call EQCYL
   call RHSEQ

CASE(-1)  ! Take metrics from exp/data_file
   ROC = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
   FTO = GP*BTOR*ROC**2
   call SETEXM   ! Main grid: (jj-0.5)*h
   call RHSEQ  !this computes ffprime and pprime

CASE(0) ! No equilibrium solver (NEQUIL=0) .or. data initiation @ 1st entry
   call EQGUESS

CASE(1)  ! EMEQ
   if (TIME == TSTART) NDTEQUILMY=0
   if (TIME > TSTART) NDTEQUILMY=1
   if (TIME-TIMEQL >= NDTEQUILMY*DTEQL) then
      call RHSEQ   ! Define p', FF', j_tor=CUTOR
      call A2EMEQ(jexit)
      if (jexit /= 0) then
         err_msg = 'Equilibrium problem at the initial iterations'
         if (IPART == 1) call astra_stop(err_msg)
      endif
      call ADDTIME(CPTEQL)
      TIMEQL = TIME
   endif

CASE(3)  ! SPIDER iterations
   if (TIME == TSTART) NDTEQUILMY=0
   if (TIME > TSTART) NDTEQUILMY=1
   if (TIME-TIMEQL >= NDTEQUILMY*DTEQL) then
      call RHSEQ
      call ADDTIME(CPTEQL)
      TIMEQL = TIME
   endif

CASE(4: 5)  ! SPIDER, FEQIS
   if (LEQ(5) == 4) then
      equil_solver = 3
   else
      equil_solver = 101
   endif
   if (TIME == TSTART) NDTEQUILMY=0
   if (TIME > TSTART) NDTEQUILMY=1
   if (TIME-TIMEQL >= NDTEQUILMY*DTEQL) then
      call RHSEQ   ! Define p', FF', j_tor=CUTOR
      if (NEQUIL > 0) then
         NEQUIL = -NEQUIL
         MEQUIL = -MEQUIL
      endif
      call A2GSSOLVER_EF(equil_solver)
      call ADDTIME(CPTEQL)
      TIMEQL = TIME
   endif

END SELECT

call ADDTIME(CPT)

end subroutine METRIC

!======================================================================|
subroutine EQCYL
!----------------------------------------------------------------------|
! Quasi-cylindrical assignment: Called if NEQUIL==LEQ(5)==-2
!----------------------------------------------------------------------|
!  In: RTOR, SHIFT, ABC, ELONG, TRIAN, NA1, NB1
! Out: NA, HRO, ROC, RHO(j), DRODA, VOLUM, IPOL, G33, GRADRO, G11, G22, SLAT
!  BDB02, B0DB2, BDB0, FOFB, BMAXT, BMINT, DRODA, GRADRO
!----------------------------------------------------------------------|

use status_inc, only: RHO, XRHO, VR, VRS, AMETR, SHIF, SHIV, &
   ELON, TRIA, SLAT, G11, G22, G33, G41, G42, G43, G44, G45, &
   BDB0, BDB02, B0DB2, IPOL, MU, &
   FOFB, BMAXT, BMINT, DRODA, GRADRO, VOLUM
use const_inc, only: VOLUME, GP, GP2, RTOR, BTOR, DELOUT, &
   ABC, HRO, ROC, FTO, ROWALL, NA, NA1, NB1
use debugger, only: markloc, debug

implicit none

integer :: j

call markloc('EQCYL', debug_lev=3*debug)

VOLUME = GP2*GP*RTOR*ABC*ABC
ROC = ABC
FTO = GP*BTOR*ROC**2
HRO = ABC/(NA1-0.5)
HRO = HRO
do j=1, NB1
   RHO(J) = XRHO(J)*ROC
   if (RHO(j) <= ROWALL) NA = j
   G22(J) = J*HRO
   VR(J)  = GP2*GP2*RTOR*RHO(J)
   VRS(j) = GP2*GP2*RTOR*G22(J)
   AMETR(J) = RHO(J)
   SHIF(J) = 0.
   SHIV(J) = 0.
   ELON(J) = 1.
   TRIA(J) = 0.
   IPOL(J) = 1.
   G33(J) = 1.
   G11(J) = VRS(j) 
   SLAT(J) = VRS(j)
   BDB02(j) = 1.+(RHO(j)*MU(j)/RTOR)**2
   B0DB2(j) = 1./BDB02(j)
   BDB0(j) = sqrt(BDB02(j))
   FOFB(j) = 1.
   BMAXT(j) = BTOR*BDB0(j)
   BMINT(j) = BMAXT(j)
   DRODA(j) = 1.
   GRADRO(j) = 1.
   G41(J) = 1.0
   G42(J) = GRADRO(J)
   G43(J) = GRADRO(J)
   G44(J) = G11(J)/VRS(J)
   G45(J) = G11(J)/VRS(J)
enddo
if (NA < NB1) then
   NB1 = NA
   DELOUT(13) = NB1
endif
NA = NA1 - 1
call INTEGR_EF(RHO, 1, VR, VOLUM, NB1)
VOLUME = VOLUM(NA1)

end subroutine EQCYL

!======================================================================|
subroutine EQGUESS
!----------------------------------------------------------------------|
! Guessed equibrium: Called if LEQ(5)==0 or data_initiation @ 1st_entry
!----------------------------------------------------------------------|
!  In: RTOR, SHIFT, ABC, ELONG, TRIAN, NA1, NB1, HRO
! Out: NA, RHO(j), DRODA, VOLUM, IPOL, G33, GRADRO, G11, G22, SLAT
!Efable: OUT: G41 G42 G43 G44 G45
!
!alling:
!       SETGEO -> SHIF, SHIV, ELONG, TRIA, DRODA, AMETR
! NEWGRD -> Takes ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC, AMETR(NA1)
!    Returns:
!    NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, AMETR(j>NA1)
! EDCELL -> 
!----------------------------------------------------------------------|

use const_inc, only: ROC, RTOR, SHIFT, ABC, ELONG, TRIAN, &
   FTO, GP, GP2, BTOR, NA, NA1, NB1, NAB, HRO, IPL, VOLUME
use status_inc, only: RHO, VR, VRS, AMETR, SHIF, &
   ELON, TRIA, SLAT, G11, G22, G33, G41, G42, G43, G44, G45, &
   BDB0, BDB02, B0DB2, IPOL, MU, FP, FV, SHEAR, &
   FOFB, BMAXT, BMINT, DRODA, GRADRO, VOLUM
use debugger, only: markloc, debug

implicit none

integer :: j, j1, JNA1O
double precision :: ROC3A, YA, YAS, YES, YDS, YDV, DFPDR, YR1
!----------------------------------------------------------------------|

call markloc('EQGUESS', debug_lev=3*debug)

ROC = ROC3A(RTOR, SHIFT, ABC, ELONG, TRIAN)
FTO = GP*BTOR*ROC**2
NA = NA1-1

! Compute NB1 
call NEWGRD

! Define AMETR, SHIF, ELON, TRIA, SHIV

call SETGEO(0)

YDV = 0.

do J=1, NB1
   if(j < nb1) then
      j1 = j + 1
   else
      j1 = j
   endif
   YA  = 0.5*(AMETR(J) + AMETR(J1))
   YES = 0.5*(ELON(J)  + ELON(J1))
   YDS = 0.5*(TRIA(J)  + TRIA(J1))
   YAS = 0.5*(SHIF(J)  + SHIF(J1))
   VOLUM(J) = GP*GP2*YA**2*YES*(RTOR + YAS - 0.25*YA*YDS)
   VR(J) = (VOLUM(J) - YDV)/HRO
   YDV = VOLUM(J)
   YR1 = min(AMETR(j)/2., 1.d-2)
   YAS = ROC3A(RTOR, SHIF(j), AMETR(j),       ELON(j), TRIA(j))
   YDS = ROC3A(RTOR, SHIF(j), AMETR(j) - YR1, ELON(j), TRIA(j))
   DRODA(J) = (YAS - YDS)/YR1
enddo
VOLUM(NA1) = GP*GP2*ABC**2*ELONG*(RTOR + SHIFT - 0.25*ABC*TRIAN)
DFPDR = (FP(NA1) - FP(NA) - (FV(NA1) - FV(NA)))/HRO
JNA1O = NA1
! Input:  ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC, AMETR(NA1)
call NEWGRD
! Output: NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, AMETR(j>NA1)

if (abs(DFPDR-IPL) < .2*IPL) then
   call EDCELL(JNA1O)
endif

! If NA1 does not change then EDCELL corrects FP(NA1) only.

!====== Definition --------------------------- Approximation ----------|
!      gradRHO = DRODA
!     <|grad(a)|> = sqrt(G1)/DRODA 
! G1=<(gradRHO)**2>;    G1=DRODA**2
! G11=VR*G1=VR*<(gradRHO)**2>   G11=VR*G1
! G2=<(gradRHO/R)**2>*VR/4/pi**2;   G22=R*G2/J;
! G22=VR*R*<(gradRHO/R)**2>/(GP2)**2/IPOL; C22=G11/R/(GP2)**2/IPOL
! SLAT=VR*DRODA*<|gradA|>;   SLAT=VR*sqrt(G1)
do J=1, NB1
   IPOL(J) = 1.
   G33(J) = (RTOR/(RTOR + SHIF(J)))**2
   GRADRO(J) = DRODA(J)
   if (j < NB1) VRS(j) = 0.5*(VR(J + 1) + VR(j))
   SLAT(J) = VRS(J)*DRODA(J)  
   G11(J) = VRS(J)*DRODA(J)**2
   G22(J) = G11(J)/GP2**2/(RTOR + SHIF(J))
   G41(J) = 1.0
   G42(J) = GRADRO(J)
   G43(J) = GRADRO(J)
   G44(J) = G11(J)/VRS(J)
   G45(J) = G11(J)/VRS(J)
enddo

!   BDB02 - <B**2/B0**2>
!   B0DB2 - <B0**2/B**2>   <(R/R0)^2>
!   BMAXT - BMAXT
!   BMINT - BMINT
!   BDB0  - <B/BTOR>
!   FOFB  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
!   SHEAR -  d[ln(q)]/d[ln(rho)] (to replace fml/shear)
! Alternative definition for BDB02 (G.Pereverzev 10.02.99)
do J=1, NA1
   if (j == 1) then
      SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
   elseif (j < NA) then
      SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
   else
      SHEAR(J) = HRO*(FP(NA1) - FP(NA))/HRO - FP(NA) + FP(NA-1)
      SHEAR(J) = 2.*SHEAR(J)*(2.*HRO - HRO)
      SHEAR(J) = SHEAR(J)/(2.*HRO*MU(NA1) - HRO*MU(NA))
   endif
   SHEAR(J) = 1. -SHEAR(J)/(GP*BTOR*HRO**2)
   YR1      = RHO(j)*G22(j)*(MU(j)/RTOR)**2
   BDB02(j) = (1. + YR1)*G33(J)*IPOL(J)**2
   B0DB2(j) = ((RTOR + SHIF(J))/RTOR)**2 + .75*(AMETR(j)/RTOR)**2
   BMAXT(j) = BTOR*RTOR/(RTOR + SHIF(j) - AMETR(j))
   BMINT(j) = BTOR*RTOR/(RTOR + SHIF(j) + AMETR(j))
   BDB0(j)  = (RTOR/(RTOR + SHIF(j)))
!  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
   YR1     = (RTOR + SHIF(j) - AMETR(j))/(RTOR + SHIF(j))
   FOFB(j) = 1. - sqrt(YR1)*(1. + 0.5*YR1)
enddo
do J=NA1, NAB
   SHEAR(j) = SHEAR(NA1)
   BDB02(j) = BDB02(NA1)
   B0DB2(j) = B0DB2(NA1)
   BMAXT(j) = BMAXT(NA1)
   BMINT(j) = BMINT(NA1)
   BDB0(j)  = BDB0(NA1)
   FOFB(j)  = FOFB(NA1)
enddo

call INTEGR_EF(RHO, 1, VR, VOLUM, NA1)

VOLUME = VOLUM(NA1)

return
end subroutine EQGUESS

!======================================================================|
subroutine SETEXM
!----------------------------------------------------------------------|
! Set external metrics 
! This subroutine is adjusted for external metric calculated for
! TJ-II   (Pereverzev 10.02.2005)
!----------------------------------------------------------------------|

use const_inc, only: RTOR, BTOR, ABC, ROC, HRO, HROX, &
   SHIFT, ELONG, TRIAN, VOLUME, GP, GP2, NA, NA1, NAB
use status_inc, only: SHIF, ELON, TRIA, SHX, ELX, TRX, &
   G11, G22, G33, G11X, G22X, G33X, GRADRO, DRODA, DRODAX, &
   IPOL, IPOLX, VR, VRS, VRX, RHO, XRHO, AMETR, SLAT, SLATX, &
   BDB0, BDB02, B0DB2, BMINT, BMAXT, FOFB, VOLUM, SHEAR, FP, MU
use debugger, only: markloc, debug
use parse_utils, only: ifdefx2

implicit none

integer :: j
double precision :: YNF, YR1
!----------------------------------------------------------------------|

call markloc('SETEXM', debug_lev=3*debug)

YNF = RTOR*GP2**2
do J=1, NA1
   if (IFDEFX2('SHX   ')) then
      SHIF(J) = SHX(j)
   else
      SHIF(J) = SHIFT
   endif
   if (IFDEFX2('ELX   ')) then
      ELON(J) = ELX(j)
   else
      ELON(J) = 1.
   endif
   if (IFDEFX2('TRX   ')) then
      TRIA(J) = TRX(j)
   else
      TRIA(J) = 0.
   endif
   if (IFDEFX2('G33X  ')) then
      G33(J) = G33X(j)
   else
      G33(J) = (RTOR/(RTOR + SHIFT))**2
   endif
   if (IFDEFX2('IPOLX ')) then
      IPOL(J) = IPOLX(j)
   else
      IPOL(J) = 1.
   endif
   if (IFDEFX2('VRX   ')) then
      VR(J) = VRX(j)
   else
      VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
   endif
   AMETR(J) = RHO(J)
enddo

!computes new roc

ROC = VR(NA1)/GP2**2 * G33(NA1)/RTOR
RHO(1: NA1) = XRHO(1: NA1)*ROC

if (debug > 0) then
   write(*, *) 'metric', RHO(1:10)
   write(*, *) VR(1:10)
   write(*, *) G33(1:10)
   write(*, *) IPOL(1:10)
endif
HRO = RHO(2) - RHO(1)
HROX = (RHO(2) - RHO(1))/ROC

!boundary
ELONG = ELON(NA1)
TRIAN = TRIA(NA1)
SHIFT = SHIF(NA1)

! Flux grid: j*h
do J=1, NA
   VRS(j) = 0.5*(VR(J+1) + VR(j))
   if (IFDEFX2('SLATX ')) then
      SLAT(J) = 0.5*(SLATX(J+1) + SLATX(j))
   else
      SLAT(J) = VRS(j) 
   endif
   if (IFDEFX2('G11X  ')) then
      G11(J) = 0.5*(G11X(j) + G11X(j+1))
   else
      G11(J) = VRS(j) 
   endif
   if (IFDEFX2('G22X  ')) then
      G22(J) = 0.5*(G22X(j) + G22X(j+1))
   else
      G22(J) = RTOR*VRS(j)/(GP2*(RTOR + SHIFT))**2
   endif
   if (IFDEFX2('DRODAX')) then
      DRODA(J) = 0.5*(DRODAX(j) + DRODAX(j+1))
   else
      DRODA(J) = 1.
   endif
enddo

! Linear extrapolation
SLAT(NA1)  = 1.5*SLAT(NA)  - 0.5*SLAT(NA-1)
VRS(NA1)   = 1.5*VRS(NA)   - 0.5*VRS(NA-1)
G11(NA1)   = 1.5*G11(NA)   - 0.5*G11(NA-1)
G22(NA1)   = 1.5*G22(NA)   - 0.5*G22(NA-1)
DRODA(NA1) = 1.5*DRODA(NA) - 0.5*DRODA(NA-1)
G11(NA1)   = 1.5*G11(NA)   - 0.5*G11(NA-1)

! Compute new minor radius
call INTEGR_EF(RHO(1: NA1), 1, 1./DRODA(1: NA1), AMETR(1: NA1), NA1)
ABC = AMETR(NA1)

do J=1, NA1
   if (j == 1) then
      SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
   elseif (j < NA) then
      SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
   else
      SHEAR(J) = HRO*(FP(NA1) - FP(NA))/HRO - FP(NA) + FP(NA-1)
      SHEAR(J) = 2.*SHEAR(J)*(2.*HRO - HRO)
      SHEAR(J) = SHEAR(J)/(2.*HRO*MU(NA1) - HRO*MU(NA))
   endif
   SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
   GRADRO(j) = DRODA(J)
   YR1 = RHO(j)*G22(j)*(MU(j)/RTOR)**2
   BDB02(j) = (1. + YR1)*G33(J)*IPOL(J)**2
   B0DB2(j) = ((RTOR + SHIF(J))/RTOR)**2 + 0.75*(AMETR(j)/RTOR)**2
   BMAXT(j) = BTOR*RTOR/(RTOR + SHIF(j) - AMETR(j))
   BMINT(j) = BTOR*RTOR/(RTOR + SHIF(j) + AMETR(j))
   BDB0(j) = (RTOR/(RTOR + SHIF(j)))
   YR1     = (RTOR + SHIF(j) - AMETR(j))/(RTOR + SHIF(j))
   FOFB(j) = 1. - sqrt(YR1)*(1. + 0.5*YR1)
enddo

if (NA1 < NAB) then

   do J=NA1+1, NAB
      SHIF(J) = SHIFT
      ELON(J) = 1.
      TRIA(J) = 0.
      G33(J) = (RTOR/(RTOR + SHIFT))**2
      IPOL(J) = 1.
      VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
      AMETR(J) = RHO(J)
      VRS(j) = 0.5*(VR(J+1) + VR(j))
      G11(J) = VRS(j)
      G22(J) = RTOR*VRS(j)/(RTOR + SHIFT)**2
      DRODA(J) = 1.
      SLAT(J) = VRS(J)*DRODA(J)
      SHEAR(j) = SHEAR(NA1)
      BDB02(j) = BDB02(NA1)
      B0DB2(j) = B0DB2(NA1)
      BMAXT(j) = BMAXT(NA1)
      BMINT(j) = BMINT(NA1)
      BDB0(j)  = BDB0(NA1)
      FOFB(j)  = FOFB(NA1)
   enddo

! Volume (on the shifted grid) is calculated using VR:
   call INTEGR_EF(RHO, 1, VR, VOLUM, NA1)

   call NEWGRD ! The RHO-grid and NA, NA1 are updated

   VOLUME = VOLUM(NA1)
endif

return
end subroutine SETEXM

!======================================================================|
! ROC3A [m]: Analytical formula for the "rho_tor" in vacuum
!     (Pereverzev 24.03.00)
double precision function ROC3A(Rmaj, shaf_shift, a_min, elongation, triangularity)
!----------------------------------------------------------------------|
! Input: Rmaj - Major radius [m]
!  shaf_shift - Shafranov shift [m]
!  a_min - Minor radius [m]
!  elongation - Elongation [d/l]
!  triangularity - Triangularity [d/l]
! Output:
!  ROC3A - Dimensional toroidal "rho" [m]
!----------------------------------------------------------------------|

use debugger, only: markloc, debug, astra_stop

implicit none

double precision, intent(in) :: Rmaj, shaf_shift, a_min, elongation, triangularity

double precision :: YD1, YT2, YGE, Y1, Y2, Y3, Y4, Y5, Y6

call markloc('ROC3A', debug_lev=2*debug)

YT2 = 2.*triangularity
YGE = a_min/(Rmaj + shaf_shift)
if (abs(YGE) > 1.d0) then
   write(*, *) " >>> Error >>> ROC3A >>> Illegal input: R+Delta < a"
   write(*, '(1P, 16X, 2(A, E10.3))') "R+Delta =", Rmaj + shaf_shift, ",   a =", a_min
   write(*, *) Rmaj, shaf_shift, a_min, elongation, triangularity
   call astra_stop
endif

if (YGE < 0) then
   write(*, *) " >>> Error >>> ROC3A: Illegal input"
   write(*, *) Rmaj, shaf_shift, a_min, elongation, triangularity
   call astra_stop
else if (YGE == 0) then
   ROC3A = 0.
else
   YD1 = 1.+YT2*(YT2-2./YGE)
   if (YD1 >= 0) then 
      YD1 = sqrt(YD1)
      Y1 = (2. - YT2*YGE)*YD1/(1. + YD1)
      Y2 = YGE - YT2
      Y3 = sqrt(1. - YT2*YGE + YGE*YD1)
      Y4 = sqrt(1. - YT2*YGE - YGE*YD1)
      Y5 = sqrt(1. + YGE)
      Y6 = sqrt(1. - YGE)
      Y1 = (Y1 + Y2)/(Y3 + (1. - YT2)*Y5) - (Y1 - Y2)/(Y4 + (1. + YT2)*Y6)
      Y1 = Y1/sqrt(2.*(1. - YT2*YGE + Y5*Y6))
   else
      Y5 = sqrt(1. + YGE)
      Y6 = sqrt(1. - YGE)
      Y1 = (1. - YT2)*Y5 + (1. + YT2)*Y6
      Y1 = (1. - Y1/sqrt(2.*(1. - YT2*YGE + Y5*Y6)))/YT2
   endif
   ROC3A = 2.*sqrt(Rmaj*a_min*elongation*Y1)
endif

return
end function ROC3A

!======================================================================|
! Module call sequence in A2EMEQ
!   
!   EQ -> mapping -> EMEQ -> NEWGRD -> mapping -> SETEDGE -> mapping
!
!                               | -> EQGB3
!                               |             | -> EQLVU3
!   inside EMEC:        EMEQ -> | -> EQAB3 -> | -> EQK3
!                               |             | -> EQC1
!                               | -> EQPPAB
!----------------------------------------------------------------------|
! This solver is called if LEQ(5) == 3 or 
!      if (NEQUIL != 42 and NEQUIL != 64.64 and NEQUIL > 0)
! EMEQ grid is defined as min(NA1, nint(NEQIL), NP)
! For (NEQUIL = 0 ) metric is prescribed by a simple formula
! For (NEQUIL = -1) metric is taken from a data file
! For (NEQUIL = 1 ) metric is frozen (can be used interactively)
!======================================================================|
subroutine A2EMEQ(jexit)

use parameter_inc, only: NRD
use emeq_mod, only: NP, emeq
use const_inc, only: HRO, ABC, ROC, RTOR, BTOR, IPL, & 
    TIME, ELONG, TRIAN, SHIFT, UPDWN, VOLUME, & 
    NA1, NA, NAB, NEQUIL, GP, GP2
use status_inc, only: TE, TI, CU, CUTOR, SHEAR, SHIV, & 
    RHO, AMETR, EQPF, EQFF, IPOL, MU, FP, SXHO, & 
    SLAT, VOLUM, SHIF, ELON, TRIA, DRODA, GRADRO, VR, VRS, XRHO, & 
    G11, G22, G33, G41, G42, G43, G44, G45, & 
    BDB0, BDB02, B0DB2, BMAXT, BMINT, FOFB
use debugger, only: markloc

implicit none

double precision, parameter :: ACEQLB=1.d-6

integer, intent(out) :: jexit

integer :: N3EQL, j, jp, jt, jna1, jcall, waitevent
double precision :: &
   ALFA, Y1, Y2, YY, YDA, GPP4, YRO, YCB, TRIABC, IINT, BTOOO
double precision, dimension(NRD) :: YA, YB, BA, BB, GR, GBD, GL, GSD, &
   A, B, C, D, BC, BD, XTR, BMODEQ, FOFBEQ, GRDAEQ, &
   XEQ, B2B0EQ, B0B2EQ, BMAXEQ, BMINEQ
character(len=80) :: STRI

external IINT, waitevent

save jcall, N3EQL
data jcall/0/

!----------------------------------------------------------------------|
! Prepare input data for the 3M equilibrium solver:

call markloc('A2EMEQ')

jexit = 0

if (jcall == -1) return

! Define N3EQL as min(NA1, nint(NEQIL), NP)
if (jcall == 0) then
   N3EQL = nint(NEQUIL)
   if (N3EQL > NP) then
      write(*, '(A, I4)') " >>> Warning >>> Maximum size of the equilibrium grid is", NP
      write(*, '(17X, A, I4)') "The grid will be reduced to NP =", NP
      N3EQL = NP
   endif
endif

! Rise triangularity (initiation stage only)
jt = 5
TRIABC = TRIAN*ABC

! Rise pressure (initiation stage only)
jp = 5
if (jcall <= jp) then
   Y1 = jcall   ! 0 <= jcall <= jp
   Y1 = Y1/jp
   do J=1, NA1
      A(J) = TE(J)
      B(J) = TI(J)
      C(J) = CU(J)
      TE(J) = Y1*A(J)
      TI(J) = Y1*B(J)
      if (jcall == 0) then
         CU(J) = (1 - (RHO(J)/ROC)**2)
      else
         CU(J) = Y1*C(J) + (1.d0 - Y1)*(1 - (RHO(J)/ROC)**2)
      endif
   enddo
   YCB = IINT(CU, ROC)
   do J=1, NA1
      CU(J) = CU(J)*IPL/YCB
   enddo
   call RHSEQ
   do J=1, NA1
      TE(J) = A(J)
      TI(J) = B(J)
      CU(J) = C(J)
   enddo
   jcall = jcall + 1
endif

YCB = RTOR/(RTOR + SHIFT)
do J=1, NA1
   XTR(J) = AMETR(J)/AMETR(NA1)
   B(J) = EQPF(J)/YCB
   A(J) = B(J) + YCB*EQFF(J)
enddo

do J=1, N3EQL
   XEQ(J) = (j - 1.)/(N3EQL - 1.)
enddo
! From transport grid in "a" to equidistant grid in "a"
ALFA = .001
call SMOOTH(ALFA, NA1, A, XTR, N3EQL, BA, XEQ)
call SMOOTH(ALFA, NA1, B, XTR, N3EQL, BB, XEQ)

if (TIME > 0.2453d16) then
   write(*, *) NA1, TIME, N3EQL, "    EQPF  EQFF  CU  CUTOR  A  B"
   write(*, '(3(2F10.5, 2X))') RTOR + SHIFT, ABC, ELONG, TRIABC, & 
       BTOR*RTOR/(RTOR + SHIFT), IPL
   do j=1, na1
      write(*, '(I4, 7F10.5)') j, A(j), B(j), CU(j), EQPF(j), EQFF(j), CU(j)
   enddo
   do j=1, N3EQL
      write(*, '(I4, 2F10.5)') j, BA(j), BB(j)
   enddo

   write(*, '(3(2F10.5, 2X))') (EQPF(j), j=NA1-5, NA1)
   write(*, '(3(2F10.5, 2X))') (EQFF(j), j=NA1-5, NA1)
   write(*, '(3(2F10.5, 2X))') (CU(j)/IPOL(j), j=NA1-5, NA1)
   write(*, '(3(2F10.5, 2X))') (CUTOR(j), j=NA1-5, NA1)
   write(*, '(3(2F10.5, 2X))') (2.*(CU(j)/IPOL(j) - CUTOR(j))/ &
      (CU(j)/IPOL(j) + CUTOR(j)), j=1, NA1)
   write(*, '(3(2F10.5, 2X))') ((EQPF(j) - CUTOR(j))* &
      RHO(j)*G22(J)*(MU(J)/RTOR)**2, j=NA1-5, NA1)
   write(*, '(3(2F10.5, 2X))') (BA(j), j=1, N3EQL)
   write(*, '(3(2F10.5, 2X))') (BB(j), j=1, N3EQL)

endif

call EMEQ &
! Input:
   (BA, BB, RTOR + SHIFT, ABC, ELONG, TRIABC, N3EQL, ACEQLB, &  ! relative accuracy
   BTOR*RTOR/(RTOR + SHIFT), IPL, &  ! Total plasma current
! output
   GR, GBD, GL, GSD, A, BD, B, BA, BB, BC, C, D, &
   B2B0EQ, B0B2EQ, BMAXEQ, BMINEQ, BMODEQ, FOFBEQ, GRDAEQ, &
! input
   TIME)
BTOOO = BTOR*RTOR/(RTOR + SHIFT)

!----------------------------------------------------------------------|
! Output: GR  - rho
!  GBD - shift
!  GL  - elongation
!  GSD - triangularity = \delta_Zakh = a*\delta_Astra
!  A   - <g^{11}>  = <[nabla(a)]^2>
!  BD  - <sqrt[g^{11}]> = <|nabla(a)|> flux surface area
!  B   - <g^{11}g^{33}> = <[nabla(a)/r]^2>
!  BA  - <g^{33}>  = <1/r^2>=G33/RTOR**2
!  BB  - IPOL*RTOR*BTOR = I
!  BC  - d\rho/da
!  C   - (dV/da)/(4\pi^2)
!  D   - V(a)
!               B2B0EQ - <B**2/B0**2>
!               B0B2EQ - <B0**2/B**2>
!               BMAXEQ - BMAXT
!               BMINEQ - BMINT
!               BMODEQ - <B/BTOR>
!               FOFBEQ - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
!               GRDAEQ - <grad a>
!----------------------------------------------------------------------|

call markloc('3-moment solver')

! N3EQL <= 1 can be returned by EQAB3 via EMEQ
j = 1
if (N3EQL > 10) then
!  Now check if GR(N3EQL) is a regular number
   write(STRI, *) GR(N3EQL)
   j = index(TRIM(STRI), "NaN") ! Not_a_Number
   if (j > 0) then !  Equilibrium crash
      j = 0
      jexit = 2
      return
   endif
endif

!----------------------------------------------------------------------|
! Define a new RHO-grid:
YRO = sqrt(RTOR/(RTOR + SHIFT))
ROC = YRO*GR(N3EQL)  ! Define a new RHO_edge
! FTN = GP*BTN*ROC*ROC
JNA1 = NA1 ! Save NA1  that can be changed by NEWGRD
YY = HRO ! Save HRO that can be changed by NEWGRD
call NEWGRD ! The RHO-grid and NA, NA1, HRO are updated
!----------------------------------------------------------------------|
! Define a new auxiliary (shifted) grid:
Y2 = 0.5d0/ROC
do J=1, NA1
   XTR(J) = SXHO(J)
enddo

GPP4 = GP2*GP2
do J=1, N3EQL
   DRODA(J) = YRO*BC(J)
   XEQ(J) = GR(J)/GR(N3EQL)
   G11(J) = A(J)*DRODA(J)**2
   G22(J) = B(J)*DRODA(J)**2
   G33(J) = BA(J)*RTOR*RTOR
   VRS(J) = GPP4*C(J)/DRODA(J)
   IPOL(J) = BB(J)/RTOR/BTOR
   GRADRO(J) = BD(J)*DRODA(J)
   VR(J) = VRS(j)
   YA(j) = G33(j)
   YB(j) = IPOL(j)
enddo
call QMAP(      N3EQL, XEQ, NA1, XTR, VRS)
call SMAP(ALFA, N3EQL, XEQ, NA1, XTR, G11)
call SMAP(ALFA, N3EQL, XEQ, NA1, XTR, G22)
call SMAP(ALFA, N3EQL, XEQ, NA1, XTR, YA) ! G33 @ aux. grid
call SMAP(ALFA, N3EQL, XEQ, NA1, XTR, YB) ! IPOL @ aux. grid
call SMAP(ALFA, N3EQL, XEQ, NA1, XTR, DRODA)
call SMAP(ALFA, N3EQL, XEQ, NA1, XTR, GRADRO)
! Multiply above quantities by linear in rho factors (i.e. f(0)=0)
do J=1, NA1
   G11(J) = G11(J)*VRS(j)
   G22(J) = G22(J)/YA(j)*(RTOR/YB(j))**2
   if (j < NA1) then
      G22(J) = G22(J)*XTR(j)*ROC
   else
      G22(J) = G22(J)*(NA*HRO + 0.5*HRO)
   endif
   SLAT(J) = GRADRO(J)*VRS(J)
enddo

if (NA1 /= JNA1) then
   call SETEDGE(JNA1, YY)
endif

!----------------------------------------------------------------------|
! Define the main transport grid:
do J = 1, NA1
   XTR(J) = XRHO(J) 
enddo
! Define VR, G33 and IPOL   on the main transport grid:
call QMAP  (      N3EQL, XEQ, NA1, XTR, VR)
call SMAP  (ALFA, N3EQL, XEQ, NA1, XTR, G33)
call SMAP  (ALFA, N3EQL, XEQ, NA1, XTR, IPOL)
call SMOOTH(ALFA, N3EQL, GBD, XEQ, NA1, SHIF, XTR)
call SMOOTH(ALFA, N3EQL, GL , XEQ, NA1, ELON, XTR)

! GL -> a
! GSD -> \delta^{ASTRA} (dimensionless)
YDA = ABC/(N3EQL - 1.)
do j=2, N3EQL
   A(j) = YDA*(J - 1.)
   B(J) = GSD(J)/A(J)
enddo
A(1) = 0.
B(1) = 0.
call TRANSF(N3EQL, B, XEQ, NA1,  TRIA, XTR)
call TRANSF(N3EQL, A, XEQ, NA1, AMETR, XTR)

do J=1, NA1
   SHIF(J) = SHIFT + SHIF(J)
   if (j == 1) then
      SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
   elseif (j < NA) then
      SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
   elseif (j == NA) then
      SHEAR(J) = HRO*(FP(j+1) - FP(j))/HRO - FP(j) + FP(j-1)
      SHEAR(J) = SHEAR(J)/(MU(j+1) + MU(j))
   endif
   SHEAR(J) = 1. - SHEAR(J)/(GP*BTOR*HRO**2)
enddo
SHEAR(NA1) = SHEAR(NA)

do J = 1, NAB
   SHIV(J) = UPDWN 
enddo

call TRANSF(N3EQL, B2B0EQ, XEQ, NA1, BDB02, XTR)
call TRANSF(N3EQL, B0B2EQ, XEQ, NA1, B0DB2, XTR)
call TRANSF(N3EQL, BMAXEQ, XEQ, NA1, BMAXT, XTR)
call TRANSF(N3EQL, BMINEQ, XEQ, NA1, BMINT, XTR)
call TRANSF(N3EQL, BMODEQ, XEQ, NA1, BDB0 , XTR)
call TRANSF(N3EQL, FOFBEQ, XEQ, NA1, FOFB , XTR)

BDB02 = BDB02*BTOOO**2./BTOR**2.
BDB0  = BDB0*BTOOO/BTOR
B0DB2 = B0DB2/BTOOO**2.*BTOR**2.

do J=NA1, NAB
   SHEAR(j) = SHEAR(NA1)
   BDB02(j) = BDB02(NA1)
   B0DB2(j) = B0DB2(NA1)
   BMAXT(j) = BMAXT(NA1)
   BMINT(j) = BMINT(NA1)
   BDB0(j)  = BDB0(NA1)
   FOFB(j)  = FOFB(NA1)
enddo

call INTEGR_EF(RHO, 1, VR, VOLUM, NA1)
VOLUME = VOLUM(NA1)

!Efable add G41-G45 : these are on the shifted grid.
do J=1, NA1
   G41(J) = 1.0
   G42(J) = GRADRO(J)
   G43(J) = GRADRO(J)
   G44(J) = G11(J)/VRS(J)
   G45(J) = G11(J)/VRS(J)
enddo

return
end subroutine A2EMEQ

!======================================================================|
subroutine A2GSSOLVER_EF(equil_solver)

use parameter_inc, only: NRD
use outcmn_inc, only: TASK, machine, DXLET, CCOIL, VCOIL, DUMCTP, DUMCT, NBNT
use const_inc
use status_inc
use plasma_state
use debugger, only: markloc
use parameters_a2spider, only: equil_now

implicit none

integer, parameter :: itfbe_ctrl=0

integer, intent(in) :: equil_solver

integer :: i, j, jneql, jnteta, jnbnd, jna1, jnstep, jstepp, j_save_bound
double precision :: rbnd(1000), zbnd(1000), tau_resistive, dampfacpsplex, &
    yrocnew, iplnew, ychipfp, dfpdrb12, yiplout, yipl, yupdwn
double precision, dimension(NRD) :: yg11, yg22, yg33, yvr, yvrs, yslat, yg41, &
   ygradro, yipol, ydroda, ypres, ybmaxt, ybmint, yfp, &
   ybdb02, ybdb0, yb0db2, yvolum, yametr, yshif, yelon, &
   ytria, yfofb, yeqpf, yeqff
double precision, dimension(NCNB) :: yccoil, yvcoil
double precision, dimension(NCTP) :: ydumctp, ydumct

save jnstep, j_save_bound, yiplout, iplnew

data jnstep/0/  
data j_save_bound/0/
data yiplout/0./

call markloc('A2GSSOLVER_EF')

jstepp = 10

if (jnstep == 0 .and. TASK(1:3) /= 'BGD') then
   call colovm(14) ! Iteration # in blue
   call textvm(62*DXLET, 2, "SPIDER iterations", 17)
   call redraw(0)
endif

jneql  = abs(nint(NEQUIL))
jnteta = abs(nint(MEQUIL))

! NBND=51 <-> ABC, ELONG, TRIAN setting for boundary

if (NBND == 0) NBND = 41 ! "NAMEXP BND" not found

! provide grid for t=TIME+TAU
if (j_save_bound == 0 .or. IPART == 1) then
   call BNDRY(rbnd(1:NBND), zbnd(1:NBND))
endif

if (nint(IPCTRL) == 50 .or. nint(IPCTRL) == 51) then
   IPCTRL = -50. + nint(IPCTRL)
   j_save_bound = 1
endif

jnbnd = NBND
if (j_save_bound == 1 .and. ipart /= 1) then
   rbnd(1: jnteta-1) = equil_now%coord_sys%position%r(jneql, 1: jnteta-1)
   zbnd(1: jnteta-1) = equil_now%coord_sys%position%z(jneql, 1: jnteta-1)
   jnbnd = jnteta - 1
endif

do j=1, NCNB
   yccoil(j) = CCOIL(j)
   yvcoil(j) = VCOIL(j)
enddo
do j=1, NCTP
   ydumctp(j) = DUMCTP(j)
   ydumct(j) = DUMCT(j)
enddo

tau_resistive = 0.01 ! this is just for AUG
dampfacpsplex = tau_resistive/TAU
dampfacpsplex = 0.

iplnew = (dampfacpsplex*iplnew + G22(NA)/RTOR/0.4/GP * &
   (FP(NA1) - FP(NA))/HRO * IPOL(NA1))/(1. + dampfacpsplex)

if (nint(inume3) == -1) iplnew = IPL !if fixed pprim ffprim, use IPL as total current for spider

if (ITFBP /= 0.) IPLFBE = iplnew     ! current for free boundary equilibrium
if (IPART == 1) then
   iplnew = IPL         ! if in initialization mode, use plasma current
   jstepp = 0           ! if in initialization mode, use plasma current
else
   jstepp = 1           ! if in initialization mode, use plasma current
endif

if (nint(inume3) == -1) then
   if (jstepp == 1) iplnew = yiplout ! if fixed ffp, pp, if in time advance, uses output current from spider if also keyplc = 0
endif

dfpdrb12 = (FP(NA1) - FP(NA))/HRO

if (plasma_up == 0 .or. plasma_trig == 1) then
   iplnew = IPL
   IPLFBE = IPL
endif

if (LEQ(4) <= 0) iplnew = IPL !if CU:AS, current is assigned from model file

yipl  = iplnew
jna1  = NA1

do j=1, na1
   yfp(j) = FP(j)
   yametr(j) = AMETR(j)

!Efable Reput values since now they are used for ff' computation
   yg11(j)   = G11(j)*VRS(j)                           !g11 = <(grad(V)^2)> 
   yg22(j)   = G22(j)*(GP2**2.0)*IPOL(j)/RTOR*VRS(j)   !g22 = <(grad(V)/R)^2>
   yg33(j)   = G33(j)/(RTOR**2.0)                      !g33 = <1/R^2>
   yvr(j)    = VR(j)
   yvrs(j)   = VRS(j)
   yslat(j)  = SLAT(j)
   ygradro(j)= GRADRO(j)*VRS(j)            !gradro = <(grad(V))> 
   yipol(j)  = IPOL(j)*RTOR*BTOR             !ipol = R*Bphi
   ypres(j)  = 1.60218E-3*((NE(j)*TE(j) + NI(j)*TI(j)) + &
   NB2EQL*(0.5*PBLON(j) + 0.5*PBPER(j)) + PFAST(j))  !thermal + fast ions  from NBI + fast alpha in keV/m^3 *1e19 to MJ/m^3
   yeqpf(j)  = EQPF(j)
   yeqff(j)  = EQFF(j)
   yvolum(j) = VOLUM(j)
   yshif(j)  = SHIF(j)
   yelon(j)  = ELON(j)
   ytria(j)  = TRIA(j)
enddo
ychipfp = FP(NA1)

! EFable call to SPIDER
i = 1   !fbe is off
if (IFBEY >= 1.) i = 2   !fbe is on
if (IPART == 1 ) i = 1   !fbe is off

if (ifbey > 0..and.plasma_up == 0) then
   jnstep = jnstep + 1
!This is to call SPIDER again...
   NEQUIL = -NEQUIL
   MEQUIL = -MEQUIL
   return
endif

! In COUPLING_SCHEME
call GS_SOLVER( &
! Input:
   equil_solver, &
   jneql, jnteta, jnbnd, jna1, &
   rbnd, zbnd, & 
   XRHO(1: jna1), RTOR, BTOR, ROC, yfp(1: jna1), ypres(1: jna1), &
   yeqpf(1: jna1), yeqff(1: jna1), &
   VOLUME, NCNB, yccoil(1: NCNB), yvcoil(1: NCNB), i, IPART, ITREQ, &
   nint(INUME3), TAU, nint(ITFBP), nint(ICIRCQ), nint(IPCTRL), nint(IFBEY), &
   TIME, ychipfp, machine, &
   PSIFB, PSIEXT, PSPLEX, &
! Output: 
   yrocnew, yipl, yg11(1: jna1), yg41(1: jna1), yg22(1: jna1), &
   yg33(1: jna1), G22E(1: jneql), G33E(1: jneql), yvr(1: jna1), yvrs(1: jna1), &
   yslat(1: jna1), ygradro(1: jna1), yipol(1: jna1), ybmaxt(1: jna1), ybmint(1: jna1), &
   ybdb02(1: jna1), ybdb0(1: jna1), yb0db2(1: jna1), ydroda(1: jna1), yvolum(1: jna1), &
   yametr(1: jna1), yupdwn, yshif(1: jna1), yelon(1: jna1), ytria(1: jna1), &
   yfofb(1: jna1), AREAT(1: jna1), PERIM(1: jna1) ) 

ROC = YROCNEW  ! Define a new RHO_edge
JNA1 = NA1     ! Save NA1  that can be changed by NEWGRD
yiplout = yipl ! new current in case

!Update control quantities
if (nint(IPCTRL) /= 0 .and. nint(IPCTRL) > -2) then
   do j=1, NCTP
      DUMCTP(j)=ydumctp(j) 
   enddo
endif   

do j=1, na1

! Meaningful lines (change evolutuion):
! Enabling G22 leads to a divergence
!Efable try removing G22, original has it, now put it back
   VR(j)    = yvr(j)
! Auxiliary lines (do not change evolutuion):
   VRS(j)   = yvrs(j)
   SLAT(j)  = yslat(j)
   BMAXT(j) = ybmaxt(j)
   BMINT(j) = ybmint(j)
   BDB02(j) = ybdb02(j)
   BDB0(j)  = ybdb0(j)
   B0DB2(j) = yb0db2(j)
   DRODA(j) = ydroda(j)
   IPOL(j)  = yipol(j)
   G11(j)   = yg11(j)
   G22(j)   = yg22(j)
   G33(j)   = yg33(j)
   GRADRO(j)= ygradro(j)
   VOLUM(j) = yvolum(j)
   AMETR(J) = yametr(J)
   SHIF(J)  = yshif(J)
   ELON(J)  = yelon(J)
   TRIA(J)  = ytria(J)
   G41(J)   = yg41(J)
   G42(J)   = GRADRO(J)
   G43(J)   = GRADRO(J)
   G44(J)   = G11(J)/VRS(J)
   G45(J)   = G11(J)/VRS(J)
   FOFB(J)  = yfofb(J)
   FP(J)    = yfp(J)     ! due to adiabatic compression done in the code
   EQPF(J)  = yeqpf(J)   ! due to adiabatic compression done in the code
   EQFF(J)  = yeqff(J)   ! due to adiabatic compression done in the code
enddo

if (NBNT > 0 .or. TIME > ITFBE) then
   UPDWN = yupdwn
   ABC   = yametr(NA1) 
   ELONG = ELON(NA1)
   TRIAN = TRIA(NA1)
   SHIFT = SHIF(NA1)
endif

if (itfbe_ctrl > 0) then
   UPDWN = yupdwn
   ABC   = yametr(NA1) 
   ELONG = ELON(NA1)
   TRIAN = TRIA(NA1)
   SHIFT = SHIF(NA1)
endif    

! Deallocate equil_out%metric_coefs%g1 & co

call NEWGRD ! The RHO-grid and NA, NA1, HRO are updated

VOLUM(NA1) = yvolum(NA1)

G22 = G22/VRS*RTOR/(GP2**2.0)/IPOL
G11 = G11/VRS
GRADRO = GRADRO/VRS
DRODA = DRODA/VRS

do J=1, NA1
   if (j == 1) then
      SHEAR(J) = (FP(2) - FP(1))/(2.*MU(1) + 0.333*(MU(1) - MU(2)))
   elseif (j < NA) then
      SHEAR(J) = (FP(j+1) - 2.*FP(j) + FP(j-1))/(MU(j+1) + MU(j))
   elseif (j == NA) then
      SHEAR(J) = HRO*(FP(j+1) - FP(j))/HRO-FP(j) + FP(j-1)
      SHEAR(J) = SHEAR(J)/(MU(j+1) + MU(j))
   endif
   SHEAR(J) = 1. -SHEAR(J)/(GP*BTOR*HRO**2)
enddo
SHEAR(NA1) = SHEAR(NA)

do J = 1, NAB
   SHIV(J) = UPDWN 
enddo

do J=NA1, NAB
   SHEAR(j) = SHEAR(NA1)
   BDB02(j) = BDB02(NA1)
   B0DB2(j) = B0DB2(NA1)
   BMAXT(j) = BMAXT(NA1)
   BMINT(j) = BMINT(NA1)
   BDB0(j)  = BDB0(NA1)
   FOFB(j)  = FOFB(NA1)
enddo

VOLUME = VOLUM(NA1)

jnstep = jnstep + 1 ! Count SPIDER calls, this was outside enddo

!This is to call SPIDER again...
NEQUIL = -NEQUIL
MEQUIL = -MEQUIL

return
end subroutine A2GSSOLVER_EF

!======================================================================|
subroutine BNDRY(RPB, ZPB)
!----------------------------------------------------------------------|
! If 3M solver is used the subroutine is not called.
! Otherwise, if a general equilibrium solver, ESC or SPIDER, is called
! then 
! 1) In case of the plasma boundary defined by 3 moments, 
!    this subroutine writes 8 points on the boundary into arrays BNDR, BNDZ
!    and into arrays RPB(1:NBND), ZPB(1:NBND)
! 2) If the plasma boundary is defined by a data file then
!    this subroutine uses the arrays BNDR, BNDZ as an input and
!    produces output in [time dependent] arrays RPB, ZPB
!----------------------------------------------------------------------|
! NBND     number of points on the plasma vacuum boundary
! NBNT     number of times for the plasma boundary evolution
!  call from ESC:
!  call BNDRY(RPB, ZPB)
!  call from SPIDER:
!  call BNDRY(RZPB, RZPB(NBND+1))
!----------------------------------------------------------------------|

use outcmn_inc, only: NBNT, BNDTIM, BNDR, BNDZ
use const_inc, only: NBND, GP2, TIME, RTOR, SHIFT, ABC, TRIAN, UPDWN, ELONG

implicit none

double precision, intent(out) :: RPB(*), ZPB(*)

integer :: j, j1, jt
double precision :: ydt, yd1, yd2, yfi

if (NBNT <= 1) then

   if (NBNT == 0) then  ! No input group "NAMEXP BND" found

      if (NBND == 0) NBND = 8   ! call from ESC

      if (NBND /= 8) then
         do j=1, NBND
            YFI = GP2*(j-1.)/(NBND)
            YD1 = sin(YFI)
            ZPB(j) = UPDWN+ABC*ELONG*YD1
            RPB(j) = RTOR+SHIFT+ABC*(cos(YFI)-TRIAN*YD1*YD1)
            BNDZ(j) = ZPB(j)
            BNDR(j) = RPB(j)
         enddo
         return
      endif

      yd1 = .75  ! sin^2(pi/3)
      yd2 = .5  ! cos(pi/3) 
      ydt = sqrt(yd1)  ! sin(pi/3)
      BNDR(1) = RTOR + SHIFT - ABC*TRIAN
      BNDZ(1) = UPDWN + ABC*ELONG
      BNDR(2) = RTOR + SHIFT - ABC*TRIAN
      BNDZ(2) = UPDWN - ABC*ELONG
      BNDR(3) = RTOR + SHIFT - ABC
      BNDZ(3) = UPDWN
      BNDR(4) = RTOR + SHIFT + ABC
      BNDZ(4) = UPDWN
      BNDR(5) = RTOR + SHIFT - ABC*(TRIAN*yd1 + yd2)
      BNDZ(5) = UPDWN + ABC*ELONG*ydt
      BNDR(6) = RTOR + SHIFT - ABC*(TRIAN*yd1 - yd2)
      BNDZ(6) = UPDWN + ABC*ELONG*ydt
      BNDR(7) = RTOR + SHIFT - ABC*(TRIAN*yd1 - yd2)
      BNDZ(7) = UPDWN - ABC*ELONG*ydt
      BNDR(8) = RTOR + SHIFT - ABC*(TRIAN*yd1 + yd2)
      BNDZ(8) = UPDWN - ABC*ELONG*ydt
   endif
   do j=1, NBND
      RPB(j) = BNDR(j)
      ZPB(j) = BNDZ(j)
   enddo
   return
endif

if (TIME <= BNDTIM(1)) then ! Take bnd at time=t1
   do j=1, NBND
      RPB(j) = BNDR(1 + (j-1)*NBNT)
      ZPB(j) = BNDZ(1 + (j-1)*NBNT)
   enddo
   return
endif
if (TIME >= BNDTIM(NBNT)) then ! Take bnd at time=t_NBNT
   do j=1, NBND
      RPB(j) = BNDR(NBNT + (j-1)*NBNT)
      ZPB(j) = BNDZ(NBNT + (j-1)*NBNT)
   enddo
   return
endif

do j=1, NBNT
   if (TIME > BNDTIM(j)) jt = j
enddo
if (jt == NBNT) write(*, *) "OGOGO"
ydt = BNDTIM(jt+1) - BNDTIM(jt)
yd1 = (TIME - BNDTIM(jt))/ydt
yd2 = (TIME - BNDTIM(jt+1))/ydt
do j=1, NBND
   RPB(j) = BNDR(jt + (j-1)*NBNT)
   ZPB(j) = BNDZ(jt + (j-1)*NBNT)
enddo
if (NBND > 12) return

!----------------------------------------------------------------------|
! The order is essential
!  Top(1), bottom(2), inward(3), outward(4), ...
YDT = -1.
do j=1, NBND
   if (ZPB(j) > YDT) then
      YDT = ZPB(j)
      j1 = j
   endif
enddo
YD1 = RPB(1)
YD2 = ZPB(1)
RPB(1) = RPB(j1)
ZPB(1) = ZPB(j1)
RPB(j1) = YD1
ZPB(j1) = YD2
! Bottom(2)
YDT = 1.
do j=2, NBND
   if (ZPB(j) < YDT) then
      YDT = ZPB(j)
      j1 = j
   endif
enddo
YD1 = RPB(2)
YD2 = ZPB(2)
RPB(2) = RPB(j1)
ZPB(2) = ZPB(j1)
RPB(j1) = YD1
ZPB(j1) = YD2
! Inward(3)
YDT = 1000.
do j=3, NBND
   if (RPB(j) < YDT) then
      YDT = RPB(j)
      j1 = j
   endif
enddo
YD1 = RPB(3)
YD2 = ZPB(3)
RPB(3) = RPB(j1)
ZPB(3) = ZPB(j1)
RPB(j1) = YD1
ZPB(j1) = YD2
! Outward(4)
YDT = -1.
do j=4, NBND
   if (RPB(j) > YDT) then
      YDT = RPB(j)
      j1 = j
   endif
enddo
YD1 = RPB(4)
YD2 = ZPB(4)
RPB(4) = RPB(j1)
ZPB(4) = ZPB(j1)
RPB(j1) = YD1
ZPB(j1) = YD2

return
end subroutine BNDRY

!======================================================================|
subroutine get_coil(tim_in, coil_arr, nt, n_coil, coil_curr)
!----------------------------------------------------------------------|
! This routine get the control quantities from the exp data at the present time
! slice.
!----------------------------------------------------------------------|

implicit none

integer, intent(in) :: nt, n_coil
double precision, intent(in) :: tim_in
double precision, intent(in), dimension(*) :: coil_arr
double precision, intent(out), dimension(n_coil) :: coil_curr

integer :: j, j1, j2, jt
double precision :: ydt, yd1, yd2

if (nt <= 1) then 
   if (nt == 0) then  
      coil_curr = 0.0
      return
   endif  
endif

! Input order:
! t_1  t_2  t_3
! coil currents(t_1)  
! coil currents(t_2) 
if (tim_in <= coil_arr(1)) then
   j1 = nt
   do j=1, n_coil
      coil_curr(j) = coil_arr(j1+j)
   enddo
   return
endif
if (tim_in >= coil_arr(nt)) then
   j1 = n_coil*(nt - 1) + nt
   do j=1, n_coil
      coil_curr(j) = coil_arr(j1+j)
   enddo
   return
endif

do j=1, nt
   if (coil_arr(j) <= tim_in) jt = j
enddo
ydt = coil_arr(jt+1) - coil_arr(jt)
yd1 = (tim_in - coil_arr(jt))/ydt  ! > 0
yd2 = (tim_in - coil_arr(jt+1))/ydt  ! < 0
do j=1, n_coil
   j1 = n_coil*(jt-1) + nt + j
   j2 = n_coil*jt + nt + j
   coil_curr(j) = yd1*coil_arr(j2) - yd2*coil_arr(j1)
enddo

return
end subroutine get_coil

!======================================================================|
subroutine GETCOILS(yvcoil, yccoil)
!----------------------------------------------------------------------|
! This routine get the coil currents from the exp data at the present time
! slice.
!----------------------------------------------------------------------|

use outcmn_inc, only: CCOIL, VCOIL, NCNBT, CCOILX, VCOILX
use const_inc, only: IPSMK, TIME, NCNB

implicit none

double precision, intent(out), dimension(NCNB) :: yvcoil, yccoil

if (nint(IPSMK) >= 1) then
   yccoil = CCOIL(1: NCNB)
   yvcoil = VCOIL(1: NCNB)
return
endif

if (NCNBT <= 1) then ! Just one time point
   if (NCNBT == 0) then
      yccoil = CCOIL(1: NCNB)
      yvcoil = VCOIL(1: NCNB)
      return
   endif  
endif

call get_coil(TIME, CCOILX, NCNBT, NCNB, yccoil)
call get_coil(TIME, VCOILX, NCNBT, NCNB, yvcoil)

return
end

!======================================================================|
subroutine GETDUMCT(ydumct)
!----------------------------------------------------------------------|
! This routine get the control quantities from the exp data at the present time
! slice.
!----------------------------------------------------------------------|

use outcmn_inc, only: DUMCTX, NCTPT
use const_inc, only: TIME, NCTP

implicit none

double precision, intent(out), dimension(NCTP) :: ydumct

call get_coil(TIME, DUMCTX, NCTPT, NCTP, ydumct)

return
end subroutine GETDUMCT

!======================================================================|
subroutine GETCTRLMS(yctrlm)   
!----------------------------------------------------------------------|
! This routine get the control quantities from the exp data at the present time
! slice.
!----------------------------------------------------------------------|

use outcmn_inc, only: NCRMT, CTRLMX
use const_inc, only: NCRM, TIME

implicit none

double precision, intent(out), dimension(NCRM) :: yctrlm

call get_coil(TIME, CTRLMX, NCRMT, NCRM, yctrlm)

return
end subroutine GETCTRLMS

!======================================================================|
subroutine RHSEQ
!-----------------------------------------------------------------------
! Input: RTOR, BTOR, NA, NA1, HRO, NB2EQL, 
!  NE, NI, TE, TI, MU, CU, AMETR, RHO, PBLON, PBPER, G22, G33, IPOL
! Output:
!       EQPF
!       EQFF
!       CUTOR
! Both quantities EQPF (~p') and EQFF (~II') are given in [MA/m^2]
! EQPF = -1.6E-3*(2*\pi*R_0)\prti{n_13*T_keV}{\psi[Vs=T*m^2]}
!      = -1.E-6*/(2*\pi*R_0)\prti{p[J/m^3=Pascal]}{\psi[Vs]}
! EQFF = -1.E-6*2*\pi/(R_0*\mu_0)*I*\prti{I}{\psi}
!      = -5./R_0*I*\prti{I}{\psi}
! Local toroidal current density j[MA/m^2] is EQPF*r/R_0+EQFF*R_0/r, i.e.
!       j(r, z) = r*(\vec j\cdot\nabla\zeta) = EQPF*r/R_0+EQFF*R_0/r , 
! ASTRA average toroidal current density is
!   R_0*<\vec j\cdot\nabla\zeta> = EQPF+EQFF*<R_0^2/r^2>
!-----------------------------------------------------------------------

use const_inc, only: INUME3, RTOR, BTOR, HRO, NA, NA1, NB2EQL
use status_inc
use debugger, only: markloc, debug

implicit none

integer :: j
double precision :: YCB, YG, YTH2

call markloc('RHSEQ', debug_lev=3*debug)

! Preparing input for the 3M equilibrium solver:
if (nint(INUME3) >= 0) then    ! if inume3 < 0 , uses eqpf, eqff from model file (user-prescribed)
   YCB = 1.6E-3*RTOR/(BTOR*HRO*HRO)
   do J=2, NA
      if (j == NA) YCB = YCB*HRO/HRO
      EQFF(J) = ( (NE(J+1)*TE(J+1) - NE(J)*TE(J)) + &
                  (NI(J+1)*TI(J+1) - NI(J)*TI(J)) )/J
      EQFF(J) = EQFF(J) + 0.5*NB2EQL * &
         (PBLON(J+1) - PBLON(J) + PBPER(J+1) - PBPER(J))/J
      EQFF(J) = EQFF(J) + (PFAST(J+1) - PFAST(J))/J
      EQFF(J) = -YCB*EQFF(J)/(MU(J))
   enddo
   EQFF(1) = EQFF(2)
   EQFF(NA1) = EQFF(NA) + (EQFF(NA) - EQFF(NA-1)) * &
      (AMETR(NA) - AMETR(NA-1))/(AMETR(NA1) - AMETR(NA))
   EQFF(NA1) = EQFF(NA) + (EQFF(NA) - EQFF(NA-1))/HRO*HRO
   do J=1, NA1
      EQPF(j) = EQFF(j)
      YTH2 = RHO(j)*G22(J)*(MU(J)/RTOR)**2
      YG = (1. + YTH2)*G33(J)
      EQFF(J) = (CU(J)/IPOL(J) - EQPF(J))/YG
      CUTOR(J) = (CU(J)/IPOL(J) + YTH2*EQPF(J))/(1. + YTH2)
   enddo
endif

return
end subroutine RHSEQ

!----------------------------------------------------------------------|
subroutine CUOFMU
!----------------------------------------------------------------------|
! Compute CU(rho) and FP(rho) from MU(rho)
!----------------------------------------------------------------------|
! Input: RHO - radial grid (m)
!        NA1 - number of grid points
!        GP2 - 2\pi
!        RTOR - R_0 [m]
!        BTOR - B_0 [T]
!        G22(1:NA1-1) - <g22/g>*............
!        G33(1:NA1-1) - 
!        IPOL(1:NA1-1)- 
!        CD(1:NA1-1)  - external (+bootstrap) current
!        MU(1:NA1)    - (1/rho)dF/d(rho) rotational transform
! Output: CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
!         FP(1:NA1) - poloidal flux [Vs]
!----------------------------------------------------------------------|

use const_inc, only: GP, GP2, RTOR, BTOR, NA1, NA, PSIBO, PSIAX, HRO
use status_inc, only: RHO, SRHO, XRHO, CU, MU, FP, G22, G33, IPOL

implicit none

integer	:: j
double precision :: YH, YM, YM1, YM2, YC, YF, HH, YAJ, YCJ, MUVAC
double precision, dimension(NA1) :: YAR
	
YC = 0.2*GP2*RTOR/BTOR
YH = RHO(2) - RHO(1)
YF = GP2*YH*YH*BTOR
YM = 0.
do J=1, NA1
   YM1 = MU(j)*j
   YM2 = YM1*G22(j)
   if (j < NA) then
      CU(j) = (YM2 - YM)*G33(j)*IPOL(j)**3/YC/RHO(j)
      FP(j+1) = FP(j) + YF*YM1
   elseif (j == NA) then
      CU(j) = (YM2 - YM)*G33(j)*IPOL(j)**3/YC/RHO(j)
      FP(j+1) = FP(j) + GP2*BTOR*YM1*YH*(RHO(NA1) - RHO(NA))
   else
! Be careful because of instability in the loop: j -> FF' -> G22 -> j ->
      CU(NA1) = CU(NA)
   endif
   YM = YM2
enddo
MUVAC = 4.*GP*0.1

do j=1, NA1
   YAR(j) = GP2*BTOR*MU(j)*SRHO(j)
enddo

call INTEGR_EF(SRHO(1:NA1), 2, YAR(1:NA1), FP(1:NA1), NA1)

HH = HRO*HRO
YAJ = 0.

do J=1, NA
   YCJ = YAJ
   if (j < NA) then
      YAJ = (FP(j+1) - FP(j))/HH
   else
      YAJ = (FP(j+1) - FP(j))/HRO**2
   endif
   YAJ = G22(j)*YAJ
   CU(j) = (YAJ - YCJ)/(HRO*(j - 0.5))
enddo

call EXTRAP_EF(XRHO(1:NA), CU(1:NA), XRHO(NA1), NA, CU(NA1), 2, NA)
YCJ = 1.25/(GP*GP*RTOR)
YAJ = 0.5/(GP*BTOR)
do J=1, NA1
   CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
enddo

PSIBO = FP(NA1)
call EXTRAP_EF(XRHO(1:NA1), FP(1:NA1), 0.0, 1, PSIAX, 1, NA1)

return
end subroutine CUOFMU

!======================================================================|
subroutine CUOFP
!----------------------------------------------------------------------|
! Compute CU(rho) and MU(rho) from FP(rho)
!----------------------------------------------------------------------|
! Input: HRO - radial step (m)
!  NA1 - number of grid points
!  G22(1:NA) - <g22/g>*............
!  G33(1:NA) - 
!  IPOL(1:NA) - 
!  CD(1:NA) - external (+bootstrap) current
!  FP(1:NA1) - poloidal flux
! Output: CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
!  MU(1:NA1) - (1/rho)dF/d(rho)     rotational transform
!----------------------------------------------------------------------|

use parameter_inc
use status_inc, only: RHO, SRHO, XRHO, FP, MU, CU, IPOL, G22, G33
use const_inc, only: GP, GP2, RTOR, HRO, BTOR, NA, NA1, PSIBO, PSIAX

implicit none

integer :: j
double precision, dimension(NA1) :: YAR, YAR1
double precision :: MUVAC, HH, YAJ, YCJ, ARRNA1

HH = HRO*HRO
YAJ = 0.

do J=1, NA
   YCJ = YAJ
   if (j < NA) then
      YAJ = (FP(j+1) - FP(j))/HH
      MU(j) = YAJ/j
      YAJ = G22(j)*YAJ
      CU(j) = (YAJ - YCJ)/HRO
   else
      YAJ = (FP(j+1) - FP(j))/HRO**2
      MU(j) = YAJ/j
      YAJ = G22(j)*YAJ
      CU(j) = (YAJ - YCJ)/HRO
   endif
   CU(j) = CU(j)/(j - 0.5)
enddo

call EXTRAP_EF(XRHO(1: NA), CU(1: NA), XRHO(NA1), NA, CU(NA1), 2, NA)

YCJ = 1.25/(GP*GP*RTOR)
YAJ = 0.5/(GP*BTOR)
do J=1, NA1
   CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
   MU(J) = YAJ*MU(j)
enddo

MUVAC = 4*GP*0.1

! Compute PSIAX, PSIBO
PSIBO = FP(NA1)
call EXTRAP_EF(XRHO(1: NA1), FP(1: NA1), 0.0, 1, PSIAX, 1, NA1)

! Compute MU from FP, MU is on shifted grid
do j=1, NA1
   YAR(j)=RHO(j)
enddo
    
call DERIV_EF(YAR(1: NA1), SRHO(1: NA1), 1, FP(1:NA1), YAR1(1:NA1), 1, NA1, 1)

do j=1, NA1
   MU(j) = YAR1(j)/(GP2*BTOR*SRHO(j))
enddo

HH = HRO*HRO
YAJ = 0.
CU = 0.
MU = 0.

do J=1, NA
   YCJ = YAJ
   if (j < NA) then
      YAJ = (FP(j+1) - FP(j))/HH
      MU(j) = YAJ/j
      YAJ = G22(j)*YAJ
      CU(j) = (YAJ - YCJ)/HRO
   else
      YAJ = (FP(j+1) - FP(j))/HRO**2
      MU(j) = YAJ/j
      YAJ = G22(j)*YAJ
      CU(j) = (YAJ - YCJ)/HRO
   endif
   CU(j) = CU(j)/(j-0.5)
enddo
MU(NA1) = ARRNA1(MU(NA), 1.)  ! See DEFARR
YCJ = 1.25/(GP*GP*RTOR)
YAJ = 0.5/(GP*BTOR)
do J=1, NA1
   CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
   MU(J) = YAJ*MU(j)
enddo

end subroutine CUOFP

!======================================================================|
subroutine NEWGRD
!----------------------------------------------------------------------|
! input:  XRHO, SXHO, HROX, ROC, NA1, AB, ABC, AMETR(NA1)
! Output: NB1, RHO, SRHO, HRO, AMETR(j>NA1)
!----------------------------------------------------------------------|
!----------------------------------------------------------------------|

use parameter_inc, only: NRD
use status_inc, only: RHO, XRHO, SRHO, SXHO, AMETR
use const_inc, only: HRO, HROX, AB, ABC, ROC, ROB, ROWALL, &
   FTO, BTN, GP, NA, NA1, NB1, NAB

implicit none

integer :: j
double precision :: YDA

HRO = HROX*ROC
do j=1, NRD
  RHO(j)  = XRHO(j)*ROC
  SRHO(j) = SXHO(j)*ROC
enddo
FTO = GP*BTN*ROC**2.0
!Efable normalized grid stuff

do j=1, NRD
   if (RHO(j) >= ROWALL) EXIT
ENDDO
NB1 = j

AMETR(NA1) = ABC

NAB = NA1
ROB = RHO(NAB)

if ( 2.*abs(AB - ABC) < AMETR(NA1) - AMETR(NA) ) then
   AB = ABC
   return
endif

NAB = NA1/ABC*AB
ROB = RHO(NAB)
AMETR(NAB) = AB

if (NA1 + 1 > NB1 .or. NA1 == NAB) return

YDA = (AB - ABC)/(NAB - NA1)
do j=NA1 + 1, NB1
   AMETR(j) = ABC + YDA*(j - NA1)
enddo

if (NA1 < NB1) then
   do j=NA1 + 1, NB1
      if (AMETR(j) < AB) NAB = j
   enddo
   if (NAB < NB1) NAB = NAB + 1
endif
ROB = RHO(NAB)
AMETR(NAB) = AB

return
end subroutine NEWGRD

!======================================================================|
subroutine SETEDGE(JNA1O, YHAO)
!----------------------------------------------------------------------|
! Input:
! YHAO - old edge-cell size
! JNA1O - old grid dimension
! Output:
! if (NA1 == JNA1O) then FP(NA1) 
!        all main arrays(NA1) otherwise
! Note! FP defined in this subroutine will be assigned to FPO in OLDNEW
!----------------------------------------------------------------------|

use const_inc, only: GP, GP2, NA, NA1, HRO, ROC, RTOR, BTOR, IPL
use status_inc

implicit none

integer, intent(in) :: JNA1O
double precision, intent(in) :: YHAO
integer :: JNAO, j
double precision :: YL, YR, YD, YM, YM1, YM2

if (NA1 > JNA1O) then ! NA1 increases: JNA1O >= NA
   JNAO = JNA1O - 1
   TE(NA1)    = TE(JNA1O)
   TI(NA1)    = TI(JNA1O)
   NE(NA1)    = NE(JNA1O)
   NI(NA1)    = NI(JNA1O)
   NIZ1(NA1)  = NIZ1(JNA1O)
   NIZ2(NA1)  = NIZ2(JNA1O)
   NIZ3(NA1)  = NIZ3(JNA1O)
   NALF(NA1)  = NALF(JNA1O)
   NHE3(NA1)  = NHE3(JNA1O)
   NHYDR(NA1) = NHYDR(JNA1O)
   NDEUT(NA1) = NDEUT(JNA1O)
   NTRIT(NA1) = NTRIT(JNA1O)
   NMAIN(NA1) = NMAIN(JNA1O)
   F0(NA1) = F0(JNA1O)
   F1(NA1) = F1(JNA1O)
   F2(NA1) = F2(JNA1O)
   F3(NA1) = F3(JNA1O)
   F4(NA1) = F4(JNA1O)
   F5(NA1) = F5(JNA1O)
   F6(NA1) = F6(JNA1O)
   F7(NA1) = F7(JNA1O)
   F8(NA1) = F8(JNA1O)
   F9(NA1) = F9(JNA1O)

   YM  = MU(JNAO)
   YM1 = MU(JNAO-1)
   YM2 = MU(JNAO-2)
   YR = YHAO/HRO
   YM = YM*8./(1. + YR)/(3. + YR) - YM2*(YR - 1.)/(3. + YR) + &
        YM1*2.*(YR - 1.)/(1. + YR)
   YD = GP2*BTOR*JNAO*HRO*YM

! find MU_old
   YM  = (FPO(JNA1O)    - FPO(JNAO)    )/(YHAO*JNAO*HRO)
   YM1 = (FPO(JNAO)     - FPO(JNAO - 1))/(HRO*(JNAO - 1)*HRO)
   YM2 = (FPO(JNAO - 1) - FPO(JNAO - 2))/(HRO*(JNAO - 2)*HRO)
   YR = YHAO/HRO
   YM = YM*8./(1. + YR)/(3. + YR) - YM2*(YR - 1.)/(3. + YR) + &
        YM1*2.*(YR - 1.)/(1. + YR)
   YD = JNAO*HRO*YM
! Extend FPO to the new grid (points NA and NA1)
   FPO(NA)  = FPO(NA-1) + (FV(NA) - FV(NA-1)) + YD*HRO
   FPO(NA1) = FPO(NA) + IPL*HRO*0.4*GP*RTOR/G22(NA)/IPOL(NA1) + &
              FV(NA) - FV(NA-1)
   YD = RHO(NA1) - RHO(JNAO)
   do j=JNA1O, NA ! Linear interpolation on the 2nd last node
      YR = (RHO(j) - RHO(JNAO))/YD
      YL = (RHO(NA1) - RHO(j))/YD
      TE(j)    = TE(NA1)   *YR + TE(JNAO)   *YL
      TI(j)    = TI(NA1)   *YR + TI(JNAO)   *YL
      NE(j)    = NE(NA1)   *YR + NE(JNAO)   *YL
      F0(j)    = F0(NA1)   *YR + F0(JNAO)   *YL
      F1(j)    = F1(NA1)   *YR + F1(JNAO)   *YL
      F2(j)    = F2(NA1)   *YR + F2(JNAO)   *YL
      F3(j)    = F3(NA1)   *YR + F3(JNAO)   *YL
      F4(j)    = F4(NA1)   *YR + F4(JNAO)   *YL
      F5(j)    = F5(NA1)   *YR + F5(JNAO)   *YL
      F6(j)    = F6(NA1)   *YR + F6(JNAO)   *YL
      F7(j)    = F7(NA1)   *YR + F7(JNAO)   *YL
      F8(j)    = F8(NA1)   *YR + F8(JNAO)   *YL
      F9(j)    = F9(NA1)   *YR + F9(JNAO)   *YL
      NI(j)    = NI(NA1)   *YR + NI(JNAO)   *YL
      VRO(j)   = VRO(NA1)  *YR + VRO(JNAO)  *YL
      NIZ1(j)  = NIZ1(NA1) *YR + NIZ1(JNAO) *YL
      NIZ2(j)  = NIZ2(NA1) *YR + NIZ2(JNAO) *YL
      NIZ3(j)  = NIZ3(NA1) *YR + NIZ3(JNAO) *YL
      NALF(j)  = NALF(NA1) *YR + NALF(JNAO) *YL
      NHE3(j)  = NHE3(NA1) *YR + NHE3(JNAO) *YL
      NHYDR(j) = NHYDR(NA1)*YR + NHYDR(JNAO)*YL
      NDEUT(j) = NDEUT(NA1)*YR + NDEUT(JNAO)*YL
      NTRIT(j) = NTRIT(NA1)*YR + NTRIT(JNAO)*YL
      NMAIN(j) = NMAIN(NA1)*YR + NMAIN(JNAO)*YL
   enddo
   call CUOFP  ! new metric is used for CU, MU EFable to be restored after VR is corrected
elseif (NA1 < JNA1O) then  ! NA1 decreases
   TE(NA1)    = TE(JNA1O)
   TI(NA1)    = TI(JNA1O)
   NE(NA1)    = NE(JNA1O)
   NI(NA1)    = NI(JNA1O)
   NIZ1(NA1)  = NIZ1(JNA1O)
   NIZ2(NA1)  = NIZ2(JNA1O)
   NIZ3(NA1)  = NIZ3(JNA1O)
   NALF(NA1)  = NALF(JNA1O)
   NHE3(NA1)  = NHE3(JNA1O)
   NHYDR(NA1) = NHYDR(JNA1O)
   NDEUT(NA1) = NDEUT(JNA1O)
   NTRIT(NA1) = NTRIT(JNA1O)
   NMAIN(NA1) = NMAIN(JNA1O)
   F0(NA1) = F0(JNA1O)
   F1(NA1) = F1(JNA1O)
   F2(NA1) = F2(JNA1O)
   F3(NA1) = F3(JNA1O)
   F4(NA1) = F4(JNA1O)
   F5(NA1) = F5(JNA1O)
   F6(NA1) = F6(JNA1O)
   F7(NA1) = F7(JNA1O)
   F8(NA1) = F8(JNA1O)
   F9(NA1) = F9(JNA1O)
   MU(NA1) = RTOR*IPL/(5.*BTOR*ROC*G22(NA)*IPOL(NA1))
endif

return
end subroutine SETEDGE

!======================================================================|
subroutine EDCELL(JNA1O)
!----------------------------------------------------------------------|
! Input:
! JNA1O - old grid size
! DFPDR = (FP(NA1)-FP(NA))/HRO
! Output:
! if (NA1 == JNA1O) then FP(NA1) 
!        all main arrays otherwise
! Note! FP defined in this subroutine will be assigned to FPO in OLDNEW
!----------------------------------------------------------------------|

use parameter_inc
use const_inc, only: NA, NA1, ROC, RTOR, BTOR, IPL, TIME
use status_inc

implicit none

integer, intent(in) :: JNA1O

integer :: JNAO, j
double precision :: YL, YR, YD

if (NA1 > JNA1O) then ! NA1 increases: JNA1O >= NA
   JNAO = JNA1O - 1
   TE(NA1)    = TE(JNA1O)
   TI(NA1)    = TI(JNA1O)
   NE(NA1)    = NE(JNA1O)
   NI(NA1)    = NI(JNA1O)
   NIZ1(NA1)  = NIZ1(JNA1O)
   NIZ2(NA1)  = NIZ2(JNA1O)
   NIZ3(NA1)  = NIZ3(JNA1O)
   NALF(NA1)  = NALF(JNA1O)
   NHE3(NA1)  = NHE3(JNA1O)
   NHYDR(NA1) = NHYDR(JNA1O)
   NDEUT(NA1) = NDEUT(JNA1O)
   NTRIT(NA1) = NTRIT(JNA1O)
   NMAIN(NA1) = NMAIN(JNA1O)
   F0(NA1) = F0(JNA1O)
   F1(NA1) = F1(JNA1O)
   F2(NA1) = F2(JNA1O)
   F3(NA1) = F3(JNA1O)
   F4(NA1) = F4(JNA1O)
   F5(NA1) = F5(JNA1O)
   F6(NA1) = F6(JNA1O)
   F7(NA1) = F7(JNA1O)
   F8(NA1) = F8(JNA1O)
   F9(NA1) = F9(JNA1O)
! (RHO(JNA1O)-RHO(JNAO)) - new value (after NEWGRD) is HRO, 

   write(*, '(/A, 2I4, A, 2I4, 2F10.6)') 'EDCELL:', JNAO, NA, ' ->', JNA1O, NA1, TIME, IPL
! Set FP(NA1) keeping I_pl

   YD = RHO(NA1) - RHO(JNAO)
   do j=JNA1O, NA ! Linear interpolation on the 2nd last node
      YR = (RHO(j) - RHO(JNAO))/YD
      YL = (RHO(NA1) - RHO(j))/YD
      TE(j) = TE(NA1)*YR + TE(JNAO)*YL
      TI(j) = TI(NA1)*YR + TI(JNAO)*YL
      NE(j) = NE(NA1)*YR + NE(JNAO)*YL
      NI(j) = NI(NA1)*YR + NI(JNAO)*YL
      F1(j) = F1(NA1)*YR + F1(JNAO)*YL
      F2(j) = F2(NA1)*YR + F2(JNAO)*YL
      F3(j) = F3(NA1)*YR + F3(JNAO)*YL
   enddo
! Define FP(NA) from sigma*E=j_OH at j=NA

elseif (NA1 < JNA1O) then  ! NA1 decreases
   TE(NA1)    = TE(JNA1O)
   TI(NA1)    = TI(JNA1O)
   NE(NA1)    = NE(JNA1O)
   NI(NA1)    = NI(JNA1O)
   NIZ1(NA1)  = NIZ1(JNA1O)
   NIZ2(NA1)  = NIZ2(JNA1O)
   NIZ3(NA1)  = NIZ3(JNA1O)
   NALF(NA1)  = NALF(JNA1O)
   NHE3(NA1)  = NHE3(JNA1O)
   NHYDR(NA1) = NHYDR(JNA1O)
   NDEUT(NA1) = NDEUT(JNA1O)
   NTRIT(NA1) = NTRIT(JNA1O)
   NMAIN(NA1) = NMAIN(JNA1O)
   F0(NA1) = F0(JNA1O)
   F1(NA1) = F1(JNA1O)
   F2(NA1) = F2(JNA1O)
   F3(NA1) = F3(JNA1O)
   F4(NA1) = F4(JNA1O)
   F5(NA1) = F5(JNA1O)
   F6(NA1) = F6(JNA1O)
   F7(NA1) = F7(JNA1O)
   F8(NA1) = F8(JNA1O)
   F9(NA1) = F9(JNA1O)
   MU(NA1) = RTOR*IPL/(5.*BTOR*ROC*G22(NA)*IPOL(NA1))
endif

end subroutine EDCELL

!======================================================================|
subroutine SETGEO(jst)
!----------------------------------------------------------------------|
! Presently the subroutine is called with jst=0 only
!----------------------------------------------------------------------|
! input:  jst, NB1, HRO, ROC, AB, RHO(j)
! Output: SHIF(jst:NB1), ELON(jst:NB1), TRIA(jst:NB1), 
!  AMETR(jst:NB1), DRODA(jst:NB1)

use parameter_inc
use const_inc, only: NB1, AB, ROC, RTOR, SHIFT, UPDWN, ELONG, TRIAN
use status_inc, only: RHO, SHIF, SHIV, ELON, TRIA, AMETR, DRODA

implicit none

integer, intent(in) :: jst

integer :: j
double precision :: YDA, YA, YR1, YR2, ROC3A
!----------------------------------------------------------------------|
if (jst + 1 > NB1) return
do j=jst + 1, NB1
   YR2 = min(1.d0, (RHO(J)/ROC)**2)
   SHIF(J) = SHIFT
   SHIV(J) = UPDWN
   ELON(J) = 0.5*(1. + ELONG + (ELONG - 1.)*YR2)
   TRIA(J) = TRIAN*YR2
enddo
YDA = 0.1*AB/NB1

YR1 = .0
if (jst == 0) then
   YA = 0.
   YR2 = 0.
else
   YA = AMETR(jst+1)
   YR2 = RHO(jst+1)
endif
do j=jst + 1, NB1
2  continue
   if (YR2 > RHO(j)) goto 3
   YR1 = YR2
   YA = YA + YDA
   YR2 = ROC3A(RTOR, SHIF(j), YA, ELON(j), TRIA(j))
   if (YR2 <= RHO(j)) goto 2
   DRODA(j) = YDA/(YR2 - YR1)
3  continue
   AMETR(j) = YA - YDA + (RHO(j) - YR1)*DRODA(j)
enddo

return
end subroutine SETGEO

!======================================================================|
subroutine ADCMP(boundary_cond, dfpdrbm12)

use const_inc
use status_inc
use debugger, only: markloc

implicit none

integer, intent(in) :: boundary_cond
double precision, intent(in) :: dfpdrbm12

integer :: J
double precision, dimension(NRD) :: XST, XHH

call markloc('ADCMP')

!Temperature e 
if (ADCMPF == 1.) then 
   call DERIV_EF(XRHO(1: NA1), SXHO(1: NA1), 1, TE(1: NA1), XHH(1: NA1), 1, NA1, 1)
   call GRID2GRID(2, SXHO(1: NA1), XHH(1: NA1), XST(1: NA1), NA1, 0)
   do J=1, NA1
      XHH(J) = TAU*RBDOT*XRHO(J)*XST(J)
      TE(J)  = TE(J) + XHH(J)
   enddo
endif

if (NA1E < NA1) then
   do J=NA1E + 1, NA1
      TE(J) = TEO(J)
   enddo
endif

! Temperature i 
if (ADCMPF == 1.) then 
   call DERIV_EF(XRHO(1: NA1), SXHO(1: NA1), 1, TI(1: NA1), XHH(1: NA1), 1, NA1, 1)
   call GRID2GRID(2, SXHO(1: NA1), XHH(1: NA1), XST(1: NA1), NA1, 0)
   do J=1, NA1
      XHH(J) = TAU*RBDOT*XRHO(J)*XST(J)
      TI(J) = TI(J) + XHH(J)
   enddo
endif
if (NA1I < NA1) then
   do J=NA1I + 1, NA1
      TI(J) = TIO(J)
   enddo
endif

!Density 
if (ADCMPF == 1.) then 
   call DERIV_EF(XRHO(1: NA1), SXHO(1: NA1), 1, NE(1: NA1), XHH(1: NA1), 1, NA1, 1)
   call GRID2GRID(2, SXHO(1: NA1), XHH(1: NA1), XST(1: NA1), NA1, 0)
   do J=1, NA1
      XHH(J) = TAU*RBDOT*XRHO(J)*XST(J)
      NE(J)  = NE(J) + XHH(J)
   enddo
endif
if (NA1N < NA1) then
   do J=NA1N+1, NA1
      NE(J) = NEO(J)
   enddo
endif

if (LEQ(4) > 0) then
   call DERIV_EF(XRHO(1: NA1), SXHO(1: NA1), 1, FP(1: NA1), XHH(1: NA1), 1, NA1, 1)
   call GRID2GRID(2, SXHO(1: NA1), XHH(1: NA1), XST(1: NA1), NA1, 0)
   do J=1, NA
      FP(J) = FP(J) + TAU*RBDOT*XRHO(J)*XST(J)
   enddo

   SELECT CASE(boundary_cond)

      CASE(1) !  Prescribed plasma current
         J = NA1
         FP(J) = FP(J-1) + HRO*dfpdrbm12
      CASE(2) ! Prescribed psi_b
         FP(NA1) = FP(NA1)
      CASE(3) ! psi_n+g dpsidrb=psiext
         FP(NA1) = (HRO*PSIEXT + FP(NA)*ROC*PSPLEX)/(HRO + ROC*PSPLEX)
      CASE(4)
         FP(NA1) = (HRO*PSIEXT + FP(NA)*PSPLEX)/(HRO + PSPLEX)  
      CASE DEFAULT
         J = NA1
         FP(J) = FP(J) + TAU*RBDOT*XRHO(J)*XST(J)

   END SELECT

! recompute current and mu
   call CUOFP

endif

return
end subroutine ADCMP

!======================================================================|

subroutine yrjkdr(YR, JK, YDR)
!computes index position JK, and volume differential dV/HRO at position JK
!input YR in units of RHO (so meters). 

use const_inc, only: HRO, ROC, NA1
use status_inc, only: VR

implicit none

double precision, intent(in) :: YR
integer, intent(out) :: JK
double precision, intent(out) :: YDR

if(YR <= 0.) then
   JK = 1
   YDR = 0.
else if(YR > ROC) then
   JK = NA1
   YDR = VR(JK)
else
   JK = int(YR/HRO) + 1      ! Next (outer) grid point label
   YDR = (JK - YR/HRO - 0.5)*VR(JK) ! To be subtracted
endif

return
end subroutine yrjkdr
