class FNC:

    header  = \
"""
double precision :: VINT, IINT, LININT, GRAD, GRADS, FRMAX, FRMIN, RFMIN, RFMAX, RFVAL, AFVAL, RFVEX, AFVEX, RFVIN, AFVIN, RFA, RFAN, XFA, XFAN, AFR, AFX, RECR, ATR, ATX, TIMINT, TIMDER, TIMAVG, GAUSS, RADIAL, RADINT, ASTEP, RSTEP, XSTEP, STEP, CUT, FTBOX, FXBOX, FABOX, FIXVAL, FTAV, FTMIN, FTMAX, FRAMP, FJUMP
external IINT
"""


class SUBPROC:

    header = \
'''subroutine SUBPROC

use const_inc
use outcmn_inc, only: IFSBP

implicit none
'''

    sbp_init = \
"""! **** Fill shared memory segments
call markloc("setvars")
call setvars(DEVAR, NA1, NB1, n_bouncon, NRD)
call markloc("setarrs")
call setarrs(plasma_profs, NRD)
"""


class POSTEP:

    header = \
"""subroutine POSTEP

use const_inc
use status_inc
use ipc_mod
use nclass_mod
use strahl_mod
use outcmn_inc, only: CPT, CPTSBR
use debugger, only: markloc

implicit none

integer :: IFSUB

"""

ininam_header = \
"""subroutine ININAM(LISTSB)

use parameter_inc, only: NSBMX, NRD
use outcmn_inc
use const_inc
use status_inc
use debugger, only: markloc

implicit none

character(len=64), dimension(NSBMX), intent(out) :: LISTSB

integer :: j

call markloc("xar_usage")
"""


class SETVAR:

    header = \
"""subroutine SETVAR

use const_inc
use status_inc
use debugger, only: markloc

implicit none

integer :: j

call markloc("setvar.tmp")

do j=1, NA1
ZEF(J)= max(1.d0, ZEFX(J))
ZMAIN(J) = ZMJ
AMAIN(J) = AMJ
NI(J) = NE(J)/ZMJ
enddo
if (ABC+abs(SHIFT) > AB) then
write(*, *) char(7), ">>> Warning >>> Inconsistent boundary setting."
if (LEQ(5) == 3) then
write(*, *) "    Plasma beyond the vacuum vessel has been cut off"
else
write(*, *) "    Plasma boundary intersects the vacuum vessel"
endif
endif
"""

class CUAS:

    header = \
'''! **** Current profile adjustment
call markloc("CU adjustment")
YB = 0.4*GP
YC = YB*RTOR/BTOR
YD = -0.8*GP*GP*RTOR
YA = 2./(HRO**2*YD)
'''

    mv1 = \
'''YFV = GP2 * HRO**2 * BTOR
FV(1) = 0.
YM = 0.
do j=1, NA1
'''

    mv2 =\
'''YM1 = MV(j)*j
YM2 = YM1*G22(j)
CV(j)=(YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)
if (j < NA) then
YDF = YFV*YM1
elseif (j == NA) then
YDF = GP2*HRO*BTOR*YM1*HRO
endif
if (j <= NA) then
FV(j+1) = FV(j) + YDF
else
CV(NA1) = CV(NA)
endif
YM = YM2'
enddo ! j (radial loop)
'''

    uloop_1 = \
'''do  j=1, NA1
FPO(j) = FV(j) + 0.2*GP*RTOR*IPL*(RHO(j)/ROC)**2
enddo
do JCALL=1, 10
YF = 1.E3
'''

    uloop_2 = \
'''FP(NA) = 0.4*GP*HRO*RTOR/(G22(NA)*IPOL(NA1))
FP(NA1) = (HRO - 0.5*HRO)*ROC*CC(NA1)/RTOR/IPOL(NA1)/TAU
FP(NA-1) = 1. + FP(NA)*FP(NA1)
'''

    uloop_3 = \
'''FP(NA) = -1.
YWA(1: NA1) = G22(1: NA1)
do j=1, NA1
YWGN(j) = 1.
YWGO(j) = 1.
YWHN(j) = 1.
YWHO(j) = 1.
YWB(j) = 0.
YWS(j) = 0.
YWR(j) = 0.
YWQ(j) = 0.
YWG11(j) = 1.
YWNB(j) = 1./RHO(j)
YWWB(j) = 1./RHO(j)
YVR(j) = CC(j)*.4*GP*RHO(j)/IPOL(j)**2.0
YWM(j) = 1./YVR(j)
YWD(j) = -(VR(j)/(GP2*RHO(j)*CC(j))) *(CUBS(j)+CD(j))
enddo
imethod = int(INUME3)
RABDOT = (abs(ADCMPF))*RBDOT
BABDOT = (abs(ADCMPF))*BBDOT

call RUNEQ_EF( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), FPO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWS(1: NA1), YWD(1: NA1), RABDOT, BABDOT, NA1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1:7), FP(1: NA1), YWQ(1: NA1), YQDCMF(1: NA1), ADCMPF,MPHIT(1: NA1) )

DFPDRB    = -YWQ(NA)/G22(NA)
dfpdrbm12 = -YWQ(NA)/G22(NA)

do J=1, NA1
UPL(J)  = (FP(J) - FPO(J))/TAU/YF
ULON(J) = IPOL(J)*G33(J)*UPL(J)
FP(J)  = FP(J) - (9.*FP(1) - FP(2))/8.
FPO(J) = FP(J)
enddo
UPL(NA1)  = UPL(NA-2)
ULON(NA1) = ULON(NA)
call CUOFP
'''

    cubs = \
'''if (j == NA1) then
CUBS(NA1) = CUBS(NA) + (CUBS(NA) - CUBS(NA-1))*HRO/HRO
else
'''

    beta = \
'''if (BETAJR(0.1*ROC) > 100.) then
write(*, *) ">>>  ERROR  >>> The initial plasma pressure is too high"
stop
endif
'''

    mu = \
'''enddo
call CUOFMU
YU = GP2*RTOR
do j=1, NA1
'''

    cu1 = \
'''YF = GP2*HRO**2 * BTOR
YC = 0.4*GP*RTOR/BTOR
YM = 0.
YMCD = 0.
do J=1, NA1
'''

    cu2 = \
'''
if (j == NA1) CYCLE
YM   = YM   +  CU(J)*RHO(J)/(G33(J)*IPOL(J)**3)
YMCD = YMCD + YWA(J)*RHO(J)/(G33(J)*IPOL(J)**3)
enddo
YIOH = GP2*YM*HRO*IPOL(NA1)
YICD = GP2*YMCD*HRO*IPOL(NA1)
CU(1: NA1) = CU(1: NA1)*(IPL - YICD)/YIOH
YM = 0.
do j=1, NA
YM = YM + YC*CU(J)*RHO(J)/(G33(J)*IPOL(J)**3)
YM1 = YM/G22(J)
MU(J) = YM1/J
FP(J+1) = FP(J) + YF*YM1
enddo
MU(NA1) = MU(NA-2)
YU = GP2*RTOR
YJ_CU = (FP(NA1) - FP(NA))/HRO * IPOL(NA1) * G22(NA)/(0.4*GP*RTOR)
YJ_CU = IPL/YJ_CU
do j=1, NA1
CU(J) = YJ_CU*CU(J)
'''

    cu_mu = \
'''if (CC(J) > 0.0) then
ULON(J) = YU*(CU(J)-YWA(J))/CC(J)
else
ULON(J) = 0.0
endif
UPL(J) = ULON(J)/(IPOL(J)*G33(J))
enddo ! j (radial loop)
ULON(NA1) = ULON(NA-2)
UPL(NA1)  = ULON(NA1)/(IPOL(NA1)*G33(NA1))
DFPDRB = 0.4*IPL/IPOL(NA1)/G22(NA1)*GP*RTOR
'''


class INIVAR:

    header = \
'''subroutine INIVAR

use outcmn_inc
use const_inc
use nclass_mod
use status_inc
use debugger, only: markloc

implicit none

include 'tmp/declar.fml'
include 'tmp/declar.fnc'

integer :: j1

call markloc("inivar")

TAU = TAUMIN
'''

    ne  = \
'''if (NA1N == NB1) then
j1 = 0
do j=1, NARRX
if(EXARNM(j) == "NEX   " .and. IFDFAX(j) < 0) j1 = j
enddo
if (j1 /= 0) write(*, *) " >>> Warning: NEX is not defined"
endif
NE(1: NA1) = NEX(1: NA1)
'''

    ni = \
"""do J=1, NA1
NI(j) = min(NI(j), NE(j)*ZEF(j)/ZMAIN(j)**2)
enddo
"""

inam_sb = \
"""do j=1, NSBMX
IFSBP(j) = 0
LISTSB(j) = char(0)
IFSBX(j) = 0
enddo
"""

class DETVAR:

    header = \
'''subroutine DETVAR

use const_inc
use status_inc
use ipc_mod
use nclass_mod
use strahl_mod
use outcmn_inc
use debugger, only: markloc

implicit none

include 'tmp/declar.fml'
include 'tmp/declar.fnc'

integer :: jdetv, ifsub
integer, external :: ifipc

call markloc("detvar.tmp (time signals)")
'''

    rad_tail  = \
"""
if (LEQ(5) == 3) then
SHIFT = min(SHIFT, 0.9*AB)
ABC = min(ABC, AB - abs(SHIFT))
ABC = max(ABC, 0.1*AB)
endif
"""


class CUEQ:

    ipl = \
"""
! Prescribed plasma current
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = 1.
YWC(6) = -1.
YWC(7) = 1./G22(NA)*IPL*RTOR/IPOL(NA1)
YWC(7) = HRO*0.4*GP*YWC(7)
bc_type_for_fp = 1
if (ITFBP /= 0.0 .and. ITFBE < TIME) then
if (ibcpsi_fb >= 0) then
if (ITFBP < 0.0 .and. ibcpsi_fb >= 2.) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX*ROC
YWC(6) = -PSPLEX*ROC
YWC(7) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (ibcpsi_fb <= 1) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = 1.
YWC(6) = -1.
YWC(7) = 1./G22(NA)*IPL*RTOR/IPOL(NA1)
YWC(7) = HRO*0.4*GP*YWC(7)
bc_type_for_fp = 1
endif
endif
"""

    lext = \
"""
! Prescribed loop voltage:
FP(NA1) = FPO(NA1) + TAU*UEXT
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 1.
bc_type_for_fp = 2
if (ITFBP /= 0.0 .and. ITFBE < TIME) then
if (ibcpsi_fb > 0) then
if (ITFBP < 0.0 .and. ibcpsi_fb >= 2.) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX*ROC
YWC(6) = -PSPLEX*ROC
YWC(7) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (ibcpsi_fb <= 1) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 1.
bc_type_for_fp = 2
endif
endif
"""

    circ = \
"""
! Circuit equation:
if (TIME-TSTART <= TAU) then
PSIEXT = FP(NA1) + LEXT*IPL
PSPLEX = (FP(NA1) - FP(NA))/HRO
PSPLEX = LEXT*IPL/PSPLEX
endif
PSIEXT = PSIEXT + TAU*UEXT
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX
YWC(6) = -PSPLEX
YWC(7) = PSIEXT*HRO
bc_type_for_fp = 4
if (ITFBP /= 0.0 .and. ITFBE < TIME) then
if (ibcpsi_fb > 0) then
if (ITFBP < 0.0 .and. ibcpsi_fb >= 2.) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX*ROC
YWC(6) = -PSPLEX*ROC
YWC(7) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (ibcpsi_fb <= 1) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 1.
bc_type_for_fp = 2
endif
endif
"""

    eqn = \
"""
do j=1, NA1
YWGN(j) = 1.
YWGO(j) = 1.
YWHN(j) = 1.
YWHO(j) = 1.
YWB(j) = 0.
YWS(j) = 0.
YWR(j) = 0.
YWQ(j) = 0.
YWG11(j) = 1.
YWNB(j) = 1./RHO(j)
YWWB(j) = 1./RHO(j)
YVR(j) = CC(j)*.4*GP*RHO(j)/IPOL(j)**2.0
YWM(j) = 1./YVR(j)
YWD(j) = -(VR(j)/(GP2*RHO(j)*CC(j))) * (CUBS(j) + CD(j))
enddo
imethod = int(INUME3)
RABDOT = (abs(ADCMPF))*RBDOT
BABDOT = (abs(ADCMPF))*BBDOT

call RUNEQ_EF( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), FPO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWS(1: NA1), YWD(1: NA1), RABDOT, BABDOT, NA1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1:7), FP(1: NA1), YWQ(1: NA1), YQDCMF(1: NA1), ADCMPF, MPHIT(1: NA1) )

DFPDRB    = -YWQ(NA)/G22(NA)
dfpdrbm12 = -YWQ(NA)/G22(NA)

YWR(1) = MU(NA1)*GP2*ROC**2 * BTOR*BABDOT
"""


class TEEQN:

    eqn = \
"""
do j=1, NA1
YWH(j) = 1.
YWR(j) = 0.
YVR(j) = VR(j)
YWGN(j) = VR(j)**(5./3.)
YWGO(j) = VRO(j)**(5./3.)
YWHN(j) = 3./2.*NE(j)
YWHO(j) = 3./2.*NEO(j)
YWM(j) = 625./VR(j)
YWWB(j) = VR(j)**(5./3.)
YWNB(j) = 3./2.*NE(j)*YWWB(j)/YVR(j)/YWM(j)
enddo
imethod = int(INUME2)
RABDOT = abs(ADCMPF)*RBDOT
BABDOT = abs(ADCMPF)*BBDOT

call RUNEQ_EF( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), TEO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), G11(1: NA1)/625, YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), 625*PET(1: NA1), 625*PETOT(1: NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1:7), TE(1: NA1), QE(1: NA1), YQDCM(1: NA1), ADCMPF, MPHIT(1: NA1) )

do j=1, NA1
te(j) = max(te(j), 0.001)
enddo
"""

    assigned = \
'''do j=1, NA
QE(J) = -G11(J)*(YWA(J)*(TE(J+1) - TE(J))/HRO + 0.5*YWB(J)*(TE(J+1) + TE(J)))*0.0016
enddo
QE(NA1) = '''


class TIEQN:

    eqn = \
"""
do j=1, NA1
YWH(j) = 1.
YWR(j) = 0.
YVR(j) = VR(j)
YWGN(j) = VR(j)**(5./3.)
YWGO(j) = VRO(j)**(5./3.)
YWHN(j) = 3./2.*NI(j)
YWHO(j) = 3./2.*NIO(j)
YWM(j) = 625./VR(j)
YWWB(j) = VR(j)**(5./3.)
YWNB(j) = 3./2.*NI(j)*YWWB(j)/YVR(j)/YWM(j)
enddo
imethod = int(INUME2)
RABDOT = abs(ADCMPF)*RBDOT
BABDOT = abs(ADCMPF)*BBDOT

call RUNEQ_EF( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), TIO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), G11(1: NA1)/625, YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), 625*PIT(1: NA1), 625*PITOT(1: NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1:7), TI(1: NA1), QI(1: NA1), YQDCM(1: NA1), ADCMPF, MPHIT(1: NA1) )

do j=1, NA1
ti(j) = max(ti(j), 0.001)
enddo
"""

    assigned = \
'''do J=1, NA
QI(J) = -G11(J)*(YWA(J)*(TI(J+1) - TI(J))/HRO + 0.5*YWB(J)*(TI(J+1) + TI(J)))*0.0016
enddo
QI(NA1) = '''


class NEEQN:

    eqn = \
"""
do j=1, NA1
YWHN(j) = 1.
YWHO(j) = 1.
YWGN(j) = VR(j)
YWGO(j) = VRO(j)
YWR(j)  = 0.
YVR(j)  = VR(j)
YWM(j)  = 1./VR(j)
YWNB(j) = VR(j)
YWWB(j) = VR(j)
enddo
imethod = int(INUME1)
RABDOT = abs(ADCMPF)*RBDOT
BABDOT = abs(ADCMPF)*BBDOT

call RUNEQ_EF( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), NEO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), G11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), SNN(1: NA1), SN(1: NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1: 7), NE(1: NA1), QN(1: NA1), YQDCM(1: NA1), ADCMPF, MPHIT(1: NA1) )
"""

    assigned = \
"""
do J=1, NA
QN(J) = G11(J)*(-YWA(J)*(NE(J+1) - NE(J))/HRO - 0.5*YWB(J)*(NE(J+1)+NE(J))) + SLAT(J)*GNX(J)
GN(J) = QN(J)/SLAT(J)
SNTOT(J) = SNTOT(J) + SNN(J)*NE(J)
enddo
"""

class NIAS:

    iondensassign = \
'''! **** Ion density assignment
call markloc("NI assignment")

if (trim(MACHINE) == 'aug') then
do J=1, NA1
NI(J) = F1(J) + F2(J) + F3(J) + F4(J) + F5(J) + F6(J) + F7(J) + F8(J) + F9(J)  ! complete AUG
enddo
endif

if (trim(MACHINE) == 'demo' .or. trim(MACHINE) == 'iter') then
do J=1, NA1
NI(J) = F1(J) + F2(J) + F3(J) + F4(J) + F6(J) + F7(J) + F8(J)  ! DEMO, ITER: D, T, H, He, Be, W, Ne
enddo
endif
'''


class UPEQN:

    eqn = \
'''YWD(ND1) = 0.
do j=1, NA1
YWH(j) = 1.
YWD(j)  = TTRQI(j)
YWGN(j) = VR(j)
YWGO(j) = VRO(j)
YWHN(j) = UPS0(j)
YWHO(j) = UPS0O(j)
YVR(j)  = VR(j)
MPHIT(j) = UPARO(j)
YWM(j) = 1./VR(j)
YWG11(j) = G11(j)
YWNB(j) = VR(j)*UPS0(j)   ! at the moment adiabatic compression is done only on UPS0
! missing UPS1 and UPS2 terms in <M_phi>
YWWB(j) = VR(j)
if (IPROT >= 1.) then
! contributions to net torque
MPHIT(j) = MPHIT(j) + UPS1O(j)/UPS0O(j)
MPHIT(j) = MPHIT(j) + UPS2O(j)/UPS0O(j)
TTRQ(j) = TTRQ(j) - ( VR(J)*(UPS1(J) + UPS2(J)) - VRO(J)*(UPS1O(J) + UPS2O(J)) )*2./ (VR(J) + VRO(J))/TAU
! contributions to stress tensor
YWgradF(J)  = 2.*(IPOL(J+1)  - IPOL(J)) /(IPOL(J+1)  + IPOL(J)) /HRO  !d log I / drho
YWgradb2(J) = 2.*(BDB02(J+1) - BDB02(J))/(BDB02(J+1) + BDB02(J))/HRO  !d log <B**2>/drho
YWR(J) = YWR(J) + RTOR/IPOL(J)*XUPAR(J)*DLNEOD(J) - RTOR/IPOL(J)*(CNPAR(J) + XUPAR(J)*YWgradF(J))*DLNEO(J) + (XUPAR(J) - XUPAP(J))*RTOR/IPOL(J)*BDB02(J)*BTOR*SGNEOD(J) + RTOR/IPOL(J)*BDB02(J)*BTOR*(XUPAR(J)*YWgradb2(J) - XUPAR(J)*YWgradF(J) + CNPAP(J) - CNPAR(J))*SGNEO(J) + RTOR/IPOL(J)*(CNPAD(J) - CNPAR(J) - XUPAD(J)*YWgradF(J))*DDNEO(J) + RTOR/IPOL(J)*(XUPAR(J) - XUPAD(J))*DDNEOD(J)
endif
enddo
imethod = int(INUME4)
RABDOT = (abs(ADCMPF))*RBDOT
BABDOT = (abs(ADCMPF))*BBDOT

call RUNEQ_EF(YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), UPARO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWD(1: NA1), TTRQ(1: NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1:7), UPAR(1: NA1), QU(1: NA1), YQDCM(1: NA1), ADCMPF, MPHIT(1: NA1))

MPHIT = 0.
'''

    assigned = \
'''do j=1, NA
QU(J) = -G11(J)*(YWA(J)*(UPAR(J+1) - UPAR(J))/HRO + 0.5*YWB(J)*(UPAR(J+1) + UPAR(J)))*0.0016
enddo
'''


    uparo = \
'''UPARO(ND1: NA1) = UPAR(ND1: NA1)
QU(4)  = 1.
YWC(4) = 1.
'''


class TETIEQN:

    eqn = \
'''do j=1, NA1
YWH(j) = 1.
YWGN(j) = VR(j)**(5./3.)
YWGO(j) = VRO(j)**(5./3.)
YWM(j) = 625./VR(j)
YWR(j) = 0.
YVR(j) = VR(j)
YWW1B(j) = VR(j)**(5./3.)
YWN1B(j) = 3./2.*NE(j)*YWW1B(j)/YVR(j)/YWM(j)
YWW2B(j) = VR(j)**(5./3.)
YWN2B(j) = 3./2.*NI(j)*YWW2B(j)/YVR(j)/YWM(j)
enddo
imethod = int(INUME2)
RABDOT = (abs(ADCMPF))*RBDOT
BABDOT = (abs(ADCMPF))*BBDOT
'''

    runeq = \
'''call RUNEQTIMP_EF(YWGN(1:NA1), 3./2.*NE(1:NA1), 3./2.*NI(1:NA1), YWGO(1:NA1), 3./2.*NEO(1:NA1), 3./2.*NIO(1:NA1), TEO(1:NA1), TIO(1:NA1), YWN1B(1:NA1), YWN2B(1:NA1), YWW1B(1:NA1), YWW2B(1:NA1), YVR(1:NA1), YWM(1:NA1), G11(1:NA1)/625, YWA1(1:NA1), YWA2(1:NA1), YWB1(1:NA1), YWB2(1:NA1), YWR(1:NA1), YWR(1:NA1), 625*PET(1:NA1), 625*PIT(1:NA1), 625*PETOT(1:NA1), 625*PITOT(1:NA1), 0.0*YWR(1:NA1), 0.0*YWR(1:NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ROC, RHO(1:NA1), imethod, YWC1(1:NA1), YWC2(1:NA1), TE(1:NA1), TI(1:NA1), QE(1:NA1), QI(1:NA1), ADCMPF)
if (ND1 < NA1) then
do j=ND1+1, NA1
QE(j) = QE(ND1)
QI(j) = QI(ND1)
enddo
endif
PETOT=PE+PET*TE
'''

    teold = \
'''TEO(ND1: NA1) = TE(ND1: NA1)
QE(4)   = 1.
YWC1(4) = 1.
'''

    tiold = \
'''TIO(ND1: NA1) = TI(ND1: NA1)
QI(4)   = 1.
YWC2(4) = 1.
'''


class RADOUT:

    header = \
"""subroutine RADOUT

!------------------------------------------------------------
! Radial profile plotting
!------------------------------------------------------------

use parameter_inc
use const_inc
use status_inc
use outcmn_inc
use debugger, only: markloc, debug

implicit none

integer irado

include 'tmp/declar.fml'
include 'tmp/declar.fnc'

call markloc('radout.tmp')
do irado=1, NAB
J = irado
"""


class TIMOUT:

    header = \
"""
enddo

return
end subroutine RADOUT

!------------------------------------------------------------
subroutine TIMOUT

!------------------------------------------------------------
! Time traces plotting
!------------------------------------------------------------

use parameter_inc
use const_inc
use status_inc
use outcmn_inc
use timeoutput_inc
use debugger, only: markloc, debug

implicit none

include 'tmp/declar.fml'
include 'tmp/declar.fnc'

call markloc('timout.tmp')
"""

class FJEQN:

    eqn = \
'''do j=1, NA1
YWHN(j) = 1.
YWHO(j) = 1.
YWGN(j) = VR(j)
YWGO(j) = VRO(j)
YWR(j)  = 0.
YVR(j)  = VR(j)
YWM(j)  = 1./VR(j)
YWNB(j) = VR(j)
YWWB(j) = VR(j)
enddo
imethod = int(INUME1)
RABDOT = (abs(ADCMPF))*RBDOT
BABDOT = (abs(ADCMPF))*BBDOT
'''

class CUEQN:

    header = \
'''! **** Current equation
call markloc("Current equation")
YC =  0.4*GP*RTOR/BTOR
YD = -0.8*GP**2 * RTOR
YA = 2./(HRO**2 * YD)
do J=1, NA1
'''

    mv = \
'''YFV = GP2*HRO**2 * BTOR
FV(1) = 0.
YM = 0.
do j=1, NA1
YM1 = MV(j)*j
YM2 = YM1*G22(j)
CV(j) = (YM2 - YM)*G33(j)*IPOL(j)**3/YC/RHO(j)
if (j < NA) then
YDF = YFV*YM1
elseif(j <= NA) then
YDF = GP2*HRO*BTOR*YM1*HRO
endif
if (j <= NA) then
FV(j+1) = FV(j) + YDF
else
CV(NA1) = CV(NA)
endif
YM = YM2
enddo
'''

    prescribed_ipl = \
'''! Prescribed plasma current:
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = 1.
YWC(6) = -1.
YWC(7) = 1./G22(NA)*IPL*RTOR/IPOL(NA1)
YWC(7) = HRO*0.4*GP*YWC(7)
bc_type_for_fp=1
! For psifb
! when using the free boundary circuit equations with free current,
! then use mixed b.c.
if (ITFBP /= 0.0 .and. ITFBE < TIME) then
if (ibcpsi_fb > 0) then
!case implicit
if (ITFBP < 0.0 .and. ibcpsi_fb >= 2.) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX*ROC
YWC(6) = -PSPLEX*ROC
YWC(7) =  PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (ibcpsi_fb <= 1) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = 1.
YWC(6) = -1.
YWC(7) = 1./G22(NA)*IPL*RTOR/IPOL(NA1)
YWC(7) = HRO*0.4*GP*YWC(7)
bc_type_for_fp = 1
endif
endif
'''

    prescribed_uloop = \
'''! Prescribed loop voltage:
FP(NA1) = FPO(NA1) + TAU*UEXT
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 1.
bc_type_for_fp = 2
!For psifb
! when using the free boundary circuit equations with free current,
! then use mixed b.c.
if (ITFBP /= 0.0 .and. ITFBE < TIME) then
if (ibcpsi_fb > 0) then
!case implicit
if (ITFBP < 0.0 .and. ibcpsi_fb >= 2.) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX*ROC
YWC(6) = -PSPLEX*ROC
YWC(7) =  PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (ibcpsi_fb < =  1) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 1.
bc_type_for_fp = 2
endif
endif
'''

    circuit_eqn = \
'''! Circuit equation:
if (TIME - TSTART <= TAU) then
PSIEXT = FP(NA1) + LEXT*IPL
PSPLEX = (FP(NA1) - FP(NA))/HRO
PSPLEX = LEXT*IPL/PSPLEX
endif
PSIEXT = PSIEXT + TAU*UEXT
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX
YWC(6) =  -PSPLEX
YWC(7) = PSIEXT*HRO
bc_type_for_fp = 4
!For psifb
! when using the free boundary circuit equations with free current,
! then use mixed b.c.
if (ITFBP /= 0.0 .and. ITFBE < TIME) then
if (ibcpsi_fb > 0) then
! case implicit
if (ITFBP < 0.0 .and. ibcpsi_fb >= 2.) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 2.
YWC(5) = HRO + PSPLEX*ROC
YWC(6) = -PSPLEX*ROC
YWC(7) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (ibcpsi_fb <= 1) then
YWC(1) = HRO
YWC(2) = 0.0
YWC(3) = 0.0
YWC(4) = 1.
bc_type_for_fp = 2
endif
endif
'''

    eqn = \
'''YWA(1: NA1) = G22(1: NA1)
do j=1, NA1
YWGN(j) = 1.
YWGO(j) = 1.
YWHN(j) = 1.
YWHO(j) = 1.
YWB(j)  = 0.
YWS(j)  = 0.
YWR(j)  = 0.
YWQ(j)  = 0.
YWG11(j) = 1.
YWNB(j) = 1./RHO(j)
YWWB(j) = 1./RHO(j)
YVR(j) = CC(j)*0.4*GP*RHO(j)/IPOL(j)**2.0
YWM(j) = 1./YVR(j)
YWD(j) = -(VR(j)/(GP2*RHO(j)*CC(j))) * (CUBS(j) + CD(j))
enddo
imethod = int(INUME3)
RABDOT = (abs(ADCMPF))*RBDOT
BABDOT = (abs(ADCMPF))*BBDOT

call RUNEQ_EF( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), FPO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWS(1: NA1), YWD(1: NA1), RABDOT, BABDOT, NA1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1: 7), FP(1: NA1), YWQ(1: NA1), YQDCMF(1: NA1), ADCMPF, MPHIT(1: NA1) )

DFPDRB    = -YWQ(NA)/G22(NA)
dfpdrbm12 = -YWQ(NA)/G22(NA)
YWR(1) = MU(NA1)*GP2*ROC**2 * BTOR*BABDOT
'''

class INIT:

    end = \
'''ROCO = ROC
!FTO = FTN
do j=1, NA1
NEO(j) = NE(j)
NIO(j) = NI(j)
TEO(j) = TE(j)
TIO(j) = TI(j)
UPARO(j) = UPAR(j)
VRO(j) = VR(j)
UPS0O(j) = UPS0(j)
UPS1O(j) = UPS1(j)
UPS2O(j) = UPS2(j)
do jj=0, 9
FJO(j, jj) = FJ(j, jj)
enddo
enddo
call markloc("init done")
'''

class INIT_CONVERGE_STEP:

    header = \
'''subroutine init_converge_step(LISTSB)

use parameter_inc, only: NSBMX, NRD
use outcmn_inc
use const_inc
use status_inc
use nclass_mod
use debugger, only: markloc

implicit none

include 'tmp/declar.fml'
include 'tmp/declar.fnc'

character(len=64), intent(in) :: LISTSB(NSBMX)

integer :: J1, IFIPC, IFSUB, JDETV, ND, ND1, &
    imethod, jcall, IFSTEP

double precision :: dfpdrbm12, ARRNA1, YHRO, YB, YC, YJ_CU, &
    YM, YMCD, YIOH, YICD, YM1, YU, RABDOT, BABDOT

double precision, dimension(NRD) :: YWA, YWB, YWC, YWD, &
    YWGN, YWHN, YWGO, YWHO, YWR, YWH, YVR, YWM,&
    YWA1, YWA2, YWB1, YWB2, YWAA, YWNB, YWWB, YWN1B, YWW1B,&
    YWN2B, YWW2B, YWC1, YWC2, YWS, YQDCM, YQDCMF,&
    YWQ, YWG11, YWgradF, YWgradb2, MPHIT
'''

    tail = \
'''

return
end subroutine init_converge_step
'''


class EQNS_INC:

    header = \
'''subroutine EQNS_INC(ibcpsi_fb, bc_type_for_fp, dfpdrbm12)
!-------------------------------------------------------------------
! Perform one time step
! Note that now time step is updated at the end of a full time cycle
!-------------------------------------------------------------------

use parameter_inc, only: NRD, NSBMX
use const_inc
use status_inc
use outcmn_inc
use nclass_mod
use strahl_mod
use plasma_state
use debugger, only: markloc

implicit none

include 'tmp/declar.fml'
include 'tmp/declar.fnc'

integer, intent(in) :: ibcpsi_fb
integer, intent(out) :: bc_type_for_fp
double precision, intent(out) :: dfpdrbm12

integer :: IFSUB, imethod, ND, ND1, NODE, JCALL

double precision :: ARRNA1, RABDOT, BABDOT, YHRO, YM1, YM2, YB, YC, YJ_CU, YM, YU, YIOH, YICD, YMCD

double precision, dimension(NRD) :: YWA, YWB, YWC, YWD, YWGN, &
    YWHN, YWGO, YWHO, YWR, YWH, YVR, YWM, YWA1, YWA2, YWB1, YWB2, &
    YWAA, YWNB, YWWB, YWN1B, YWW1B, YWN2B, YWW2B, & 
    YWC1, YWC2, YWS, YQDCM, MPHIT, YQDCMF, YWQ, YWG11, YWgradF, YWgradb2

MPHIT = 0.
'''

    tail = \
'''
return
end subroutine EQNS_INC
'''
