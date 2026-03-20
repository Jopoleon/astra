pol_flux = \
'''PSIAX = EXTRAP(XRHO(1:NA1), FP(1:NA1), 0.d0, NA1, 1, .true.)
PSIBO = FP(NA1)
FP_NORM(1:NA1) = (FP(1:NA1) - PSIAX)/(PSIBO - PSIAX)
rho_pol(1:NA1) = SQRT(FP_NORM(1: NA1))
'''


class POSTEP:

    header = \
"""subroutine POSTEP

use scalars
use status
use nclass_mod
use strahl_mod
use a2tglf, only: tglf_ipc
use a2qlk, only: qlk_ipc
use a2neo, only: neo_ipc
use a2rabbit, only: rabbit
use a2torbeam, only: torba
use cpu_usage, only: wallTime_sbr, cpuTime_sbr
use debugger, only: markloc

implicit none

integer :: IFSUB, t_wall1, t_wall2, rate
double precision :: t_cpu1, t_cpu2
"""


class ININAM:
    
    header = \
"""subroutine ININAM

use read_input, only: sbr_name, n_sbr, awd
use graph_utils
use scalars
use status
use debugger, only: markloc
use json_vars, only: profxNames, n_intern

implicit none

integer :: j

call markloc("xar_usage")

allocate(DTNAME(n_intern+4*n_sbr))
"""

class SETVAR:

    header = \
"""subroutine SETVAR

use scalars
use status
use debugger, only: markloc

implicit none

integer :: j

call markloc("setvar")

do j=1, NA1
ZEF(J)= max(1.d0, ZEFX(J))
ZMAIN(J) = ZMJ
AMAIN(J) = AMJ
NI(J) = NE(J)/ZMJ
enddo
if (ABC+abs(SHIFT) > AB) then
write(*, *) char(7), ">>> Warning >>> Inconsistent boundary setting."
if (IPEQL == 3) then
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
YB = mu0
YC = YB*RTOR/BTOR
YD = -0.8*GP**2 * RTOR
YA = 2./(HRO**2 * YD)
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
if (j < NA) then
YDF = YFV*YM1
elseif (j == NA) then
YDF = GP2*HRO*BTOR*YM1*HRO
endif
if (j <= NA) then
FV(j+1) = FV(j) + YDF
endif
YM = YM2
enddo ! j (radial loop)
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
YC = mu0*RTOR/BTOR
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
MU(NA1) = EXTRAP(SXHO(1:NA), MU(1:NA), SXHO(NA1), NA, 2, .false.)
YU = GP2*RTOR
YJ_CU = (FP(NA1) - FP(NA))/HRO * IPOL(NA1) * G22(NA)/(mu0*RTOR)
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
'''


class INIVAR:

    header = \
'''subroutine INIVAR

use read_input, only: IFDFAX
use scalars
use status
use debugger, only: markloc
use json_vars, only: profxNames, n_profx

implicit none

include 'src/tmp/declar.fml'
include 'src/tmp/declar.fnc'

integer :: j1

call markloc("inivar")

TAU = TAUMIN
'''

    ne  = \
'''NE(1: NA1) = NEX(1: NA1)
'''

    ni = \
"""do J=1, NA1
NI(j) = min(NI(j), NE(j)*ZEF(j)/ZMAIN(j)**2)
enddo
"""

    
class DETVAR:

    header = \
'''subroutine DETVAR

use scalars
use status
use nclass_mod
use strahl_mod
use standard_functions
use read_input, only: IFDFVX
use a2tglf, only: tglf_alloc, tglf_out
use a2qlk, only: qlk_alloc, qlk_out
use a2neo, only: neo_alloc, neo_out 
use cpu_usage, only: wallTime_sbr, cpuTime_sbr
use json_vars, only: profxNames
use debugger, only: markloc

implicit none

include 'src/tmp/declar.fml'
include 'src/tmp/declar.fnc'

integer :: jdetv, ifsub
integer :: t_wall1, t_wall2, rate
double precision :: t_cpu1, t_cpu2

call markloc("detvar (time signals)")
call tglf_alloc
call qlk_alloc
call neo_alloc
'''

    rad_tail  = \
"""
if (IPEQL == 3) then
SHIFT = min(SHIFT, 0.9*AB)
ABC = min(ABC, AB - abs(SHIFT))
ABC = max(ABC, 0.1*AB)
endif
"""


class CUEQ:

    ipl = \
"""
! Prescribed plasma current
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = 1.
bc_values(4) = -1.
bc_values(5) = HRO*mu0/G22(NA)*IPL*RTOR/IPOL(NA1)
bc_type_for_fp = 1
if (ITFBP /= 0 .and. ITFBE < TIME) then
if (IBCPSI >= 0) then
if (ITFBP < 0 .and. IBCPSI >= 2) then
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (IBCPSI <= 1) then
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = 1.
bc_values(4) = -1.
bc_values(5) = HRO*mu0/G22(NA)*IPL*RTOR/IPOL(NA1)
bc_type_for_fp = 1
endif
endif
"""

    lext = \
"""
! Prescribed loop voltage:
FP(NA1) = FPO(NA1) + TAU*UEXT
bctype = 1
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_type_for_fp = 2
if (ITFBP /= 0 .and. ITFBE < TIME) then
if (IBCPSI > 0) then
if (ITFBP < 0 .and. IBCPSI >= 2) then
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (IBCPSI <= 1) then
bctype = 1
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_type_for_fp = 2
endif
endif
"""

    circ = \
"""
! Circuit equation:
if (TIME-TSTART <= TAU) then
PSIEXT = FP(NA1) + LEXT*IPL
PSPLEX = LEXT/ROC*5.*IPOL(NA1)*G22(NA)/GP2/RTOR
endif
PSPLEX = LEXT/ROC*5.*IPOL(NA1)*G22(NA)/GP2/RTOR
PSIEXT = PSIEXT + TAU*UEXT
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) = PSIEXT*HRO
bc_type_for_fp = 3
if (ITFBP /= 0 .and. ITFBE < TIME) then
if (IBCPSI > 0) then
if (ITFBP < 0 .and. IBCPSI >= 2) then
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (IBCPSI <= 1) then
bctype = 1
bc_values(1) = 0.0
bc_values(2) = 0.0
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
YWWB(j) = 1./RHO(j)
YVR(j) = CC(j)*mu0*RHO(j)/IPOL(j)**2
unit_coeff = 1.
YWD(j) = -(VR(j)/(GP2*RHO(j)*CC(j))) * (CUBS(j) + CD(j))
enddo
imethod = INUME3

call RUNEQ( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), FPO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWS(1: NA1), YWD(1: NA1), RBDOT, BBDOT, NA1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, FP(1: NA1), YWQ(1: NA1), YQDCMF(1: NA1), MPHIT(1: NA1) )

dfpdrbm12 = -YWQ(NA)/G22(NA)

"""


class TEEQN:

    eqn = \
"""
if (TIME >= tbeg_eq .and. TIME <= tend_eq) then
do j=1, NA1
YWR(j) = 0.
YVR(j) = VR(j)
YWGN(j) = VR(j)**(5./3.)
YWGO(j) = VRO(j)**(5./3.)
YWHN(j) = 3./2.*NE(j)
YWHO(j) = 3./2.*NEO(j)
unit_coeff = 625.
YWWB(j) = VR(j)**(5./3.)
enddo
imethod = INUME2

call RUNEQ( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), TEO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, G11(1: NA1)/unit_coeff, YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), unit_coeff*PET(1: NA1), unit_coeff*PETOT(1: NA1), RBDOT, BBDOT, ND1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, TE(1: NA1), QE(1: NA1), YQDCM(1: NA1), MPHIT(1: NA1) )
do j=1, NA1
te(j) = max(te(j), 0.001)
enddo
NA1E = ND1
else
TE(1: NA1) = TEO(1: NA1)
do j=1, NA
QE(J) = -G11(J)*(YWA(J)*(TE(J+1) - TE(J))/HRO + 0.5*YWB(J)*(TE(J+1) + TE(J)))*0.0016
enddo
QE(NA1) = QE(NA)
NA1E = NA1
endif
PETOT=PE+PET*TE
"""

    assigned = \
'''do j=1, NA
QE(J) = -G11(J)*(YWA(J)*(TE(J+1) - TE(J))/HRO + 0.5*YWB(J)*(TE(J+1) + TE(J)))*0.0016
enddo
QE(NA1) = QE(NA)
NA1E = NA1
PETOT=PE+PET*TE
'''


class TIEQN:

    eqn = \
"""
if (TIME >= tbeg_eq .and. TIME <= tend_eq) then
do j=1, NA1
YWR(j) = 0.
YVR(j) = VR(j)
YWGN(j) = VR(j)**(5./3.)
YWGO(j) = VRO(j)**(5./3.)
YWHN(j) = 3./2.*NI(j)
YWHO(j) = 3./2.*NIO(j)
unit_coeff = 625.
YWWB(j) = VR(j)**(5./3.)
enddo
imethod = INUME2

call RUNEQ( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), TIO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, G11(1: NA1)/unit_coeff, YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), unit_coeff*PIT(1: NA1), unit_coeff*PITOT(1: NA1), RBDOT, BBDOT, ND1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, TI(1: NA1), QI(1: NA1), YQDCM(1: NA1), MPHIT(1: NA1) )
do j=1, NA1
ti(j) = max(ti(j), 0.001)
enddo
NA1I = ND1
else
TI(1: NA1) = TIO(1: NA1)
do J=1, NA
QI(J) = -G11(J)*(YWA(J)*(TI(J+1) - TI(J))/HRO + 0.5*YWB(J)*(TI(J+1) + TI(J)))*0.0016
enddo
QI(NA1) = QI(NA)
NA1I = NA1
endif
PITOT=PI+PIT*TI
"""

    assigned = \
'''do J=1, NA
QI(J) = -G11(J)*(YWA(J)*(TI(J+1) - TI(J))/HRO + 0.5*YWB(J)*(TI(J+1) + TI(J)))*0.0016
enddo
QI(NA1) = QI(NA)
NA1I = NA1
PITOT=PI+PIT*TI
'''


class NEEQN:

    eqn = \
"""
if (TIME >= tbeg_eq .and. TIME <= tend_eq) then
do j=1, NA1
YWHN(j) = 1.
YWHO(j) = 1.
YWGN(j) = VR(j)
YWGO(j) = VRO(j)
YWR(j)  = 0.
YVR(j)  = VR(j)
unit_coeff  = 1.
YWWB(j) = VR(j)
enddo
imethod = INUME1

call RUNEQ( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), NEO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, G11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), SNN(1: NA1), SN(1: NA1), RBDOT, BBDOT, ND1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, NE(1: NA1), QN(1: NA1), YQDCM(1: NA1), MPHIT(1: NA1) )
do J=1, NA
QN(J) = QN(J) + SLAT(J)*GNX(J)
GN(J) = QN(J)/SLAT(J)
enddo
NA1N = ND1
else
NE(1: NA1) = NEO(1: NA1)
do J=1, NA
QN(J) = G11(J)*(-YWA(J)*(NE(J+1) - NE(J))/HRO - 0.5*YWB(J)*(NE(J+1)+NE(J))) + SLAT(J)*GNX(J)
GN(J) = QN(J)/SLAT(J)
enddo
QN(NA1) = QN(NA)
NA1N = NA1
endif
"""

    assigned = \
"""
do J=1, NA
QN(J) = G11(J)*(-YWA(J)*(NE(J+1) - NE(J))/HRO - 0.5*YWB(J)*(NE(J+1)+NE(J))) + SLAT(J)*GNX(J)
GN(J) = QN(J)/SLAT(J)
enddo
QN(NA1) = QN(NA)
NA1N = NA1
"""


class UPEQN:

    eqn = \
'''QU(1) = HRO
YWD(ND1) = 0.
do j=1, NA1
YWD(j)  = TTRQI(j)
YWGN(j) = VR(j)
YWGO(j) = VRO(j)
YWHN(j) = UPS0(j)
YWHO(j) = UPS0O(j)
YVR(j)  = VR(j)
MPHIT(j) = UPARO(j)
unit_coeff = 1.
YWG11(j) = G11(j)
! missing UPS1 and UPS2 terms in <M_phi>
YWWB(j) = VR(j)
if (IPROT >= 1) then
! contributions to net torque
MPHIT(j) = MPHIT(j) + UPS1O(j)/UPS0O(j)
MPHIT(j) = MPHIT(j) + UPS2O(j)/UPS0O(j)
TTRQ(j) = TTRQ(j) - ( VR(J)*(UPS1(J) + UPS2(J)) - VRO(J)*(UPS1O(J) + UPS2O(J)) )*2./ (VR(J) + VRO(J))/TAU
! contributions to stress tensor
YWgradF(J)  = 2.*(IPOL(J+1)  - IPOL(J)) /(IPOL(J+1)  + IPOL(J)) /HRO  !d log I / drho
YWgradb2(J) = 2.*(BDB02(J+1) - BDB02(J))/(BDB02(J+1) + BDB02(J))/HRO  !d log <B**2>/drho
YWR(J) = RTOR/IPOL(J)*XUPAR(J)*DLNEOD(J) - RTOR/IPOL(J)*(CNPAR(J) + XUPAR(J)*YWgradF(J))*DLNEO(J) + (XUPAR(J) - XUPAP(J))*RTOR/IPOL(J)*BDB02(J)*BTOR*SGNEOD(J) + RTOR/IPOL(J)*BDB02(J)*BTOR*(XUPAR(J)*YWgradb2(J) - XUPAR(J)*YWgradF(J) + CNPAP(J) - CNPAR(J))*SGNEO(J) + RTOR/IPOL(J)*(CNPAD(J) - CNPAR(J) - XUPAD(J)*YWgradF(J))*DDNEO(J) + RTOR/IPOL(J)*(XUPAR(J) - XUPAD(J))*DDNEOD(J)
endif
enddo
imethod = INUME4

call RUNEQ(YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), UPARO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWD(1: NA1), TTRQ(1: NA1), RBDOT, BBDOT, ND1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, UPAR(1: NA1), QU(1: NA1), YQDCM(1: NA1), MPHIT(1: NA1))

MPHIT = 0.
NA1U = ND1
'''

    assigned = \
'''do j=1, NA
QU(J) = -G11(J)*(YWA(J)*(UPAR(J+1) - UPAR(J))/HRO + 0.5*YWB(J)*(UPAR(J+1) + UPAR(J)))*0.0016
enddo
QU(NA1) = QU(NA)
'''

    uparo = \
'''UPARO(ND1: NA1) = UPAR(ND1: NA1)
QU(4)  = 1.
bctype = 1
'''


class TETIEQN:

    eqn = \
'''
QE(1)  = HRO
NA1E = ND1

do j=1, NA1
YWGN(j) = VR(j)**(5./3.)
YWGO(j) = VRO(j)**(5./3.)
unit_coeff = 625.
YVR(j) = VR(j)
enddo
imethod = INUME2
'''

    runeq = \
'''NA1I = ND1
call RUNEQ_TETI(YWGN(1:NA1), 3./2.*NE(1:NA1), 3./2.*NI(1:NA1), YWGO(1:NA1), 3./2.*NEO(1:NA1), 3./2.*NIO(1:NA1), TEO(1:NA1), TIO(1:NA1), YVR(1:NA1), unit_coeff, G11(1:NA1)/unit_coeff, YWA1(1:NA1), YWA2(1:NA1), YWB1(1:NA1), YWB2(1:NA1), unit_coeff*PET(1:NA1), unit_coeff*PIT(1:NA1), unit_coeff*PETOT(1:NA1), unit_coeff*PITOT(1:NA1), ND1, NA1, HRO, TAU, RHO(1:NA1), imethod, TE(1:NA1), TI(1:NA1), QE(1:NA1), QI(1:NA1))
if (ND1 < NA1) then
do j=ND1+1, NA1
QE(j) = QE(ND1)
QI(j) = QI(ND1)
enddo
endif
PETOT=PE+PET*TE
PITOT=PI+PIT*TI
'''

    teold = \
'''TEO(ND1: NA1) = TE(ND1: NA1)
QE(4)   = 1.
bc_type_imp(1) = 1
'''

    tiold = \
'''TIO(ND1: NA1) = TI(ND1: NA1)
QI(4)   = 1.
bc_type_imp(2) = 1
'''


class RADOUT:

    header = \
"""subroutine RADOUT

!------------------------------------------------------------
! Radial profile plotting
!------------------------------------------------------------

use pi_const, only: GP, GP2, mu0
use scalars
use status
use graph_utils
use standard_functions
use debugger, only: markloc, debug

implicit none

include 'src/tmp/declar.fml'
include 'src/tmp/declar.fnc'

call markloc('radout')
do J=1, NAB
"""


class TIMOUT:

    header = \
"""
enddo

end subroutine RADOUT

!------------------------------------------------------------
subroutine TIMOUT

!------------------------------------------------------------
! Time traces plotting
!------------------------------------------------------------

use pi_const, only: GP, GP2, mu0
use scalars
use status
use graph_utils
use standard_functions
use debugger, only: markloc, debug

implicit none

include 'src/tmp/declar.fml'
include 'src/tmp/declar.fnc'

call markloc('timout')
"""

class FJEQN:

    eqn = \
'''
do j=1, NA1
YWHN(j) = 1.
YWHO(j) = 1.
YWGN(j) = VR(j)
YWGO(j) = VRO(j)
YWR(j)  = 0.
YVR(j)  = VR(j)
unit_coeff  = 1.
YWWB(j) = VR(j)
enddo
imethod = INUME1
'''

class CUEQN:

    header = \
'''! **** Current equation
call markloc("Current equation")
YC =  mu0*RTOR/BTOR
YD = -0.8*GP**2 * RTOR
YA = 2./(HRO**2 * YD)
do J=1, NA1
'''

    mv = \
'''
'''

    prescribed_ipl = \
'''! Prescribed plasma current:
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = 1.
bc_values(4) = -1.
bc_values(5) = HRO*mu0/G22(NA)*IPL*RTOR/IPOL(NA1)
bc_type_for_fp=1
! For psifb
! when using the free boundary circuit equations with free current,
! then use mixed b.c.
if (ITFBP /= 0 .and. ITFBE < TIME) then
if (IBCPSI > 0) then
!case implicit
if (ITFBP < 0 .and. IBCPSI >= 2) then
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) =  PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (IBCPSI <= 1) then
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = 1.
bc_values(4) = -1.
bc_values(5) = HRO*mu0/G22(NA)*IPL*RTOR/IPOL(NA1)
bc_type_for_fp = 1
endif
endif
'''

    prescribed_uloop = \
'''! Prescribed loop voltage:
FP(NA1) = FPO(NA1) + TAU*UEXT
bctype = 1
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_type_for_fp = 2
!For psifb
! when using the free boundary circuit equations with free current,
! then use mixed b.c.
if (ITFBP /= 0 .and. ITFBE < TIME) then
if (IBCPSI > 0) then
!case implicit
if (ITFBP < 0 .and. IBCPSI >= 2) then
bc_type = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) =  PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (IBCPSI <=  1) then
bctype = 1
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_type_for_fp = 2
endif
endif
'''

    circuit_eqn = \
'''! Circuit equation:
if (TIME - TSTART <= TAU) then
PSIEXT = FP(NA1) + LEXT*IPL
PSPLEX = LEXT/ROC*5.*IPOL(NA1)*G22(NA)/GP2/RTOR
endif
PSPLEX = LEXT/ROC*5.*IPOL(NA1)*G22(NA)/GP2/RTOR
PSIEXT = PSIEXT + TAU*UEXT
bctype = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) =  -PSPLEX*ROC
bc_values(5) = PSIEXT*HRO
bc_type_for_fp = 3
!For psifb
! when using the free boundary circuit equations with free current,
! then use mixed b.c.
if (ITFBP /= 0 .and. ITFBE < TIME) then
if (IBCPSI > 0) then
! case implicit
if (ITFBP < 0 .and. IBCPSI >= 2) then
bc_type = 3
bc_values(1) = 0.0
bc_values(2) = 0.0
bc_values(3) = HRO + PSPLEX*ROC
bc_values(4) = -PSPLEX*ROC
bc_values(5) = PSIEXT*HRO
bc_type_for_fp = 3
endif
endif
if (IBCPSI <= 1) then
bc_type = 1
bc_values(1) = 0.0
bc_values(2) = 0.0
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
YWWB(j) = 1./RHO(j)
YVR(j) = CC(j)*mu0*RHO(j)/IPOL(j)**2
unit_coeff = 1.
YWD(j) = -(VR(j)/(GP2*RHO(j)*CC(j))) * (CUBS(j) + CD(j))
enddo
imethod = INUME3

call RUNEQ( YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), FPO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, YWG11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), YWS(1: NA1), YWD(1: NA1), RBDOT, BBDOT, NA1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, FP(1: NA1), YWQ(1: NA1), YQDCMF(1: NA1), MPHIT(1: NA1) )

dfpdrbm12 = -YWQ(NA)/G22(NA)
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
call markloc("init done")
'''

class INIT_CONVERGE_STEP:

    header = \
'''subroutine init_converge_step

use pi_const, only: GP, GP2, mu0
use read_input, only: equ_file, exp_file
use cpu_usage, only: wallTime_sbr, cpuTime_sbr
use scalars
use status
use nclass_mod
use a2tglf, only: tglf_ipc, tglf_out
use a2qlk, only: qlk_ipc, qlk_out
use debugger, only: markloc
use numerical_tools, only: extrap
use metrics, only: cuofmu

implicit none

include 'src/tmp/declar.fml'
include 'src/tmp/declar.fnc'

integer :: IFSUB, t_wall1, t_wall2, rate
double precision :: YB, YC, YU, YJ_CU, YM, YMCD, YIOH, YICD, YM1, t_cpu1, t_cpu2
double precision, dimension(NRD) :: YWA
'''

class EQNS_INC:

    header = \
'''subroutine EQNS_INC(ibcpsi_fb, bc_type_for_fp, dfpdrbm12)
!-------------------------------------------------------------------
! Perform one time step
! Note that now time step is updated at the end of a full time cycle
!-------------------------------------------------------------------

use pi_const, only: GP, GP2, mu0
use scalars
use status
use a2tglf, only: tglf_ipc, tglf_out
use a2qlk, only: qlk_ipc, qlk_out
use cpu_usage, only: wallTime_sbr, cpuTime_sbr
use nclass_mod
use strahl_mod
use standard_functions
use debugger, only: markloc
use numerical_tools, only: extrap
use transport_solver
use metrics, only: cuofp, cuofmu

implicit none

include 'src/tmp/declar.fml'
include 'src/tmp/declar.fnc'

integer, intent(in) :: ibcpsi_fb
integer, intent(out) :: bc_type_for_fp
double precision, intent(out) :: dfpdrbm12

integer :: IFSUB, imethod, ND, ND1, JCALL, bctype, bc_type_imp(2), t_wall1, t_wall2, rate

double precision :: YHRO, YM1, YM2, YB, YC, YJ_CU, YM, YU, YIOH, YICD, YMCD, bc_value_imp(2), t_cpu1, t_cpu2, unit_coeff
double precision, dimension(5) :: bc_values
double precision, dimension(NRD) :: YWA, YWB, YWC, YWD, YWGN, &
    YWHN, YWGO, YWHO, YWR, YVR, YWA1, YWA2, YWB1, YWB2, YWAA, YWWB, &
    YWC1, YWC2, YWS, YQDCM, MPHIT, YQDCMF, YWQ, YWG11, YWgradF, YWgradb2

MPHIT = 0.
'''
