! LI3 []: Internal inductance li(r) [4]
!   Pereverzev 03-APR-08
!   /
!       c*c |
!     Li = ------------*|(Bpol**2)dV   li = Li/(2*pi*R0)
!    4*pi*Ipl*Ipl |
!   /
!   V
!
!   SI: Li[H] = 2*Wi/Ipl**2 ; li[dim.less] = Li[mkHn]/(0.2*pi*R0[m])
!          2*pi*R0[m]*li[dim.less] = 10*Li[mkHn]
!----------------------------------------------------------------------|
! The same as LINT except for calculation of the denominator 
!     as an integral of the current density
!----------------------------------------------------------------------|

double precision function LI3R(YRO)

use pi_const, only: GP
use scalars, only: ROC, HRO, BTOR, NA
use status, only: IPOL, G22, CU, VR, MU

implicit none

double precision, intent(in) :: YRO
double precision :: YK, YR, YIPL, YWBP
integer :: J, J1, JK

if (YRO <= HRO) then
   YK = 1.
elseif (YRO > ROC) then
   YK = ROC/HRO
else
   YK = YRO/HRO
endif
JK = YK
YR = JK + 1.
if (JK >= NA) then
   JK = NA
   YR = JK + 1.
endif
YWBP = (YK*YK - JK*JK)*(YK*YK + JK*JK)*IPOL(JK+1)*MU(JK+1)**2  *G22(JK+1)/(JK+1)
YIPL = 0.
do J=1, JK
   J1  = 2*J - 1
   YWBP = YWBP + J1*(2.*J*J - J1)*IPOL(J)*G22(J)/J*MU(J)**2
   YIPL = YIPL + 0.5*(CU(j)*VR(j) + CU(j+1)*VR(j+1))/IPOL(j)**2
enddo

LI3R = 2.*HRO*YWBP*(5.*GP*BTOR/IPOL(JK+1)/YIPL)**2

end function LI3R
