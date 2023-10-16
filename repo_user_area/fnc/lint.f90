! LINT []: Internal inductance li(r) []
!      Pereverzev 27-DEC-90
!
!     Li = ------------*|(Bpol**2)dV
!     li = Li/(2*pi*R0)
!          4*pi*Ipl*Ipl |
!  SI: Li[H] = 2*Wi/Ipl**2 ;  li[dim.less] = Li[mkHn]/(0.2*pi*R0[m])
!       2*pi*R0[m]*li[dim.less] = 10*Li[mkHn]

double precision function LINTR(YRO)

use const_inc, only: NA, HRO, ROC
use status_inc, only: IPOL, G22, MU

implicit none

double precision, intent(in) :: YRO
integer :: J, J1, JK
double precision :: YK, YR

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
LINTR = 0.
do J=1, JK
    J1  = J + J - 1
    LINTR = LINTR + J1*(2.*J**2 - J1)*IPOL(J)*G22(J)/J*MU(J)**2
enddo
LINTR = LINTR + (YK**2 - JK**2)*(YK**2 + JK**2)*MU(JK+1)**2 * IPOL(JK+1)*G22(JK+1)/YR
LINTR = 0.5*LINTR*HRO*(YR/(YK**2 * MU(JK+1)*IPOL(JK+1)*G22(JK+1)))**2

return
end function LINTR
