! BETA % Beta toroidal  % 
!
! SI:   BETA=2\mu_0*p/B^2 = 12.8e-3*pi*(n_20)*(T_keV)/(B_T)^2 * 100[%]
! CGS:  BETA=8\pi*p/B^2   = 12.8e-3*pi*[0.1*NE*TE+...]/BTOR^2 * 100[%]
!
!    NOTE! NBI Pressure is included as in equilibrium
!            (Pereverzev 02-MAY-2006)

double precision function BETAR(YR)

use const_inc, only: ROC, HRO, NA1, BTOR
use status_inc, only: VR, TE, TI, NE, NI, PBLON, PBPER

implicit  none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: YV, YDR

call yrjkdr(YR, JK, YDR)

YV = 0.
BETAR = 0.
do J=1, JK
   YV = YV + VR(J)
   BETAR = BETAR + VR(j)*(TE(J)*NE(J) + TI(J)*NI(J) + .5*(PBLON(J) + PBPER(J)))
enddo
YV = YV - YDR
BETAR = BETAR - YDR*(TE(JK)*NE(JK) + TI(JK)*NI(JK) + .5*(PBLON(JK) + PBPER(JK)))
BETAR = 0.402*BETAR/(YV*BTOR*BTOR)

return
end function betar
