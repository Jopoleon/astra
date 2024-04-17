! BETA N: Beta toroidal normalized thermal
!
! SI: BETA = 2\mu_0*p/B^2 = 12.8e - 3*pi*(n_20)*(T_keV)/(B_T)^2 * 100[%]
! CGS: BETA = 8\pi*p/B^2   = 12.8e - 3*pi*[0.1*NE*TE + ...]/BTOR^2 * 100[%]
!
!    NOTE! NBI Pressure is included as in equilibrium
!     (Pereverzev 02 - MAY - 2006) (Fable, Apr 2024)
double precision function BETATR(YR)

use const_inc, only: BTOR, IPL, ABC, NA1, HRO, ROC
use status_inc, only: TE, TI, NE, NI, PBLON, PBPER, PFAST, VR

implicit  none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: YV, YDR

call yrjkdr(YR, JK, YDR)

YV = 0.
BETATR = 0.
do J = 1, JK
   YV = YV + VR(J)
   BETATR = BETATR + VR(j)*(TE(J)*NE(J) + TI(J)*NI(J))
enddo
YV = YV - YDR
BETATR = BETATR - YDR*(TE(JK)*NE(JK) + TI(JK)*NI(JK))
BETATR = 0.402*BETATR*ABC/(YV*BTOR*IPL)

return
end function BETATR
