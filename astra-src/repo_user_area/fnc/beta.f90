! BETA % Beta toroidal  % 
!
! SI:   BETA=2\mu_0*p/B^2 = 12.8e-3*pi*(n_20)*(T_keV)/(B_T)^2 * 100[%]
! CGS:  BETA=8\pi*p/B^2   = 12.8e-3*pi*[0.1*NE*TE+...]/BTOR^2 * 100[%]
!
!    NOTE! NBI Pressure is included as in equilibrium
!            (Pereverzev 02-MAY-2006)

double precision function BETAR(YR)

use scalars, only: ROC, HRO, NA1, BTOR
use status, only: VR, TE, TI, NE, NI, PBLON, PBPER, PFAST
use standard_functions, only: jrho_drho

implicit  none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: YV, YDR

call jrho_drho(YR, JK, YDR)

YV = 0.
BETAR = 0.
do J=1, JK
    YV = YV + VR(J)
    BETAR = BETAR + VR(j)*(TE(J)*NE(J) + TI(J)*NI(J) + 0.5*(PBLON(J) + PBPER(J)) + PFAST(J))
enddo
YV = YV - YDR
BETAR = BETAR - YDR*(TE(JK)*NE(JK) + TI(JK)*NI(JK) + 0.5*(PBLON(JK) + PBPER(JK)) + PFAST(J))
BETAR = 0.402*BETAR/(YV*BTOR*BTOR)

end function betar
