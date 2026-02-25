! QOH [MW]:  Integral {0, R} (POH) dV
!   (Pereverzev 12-FEB-90)
double precision FUNCTION QOHR(YR)

use const_inc, only: RTOR, GP2, NA1
use status_inc, only: ULON, IPOL, CC, G33
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: POH(NA1)


do J=1, NA1
    POH(J) = CC(J)*(ULON(J)/(GP2*RTOR*IPOL(J)))**2/G33(J)
enddo

QOHR = VINT(POH, YR)

end function QOHR
