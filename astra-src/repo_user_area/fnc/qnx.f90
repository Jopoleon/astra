! QNX [10#19/s]:  Integral {0, R} (SNX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QNXR(YR)

use status, only: SNX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QNXR = VINT(SNX, YR)

end function QNXR
