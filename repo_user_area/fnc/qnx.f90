! QNX [10#19/s]:  Integral {0, R} (SNX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QNXR(YR)

use status_inc, only: SNX

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QNXR = VINT(SNX, YR)

end function QNXR
