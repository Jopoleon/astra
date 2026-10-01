! QNTOT [10#19/s]:  Integral {0, R} (SNTOT) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QNTOTR(YR)

use status, only: SNTOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QNTOTR = VINT(SNTOT, YR)

end function QNTOTR
