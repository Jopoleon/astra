! QNTOT [10#19/s]:  Integral {0, R} (SNTOT) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QNTOTR(YR)

use status_inc, only: SNTOT

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QNTOTR = VINT(SNTOT, YR)

end function QNTOTR
