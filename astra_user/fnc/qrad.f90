! QRAD [MW]:  Integral {0, R} (PRAD) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QRADR(YR)

use status_inc, only: PRAD

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QRADR = VINT(PRAD, YR)

end function QRADR
