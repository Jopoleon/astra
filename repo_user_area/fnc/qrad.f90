! QRAD [MW]:  Integral {0, R} (PRAD) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QRADR(YR)

use status, only: PRAD
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QRADR = VINT(PRAD, YR)

end function QRADR
