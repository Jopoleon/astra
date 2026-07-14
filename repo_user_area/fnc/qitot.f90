! QITOT [MW]:  Integral {0, R} (PITOT) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QITOTR(YR)

use status, only: PITOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QITOTR = VINT(PITOT, YR)

end function QITOTR
