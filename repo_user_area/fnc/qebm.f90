! QEBM [MW]:  Integral {0, R} (PEBM) dV
!   (Pereverzev 20-MAY-08)
double precision FUNCTION QEBMR(YR)

use status, only: PEBM
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QEBMR = VINT(PEBM, YR)

end function QEBMR
