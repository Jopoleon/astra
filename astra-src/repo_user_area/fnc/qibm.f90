! QIBM [MW]:  Integral {0, R} (PIBM) dV
!   (Pereverzev 20-MAY-08)
double precision FUNCTION QIBMR(YR)

use status, only: PIBM
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QIBMR = VINT(PIBM, YR)

end function QIBMR
