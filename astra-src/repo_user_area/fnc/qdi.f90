! QDI [MW]:  Integral {0, R} ( PDI ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QDIR(YR)

use status, only: PDI
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QDIR = VINT(PDI, YR)

end function QDIR
