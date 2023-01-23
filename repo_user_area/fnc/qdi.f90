! QDI [MW]:  Integral {0, R} ( PDI ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QDIR(YR)

use status_inc, only: PDI

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QDIR = VINT(PDI, YR)

end function QDIR
