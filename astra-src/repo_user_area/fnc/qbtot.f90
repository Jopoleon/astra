! QBTOT [MW]:  Integral {0, R} (PBEAM) dV
!   (Polevoy 28.09.89)
double precision FUNCTION QBTOTR(YR)

use status, only: PBEAM
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QBTOTR = VINT(PBEAM, YR)

end function QBTOTR
