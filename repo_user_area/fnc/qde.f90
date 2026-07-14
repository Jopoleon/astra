! QDE [MW]:  Integral {0, R} ( PDE ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QDER(YR)

use status, only: PDE
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QDER = VINT(PDE, YR)

end function QDER
