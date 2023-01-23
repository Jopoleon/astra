! QDE [MW]:  Integral {0, R} ( PDE ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QDER(YR)

use status_inc, only: PDE

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QDER = VINT(PDE, YR)

end function QDER
