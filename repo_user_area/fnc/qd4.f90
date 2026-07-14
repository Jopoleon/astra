! QD4 [F9/s]:  Integral {0, R} ( SD4 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD4R(YR)

use status, only: SD4
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD4R = VINT(SD4, YR)

end function QD4R
