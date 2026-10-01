! QD9 [F9/s]:  Integral {0, R} ( SD9 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD9R(YR)

use status, only: SD9
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD9R = VINT(SD9, YR)

end function QD9R
