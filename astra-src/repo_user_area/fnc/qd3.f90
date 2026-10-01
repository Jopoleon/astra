! QD3 [F9/s]:  Integral {0, R} ( SD3 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD3R(YR)

use status, only: SD3
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD3R = VINT(SD3, YR)

end function QD3R
