! QD5 [F9/s]:  Integral {0, R} ( SD5 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD5R(YR)

use status, only: SD5
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD5R = VINT(SD5, YR)

end function QD5R
