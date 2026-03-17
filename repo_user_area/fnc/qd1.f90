! QD1 [F9/s]:  Integral {0, R} ( SD1 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD1R(YR)

use status, only: SD1
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD1R = VINT(SD1, YR)

end function QD1R
