! QD2 [F9/s]:  Integral {0, R} ( SD2 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD2R(YR)

use status, only: SD2
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD2R = VINT(SD2, YR)

end function QD2R
