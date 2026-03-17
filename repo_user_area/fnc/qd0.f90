! QD0 [F9/s]:  Integral {0, R} ( SD0 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD0R(YR)

use status, only: SD0
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QD0R = VINT(SD0, YR)

end function QD0R
