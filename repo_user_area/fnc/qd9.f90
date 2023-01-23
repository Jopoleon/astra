! QD9 [F9/s]:  Integral {0, R} ( SD9 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD9R(YR)

use status_inc, only: SD9

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD9R = VINT(SD9, YR)

end function QD9R
