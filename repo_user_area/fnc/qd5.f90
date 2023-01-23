! QD5 [F9/s]:  Integral {0, R} ( SD5 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD5R(YR)

use status_inc, only: SD5

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD5R = VINT(SD5, YR)

end function QD5R
