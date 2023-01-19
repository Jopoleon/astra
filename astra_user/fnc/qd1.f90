! QD1 [F9/s]:  Integral {0, R} ( SD1 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD1R(YR)

use status_inc, only: SD1

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD1R = VINT(SD1, YR)

end function QD1R
