! QD8 [F9/s]:  Integral {0, R} ( SD8 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD8R(YR)

use status_inc, only: SD8

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD8R = VINT(SD8, YR)

end function QD8R
