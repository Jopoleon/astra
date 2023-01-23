! QD4 [F9/s]:  Integral {0, R} ( SD4 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD4R(YR)

use status_inc, only: SD4

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD4R = VINT(SD4, YR)

end function QD4R
