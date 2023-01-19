! QD6 [F9/s]:  Integral {0, R} ( SD6 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD6R(YR)

use status_inc, only: SD6

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD6R = VINT(SD6, YR)

end function QD6R
