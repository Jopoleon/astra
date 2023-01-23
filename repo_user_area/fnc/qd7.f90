! QD7 [F9/s]:  Integral {0, R} ( SD7 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD7R(YR)

use status_inc, only: SD7

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD7R = VINT(SD7, YR)

end function QD7R
