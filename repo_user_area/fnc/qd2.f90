! QD2 [F9/s]:  Integral {0, R} ( SD2 ) dV
!   (Pereverzev 9-OCT-07)
double precision FUNCTION QD2R(YR)

use status_inc, only: SD2

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QD2R = VINT(SD2, YR)

end function QD2R
