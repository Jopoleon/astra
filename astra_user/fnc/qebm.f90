! QEBM [MW]:  Integral {0, R} (PEBM) dV
!   (Pereverzev 20-MAY-08)
double precision FUNCTION QEBMR(YR)

use status_inc, only: PEBM

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QEBMR = VINT(PEBM, YR)

end function QEBMR
