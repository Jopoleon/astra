! TIAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TI/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TIAVR(YR)

use status, only: TI
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
double precision, external :: VOLR

TIAVR = VINT(TI, YR)/VOLR(YR)

end function TIAVR
