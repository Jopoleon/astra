! TIXAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TIX/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TIXAVR(YR)

use status, only: TIX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
double precision, external :: VOLR

TIXAVR = VINT(TIX, YR)/VOLR(YR)

end function TIXAVR
