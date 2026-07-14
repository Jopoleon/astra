! TEXAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TEX/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TEXAVR(YR)

use status, only: TEX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
double precision, external :: VOLR

TEXAVR = VINT(TEX, YR)/VOLR(YR)

end function TEXAVR
