! TEAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TE/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TEAVR(YR)

use status, only: TE
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
double precision :: VOLR

TEAVR = VINT(TE, YR)/VOLR(YR)

end function TEAVR
