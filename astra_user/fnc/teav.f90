! TEAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TE/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TEAVR(YR)

use status_inc, only: TE

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VOLR

TEAVR = VINT(TE, YR)/VOLR(YR)

return
end function TEAVR
