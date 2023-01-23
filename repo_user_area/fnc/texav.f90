! TEXAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TEX/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TEXAVR(YR)

use status_inc, only: TEX

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VOLR

TEXAVR = VINT(TEX, YR)/VOLR(YR)

return
end function TEXAVR
