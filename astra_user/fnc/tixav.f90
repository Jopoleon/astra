! TIXAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TIX/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TIXAVR(YR)

use status_inc, only: TIX

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VOLR

TIXAVR = VINT(TIX, YR)/VOLR(YR)

return
end function TIXAVR
