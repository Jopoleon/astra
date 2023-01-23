! TIAV [keV]: Volume average ion temperature (r)
! Integral {0:R} (TI/VOL) dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION TIAVR(YR)

use status_inc, only: TI

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VOLR

TIAVR = VINT(TI, YR)/VOLR(YR)

return
end function TIAVR
