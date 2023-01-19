! NEXAV [10#19/m#3]: Volume average density (r)
! Integral {0:R} (NEX/VOL) dV
!   (Yushmanov 11-MAY-87)
double precision FUNCTION NEXAVR(YR)

use status_inc, only: NEX

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VOLR

NEXAVR = VINT(NEX, YR)/VOLR(YR)

return
end function NEXAVR
