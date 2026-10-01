! NEXAV [10#19/m#3]: Volume average density (r)
! Integral {0:R} (NEX/VOL) dV
!   (Yushmanov 11-MAY-87)
double precision FUNCTION NEXAVR(YR)

use status, only: NEX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
double precision, external :: VOLR

NEXAVR = VINT(NEX, YR)/VOLR(YR)

end function NEXAVR
