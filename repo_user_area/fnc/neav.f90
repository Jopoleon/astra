! NEAV [10#19/m#3]: Volume average density (r)
! Integral {0:R} (NE/VOL) dV
!   (Yushmanov 11-MAY-87)
double precision FUNCTION NEAVR(YR)

use status, only: NE, VOLUM
use scalars, only: HRO
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR
double precision, external :: VOLR

NEAVR = VINT(NE, YR)/VOLR(YR)

end function NEAVR
