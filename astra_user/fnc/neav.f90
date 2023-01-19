! NEAV [10#19/m#3]: Volume average density (r)
! Integral {0:R} (NE/VOL) dV
!   (Yushmanov 11-MAY-87)
double precision FUNCTION NEAVR(YR)

use status_inc, only: NE

implicit none

double precision, intent(in) :: YR
double precision :: VINT, VOLR

NEAVR = VINT(NE, YR)/VOLR(YR)

return
end function NEAVR
