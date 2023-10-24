! NEAV [10#19/m#3]: Volume average density (r)
! Integral {0:R} (NE/VOL) dV
!   (Yushmanov 11-MAY-87)
double precision FUNCTION NEAVR(YR)

use status_inc, only: NE, VOLUM
use const_inc, only: HRO

implicit none

double precision, intent(in) :: YR
integer :: i
double precision, external :: VINT

i = int(YR/HRO) + 1

NEAVR = VINT(NE, YR)/VOLUM(i) !VOLR(YR)

return
end function NEAVR
