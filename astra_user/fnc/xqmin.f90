! XQMIN []: Position of minimal "q" value on 0<=a/ABC<=1
!   (Pereverzev 13-JAN-97)
! Usage: xmin_XQMINB;
!  CV1=XQMINB; HE=HE*XSTEP(CV1);
!
double precision FUNCTION XQMINR(YR)

use const_inc, only: ABC, NA1
use status_inc, only: AMETR, MU

implicit none

double precision :: YR, YMIN
integer :: J

YMIN = 0.
XQMINR = 0.
do J=1, NA1
   if (MU(j) < YMIN) CYCLE
   YMIN = MU(j)
   if (j /= 1) XQMINR = AMETR(j)/ABC
enddo

return
end function XQMINR
