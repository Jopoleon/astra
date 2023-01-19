! VOL [m**3]: Volume(r)
! Integral {0:R} dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION VOLR(YR)

use status_inc, only: VR

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: YV, YDR

call yrjkdr(YR, JK, YDR)

YV = 0.
DO J=1, JK
   YV = YV + VR(J)
enddo
YV = YV - YDR

return
end function VOLR
