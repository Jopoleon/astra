! VOL [m**3]: Volume(r)
! Integral {0:R} dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION VOLR(YR)

use const_inc, only: HRO
use status_inc, only: VR

implicit none

double precision, intent(in) :: YR
integer :: JK
double precision :: YDR

call yrjkdr(YR, JK, YDR)
VOLR = SUM(VR(1:JK))*HRO - YDR*VR(JK)

return
end function VOLR
