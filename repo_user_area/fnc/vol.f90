! VOL [m**3]: Volume(r)
! Integral {0:R} dV
!   (Perverzev 13-MARCH-91)

double precision FUNCTION VOLR(YR)

use scalars, only: HRO
use status, only: VR
use standard_functions, only: jrho_drho

implicit none

double precision, intent(in) :: YR
integer :: JK
double precision :: YDR

call jrho_drho(YR, JK, YDR)
VOLR = SUM(VR(1:JK))*HRO - YDR*VR(JK)

end function VOLR
