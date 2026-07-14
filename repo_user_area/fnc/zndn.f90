! ZNDN []: Volume-density averaged Zef (r)
!   <Zef*Ne*dV>/<Ne*dV>
!     (Polevoy 28.09.89)
double precision function ZNDNR(YR)

use status, only: TE, NE, VR, ZEF
use standard_functions, only: jrho_drho

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: YQ, YV, YDR

call jrho_drho(YR, JK, YDR)

YQ = 0.
YV = 0.
ZNDNR = ZEF(1)
DO J=1, JK
   YV = YV + NE(J)*VR(J)
   YQ = YQ + ZEF(J)*NE(J)*VR(J)
enddo

YV = YV - YDR*NE(JK)
YQ = YQ - ZEF(JK)*YDR*NE(JK)
ZNDNR = YQ/YV

end function ZNDNR
