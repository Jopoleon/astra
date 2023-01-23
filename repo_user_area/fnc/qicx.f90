! QICX [MW]:  Integral {0, R} ( PICX ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QICXR(YR)

use const_inc, only: HRO, NNWM, NNCL
use status_inc, only: VR, TE, TN, TI, NE, NN, NI

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: PICX, SVCX, YDR

call yrjkdr(YR, JK, YDR)

QICXR = 0.
DO J=1, JK
   include 'fml/picx'
   QICXR = QICXR + PICX*VR(J)
enddo
QICXR = HRO*(QICXR - PICX*YDR)

end function QICXR
