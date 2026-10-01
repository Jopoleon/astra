! QINEU [MW]:  Integral {0, R} ( PINEU ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QINEUR(YR)

use scalars, only: HRO, NNCL, NNWM
use status, only: VR, NE, NN, NI, TE, TN, TI
use standard_functions, only: jrho_drho

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: PINEU, YDR, PICX, SVCX, SVIE, SVII

call jrho_drho(YR, JK, YDR)
QINEUR = 0.
DO J=1, JK
   include 'fml/pineu'
   QINEUR = QINEUR + PINEU*VR(J)
enddo
QINEUR = HRO*(QINEUR - PINEU*YDR)

end function QINEUR
