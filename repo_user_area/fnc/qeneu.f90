! QENEU [MW]:  Integral {0, R} ( PENEU ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QENEUR(YR)

use const_inc, only: HRO, NNCL, NNWM
use status_inc, only: VR, NE, NN, NI, TE, TN, TI
use standard_functions, only: jrho_drho

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: PENEU, YDR, PECX, SVCX, SVIE, SVII

call jrho_drho(YR, JK, YDR)
QENEUR = 0.
DO J=1, JK
   include 'fml/peneu'
   QENEUR = QENEUR + PENEU*VR(J)
enddo
QENEUR = HRO*(QENEUR - PENEU*YDR)

end function QENEUR
