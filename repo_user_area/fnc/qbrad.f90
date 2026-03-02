! QBRAD [MW]:  Integral {0, R} ( PBRAD ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QBRADR(YR)

use const_inc, only: HRO
use status_inc, only: VR, NE, TE, ZEF
use standard_functions, only: jrho_drho

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision :: PBRAD, YDR

call jrho_drho(YR, JK, YDR)

QBRADR = 0.
DO J=1, JK
   include 'fml/pbrad'
   QBRADR = QBRADR + PBRAD*VR(J)
enddo
QBRADR = HRO*(QBRADR - PBRAD*YDR)

end function QBRADR
