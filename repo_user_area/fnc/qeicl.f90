! QEICL [MW]:  Integral {0, R} ( PEICL ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QEICLR(YR)

use scalars, only: HRO
use status, only: VR, TE, TI, NE, NI, AMAIN, ZMAIN
use standard_functions, only: jrho_drho

implicit none

double precision, intent(in) :: YR
integer :: J, JK
double precision PEICL, COULG, YDR

call jrho_drho(YR, JK, YDR)
QEICLR = 0.
DO J=1, JK
   include  'fml/peicl'
   QEICLR = QEICLR + PEICL*VR(J)
enddo
QEICLR = HRO*(QEICLR - PEICL*YDR)

end function QEICLR
