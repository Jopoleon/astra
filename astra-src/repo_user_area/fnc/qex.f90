! QEX [MW]:  Integral {0, R} (PEX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QEXR(YR)

use status, only: PEX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QEXR = VINT(PEX, YR)

end function QEXR
