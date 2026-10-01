! QIX [MW]:  Integral {0, R} (PIX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QIXR(YR)

use status, only: PIX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QIXR = VINT(PIX, YR)

end function QIXR
