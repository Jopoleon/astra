! QIX [MW]:  Integral {0, R} (PIX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QIXR(YR)

use status_inc, only: PIX

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QIXR = VINT(PIX, YR)

end function QIXR
