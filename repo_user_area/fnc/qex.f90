! QEX [MW]:  Integral {0, R} (PEX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QEXR(YR)

use status_inc, only: PEX

implicit none

double precision, intent(in) :: YR
double precision :: VINT

QEXR = VINT(PEX, YR)

end function QEXR
