! QTOT [MW]:  Integral {0, R} ( PETOT+PITOT ) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QTOTR(YR)

use status, only: PE, PI
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QTOTR = VINT(PE + PI, YR)

end function QTOTR
