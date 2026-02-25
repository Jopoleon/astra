! QETOT [MW]:  Integral {0, R} (PETOT) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QETOTR(YR)

use status_inc, only: PETOT
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QETOTR = VINT(PETOT, YR)

end function QETOTR
