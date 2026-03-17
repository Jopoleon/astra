! QRADX [MW]:  Integral {0, R} (PRADX) dV
!   (Yushmanov 11-JAN-89)
double precision FUNCTION QRADXR(YR)

use status, only: PRADX
use standard_functions, only: VINT

implicit none

double precision, intent(in) :: YR

QRADXR = VINT(PRADX, YR)

end function QRADXR
